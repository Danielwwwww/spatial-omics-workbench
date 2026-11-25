"""Project adult mouse CCFv3 coronal labels onto a section in H&E coordinates."""
import base64
import gzip
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from scipy.optimize import minimize


SOURCE = "https://alleninstitute.github.io/CCF-MAP/descriptions/mouse_ccf.html"
HIPPOCAMPAL = ("海马 CA1", "海马 CA2", "海马 CA3", "齿状回")


@lru_cache(maxsize=1)
def reference():
    path = Path(__file__).with_name("allen-ccfv3.json")
    payload = json.loads(path.read_text(encoding="utf-8"))
    labels = np.frombuffer(gzip.decompress(base64.b64decode(payload["labels"])),dtype=np.uint8)
    return payload["names"], tuple(payload["ap_voxels"]), labels.reshape(payload["shape"])


def resize(raster, shape):
    return np.asarray(Image.fromarray(raster).resize(shape[::-1], Image.Resampling.NEAREST))


def crop_bounds(mask):
    y,x = np.where(mask)
    return np.s_[y.min():y.max()+1,x.min():x.max()+1]


def dice(a,b):
    return 2 * np.count_nonzero(a & b) / max(1,np.count_nonzero(a)+np.count_nonzero(b))


def corrected(layer, params):
    angle, sx, sy, dy, dx = params
    c,s = np.cos(angle),np.sin(angle)
    inverse = np.linalg.inv(np.array([[c,-s],[s,c]]) @ np.diag([sy,sx]))
    center = (np.array(layer.shape)-1)/2
    return ndi.affine_transform(layer,inverse,center-inverse@(center+[dy,dx]),order=0,prefilter=False)


def project(tissue, settings, hippocampus=None):
    names, levels, planes = reference()
    orientation = settings.get("atlas_rotation", "auto")
    turns = range(4) if orientation == "auto" else [int(orientation)//90]
    if orientation != "auto" and str(orientation) not in ("0","90","180","270"):
        raise ValueError("atlas_rotation must be auto, 0, 90, 180 or 270")
    ap = settings.get("atlas_ap_mm")
    if ap is not None:
        ap = float(ap)
        if not np.isfinite(ap) or not 0 <= ap <= 13.2:
            raise ValueError("Atlas AP position must be 0–13.2 mm from the CCF volume origin")
    bbox = crop_bounds(tissue)
    target = tissue[bbox]
    aspect = target.shape[1]/target.shape[0]
    shape = tuple(max(8,round(size/max(target.shape)*150)) for size in target.shape)
    mask = resize(target,shape)
    guide = resize(hippocampus[bbox],shape) if hippocampus is not None else None
    if guide is not None and guide.sum()<15:
        guide = None
    hip_ids = [i for i,name in enumerate(names) if name.endswith(HIPPOCAMPAL)]
    ap_positions = np.array(levels)*.025
    nearest = int(np.argmin(abs(ap_positions-ap))) if ap is not None else None
    candidates = []

    def quality(layer, ratio):
        outline = dice(layer>0,mask)
        score = .85*outline+.15*ratio
        if guide is not None:
            score = .43*score+.57*dice(np.isin(layer,hip_ids),guide)
        return score,outline

    for index,plane in enumerate(planes):
        if nearest is not None and index!=nearest:
            continue
        for turn in turns:
            rotated = np.rot90(plane,turn)
            crop = rotated[crop_bounds(rotated>0)]
            ratio = min(aspect/(crop.shape[1]/crop.shape[0]),(crop.shape[1]/crop.shape[0])/aspect)
            layer = resize(crop,shape)
            score,_ = quality(layer,ratio)
            candidates.append((score,index,turn,ratio,layer,crop))
    trials = []
    start = [0,1,1,0,0]
    for initial,index,turn,ratio,layer,crop in sorted(candidates,key=lambda x:x[0],reverse=True)[:6]:
        # Refine only shortlisted planes. Keep the original if refinement cannot improve the score.
        opt = minimize(lambda p:-quality(corrected(layer,p),ratio)[0],start,method="Powell",
                       bounds=[(-.16,.16),(.88,1.12),(.88,1.12),(-7,7),(-7,7)],
                       options={"maxiter":8,"maxfev":170})
        params = opt.x if -opt.fun > initial else np.array(start)
        score,outline = quality(corrected(layer,params),ratio)
        trials.append((score,outline,index,turn,params,crop))
    score,outline,index,turn,params,crop = max(trials,key=lambda x:x[0])
    metadata = {"source":SOURCE,"reference":"https://doi.org/10.1016/j.cell.2020.04.007",
                "ap_mm_from_ccf_origin":round(float(ap_positions[index]),2),
                "rotation":turn*90,"outline_dice":outline,"guided":guide is not None,"ap_fixed":ap is not None}
    if outline < .72:
        metadata["method"] = "Allen CCFv3 已加载。H&E 与候选冠状层面轮廓匹配度偏低，请核对层面和朝向后重算。"
        return tissue.astype(np.int16), names[:2], metadata
    scaled = params.copy()
    scaled[3:] *= np.array(target.shape)/np.array(shape)
    projected = corrected(resize(crop,target.shape),scaled)
    projected = np.where(target,np.maximum(1,projected),0)
    result = np.zeros(tissue.shape,dtype=np.int16)
    result[bbox] = projected
    # Keep only labels represented on this section. ID 1 always denotes unreviewed tissue.
    active = [0,1]+[i for i in np.unique(result) if i>1]
    lut = np.zeros(len(names),dtype=np.int16)
    lut[active] = np.arange(len(active))
    metadata["named_fraction"] = float(np.count_nonzero(result>1)/tissue.sum())
    metadata["method"] = ("Allen CCFv3 成年小鼠冠状图谱，按 H&E 轮廓"
                          + ("及海马参考区域" if guide is not None else "")
                          + "匹配。边界和朝向请结合组织结构复核。")
    return lut[result], [names[i] for i in active], metadata
