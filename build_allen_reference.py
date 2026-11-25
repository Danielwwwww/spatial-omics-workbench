"""Prepare compact public coronal planes from the Allen Mouse CCFv3 annotations."""
import argparse
import base64
import gzip
import json
from pathlib import Path

import numpy as np


GROUPS = [
    ("MOp", "初级运动皮层"), ("MOs", "次级运动皮层"), ("SSp", "初级躯体感觉皮层"),
    ("SSs", "补充躯体感觉皮层"), ("AUD", "听觉皮层"), ("VIS", "视觉皮层"),
    ("RSP", "压后皮层"), ("PTLp", "后顶叶联络皮层"), ("Isocortex", "其他新皮层"),
    ("OLF", "嗅觉区"), ("CA1", "海马 CA1"), ("CA2", "海马 CA2"), ("CA3", "海马 CA3"),
    ("DG", "齿状回"), ("SUB", "下托"), ("HPF", "其他海马结构"),
    ("TH", "丘脑"), ("HY", "下丘脑"), ("MB", "中脑"), ("STR", "纹状体"),
    ("PAL", "苍白球区"), (("LA","BLA","BMA","PA","CEA","MEA","COA","AAA"), "杏仁核复合区"),
    ("CTXsp", "皮层下板"), ("CB", "小脑"), ("P", "脑桥"), ("MY", "延髓"),
    ("cc", "胼胝体"), ("ec", "外囊"), ("int", "内囊"), ("fi", "海马伞"),
    ("fiber tracts", "其他纤维束"), ("VS", "脑室系统"),
]


def build(annotation, structure_file, output):
    header, raw = annotation.read_bytes().split(b"\n\n", 1)
    fields = dict(line.split(": ",1) for line in header.decode().splitlines() if ": " in line)
    shape = tuple(map(int, fields["sizes"].split()))
    if shape != (528,320,456) or (fields["encoding"],fields["type"],fields["endian"]) != ("gzip","unsigned int","little"):
        raise ValueError("Expected Allen CCFv3 annotation_25.nrrd")
    volume = np.frombuffer(gzip.decompress(raw), dtype="<u4").reshape(shape,order="F")
    trees = {}

    def descend(node):
        ids = {node["id"]}
        for child in node["children"]:
            ids.update(descend(child))
        trees[node["acronym"]] = ids
        return ids

    tissues = descend(json.loads(structure_file.read_text(encoding="utf-8"))["msg"][0])
    known = max(tissues)
    # Structure IDs can be large and sparse, so remap only the distinct atlas IDs.
    full_ids = np.array(sorted(tissues),dtype=np.uint32)
    names = ["组织外","待复核区域"]
    classes = np.ones(len(full_ids),dtype=np.uint8)
    for codes,name in reversed(GROUPS):
        names.append(name)
        for code in ((codes,) if isinstance(codes,str) else codes):
            classes[np.searchsorted(full_ids,sorted(trees[code]))] = len(names)-1
    levels = np.arange(0,shape[0],8)
    planes = volume[levels,::4,::4]
    bins = np.minimum(np.searchsorted(full_ids,planes),len(full_ids)-1)
    labels = np.where(planes==0,0,np.where((planes<=known)&(full_ids[bins]==planes),classes[bins],1)).astype(np.uint8)
    for name in ("海马 CA1","海马 CA2","海马 CA3","齿状回"):
        region = labels==names.index(name)
        labels[region & (np.arange(labels.shape[2]) < labels.shape[2]/2)] = len(names)
        names.append("图谱左"+name)
        labels[region & (np.arange(labels.shape[2]) >= labels.shape[2]/2)] = len(names)
        names.append("图谱右"+name)
    valid = (labels>0).sum(axis=(1,2))>=70
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps({"names":names,"ap_voxels":levels[valid].tolist(),
        "shape":list(labels[valid].shape),
        "labels":base64.b64encode(gzip.compress(labels[valid].tobytes(),mtime=0)).decode()},
        ensure_ascii=False,separators=(",",":")),encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("annotation",type=Path)
    parser.add_argument("structures",type=Path)
    parser.add_argument("output",type=Path)
    args = parser.parse_args()
    build(args.annotation,args.structures,args.output)
