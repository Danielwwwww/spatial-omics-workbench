"""Portable H&E, mass-spectrometry imaging and spatial RNA workflow."""
import argparse
import base64
import gzip
import heapq
import io
import json
import re
import uuid
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
from scipy import ndimage as ndi, sparse

from .registration import apply, register, warp, largest, otsu, affine, sample
from .manual_regions import validate_regions
from .atlas import project as project_atlas
from .st_display import STDisplay


def table(path):
    name = path.name.lower().removesuffix('.gz')
    return pd.read_csv(path, sep=',' if name.endswith('.csv') else '\t', comment='#')


def points(table):
    xy = table[['x','y']].to_numpy(dtype=float)
    if len(xy) < 8 or not np.isfinite(xy).all() or (np.ptp(xy, axis=0) == 0).any():
        raise ValueError('At least eight finite spatial points spanning x and y are required')
    return xy


def numeric(counts, name):
    counts = np.asarray(counts, dtype=float)
    if not np.isfinite(counts).all() or (counts < 0).any():
        raise ValueError(f'{name} must contain finite nonnegative values')
    return counts


def load_st(path, dge, selected, max_genes):
    if path.suffix.lower() == '.h5ad':
        import anndata as ad

        data = ad.read_h5ad(path)
        if 'spatial' not in data.obsm:
            raise ValueError('ST h5ad requires obsm["spatial"]')
        xy = points(pd.DataFrame(np.asarray(data.obsm['spatial'])[:, :2], columns=['x','y']))
        x = data.layers['counts'] if 'counts' in data.layers else data.X
        if not data.var_names.is_unique or x is None:
            raise ValueError('ST h5ad requires unique gene names and a count matrix')
        raw = np.asarray(x.sum(axis=1)).ravel()
        counts = x.data if sparse.issparse(x) else x
        numeric(counts, 'ST counts')
        indices = choose_genes(list(data.var_names), np.asarray((x > 0).sum(axis=0)).ravel(), selected, max_genes)
        return xy, raw, list(data.var_names[indices]), x[:, indices]
    t = table(path)
    if 'CELL_BARCODE' in t:
        xy = points(t['CELL_BARCODE'].str.extract(r'^(?P<x>\d+(?:\.\d+)?)-(?P<y>\d+(?:\.\d+)?)$'))
        raw = numeric(t['NUM_TRANSCRIPTS'], 'ST UMI')
        if dge is None:
            return xy, raw, [], np.empty((len(xy), 0))
        selected = list(dict.fromkeys(selected))[:max_genes]
        keep = set(selected); expression = {}
        reader = gzip.open if dge.name.endswith('.gz') else open
        with reader(dge, 'rb') as f:
            header = f.readline().decode().rstrip('\r\n').split('\t')
            barcode_index = pd.Index(header[1:])
            if barcode_index.has_duplicates or t.CELL_BARCODE.duplicated().any():
                raise ValueError('DGE and summary barcodes must be unique')
            index = barcode_index.get_indexer(t.CELL_BARCODE)
            if header[0] != 'GENE' or (index < 0).any() or len(index) != len(header) - 1:
                raise ValueError('DGE header must match summary barcodes')
            def value_array(gene, values):
                result = numeric(np.fromstring(values, sep='\t'), gene)
                if len(result) != len(index):
                    raise ValueError(f'Bad DGE row for {gene}')
                return result[index]

            candidates = []
            limit = min(24, max_genes)
            for row_number, line in enumerate(f):
                name, _, values = line.partition(b'\t')
                gene = name.decode()
                if gene in keep:
                    expression[gene] = value_array(gene, values)
                    if len(expression) == len(keep):
                        break
                if not selected:
                    # Wide DGE stores canonical tab-separated integer counts.
                    zeros = values.count(b'\t0') + values.startswith(b'0')
                    entry = (len(index) - zeros, -row_number, gene, values)
                    if entry[0] <= 0:
                        continue
                    if len(candidates) < limit:
                        heapq.heappush(candidates, entry)
                    elif entry[:2] > candidates[0][:2]:
                        heapq.heapreplace(candidates, entry)
            if not selected:
                for _, _, gene, values in sorted(candidates, reverse=True):
                    expression[gene] = value_array(gene, values)
                selected = list(expression)
        names = [n for n in selected if n in expression][:max_genes]
        return xy, raw, names, np.column_stack([expression[n] for n in names]) if names else np.empty((len(raw),0))
    xy = points(t)
    gene_columns = [c for c in t if c.startswith('gene:')]
    x = numeric(t[gene_columns], 'gene counts')
    raw = numeric(t['umi'], 'ST UMI') if 'umi' in t else x.sum(axis=1)
    ix = choose_genes([n[5:] for n in gene_columns], (x > 0).sum(0), selected, max_genes)
    return xy, raw, [gene_columns[i][5:] for i in ix], x[:, ix]


def choose_genes(names, detected, requested, maximum):
    if requested:
        name_index = {v: i for i,v in enumerate(names)}
        return np.array([name_index[v] for v in dict.fromkeys(requested) if v in name_index][:maximum], dtype=int)
    return np.argsort(-detected, kind='stable')[:maximum]


def discover_masses(spectra, tic, ppm):
    """Suggest recurring intense m/z features from evenly sampled spectra."""
    support = {}
    step = max(2 * ppm, 12) / 1e6
    for i in np.linspace(0, len(spectra) - 1, min(512, len(spectra)), dtype=int):
        mz, intensity = spectra[i]
        if tic[i] <= 0 or not len(mz):
            continue
        top = np.argpartition(intensity, -min(32, len(mz)))[-32:]
        seen = {}
        for m, value in zip(mz[top], intensity[top]):
            if 50 <= m <= 2000 and value > 0:
                key = round(np.log(m) / step)
                weight = float(value / tic[i])
                if key not in seen or weight > seen[key][0]:
                    seen[key] = (weight, float(m))
        for key, (weight, mass) in seen.items():
            record = support.setdefault(key, [0, 0., 0.])
            record[0] += 1
            record[1] += weight
            record[2] += weight * mass
    ranked = sorted((row for row in support.values() if row[0] >= 2), key=lambda row: -row[1])
    result = []
    for _, weight, weighted_mass in ranked:
        mass = weighted_mass / weight
        if all(abs(mass - other) > 3 * ppm * mass / 1e6 for other in result):
            result.append(round(mass, 5))
        if len(result) == 8:
            break
    if not result:
        raise ValueError('No recurrent ion features found; supply target m/z values')
    return result


def read_imzml(path, ibd, masses, ppm):
    if ibd is None:
        raise ValueError('imzML requires its IBD companion')
    params = lambda e: {c.get('accession'): c.get('value') for c in e.findall('{*}cvParam')}
    groups, positions, spectra = {}, [], []
    declared_uuid = None
    for _, node in ET.iterparse(path, events=('end',)):
        name = node.tag.rsplit('}',1)[-1]
        if name == 'referenceableParamGroup':
            groups[node.get('id')] = params(node)
        elif name == 'fileContent':
            declared_uuid = params(node).get('IMS:1000080')
        elif name == 'spectrum':
            arrays = {}
            for a in node.findall('./{*}binaryDataArrayList/{*}binaryDataArray'):
                ref = a.find('{*}referenceableParamGroupRef')
                p = {**(groups.get(ref.get('ref'),{}) if ref is not None else {}), **params(a)}
                if 'MS:1000574' in p or 'MS:1002312' in p:
                    raise ValueError('Re-export imzML with uncompressed external arrays')
                dtype = next((d for k,d in [('MS:1000523','<f8'),('MS:1000521','<f4'),('MS:1000519','<i4'),('MS:1000522','<i8')] if k in p),None)
                key = 'mz' if 'MS:1000514' in p else 'intensity'
                if dtype is not None:
                    arrays[key] = (np.dtype(dtype), int(p['IMS:1000103']), int(p['IMS:1000102']))
            if 'mz' not in arrays or 'intensity' not in arrays:
                raise ValueError('imzML spectrum requires supported m/z and intensity arrays')
            if arrays['mz'][1] != arrays['intensity'][1]:
                raise ValueError('m/z and intensity arrays have different lengths')
            scan = params(node.find('./{*}scanList/{*}scan'))
            positions.append([int(scan['IMS:1000050']), int(scan['IMS:1000051'])])
            spectra.append((arrays['mz'], arrays['intensity']))
            node.clear()
    xy = points(pd.DataFrame(positions,columns=['x','y']))
    if len(np.unique(xy, axis=0)) != len(xy):
        raise ValueError('imzML contains repeated x/y positions; supply a single 2D section')
    with ibd.open('rb') as stream:
        if declared_uuid and uuid.UUID(declared_uuid).bytes != stream.read(16):
            raise ValueError('The imzML UUID does not match IBD')

        def read_array(spec):
            dtype, length, offset = spec
            if length < 0 or offset < 0:
                raise ValueError('imzML array length and offset must be nonnegative')
            stream.seek(offset)
            buffer = stream.read(length * dtype.itemsize)
            if len(buffer) != length * dtype.itemsize:
                raise ValueError('An imzML spectrum points outside its IBD file')
            return np.frombuffer(buffer, dtype=dtype)

        if not masses:
            order = np.lexsort(xy.T)
            sample_indices = np.linspace(0, len(order) - 1, min(512, len(order)), dtype=int)
            sampled = [tuple(map(read_array, spectra[order[i]])) for i in sample_indices]
            sampled_tic = [intensity.sum(dtype=np.float64) for _, intensity in sampled]
            masses = discover_masses(sampled, sampled_tic, ppm)

        tic, ions = [], []
        for mz_spec, intensity_spec in spectra:
            mz, intensity = read_array(mz_spec), read_array(intensity_spec)
            tic.append(intensity.sum(dtype=np.float64))
            ions.append([intensity[np.abs(mz-m) <= m*ppm/1e6].sum(dtype=np.float64) for m in masses])
    tic = numeric(tic,'TIC')
    return xy, tic, [f'm/z {m}' for m in masses], numeric(ions,'ion intensity')


def load_msi(path, ibd, masses, ppm):
    if path.suffix.lower() == '.imzml':
        return read_imzml(path, ibd, masses, ppm)
    t = table(path)
    names = [n for n in t if n.startswith('ion:')]
    if not names:
        raise ValueError('MSI table must include at least one ion:<feature> column')
    x = numeric(t[names], 'MSI intensity')
    tic = numeric(t['tic'], 'TIC') if 'tic' in t else x.sum(1)
    return points(t), tic, [n[4:] for n in names], x


def read_he(path, override=None, max_size=900):
    with Image.open(path) as image:
        rgb = image.convert('RGB')
    original_size = rgb.size
    rgb.thumbnail((max_size, max_size), Image.Resampling.LANCZOS)
    scale = np.array(rgb.size) / original_size
    transform = affine(np.diag(scale), (scale-1)/2)
    arr = np.asarray(rgb).astype(float)
    extent = max(1, round(max(rgb.size)/450))
    if override is not None:
        with Image.open(override) as mask:
            if mask.size != original_size:
                raise ValueError('H&E mask must have the original image dimensions')
            tissue = np.asarray(mask.convert('L').resize(rgb.size, Image.Resampling.NEAREST)) > 0
    else:
        feature = (255-arr.min(2)) * .55 + np.ptp(arr, axis=2) * .45
        tissue = ndi.binary_closing(feature > max(otsu(feature)*.68, 11), iterations=extent)
        labels, _ = ndi.label(tissue)
        sizes = np.bincount(labels.ravel()); sizes[0]=0
        threshold = max(20, int(sizes.max() * .015))
        tissue = ndi.binary_fill_holes(sizes[labels] >= threshold)
    if tissue.sum() < 100 or tissue.sum() == tissue.size:
        raise ValueError('Uninformative H&E tissue mask; upload a binary mask')
    buffer = io.BytesIO()
    rgb.save(buffer, format='JPEG', quality=86)
    feature = -np.log((arr[:,:,0]+1) / 256)
    return base64.b64encode(buffer.getvalue()).decode(), feature, tissue, transform, original_size


def raster(xy, values, kind, pitch=None):
    origin = np.min(xy, axis=0)
    factor = np.array(pitch if pitch is not None else (1,1), dtype=float)
    if factor.shape != (2,) or (factor <= 0).any() or not np.isfinite(factor).all():
        raise ValueError('Coordinate pitch must have two positive finite elements')
    scaled = (xy-origin) * factor
    spacing = max(np.ptp(scaled, axis=0).max() / 464, 1e-12)
    # Retain native spacing on small regular grids instead of upsampling sparse observations.
    steps = []
    for axis in range(2):
        d = np.diff(np.unique(scaled[:,axis])); d = d[d > spacing * .02]
        if d.size:
            steps.append(np.quantile(d,.5))
    spacing = max(spacing, min(steps))
    mat = affine(np.diag(factor/spacing), 3-origin*factor/spacing)
    x,y = np.floor(apply(mat,xy)+.5).astype(int).T
    dims = (y.max()+4,x.max()+4)
    counts = np.zeros(dims); levels = counts.copy()
    raw = np.log1p(values)
    np.add.at(counts,(y,x),1); np.add.at(levels,(y,x),raw)
    observed = counts>0
    image = np.divide(levels, counts, out=np.zeros_like(levels), where=observed)
    if kind == 'st':
        # Treat low-count rows as background for initialization only.
        threshold = max(otsu(image[observed]), 0)
        area = observed & (image >= threshold) & (image > 0)
        area = ndi.binary_fill_holes(largest(ndi.binary_closing(area,iterations=2)))
    else:
        area = ndi.binary_fill_holes(ndi.binary_closing(observed, iterations=1) | observed)
    if area.sum() < 30:
        area = ndi.binary_fill_holes(ndi.binary_closing(observed) | observed)
    # Smooth and fill gaps near measured data without pulling zero padding into the signal.
    den = ndi.gaussian_filter(observed.astype(float), .8)
    num = ndi.gaussian_filter(image, .8)
    # Keep fine anatomical contrast on measured pixels. Fill only gaps by weighted smoothing.
    image = np.where(observed, image, np.divide(num,den,out=np.zeros_like(num),where=den>1e-6))
    return image, area, mat


def convert_landmarks(rows, source_transform, target_transform):
    fitted = []
    for row in rows:
        src,dst = np.array([row['source'],row['target']],dtype=float)
        weight = float(row.get('weight',1))
        if (not np.isfinite([src,dst]).all() or row['role'] not in ('fit','check') or not np.isfinite(weight) or weight <= 0):
            raise ValueError('Landmarks require finite x/y, a fit or check role, and positive weight')
        fitted.append({**row,'weight':weight,
                       'source':apply(source_transform,src[None])[0].tolist(),
                       'target':apply(target_transform,dst[None])[0].tolist()})
    return fitted


def polygon_pixels(corners, tissue, he_to_view):
    xy = np.asarray(corners,dtype=float)
    if xy.ndim != 2 or xy.shape[1] != 2 or len(xy)<3 or not np.isfinite(xy).all():
        raise ValueError('ROI polygons require at least three H&E x/y vertices')
    bitmap = Image.new('1',tissue.shape[::-1])
    ImageDraw.Draw(bitmap).polygon([tuple(row) for row in apply(he_to_view,xy)],fill=1)
    return np.asarray(bitmap) & tissue


def region_labels(config, tissue, he_to_view, required_regions=(), base=None):
    result, names = base if base is not None else (
        tissue.astype(np.int16), ['Outside H&E', 'Other tissue' if config else 'Tissue'])
    for name, corners in config.items():
        inside = polygon_pixels(corners,tissue,he_to_view)
        if name in required_regions and not inside.any():
            raise ValueError(f"{name} 没有覆盖 H&E 组织。")
        if name not in names:
            names.append(str(name))
        result[inside] = names.index(name)
    return result, names


def display_labels(region, mask):
    # Row-wise runs keep the web artifact much smaller than a full label raster.
    boundaries = []
    for val in np.unique(region[mask]):
        if val == 0:
            continue
        active = mask & (region == val)
        bounds = active & ~ndi.binary_erosion(active)
        y,x = np.nonzero(bounds)
        boundaries.append({'id':int(val),'xy':np.column_stack((x,y))[::2].tolist()})
    return boundaries


def round_rows(a, digits=4):
    return np.round(np.asarray(a,dtype=float),digits).tolist()


def spatial_gene_layers(xy, labels, umi, matrix, shape):
    width, height = shape[::-1]
    display = STDisplay(xy, labels, umi, width, height)
    columns = matrix.tocsc() if sparse.issparse(matrix) else None
    result = []
    for j in range(matrix.shape[1]):
        if columns is not None:
            start, end = columns.indptr[j:j+2]
            rows, counts = columns.indices[start:end], columns.data[start:end]
        else:
            values = np.asarray(matrix[:, j]).ravel()
            rows = np.flatnonzero(values)
            counts = values[rows]
        result.append(display.gene(rows, counts))
    return result


def integrate(mxy, sxy, ms_data, st_data, regions, tile, mask, st_support, min_umi=1, footprint=None):
    m_label, s_label = sample(regions,mxy,0),sample(regions,sxy,0)
    sgood = (s_label > 0) & sample(mask & st_support,sxy,0) & (st_data[0]>=min_umi)
    nw = int(np.ceil(regions.shape[1]/tile))+1
    n = (int(np.ceil(regions.shape[0]/tile))+1)*nw
    def keys(xy, label, good):
        i,j = np.floor(xy[good]/tile).astype(int).T
        return ((j*nw+i)*int(regions.max()+1)+label[good]).astype(int)
    sk = keys(sxy,s_label,sgood)
    samples = mxy[:, None, :] if footprint is None else footprint
    mk, pixel_ids = [], []
    for points in samples.transpose(1, 0, 2):
        label = sample(regions,points,0)
        eligible = (label > 0) & sample(st_support,points,0) & (points >= 0).all(axis=1)
        mk.append(keys(points,label,eligible))
        pixel_ids.append(np.flatnonzero(eligible))
    mk, pixel_ids = np.concatenate(mk), np.concatenate(pixel_ids)
    size=n*(regions.max()+1)
    mmat=sparse.csr_matrix((np.full(len(mk),1/samples.shape[1]),(mk,pixel_ids)),shape=(size,len(mxy)))
    mc,sc=mmat.getnnz(axis=1),np.bincount(sk,minlength=size)
    ids=np.flatnonzero((mc>0)&(sc>0))
    if not ids.size:
        raise ValueError('No common grid cells after registration. Review coordinates and add landmarks.')
    smat=sparse.csr_matrix((np.ones(len(sk)),(sk,np.arange(len(sk)))),shape=(size,len(sk)))
    tic, ion = ms_data
    normalized_ion = np.divide(ion,tic[:,None],out=np.zeros_like(ion),where=tic[:,None]>0)
    # Exclude zero-TIC pixels from normalized means. Measured zero ion intensity remains zero.
    support = mmat[ids]
    mass = np.asarray(support @ (tic > 0)).ravel()
    normalized = support @ normalized_ion
    mz=np.divide(normalized,mass[:,None],out=np.zeros_like(normalized),where=mass[:,None]>0)
    m_norm_n = support.sign() @ (tic > 0)
    genes=smat[ids]@st_data[1][sgood]; genes=genes.toarray() if sparse.issparse(genes) else genes
    umi=smat[ids]@st_data[0][sgood]
    count=(ids//(regions.max()+1)).astype(int)
    centers = (smat[ids] @ sxy[sgood]) / sc[ids,None]
    regional = []
    for label in np.unique(s_label[sgood]):
        selected = sgood & (s_label == label)
        total = float(st_data[0][selected].sum())
        expression = np.asarray(st_data[1][selected].sum(axis=0)).ravel()
        regional.append({'id':int(label), 'stN':int(selected.sum()), 'umi':total,
                         'cp10k':round_rows(expression * (10000 / total), 6)})
    return {'xy': round_rows(centers,3), 'tile_xy': np.column_stack((count%nw,count//nw)).tolist(),
            'region': (ids%(regions.max()+1)).tolist(),'msiN':mc[ids].tolist(),'stN':sc[ids].tolist(),
            'msiTicN':m_norm_n.astype(int).tolist(),'msiWeight':round_rows(mass,8),
            'umi':round_rows(umi), 'ions':round_rows(mz.T,9),
            'genes':round_rows(np.log1p(1e4*genes/umi[:,None]).T,5)}, m_label, regional


def analyze(he, msi, st, config, output=None, ibd=None, dge=None, he_mask=None, coordinate_dir=None):
    strength = float(config.get('ppm',10)); masses=list(map(float,config.get('mz',[])))
    if not 0<strength<=1000 or any(not np.isfinite(m) or m<=0 for m in masses) or len(masses)>75:
        raise ValueError('Use 1–75 positive m/z targets and a positive tolerance of at most 1000 ppm')
    ne = int(config.get('max_genes',48))
    if not 0 < ne <= 2000:
        raise ValueError('max_genes must be in 1–2000')
    stxy, umi, genes, x = load_st(st,dge,config.get('genes',[]),ne)
    mxy, tic, ions, inten = load_msi(msi,ibd,masses,strength)
    jpeg, target, tmask, hmat, original = read_he(he,he_mask)
    scales = config.get('he_um_per_px')
    scales = np.asarray(scales if scales is not None else (1,1),dtype=float)
    if scales.shape != (2,) or not np.isfinite(scales).all() or (scales<=0).any():
        raise ValueError('he_um_per_px must be [x, y] and positive')
    ev = scales / np.diag(hmat[:2,:2]); units = 'µm' if 'he_um_per_px' in config else 'H&E original px'
    alignments, layers, masks = {}, {}, {}
    footprint = None
    pitch = config.get('msi_pixel_pitch_xy', (1,1) if msi.suffix.lower()=='.imzml' else None)
    if pitch is not None:
        pitch = np.asarray(pitch,dtype=float)
        if pitch.shape != (2,) or not np.isfinite(pitch).all() or (pitch <= 0).any():
            raise ValueError('msi_pixel_pitch_xy must contain two positive finite widths')
    for name, xy, signal in [('msi',mxy,tic),('st',stxy,umi)]:
        image, mask, smat = raster(xy,signal,name,config.get(f'{name}_pitch_xy'))
        lm = convert_landmarks(config.get('landmarks',{}).get(name,[]), smat,hmat)
        matrix, field, qa, candidate = register(image,mask,target,tmask,lm,
                                               refine=config.get('refine',True),
                                               allow_refined=config.get('use_refined',{}).get(name,True),
                                               evaluation_scale=ev, use_stalign=config.get('stalign',False),
                                               tre_priority=config.get('tre_priority','balanced'))
        to_view = matrix @ smat
        q = apply(to_view,xy)
        if field is not None:
            q = field.forward(q)
        if name == 'msi' and pitch is not None and config.get('area_weighted',True):
            axis = np.array([-1,0,1])/3
            offsets = np.array(np.meshgrid(axis,axis)).reshape(2,-1).T * pitch
            offsets_xy = (xy[:,None,:] + offsets).reshape(-1,2)
            subpoints = apply(to_view, offsets_xy)
            if field is not None:
                subpoints = field.forward(subpoints)
            footprint = subpoints.reshape(len(xy),len(offsets),2)
        layers[name]=q; masks[name]=warp(mask,matrix,tmask.shape,field,order=0)
        qa.update({'source_to_raster':smat.tolist(),'source_to_he_view_affine':to_view.tolist(),'nonlinear_inverse_field':field.to_json() if field else None,
                   'candidate_inverse_field':candidate.to_json() if candidate else None,'tre_units':units})
        alignments[name]=qa
    atlas = None
    base = None
    if config.get('atlas') == 'allen-ccfv3':
        guides = config.get('atlas_hippocampus', [])
        guide = np.zeros_like(tmask)
        for corners in guides:
            guide |= polygon_pixels(corners,tmask,hmat)
        atlas_labels, atlas_names, atlas = project_atlas(tmask,config,guide if guides else None)
        base = (atlas_labels,atlas_names)
    base_regions = config.get('regions', {})
    reserved = base[1] if base is not None else ("Outside H&E","Other tissue","Tissue")
    manual = validate_regions(config.get('manual_regions', []), original,
                              (*reserved, *base_regions))
    rois = dict(base_regions)
    rois.update((r["name"], r["points"]) for r in manual)
    region, region_names = region_labels(rois,tmask,hmat,{r["name"] for r in manual},base)
    tile=float(config.get('tile_size_he_px',max(target.shape)/22 / np.diag(hmat[:2,:2]).mean()))
    tile *= np.diag(hmat[:2,:2]).mean()
    if not np.isfinite(tile) or tile < 3 or tile > max(target.shape):
        raise ValueError('Tile size must be between 3 preview pixels and the full image width or height')
    combined = masks['msi'] & masks['st'] & tmask
    min_umi = int(config.get('min_umi', 3))
    if min_umi < 1:
        raise ValueError('min_umi must be a positive integer')
    st_labels = sample(region, layers['st'], 0)
    st_labels[umi < min_umi] = 0
    if region.max() >= 256:
        raise ValueError('Supply at most 254 named brain regions')
    st_maps = spatial_gene_layers(layers['st'], st_labels, umi, x, region.shape)
    grid, ml, regional = integrate(layers['msi'], layers['st'], (tic,inten), (umi,x), region, tile,
                                  masks['msi'], masks['st'], min_umi, footprint)
    overall = grid
    if len(region_names) > 2:
        overall, _, _ = integrate(layers['msi'], layers['st'], (tic,inten), (umi,x), tmask.astype(np.int16),
                                  tile, masks['msi'], masks['st'], min_umi, footprint)
    ratio = np.divide(inten,tic[:,None],out=np.zeros_like(inten),where=tic[:,None]>0)
    if coordinate_dir is not None:
        coordinate_dir.mkdir(parents=True, exist_ok=True)
        for kind, original_xy in [('msi', mxy), ('st', stxy)]:
            mapped = layers[kind]
            frame = pd.DataFrame(original_xy, columns=['source_x','source_y'])
            frame.insert(0, 'input_row_0based', np.arange(len(frame)))
            frame[['he_x_px','he_y_px']] = apply(np.linalg.inv(hmat), mapped)
            frame['he_region'] = np.array(region_names)[sample(region, mapped, 0)]
            frame['in_he_tissue'] = sample(tmask, mapped, 0)
            frame['in_three_way_support'] = sample(combined, mapped, 0)
            frame.to_csv(coordinate_dir/f'{kind}_coordinates.csv.gz', index=False, float_format='%.7g')
    bundle={'schema':1,'name':str(config.get('name','New dataset')),'origin':str(config.get('origin','uploaded data')),
            'image':{'width':target.shape[1],'height':target.shape[0],'originalSize':list(original),'jpeg':jpeg},
            'ionNames':ions,'geneNames':genes,'regions':region_names,'atlas':atlas,'boundaries':display_labels(region,tmask),
            'regionExpression':regional,'stMap':st_maps,'manualRegions':manual,'globalGrid':overall,
            'msi':{'xy':round_rows(layers['msi'],3),'region':ml.tolist(),'tic':round_rows(tic),'raw':round_rows(inten.T),
                   'norm':round_rows(ratio.T,9)},'grid':grid, 'tileSize':tile,
            'qa':{'alignment':alignments,'gridMapping':'footprint' if footprint is not None else 'centroid',
                  'tissueCommonFraction':float(combined.sum()/tmask.sum()),
                  'msiPixels':len(mxy),'stPoints':len(stxy),'overlapGridCells':len(set(map(tuple,overall['tile_xy']))),
                  'regionGridCells':len(grid['xy']),'minUmi':min_umi,'autoMz':msi.suffix.lower()=='.imzml' and not masses,
                  'ticMode':'per-spectrum TIC' if msi.suffix.lower()=='.imzml' else ('CSV TIC if supplied; otherwise sum of ion columns'),
                  'ppm':strength,'heOriginalToView':hmat.tolist()}}
    if output:
        output.parent.mkdir(parents=True,exist_ok=True)
        temp=output.with_suffix(output.suffix+'.tmp')
        temp.write_text(json.dumps(bundle,ensure_ascii=False,separators=(',',':'),allow_nan=False),encoding='utf-8')
        temp.replace(output)
    return bundle


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('he','msi','st','config','out'):
        p.add_argument('--'+name,type=Path,required=True)
    for name in ('ibd','dge','he-mask','coordinate-dir'):
        p.add_argument('--'+name,type=Path)
    args=p.parse_args()
    results=analyze(args.he,args.msi,args.st,json.loads(args.config.read_text(encoding='utf-8')),args.out,args.ibd,args.dge,
                    args.he_mask,args.coordinate_dir)
    print(f"Wrote {args.out} · {results['qa']['overlapGridCells']} matched grid cells")


if __name__=='__main__':
    main()
