"""Local upload server for the same browser viewer served by GitHub Pages."""
import json
import shutil
from importlib.util import find_spec
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Optional

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles

from pipeline import analyze


ROOT = Path(__file__).resolve().parent
app = FastAPI(title='Spatial Omics Workbench')


@app.get('/api/health')
def health():
    return {'upload': True, 'stalign': find_spec('STalign') is not None}


def save(upload, folder, name, suffixes):
    if upload is None:
        return None
    suffix = next((s for s in suffixes if (upload.filename or '').lower().endswith(s)), None)
    if suffix is None:
        raise ValueError(f'{name} requires one of {", ".join(suffixes)}')
    path = folder / (name + suffix)
    with path.open('wb') as stream:
        shutil.copyfileobj(upload.file, stream)
    return path


@app.post('/api/analyze')
def analyze_upload(he:UploadFile=File(...), msi:UploadFile=File(...), st:UploadFile=File(...),
                   ibd:Optional[UploadFile]=File(None), dge:Optional[UploadFile]=File(None),
                   he_mask:Optional[UploadFile]=File(None), config: str=Form('{}')):
    try:
        settings=json.loads(config)
        if not isinstance(settings,dict):
            raise ValueError('Configuration must be a JSON object')
        with TemporaryDirectory(prefix='spatial-workbench-') as temp:
            root=Path(temp)
            he_path=save(he,root,'he',['.png','.jpg','.jpeg','.tif','.tiff'])
            msi_path=save(msi,root,'msi',['.imzml','.csv.gz','.tsv.gz','.csv','.tsv'])
            st_path=save(st,root,'st',['.h5ad','.csv.gz','.tsv.gz','.txt.gz','.csv','.tsv','.txt'])
            companion=save(ibd,root,'msi',['.ibd'])
            dge_path=save(dge,root,'dge',['.txt.gz','.tsv.gz','.txt','.tsv'])
            he_mask_path=save(he_mask,root,'mask',['.png','.tif','.tiff'])
            result=analyze(he_path,msi_path,st_path,settings,ibd=companion,dge=dge_path,he_mask=he_mask_path)
        return JSONResponse(result)
    except (ValueError, KeyError, TypeError, IndexError, OSError) as exc:
        raise HTTPException(status_code=422,detail=str(exc)) from exc


app.mount('/',StaticFiles(directory=ROOT/'docs',html=True),name='viewer')
