from __future__ import annotations
import os
for task_option in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:
    os.environ[task_option]='1'
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import csv,hashlib,json,sys
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
PROJECT=ROOT.parent
STAGES={i:next(p for p in PROJECT.iterdir() if p.is_dir() and p.name.startswith(str(i))) for i in range(1,9)}
S1,S4,S5,S6,S7,S8=[STAGES[i] for i in [1,4,5,6,7,8]]
RESULTS=ROOT/'results';PRIVATE=ROOT/'00_restricted';EVIDENCE=ROOT/'evidence';FIGURES=ROOT/'figures';LOGS=ROOT/'logs'
for task_path in [RESULTS,PRIVATE,EVIDENCE,FIGURES,LOGS]:task_path.mkdir(parents=True,exist_ok=True)
RAW=Path(pd.read_csv(S6/'results/LEGACY_INPUT_MANIFEST.csv').iloc[0]['path']).parents[1]
CANDIDATES={'Q15818':'NPTX1','P09619':'PDGFRB','Q92876':'KLK6'}
PEPTIDES={'NPTX1':['FQLTFPLR','LENLEQYSR'],'KLK6':['ESSQEQSSVVR','LSELIQPLPLER','YTNWIQK']}
FAMILY=json.loads((S5/'results/MAIN_ANALYSIS_SAP_FROZEN.json').read_text(encoding='utf-8'))['primary_family_accessions']
MODELS=['M0','APOE','BOTH']
def now():return datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(timespec='seconds')
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()
def clean(x):
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if isinstance(x,np.ndarray):return clean(x.tolist())
    if isinstance(x,np.generic):return clean(x.item())
    if isinstance(x,float) and not np.isfinite(x):return None
    return x
def save_json(p,x):Path(p).write_text(json.dumps(clean(x),ensure_ascii=False,indent=2),encoding='utf-8')
def csvout(p,x):
    if not isinstance(x,pd.DataFrame):x=pd.DataFrame(x)
    x.to_csv(p,index=False,encoding='utf-8-sig')
def log(s):
    line=now()+' '+s;print(line,flush=True)
    with (LOGS/'STAGE9_RUN.log').open('a',encoding='utf-8') as f:f.write(line+'\n')
def by(p):
    p=np.asarray(p,float);m=len(p);order=np.argsort(p);q=np.minimum.accumulate((p[order]*m*np.sum(1/np.arange(1,m+1))/np.arange(1,m+1))[::-1])[::-1]
    ans=np.ones(m);ans[order]=np.minimum(1,q);return ans

def assert_unsealed():
    if (ROOT/'SHA256_MANIFEST.csv').exists():
        raise SystemExit('Stage9 is sealed. Use verify_delivery.py for read-only checks, or reproduce_in_copy.py for a new sibling analysis.')
