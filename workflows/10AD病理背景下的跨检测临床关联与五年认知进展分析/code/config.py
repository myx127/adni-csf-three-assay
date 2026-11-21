from __future__ import annotations
import os
for task_thread_limit in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[task_thread_limit]='1'
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import csv,hashlib,json,re,sys
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];PROJECT=ROOT.parent
STAGES={i:next(p for p in PROJECT.iterdir() if p.is_dir() and re.match(r'^'+str(i)+r'(?!\d)',p.name)) for i in range(1,10)}
S1,S4,S5,S6,S7,S8,S9=[STAGES[i] for i in [1,4,5,6,7,8,9]]
RESULTS=ROOT/'results';PRIVATE=ROOT/'00_restricted';FIGURES=ROOT/'figures';EVIDENCE=ROOT/'evidence';LOGS=ROOT/'logs';SIM=RESULTS/'simulation_replicates'
for task_path in [RESULTS,PRIVATE,FIGURES,EVIDENCE,LOGS,SIM]:task_path.mkdir(parents=True,exist_ok=True)
manifest=next(PROJECT.rglob('CORE_INPUT_MANIFEST.csv'))
INPUTS={r['input_id']:r for r in csv.DictReader(manifest.open(encoding='utf-8-sig'))}
INPUTS.update({r['input_id']:r for r in csv.DictReader((S4/'results/ADDITIONAL_INPUT_MANIFEST.csv').open(encoding='utf-8-sig'))})
CANDIDATES={'Q15818':'NPTX1','Q92876':'KLK6','P09619':'PDGFRB'}
PRIMARY='MEM5Y'
FORMAL_BRANCHES=['MEM5Y','MEM_BASE','EF5Y','AMY_CL','MEM_AMY_INTERACTION','MEM_APOE_INTERACTION']
SENSITIVITIES=['MEM_NO_PTAU','MEM_RAW_P','MEM3Y','MEM_HIGHINFO','MEM_EQUAL_PERSON','MEM_BATCH9','MEM_BASE_ADJUSTED','MEM_QUADRATIC','MEM_PET_MAIN']
SCENARIOS=['IID_VISITS','AR1_AND_PERSON_SLOPE','HETERO_P','TECHNICAL_AGE_SLOPES','T5_TECHNICAL_SLOPES']
FORMAL_REPS=20000;SENS_REPS=6000;BATCH=256;WORKERS=4;ROOT_SEED=202610290010
ALPHAS={'NOMINAL':.05,'BY6_SUFFICIENT':.05/(6*sum(1/i for i in range(1,7))),'BY32_SUFFICIENT':.05/(32*sum(1/i for i in range(1,33)))}
def now():return datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(timespec='seconds')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
 return h.hexdigest()
def clean(v):
 if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [clean(x) for x in v]
 if isinstance(v,np.ndarray):return clean(v.tolist())
 if isinstance(v,np.generic):return clean(v.item())
 if isinstance(v,float) and not np.isfinite(v):return None
 return v
def save_json(p,v):Path(p).write_text(json.dumps(clean(v),ensure_ascii=False,indent=2),encoding='utf-8')
def csvout(p,v):
 if not isinstance(v,pd.DataFrame):v=pd.DataFrame(v)
 v.to_csv(p,index=False,encoding='utf-8-sig')
def log(s):
 line=now()+' '+s;print(line,flush=True)
 with (LOGS/'STAGE10_RUN.log').open('a',encoding='utf-8') as f:f.write(line+'\n')
def by(p):
 p=np.asarray(p,float);m=len(p);order=np.argsort(p);rank=np.arange(1,m+1);q=np.minimum.accumulate((p[order]*m*np.sum(1/rank)/rank)[::-1])[::-1]
 out=np.ones(m);out[order]=np.minimum(1,q);return out
def assert_unsealed():
 if (ROOT/'SHA256_MANIFEST.csv').exists():raise SystemExit('Sealed directory: use verify_delivery.py or reproduce in a new sibling directory')
