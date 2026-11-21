from __future__ import annotations
import os
for task_thread_limit in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[task_thread_limit]='1'
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import sys,re,json,csv,hashlib
import numpy as np
import pandas as pd
if hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8',errors='backslashreplace')
ROOT=Path(__file__).resolve().parents[1];PROJECT=ROOT.parent
STAGES={i:next(p for p in PROJECT.iterdir() if p.is_dir() and re.match(r'^'+str(i)+r'(?!\d)',p.name)) for i in range(1,12)}
S1,S5,S11=[STAGES[i] for i in [1,5,11]]
RESULTS=ROOT/'results';PRIVATE=ROOT/'00_restricted';FIGURES=ROOT/'figures';EVIDENCE=ROOT/'evidence';LOGS=ROOT/'logs';SIM=RESULTS/'simulation_replicates'
for task_path in [RESULTS,PRIVATE,FIGURES,EVIDENCE,LOGS,SIM]:task_path.mkdir(parents=True,exist_ok=True)
manifest=next(PROJECT.rglob('CORE_INPUT_MANIFEST.csv'));INPUTS={r['input_id']:r for r in csv.DictReader(manifest.open(encoding='utf-8-sig'))}
GENES={'P47972':'NPTX2','Q7Z5A7':'TAFA5','P35625':'TIMP3'}
MODELS=['NO_P','LOG_P','RAW_P','RCS_P','BG2_NO_P','BG2_LOG_P','BG2_RCS_P','BG3_NO_P','BG3_LOG_P']
COMPARISONS={'LOG_SHIFT':('NO_P','LOG_P'),'RAW_SHIFT':('NO_P','RAW_P'),'RCS_SHIFT':('NO_P','RCS_P'),'BG2_SHIFT':('BG2_NO_P','BG2_LOG_P'),'BG3_SHIFT':('BG3_NO_P','BG3_LOG_P')}
PRIMARY_COMPARISONS=['LOG_SHIFT','RAW_SHIFT','RCS_SHIFT']
SCENARIOS=['IID_VISITS','AR1_AND_PERSON_SLOPE','HETERO_P','TECHNICAL_AGE_SLOPES','T5_TECHNICAL_SLOPES']
REPS=60000;BATCH=512;WORKERS=8;ROOT_SEED=202610290012
def alpha(m):return .05/(m*sum(1/i for i in range(1,m+1)))
ALPHAS={'NOMINAL':.05,'BY9_SUFFICIENT':alpha(9),'BY33_SUFFICIENT':alpha(33),'BY42_SUFFICIENT':alpha(42)}
def now():return datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(timespec='seconds')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
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
 msg=now()+' '+s;print(msg,flush=True)
 with (LOGS/'STAGE12_RUN.log').open('a',encoding='utf-8') as f:f.write(msg+'\n')
def by(p):
 p=np.asarray(p,float);m=len(p);ix=np.argsort(p);q=np.minimum.accumulate((p[ix]*m*sum(1/i for i in range(1,m+1))/np.arange(1,m+1))[::-1])[::-1];out=np.ones(m);out[ix]=np.minimum(q,1);return out
def assert_unsealed():
 if (ROOT/'SHA256_MANIFEST.csv').exists():raise SystemExit('Sealed. Verify or reproduce in a fresh sibling directory.')

