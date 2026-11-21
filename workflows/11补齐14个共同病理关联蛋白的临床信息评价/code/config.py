from __future__ import annotations
import os
for task_thread_limit in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']: os.environ[task_thread_limit]='1'
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import csv,hashlib,json,re,sys
if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8',errors='backslashreplace')
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; PROJECT=ROOT.parent
STAGES={i:next(p for p in PROJECT.iterdir() if p.is_dir() and re.match(r'^'+str(i)+r'(?!\d)',p.name)) for i in range(1,11)}
S1,S4,S5,S8,S9,S10=[STAGES[i] for i in [1,4,5,8,9,10]]
RESULTS=ROOT/'results'; PRIVATE=ROOT/'00_restricted'; FIGURES=ROOT/'figures'; EVIDENCE=ROOT/'evidence'; LOGS=ROOT/'logs'; SIM=RESULTS/'simulation_replicates'
for task_path in [RESULTS,PRIVATE,FIGURES,EVIDENCE,LOGS,SIM]: task_path.mkdir(parents=True,exist_ok=True)
FAMILY34={'O00468':'AGRN','O94813':'SLIT2','P00441':'SOD1','P00558':'PGK1','P01034':'CST3','P01138':'NGF','P01303':'NPY','P05362':'ICAM1','P07196':'NEFL','P08069':'IGF1R','P08758':'ANXA5','P09619':'PDGFRB','P12036':'NEFH','P13500':'CCL2','P15328':'FOLR1','P21246':'PTN','P22676':'CALB2','P30041':'PRDX6','P35247':'SFTPD','P35625':'TIMP3','P35968':'KDR','P36222':'CHI3L1','P40925':'MDH1','P47972':'NPTX2','P78423':'CX3CL1','Q02246':'CNTN2','Q13231':'CHIT1','Q13421':'MSLN','Q15818':'NPTX1','Q16270':'IGFBP7','Q7Z5A7':'TAFA5','Q92686':'NRGN','Q92876':'KLK6','Q99988':'GDF15'}
FOCUS14=['AGRN','CALB2','CNTN2','CX3CL1','KLK6','MDH1','NEFL','NPTX1','NPTX2','NRGN','PGK1','SLIT2','TAFA5','TIMP3']
CANDIDATES=FAMILY34
PRIMARY='MEM5Y'; FORMAL_BRANCHES=['MEM5Y','MEM_BASE','EF5Y']
SENSITIVITIES=['MEM_NO_PTAU','MEM_HIGHINFO','MEM_EQUAL_PERSON']
SCENARIOS=['IID_VISITS','AR1_AND_PERSON_SLOPE','HETERO_P','TECHNICAL_AGE_SLOPES','T5_TECHNICAL_SLOPES']
FORMAL_REPS=100000; SENS_REPS=10000; BATCH=512; WORKERS=8; ROOT_SEED=202610290011
def alpha(m): return .05/(m*sum(1/i for i in range(1,m+1)))
ALPHAS={'NOMINAL':.05,'BY68_SUFFICIENT':alpha(68),'BY136_SUFFICIENT':alpha(136),'BY204_SUFFICIENT':alpha(204),'BY28_DIAGNOSTIC_SUFFICIENT':alpha(28)}
def now(): return datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(timespec='seconds')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(4*1024*1024),b''): h.update(block)
 return h.hexdigest()
def clean(v):
 if isinstance(v,dict): return {str(k):clean(x) for k,x in v.items()}
 if isinstance(v,(list,tuple)): return [clean(x) for x in v]
 if isinstance(v,np.ndarray): return clean(v.tolist())
 if isinstance(v,np.generic): return clean(v.item())
 if isinstance(v,float) and not np.isfinite(v): return None
 return v
def save_json(p,v): Path(p).write_text(json.dumps(clean(v),ensure_ascii=False,indent=2),encoding='utf-8')
def csvout(p,v):
 if not isinstance(v,pd.DataFrame): v=pd.DataFrame(v)
 v.to_csv(p,index=False,encoding='utf-8-sig')
def log(s):
 line=now()+' '+s; print(line,flush=True)
 with (LOGS/'STAGE11_RUN.log').open('a',encoding='utf-8') as f: f.write(line+'\n')
def by(p):
 p=np.asarray(p,float); m=len(p); order=np.argsort(p); rank=np.arange(1,m+1)
 q=np.minimum.accumulate((p[order]*m*np.sum(1/rank)/rank)[::-1])[::-1]
 out=np.ones(m); out[order]=np.minimum(1,q); return out
def assert_unsealed():
 if (ROOT/'SHA256_MANIFEST.csv').exists(): raise SystemExit('Sealed: verify delivery or reproduce in a fresh sibling directory.')

