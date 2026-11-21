from __future__ import annotations
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
from datetime import datetime
import hashlib,json,csv
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];PROJECT=ROOT.parent
S1=PROJECT/'1候选逐蛋白可评价性审计';S3=PROJECT/'3有限簇推断校准＋识别假设敏感性验证';S4=PROJECT/'4真实预处理规则决策＋独立病理和认知覆盖审计'
RESULTS=ROOT/'results';RESTRICTED=ROOT/'00_restricted';LOGS=ROOT/'logs';FIGURES=ROOT/'figures';SIM=RESULTS/'simulation_replicates'
for p in [RESULTS,RESTRICTED,LOGS,FIGURES,SIM]:p.mkdir(parents=True,exist_ok=True)
MANIFEST=PROJECT/'规划和详细数据情况以及分析方案和研究计划/01_数据清单/CORE_INPUT_MANIFEST.csv'
INPUTS={r['input_id']:r for r in csv.DictReader(MANIFEST.open(encoding='utf-8-sig'))}
INPUTS.update({r['input_id']:r for r in csv.DictReader((S4/'results/ADDITIONAL_INPUT_MANIFEST.csv').open(encoding='utf-8-sig'))})
SEEDS={'DEV':2026102801,'CONFIRM':2026102802,'POWER':2026102803,'STRESS':2026102804,'SECONDARY':2026102805}
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
 t=datetime.now().astimezone().isoformat(timespec='seconds')+' '+s
 print(t,flush=True)
 with (LOGS/'STAGE5_RUN.log').open('a',encoding='utf-8') as f:f.write(t+'\n')
