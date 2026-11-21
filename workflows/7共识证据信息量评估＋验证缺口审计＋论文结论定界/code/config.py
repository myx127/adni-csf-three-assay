from __future__ import annotations
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
from datetime import datetime
import csv,hashlib,json,sys
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];PROJECT=ROOT.parent
S1=PROJECT/'1候选逐蛋白可评价性审计';S4=PROJECT/'4真实预处理规则决策＋独立病理和认知覆盖审计';S5=PROJECT/'5真实设计下完整模型校准＋逐蛋白功效与主分析';S6=PROJECT/'6候选信号归因＋分子映射与正交验证'
RESULTS=ROOT/'results';PRIVATE=ROOT/'00_restricted';RESTRICTED=PRIVATE;LOGS=ROOT/'logs';FIGURES=ROOT/'figures';EVIDENCE=ROOT/'evidence';SIM=RESULTS/'simulation_replicates'
for p in [RESULTS,PRIVATE,LOGS,FIGURES,EVIDENCE,SIM]:p.mkdir(parents=True,exist_ok=True)
INPUTS={r['input_id']:r for r in csv.DictReader((PROJECT/'规划和详细数据情况以及分析方案和研究计划/01_数据清单/CORE_INPUT_MANIFEST.csv').open(encoding='utf-8-sig'))}
INPUTS.update({r['input_id']:r for r in csv.DictReader((S4/'results/ADDITIONAL_INPUT_MANIFEST.csv').open(encoding='utf-8-sig'))})
MARGINS=[.2,.35,.5];K_GRID=[1,2,4,8];GAMMAS=[0.,.55];SCENARIOS=['IID','SLOPES3','T5_SLOPES3'];SEED_ROOT=2026103001
FAMILY=json.loads((S5/'results/MAIN_ANALYSIS_SAP_FROZEN.json').read_text(encoding='utf-8'))['primary_family_accessions'];M=len(FAMILY);HARMONIC=sum(1/j for j in range(1,M+1));ALPHA_BY_SUFFICIENT=.05/(M*HARMONIC)
CANDIDATES={'Q15818':'NPTX1','P09619':'PDGFRB','Q92876':'KLK6'}
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
def save_json(path,x):Path(path).write_text(json.dumps(clean(x),ensure_ascii=False,indent=2),encoding='utf-8')
def csvout(path,x):
 if not isinstance(x,pd.DataFrame):x=pd.DataFrame(x)
 x.to_csv(path,index=False,encoding='utf-8-sig')
def log(message):
 s=datetime.now().astimezone().isoformat(timespec='seconds')+' '+message;print(s,flush=True)
 with (LOGS/'STAGE7_RUN.log').open('a',encoding='utf-8') as f:f.write(s+'\n')
