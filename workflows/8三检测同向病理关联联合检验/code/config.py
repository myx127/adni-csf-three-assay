from __future__ import annotations
import os
for task_thread_option in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[task_thread_option]='1'
from pathlib import Path
from datetime import datetime
import hashlib,json,csv,sys
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];PROJECT=ROOT.parent
S1=PROJECT/'1候选逐蛋白可评价性审计';S4=PROJECT/'4真实预处理规则决策＋独立病理和认知覆盖审计'
S5=PROJECT/'5真实设计下完整模型校准＋逐蛋白功效与主分析';S6=PROJECT/'6候选信号归因＋分子映射与正交验证'
S7=PROJECT/'7共识证据信息量评估＋验证缺口审计＋论文结论定界'
RESULTS=ROOT/'results';PRIVATE=ROOT/'00_restricted';RESTRICTED=PRIVATE;LOGS=ROOT/'logs';FIGURES=ROOT/'figures';EVIDENCE=ROOT/'evidence';SIM=RESULTS/'simulation_replicates'
for task_output_path in [RESULTS,PRIVATE,LOGS,FIGURES,EVIDENCE,SIM]:task_output_path.mkdir(parents=True,exist_ok=True)
FAMILY=json.loads((S5/'results/MAIN_ANALYSIS_SAP_FROZEN.json').read_text(encoding='utf-8'))['primary_family_accessions']
MODELS=['M0','APOE','BOTH'];SCENARIOS=['IID','HETERO','SLOPE_N','SLOPES3','T5_SLOPES3']
REPS=40000;BATCH=128;WORKERS=6;ROOT_SEED=202610280801
ALPHAS={'NOMINAL':.05,'BY34_SUFFICIENT':.05/(34*sum(1/j for j in range(1,35))),'BY68_SUFFICIENT':.05/(68*sum(1/j for j in range(1,69)))}
CANDIDATES={'Q15818':'NPTX1','P09619':'PDGFRB','Q92876':'KLK6'}
def patterns():
 rows=[{'pattern':'ZERO_ALL','kind':'NULL','effect':[0.,0.,0.],'weak_axis':'ALL','magnitude':0.,'sign':0}]
 for axis in range(3):
  for sign in [-1,1]:
   for group,base,value in [('ONE_ZERO',1.,0.),('TWO_ZERO',0.,1.),('MIXED_STRONG',.5,-.35),('MIXED_NEAR_ZERO',.5,-.02)]:
    v=np.full(3,base);v[axis]=value;v*=sign
    rows.append({'pattern':f'{group}_{axis}_{sign}','kind':'NULL','effect':v.tolist(),'weak_axis':['N','S','T'][axis],'magnitude':abs(value),'sign':sign})
 for effect in [.10,.20,.35,.50,.75]:
  for sign in [-1,1]:rows.append({'pattern':f'BALANCED_{effect}_{sign}','kind':'POWER','effect':[effect*sign]*3,'weak_axis':'ALL','magnitude':effect,'sign':sign})
 for axis in range(3):
  for effect in [.10,.20,.35]:
   for sign in [-1,1]:
    v=np.full(3,.5);v[axis]=effect;v*=sign
    rows.append({'pattern':f'UNBALANCED_{axis}_{effect}_{sign}','kind':'POWER','effect':v.tolist(),'weak_axis':['N','S','T'][axis],'magnitude':effect,'sign':sign})
 return rows
PATTERNS=patterns()
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for chunk in iter(lambda:f.read(4*1024*1024),b''):h.update(chunk)
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
 line=datetime.now().astimezone().isoformat(timespec='seconds')+' '+message;print(line,flush=True)
 with (LOGS/'STAGE8_RUN.log').open('a',encoding='utf-8') as f:f.write(line+'\n')
