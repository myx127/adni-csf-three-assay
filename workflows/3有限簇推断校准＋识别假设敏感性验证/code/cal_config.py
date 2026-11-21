from __future__ import annotations
import hashlib,json
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
PROJECT=ROOT.parent
STAGE1=PROJECT/'1候选逐蛋白可评价性审计'
STAGE2=PROJECT/'2候选模型的识别模拟'
RESULTS,LOGS,FIGURES=ROOT/'results',ROOT/'logs',ROOT/'figures'
for d in [RESULTS,LOGS,FIGURES]:d.mkdir(parents=True,exist_ok=True)
METHODS=['IID_CHI2','CR1_CHI2','CR1_HT_F','WILD_SCORE_T','SUM3_CHI2','SUM3_MIN_F','CGM3_CHI2','PIGEONHOLE_CENTERED']
NEW_METHODS=['CR1_HT_F','WILD_SCORE_T','SUM3_CHI2','SUM3_MIN_F','PIGEONHOLE_CENTERED']
ROOT_SEEDS={'DEVELOPMENT':2026102701,'CONFIRMATION':2026102702,'N_GRID':2026102702,'SENSITIVITY':2026102703}

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
 return h.hexdigest()

def save_json(path,value):
 def clean(x):
  if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
  if isinstance(x,(list,tuple)):return [clean(v) for v in x]
  if isinstance(x,np.ndarray):return clean(x.tolist())
  if isinstance(x,np.generic):return clean(x.item())
  if isinstance(x,float) and not np.isfinite(x):return None
  return x
 Path(path).write_text(json.dumps(clean(value),ensure_ascii=False,indent=2),encoding='utf-8')

def csv(path,value):
 if not isinstance(value,pd.DataFrame):value=pd.DataFrame(value)
 value.to_csv(path,index=False,encoding='utf-8-sig')

def log(text):
 line=datetime.now().astimezone().isoformat(timespec='seconds')+' '+str(text)
 print(line,flush=True)
 with (LOGS/'CALIBRATION_RUN.log').open('a',encoding='utf-8') as f:f.write(line+'\n')

def calibration_scenarios():
 rows=[
 ('C01','IID等精度',True,2,{'sd':[.8,.8,.8]}),
 ('C02','IID不等精度',True,2,{'sd':[.25,.8,2.]}),
 ('C03','65簇共享误差',True,2,{'batch_sd':[0,0,.6],'exposure_cluster_axis':2}),
 ('C04','36簇共享误差',True,1,{'batch_sd':[0,.6,0],'exposure_cluster_axis':1}),
 ('C05','19簇共享误差',True,0,{'batch_sd':[.6,0,0],'exposure_cluster_axis':0}),
 ('C06','三向共享误差0.45',True,2,{'batch_sd':[.45,.45,.45]}),
 ('C07','三向共享误差0.8',True,2,{'batch_sd':[.8,.8,.8],'exposure_cluster_axis':2}),
 ('C08','不平衡三组',True,2,{'group_probs':[.6,.3,.1],'batch_sd':[0,0,.6],'exposure_cluster_axis':2}),
 ('C09','MCAR20及65簇',True,2,{'missing':'mcar','batch_sd':[0,0,.6]}),
 ('C10','共同及误差方差变化',True,2,{'uvar':[.5,1.,1.5],'hetero':.3,'batch_sd':[0,0,.6]}),
 ('A01','M2相对截距0.35',False,2,{'delta':[0,.35,0],'batch_sd':[0,0,.6]}),
 ('A02','M3相对截距0.50',False,2,{'delta':[0,0,.5],'batch_sd':[0,0,.6]}),
 ('A03','M3相对载荷0.40',False,2,{'load_delta':[0,0,.4],'batch_sd':[0,0,.6]}),
 ('A04','M2截距0.50及三向批次',False,2,{'delta':[0,.5,0],'batch_sd':[.45]*3}),
 ]
 return [{'scenario_id':sid,'index':i,'label_zh':label,'compatible_null':null,'cluster_axis':axis,'settings':settings}
         for i,(sid,label,null,axis,settings) in enumerate(rows)]

def sensitivity_scenarios():
 rows=[]
 def add(sid,label,settings):rows.append({'scenario_id':sid,'index':len(rows),'label_zh':label,'compatible_null':False,'cluster_axis':2,'settings':{'batch_sd':[0,0,.45],**settings}})
 for i,rho in enumerate([-.4,-.2,0.,.2,.4]):add(f'R{i}',f'固定误差相关{rho}',{'rho':rho})
 for i,q in enumerate([0.,.1,.2,.35,.5]):add(f'Q{i}',f'固定二次响应{q}',{'quadratic':q})
 add('SAT','固定饱和响应',{'saturation':True})
 add('MCAR','MCAR20',{'missing':'mcar'})
 for i,s in enumerate([.35,.7,1.]):add(f'SEL{i}',f'选择性缺失强度{s}',{'missing':'selective','selection_strength':s})
 add('LOD30','LOD30删除',{'missing':'lod30'});add('LOD70','LOD70删除',{'missing':'lod70'})
 add('FLOOR','LOD30地板替代',{'missing':'floor30'})
 add('RAW','未处理plex混杂',{'technical':'raw','batch_sd':[0,0,0]})
 add('GIS','GIS后混杂设计',{'technical':'gis','batch_sd':[0,0,0]})
 add('WEAK','弱共同信号',{'beta':.1,'uvar':[.05]*3,'sd':[.8,1.,1.5]})
 add('OFF2','M2相对截距0.20',{'delta':[0,.2,0]});add('OFF3','M3相对截距0.10',{'delta':[0,0,.1]})
 add('RHOV','误差相关随组变化',{'rho':[-.2,0,.4]})
 assert len(rows)==24
 return rows

def grid():
 rows=[]
 def add(phase,s,n,reps,B):
  phaseid={'DEVELOPMENT':1,'CONFIRMATION':2,'N_GRID':3,'SENSITIVITY':4}[phase]
  seed=int(np.random.SeedSequence([ROOT_SEEDS[phase],phaseid,s['index'],n]).generate_state(1,dtype=np.uint64)[0])
  rows.append({'cell_id':f'{phase}_{s["scenario_id"]}_n{n}','phase':phase,'scenario_id':s['scenario_id'],
               'scenario_index':s['index'],'n':n,'repetitions':reps,'bootstrap_B':B,'seed':str(seed),'cluster_axis':s['cluster_axis']})
 for phase,reps,B in [('DEVELOPMENT',1000,199),('CONFIRMATION',2000,499)]:
  for s in calibration_scenarios():add(phase,s,571,reps,B)
 for s in calibration_scenarios():
  if s['scenario_id'] in ['C03','C05','C06','A01']:
   for n in [200,400]:add('N_GRID',s,n,1000,499)
 for s in sensitivity_scenarios():add('SENSITIVITY',s,571,1000,199)
 assert sum(r['repetitions'] for r in rows)==74000
 return rows
