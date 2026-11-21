from __future__ import annotations
import csv,hashlib,json,math,re
from dataclasses import dataclass
from datetime import datetime,date
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];PROJECT=ROOT.parent;STAGE1=PROJECT/'1候选逐蛋白可评价性审计'
RESULTS=ROOT/'results';RESTRICTED=ROOT/'00_restricted';LOGS=ROOT/'logs';FIGURES=ROOT/'figures'
for p in [RESULTS,RESTRICTED,LOGS,FIGURES]:p.mkdir(parents=True,exist_ok=True)
MANIFEST=PROJECT/'规划和详细数据情况以及分析方案和研究计划/01_数据清单/CORE_INPUT_MANIFEST.csv'
INPUTS={r['input_id']:r for r in csv.DictReader(MANIFEST.open(encoding='utf-8-sig'))}
EXTRA={r['input_id']:r for r in csv.DictReader((RESULTS/'ADDITIONAL_INPUT_MANIFEST.csv').open(encoding='utf-8-sig'))}
INPUTS.update(EXTRA)
MISSING={'','NA','N/A','NAN','NULL','NONE','NAT','MISSING','.'}
WINDOWS=[0,30,90,180,365,730]

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save_json(path,x):
 def clean(v):
  if isinstance(v,dict):return {str(k):clean(y) for k,y in v.items()}
  if isinstance(v,(list,tuple)):return [clean(y) for y in v]
  if isinstance(v,np.generic):return clean(v.item())
  if isinstance(v,(date,datetime)):return v.isoformat()
  if isinstance(v,float) and not math.isfinite(v):return None
  return v
 Path(path).write_text(json.dumps(clean(x),ensure_ascii=False,indent=2),encoding='utf-8')
def write_csv(path,x):
 if not isinstance(x,pd.DataFrame):x=pd.DataFrame(x)
 x.to_csv(path,index=False,encoding='utf-8-sig')
def log(message):
 line=datetime.now().astimezone().isoformat(timespec='seconds')+' '+message
 print(line,flush=True)
 with (LOGS/'AUDIT_RUN.log').open('a',encoding='utf-8') as f:f.write(line+'\n')
def text(x):
 s=str(x or '').strip()
 return '' if s.upper() in MISSING else s
def rid(x):
 try:
  v=float(x);return int(v) if v>0 and v.is_integer() else None
 except (ValueError,TypeError):return None
def parse_date(x):
 s=text(x)
 if not s:return None
 for fmt in ['%Y-%m-%d','%m/%d/%Y','%Y/%m/%d','%Y-%m-%d %H:%M:%S','%Y-%m-%d %H:%M:%S.%f']:
  try:return datetime.strptime(s,fmt).date()
  except ValueError:pass
 return None
def rows(input_id,fields):
 # Read only named fields into dictionaries; non-selected values are never interpreted.
 with Path(INPUTS[input_id]['path']).open(encoding='utf-8-sig',newline='') as f:
  reader=csv.reader(f);header=next(reader);positions={name:header.index(name) for name in fields if name in header}
  for index,record in enumerate(reader,2):
   yield index,{name:record[j] if j<len(record) else '' for name,j in positions.items()}

@dataclass
class Bound:
 kind:str
 low:float=float('nan')
 high:float=float('nan')
 low_open:bool=False
 high_open:bool=False
 @property
 def exact(self):return self.kind=='EXACT'
 @property
 def known(self):return self.kind in ['EXACT','LEFT_CENSORED','RIGHT_CENSORED','INTERVAL_CENSORED','OUTSIDE_TECHNICAL_LOW','OUTSIDE_TECHNICAL_HIGH']
 def category(self,cutoff,direction):
  if not self.known:return 'UNKNOWN'
  if direction=='lt':
   if self.high<cutoff or (self.high==cutoff and self.high_open):return 'POS'
   if self.low>=cutoff:return 'NEG'
  if direction=='gt':
   if self.low>cutoff or (self.low==cutoff and self.low_open):return 'POS'
   if self.high<=cutoff:return 'NEG'
  return 'UNKNOWN'
 def signature(self):return (self.kind,self.low,self.high,self.low_open,self.high_open) if self.known else (self.kind,)

NUMBER=re.compile(r'^\s*([<>≤≥]?)\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$')
def numeric(x,positive=False,limits=None):
 s=text(x)
 if not s:return Bound('MISSING')
 m=NUMBER.fullmatch(s)
 if not m:return Bound('INVALID')
 v=float(m.group(2));op=m.group(1)
 if not math.isfinite(v) or positive and v<=0:return Bound('INVALID')
 if op in ['<','≤']:return Bound('LEFT_CENSORED',0. if positive else -math.inf,v,positive,op=='<')
 if op in ['>','≥']:return Bound('RIGHT_CENSORED',v,math.inf,op=='>',False)
 if limits and v<limits[0]:return Bound('OUTSIDE_TECHNICAL_LOW',0.,limits[0],True,True)
 if limits and v>limits[1]:return Bound('OUTSIDE_TECHNICAL_HIGH',limits[1],math.inf,True,False)
 return Bound('EXACT',v,v)
def ratio(a,b):
 if a.exact and b.exact and b.low>0:return Bound('EXACT',a.low/b.low,a.low/b.low)
 if not a.known or not b.known or b.low<0:return Bound('MISSING')
 low=0. if math.isinf(b.high) else a.low/b.high
 high=math.inf if b.low<=0 or math.isinf(a.high) else a.high/b.low
 return Bound('INTERVAL_CENSORED',low,high,a.low_open or b.high_open,a.high_open or b.low_open)

def register_thresholds():
 records=[
 ('E01','ELECSYS','A42','ABETA42','lt',980.,'pg/mL','A_OPERATIONAL','M02 local methods: Table3; 980 is documented rounded cutoff','LOCAL_M02'),
 ('E01','ELECSYS','A4240','ABETA42_OVER_ABETA40','lt',.061,'ratio','A_RATIO_OPERATIONAL','M02 local methods; requires independently collected ABETA40','LOCAL_M02'),
 ('E01','ELECSYS','T27','PTAU','gt',27.,'pg/mL','PTAU_HIGH_OPERATIONAL','Primary study clinical-progression cutpoint; not interchangeable with tau-PET positivity','https://pmc.ncbi.nlm.nih.gov/articles/PMC6911086/'),
 ('E01','ELECSYS','AD_RATIO','PTAU_OVER_ABETA42','gt',.025,'ratio','AD_SIGNATURE_NOT_T','M02 Table3; AD/amyloid signature ratio, not a separate T axis','LOCAL_M02'),
 ('I10','ALZBIO3','A192','ABETA','lt',192.,'pg/mL','LEGACY_A','M01 rescaled baseline threshold, not applicable to raw or Elecsys values','LOCAL_M01'),
 ('I10','ALZBIO3','T23','PTAU','gt',23.,'pg/mL','LEGACY_PTAU_HIGH','Published AlzBio3 convention; not tau-PET threshold','https://pmc.ncbi.nlm.nih.gov/articles/PMC2696350/'),
 ('I11','PET_AMY','A_RELEASE','AMYLOID_STATUS','provided',None,'binary','A_RELEASED','Use tracer-specific released status after primary-ROI QC','LOCAL_M03'),
 ('I11','PET_AMY','A_CL18','CENTILOIDS','gt',18.,'CL','A_CENTILOID_SENSITIVITY','M03 recommends18CL; report equality separately, operational >18','LOCAL_M03'),
 ('I12','PET_TAU','T_UNFROZEN','META_TEMPORAL_SUVR','none',None,'SUVR','CONTINUOUS_ONLY','No single validated cross-tracer T cutoff assumed','LOCAL_M04')]
 d=pd.DataFrame(records,columns=['input_id','assay','rule_id','marker','direction','cutoff','units','role','basis','source'])
 write_csv(RESULTS/'THRESHOLD_REGISTRY.csv',d)
 save_json(RESULTS/'THRESHOLD_FREEZE_BEFORE_COUNTS.json',{'threshold_registry_sha256':sha(RESULTS/'THRESHOLD_REGISTRY.csv'),
      'fixed_before_clinical_value_reduction':True,'tauPET_binary_definition_frozen':False,'clinical_association_analysis':False})
 return d
