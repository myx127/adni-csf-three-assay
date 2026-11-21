from __future__ import annotations
from collections import Counter,defaultdict
import hashlib,json,re
import pandas as pd
from audit_common import *

FIELD_ACCESS=[];STRUCTURE=[];VALUE_STATUS=[];DUPLICATES=[]

def stream(key,fields):
 FIELD_ACCESS.append({'input_id':key,'fields':fields,'purpose':'Date/validity/categories only; no values exported'})
 yield from rows(key,fields)

def signature(bound):
 # Temporary equality tokens are removed before any file is written.
 return hashlib.sha256(repr(bound.signature()).encode()).hexdigest() if bound.known else ''

def consolidate(records,keyfields,source):
 if not records:return pd.DataFrame()
 frame=pd.DataFrame(records);clean=[];duplicate_groups=conflicts=0
 sigcols=[c for c in frame if c.startswith('_sig_')]
 for _,part in frame.groupby(keyfields,dropna=False,sort=False):
  if len(part)>1:duplicate_groups+=1
  # PET top QC / latest processing priority is metadata-based, not value-based.
  if '_qc_rank' in part:
   part=part.loc[part._qc_rank.eq(part._qc_rank.min())]
   latest=part._process_date.max()
   if latest:part=part.loc[part._process_date.eq(latest)]
  out=part.iloc[0].to_dict();raw_conflict={}
  for col in sigcols:
   values=set(part[col])-set(['']);name=col[5:];raw_conflict[name]=len(values)>1
   if raw_conflict[name]:conflicts+=1
   out[name+'_CONFLICT']=raw_conflict[name]
  for col in frame.columns:
   if col.endswith('_EXACT') or col.endswith('_KNOWN') or col.endswith('_VALID'):
    marker=col.rsplit('_',1)[0]
    out[col]=bool(part[col].any() and not raw_conflict.get(marker,False))
   if col.endswith('_CAT'):
    values=set(part[col])-{'UNKNOWN',''}
    out[col]=next(iter(values)) if len(values)==1 else 'UNKNOWN'
  for prefix in ['A','T','AT']:
   if prefix+'_CAT' in out and prefix+'_VALID' in out:out[prefix+'_VALID']=out[prefix+'_CAT']!='UNKNOWN'
  if 'A_CAT' in out and 'T_CAT' in out:
   out['AT_CAT']='A'+('+' if out['A_CAT']=='POS' else '-')+'T'+('+' if out['T_CAT']=='POS' else '-') if out['A_CAT']!='UNKNOWN' and out['T_CAT']!='UNKNOWN' else 'UNKNOWN'
   out['AT_VALID']=out['AT_CAT']!='UNKNOWN'
  clean.append({k:v for k,v in out.items() if not k.startswith('_sig_')})
 DUPLICATES.append({'source':source,'raw_reduced_records':len(records),'unique_date_assay_records':len(clean),
                    'duplicate_key_groups':duplicate_groups,'conflicting_marker_groups':conflicts})
 return pd.DataFrame(clean)

def load_csf(key):
 elec=key=='E01';markers=['ABETA40','ABETA42','TAU','PTAU'] if elec else ['ABETA','TAU','PTAU']
 fields=['RID','BATCH','EXAMDATE' if elec else 'DRAWDTE','RUNDATE','VISCODE2' if elec else 'VISCODE']+markers
 if elec:fields+=['COMMENT']
 result=[];status={m:Counter() for m in markers};counts=Counter();people=set()
 for line,row in stream(key,fields):
  counts['rows']+=1;person=rid(row.get('RID'));day=parse_date(row.get('EXAMDATE' if elec else 'DRAWDTE'))
  if not person:counts['invalid_RID']+=1;continue
  people.add(person)
  if not day:counts['invalid_date']+=1;continue
  batch=text(row.get('BATCH'));match=re.search(r'(\d+)$',batch)
  priority=0 if batch=='MEDIAN' else int(match.group(1)) if match else 1
  out={'RID':person,'date':day,'source':key,'assay':'ELECSYS' if elec else 'ALZBIO3','batch':batch,
       'batch_priority':priority,'source_row':line,'VISCODE':text(row.get('VISCODE2' if elec else 'VISCODE')),
       'comment_present':bool(text(row.get('COMMENT')))}
  values={}
  for marker in markers:
   limits={'ABETA42':(200.,1700.),'TAU':(80.,1300.),'PTAU':(8.,120.)}.get(marker) if elec else None
   value=numeric(row.get(marker),positive=True,limits=limits);values[marker]=value;status[marker][value.kind]+=1
   out[marker+'_EXACT']=value.exact;out[marker+'_KNOWN']=value.known;out[marker+'_KIND']=value.kind
   out['_sig_'+marker]=signature(value)
  if elec:
   pairs={'A4240':ratio(values['ABETA42'],values['ABETA40']),'PTAU_ABETA42':ratio(values['PTAU'],values['ABETA42'])}
   for name,value in pairs.items():out[name+'_EXACT']=value.exact;out[name+'_KNOWN']=value.known;out['_sig_'+name]=signature(value)
   out['A_CAT']=values['ABETA42'].category(980.,'lt');out['T_CAT']=values['PTAU'].category(27.,'gt')
   out['A4240_CAT']=pairs['A4240'].category(.061,'lt');out['AD_RATIO_CAT']=pairs['PTAU_ABETA42'].category(.025,'gt')
  else:
   out['A_CAT']=values['ABETA'].category(192.,'lt');out['T_CAT']=values['PTAU'].category(23.,'gt')
  out['AT_CAT']='A'+('+' if out['A_CAT']=='POS' else '-')+'T'+('+' if out['T_CAT']=='POS' else '-') if out['A_CAT']!='UNKNOWN' and out['T_CAT']!='UNKNOWN' else 'UNKNOWN'
  out['AT_VALID']=out['AT_CAT']!='UNKNOWN';out['A_VALID']=out['A_CAT']!='UNKNOWN';out['T_VALID']=out['T_CAT']!='UNKNOWN'
  result.append(out)
 for marker,counter in status.items():
  for state,n in counter.items():VALUE_STATUS.append({'source':key,'marker':marker,'state':state,'n_records':n})
 STRUCTURE.append({'source':key,**counts,'unique_RID':len(people)})
 return consolidate(result,['RID','date','batch'],key)

def load_pet(key):
 amy=key=='I11';markers=['CENTILOIDS','SUMMARY_SUVR'] if amy else ['META_TEMPORAL_SUVR']
 fields=['RID','PTID','SCANDATE','PROCESSDATE','VISCODE2','qc_flag','TRACER','AMYLOID_STATUS']+markers
 records=[];counts=Counter();people=set();status={m:Counter() for m in markers}
 for line,row in stream(key,fields):
  counts['rows']+=1;person=rid(row.get('RID'));day=parse_date(row.get('SCANDATE'))
  if not person:counts['invalid_RID']+=1;continue
  people.add(person)
  if not day:counts['invalid_date']+=1;continue
  qc=text(row.get('qc_flag'));accepted=qc in ['1','2'];process=parse_date(row.get('PROCESSDATE'))
  out={'RID':person,'date':day,'source':key,'tracer':text(row.get('TRACER')),'qc':qc,'qc_accepted':accepted,
       'qc_full_pass':qc=='2','_qc_rank':0 if qc=='2' else 1 if qc=='1' else 2,'_process_date':str(process or ''),
       'source_row':line,'VISCODE':text(row.get('VISCODE2')),'PTID_RID_MATCH':True}
  ptid=text(row.get('PTID'));m=re.search(r'_S_(\d+)$',ptid)
  if m and int(m.group(1))!=person:out['PTID_RID_MATCH']=False;counts['PTID_RID_conflicts']+=1;accepted=False
  for marker in markers:
   value=numeric(row.get(marker),positive=marker!='CENTILOIDS');status[marker][value.kind]+=1
   out[marker+'_EXACT']=value.exact and accepted;out['_sig_'+marker]=signature(value)
   if marker=='CENTILOIDS':
    out['CL18_CAT']=value.category(18.,'gt') if accepted else 'UNKNOWN'
    out['CL18_EQUAL']=value.exact and value.low==18 and accepted
  if amy:
   a=numeric(row.get('AMYLOID_STATUS'));category='POS' if a.exact and a.low==1 else 'NEG' if a.exact and a.low==0 else 'UNKNOWN'
   out['A_CAT']=category if accepted else 'UNKNOWN';out['A_VALID']=out['A_CAT']!='UNKNOWN';out['_sig_A']=category if category!='UNKNOWN' else ''
  records.append(out)
 for marker,counter in status.items():
  for state,n in counter.items():VALUE_STATUS.append({'source':key,'marker':marker,'state':state,'n_records':n})
 STRUCTURE.append({'source':key,**counts,'unique_RID':len(people)})
 return consolidate(records,['RID','date','tracer'],key)

def load_cognition():
 fields=['RID','EXAMDATE','VISCODE2','PHASE','ADNI_MEM','ADNI_EF','ADNI_EF2','ADNI_LAN','ADNI_VS'];records=[];counts=Counter();people=set()
 for line,row in stream('I09',fields):
  counts['rows']+=1;person=rid(row.get('RID'));day=parse_date(row.get('EXAMDATE'))
  if not person:counts['invalid_RID']+=1;continue
  people.add(person)
  if not day:counts['invalid_date']+=1;continue
  out={'RID':person,'date':day,'source_row':line,'VISCODE':text(row.get('VISCODE2')),'phase':text(row.get('PHASE'))}
  for marker in fields[4:]:
   value=numeric(row.get(marker));out[marker+'_EXACT']=value.exact;out['_sig_'+marker]=signature(value)
   if value.exact and value.low==-4:counts[marker+'_literal_minus4_finite']+=1
  records.append(out)
 STRUCTURE.append({'source':'I09',**counts,'unique_RID':len(people)})
 return consolidate(records,['RID','date'],'I09')
