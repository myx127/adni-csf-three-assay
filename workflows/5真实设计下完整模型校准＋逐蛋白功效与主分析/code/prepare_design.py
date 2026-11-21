from __future__ import annotations
from collections import defaultdict,Counter
from datetime import datetime
import csv,re
from config import *

def day(x):
 s=str(x or '').strip()
 for fmt in ['%Y-%m-%d','%m/%d/%Y','%Y/%m/%d','%Y-%m-%d %H:%M:%S','%Y-%m-%d %H:%M:%S.%f']:
  try:return datetime.strptime(s,fmt).date()
  except ValueError:pass
 return None
def records(key):
 with Path(INPUTS[key]['path']).open(encoding='utf-8-sig',newline='') as f:
  for i,row in enumerate(csv.DictReader(f),2):yield i,row
def value(x):
 try:
  v=float(x);return v if np.isfinite(v) else np.nan
 except (ValueError,TypeError):return np.nan

def verify():
 checks=[]
 for k,row in INPUTS.items():
  digest=sha(row['path']);checks.append({'id':k,'path':row['path'],'sha256':digest,'expected_sha256':row['sha256'],'passed':digest==row['sha256'] and Path(row['path']).stat().st_size==int(row['bytes'])})
 for stage,manifestname,col,rels in [
  (S1,'STAGE1_INPUT_OUTPUT_SHA256.csv','path',['00_restricted/TRIPLE_SAMPLE_REGISTRY_RESTRICTED.csv','00_restricted/TRIPLE_CANDIDATE_AVAILABILITY_RESTRICTED.parquet','results/EXTENDED_CORE_MAPPING_DRAFT.csv','results/PROTEIN_EVALUABILITY_70.csv']),
  (S4,'STAGE4_INPUT_OUTPUT_SHA256.csv','relative_path',['00_restricted/INDEX_CLINICAL_COVERAGE_FLAGS_RESTRICTED.csv','00_restricted/CSF_ELECSYS_FLAGS_RESTRICTED.csv','00_restricted/AMYLOID_PET_FLAGS_RESTRICTED.csv','00_restricted/COVARIATE_FLAGS_RESTRICTED.csv','results/PREPROCESSING_RULES_v1.0.json','results/ADDITIONAL_INPUT_MANIFEST.csv'])]:
  manifest=pd.read_csv(stage/manifestname)
  for rel in rels:
   expected=manifest.loc[manifest[col].eq(rel),'sha256'].iloc[0];digest=sha(stage/rel)
   checks.append({'id':stage.name+'/'+rel,'path':str(stage/rel),'sha256':digest,'expected_sha256':expected,'passed':digest==expected})
 csvout(RESULTS/'SOURCE_HASH_CHECKS.csv',checks);assert all(x['passed'] for x in checks),'Source hash mismatch'
 initial=RESULTS/'ACCESS_BOUNDARY_INITIAL.json'
 if initial.exists():assert json.loads(initial.read_text(encoding='utf-8'))['protocol_sha256']==sha(ROOT/'STAGE5_PROTOCOL_v1.0.md')
 else:save_json(initial,{'protocol_sha256':sha(ROOT/'STAGE5_PROTOCOL_v1.0.md'),'before_stage5_clinical_value_extraction':True,'before_stage5_protein_value_connection':True,'user_authorized_main_analysis':True})

def clinical():
 index=pd.read_csv(S1/'00_restricted/TRIPLE_SAMPLE_REGISTRY_RESTRICTED.csv')
 flags=pd.read_csv(S4/'00_restricted/INDEX_CLINICAL_COVERAGE_FLAGS_RESTRICTED.csv').set_index('RID')
 cov=pd.read_csv(S4/'00_restricted/COVARIATE_FLAGS_RESTRICTED.csv').set_index('RID')
 assert len(index)==571 and index.RID.is_unique
 ef=pd.read_csv(S4/'00_restricted/CSF_ELECSYS_FLAGS_RESTRICTED.csv');ef=ef.loc[ef.PTAU_EXACT]
 emap={(int(r.RID),str(r.date),str(r.batch)):int(r.source_row) for r in ef.itertuples()}
 p_rows={}
 for r in index.itertuples():
  if flags.loc[r.RID,'ELECSYS_PTAU_W0']:
   p_rows[int(r.RID)]=emap[(int(r.RID),str(r.EXAMDATE),str(flags.loc[r.RID,'PTAU_ANCHOR_BATCH_W0']))]
 targetlines=set(p_rows.values());pvalues={line:value(row['PTAU']) for line,row in records('E01') if line in targetlines}
 assert all(8<=p<=120 for p in pvalues.values()) and len(pvalues)==567
 targets={int(r.RID):day(r.EXAMDATE) for r in index.itertuples()}
 demo=defaultdict(list)
 for line,row in records('I13'):
  p=value(row.get('RID'));p=int(p) if np.isfinite(p) else None
  if p not in targets or row.get('HAS_QC_ERROR')=='1':continue
  m=re.fullmatch(r'(\d{1,2})/(\d{4})',row.get('PTDOB','').strip());full=day(row.get('PTDOB'));yd=day(row.get('PTDOBYY'));yn=value(row.get('PTDOBYY'))
  year=full.year if full else int(m[2]) if m else int(yn) if np.isfinite(yn) else yd.year if yd else None
  month=full.month if full else int(m[1]) if m else None
  demo[p].append({'year':year,'month':month,'sex':row.get('PTGENDER'),'education':value(row.get('PTEDUCAT')),'date':day(row.get('VISDATE')),'source_row':line})
 apo=defaultdict(set)
 for line,row in records('I14'):
  p=value(row.get('RID'));p=int(p) if np.isfinite(p) else None
  if p not in targets or row.get('APUSABLE')=='0':continue
  alleles=re.findall('[234]',row.get('GENOTYPE',''))
  if len(alleles)==2:apo[p].add(int('4' in alleles))
 amy=pd.read_csv(S4/'00_restricted/AMYLOID_PET_FLAGS_RESTRICTED.csv');amy=amy.loc[amy.CENTILOIDS_EXACT]
 ab=defaultdict(list)
 for r in amy.to_dict('records'):ab[int(r['RID'])].append(r)
 pet_selected={}
 for p,t in targets.items():
  opts=[r for r in ab[p] if abs((day(r['date'])-t).days)<=90]
  if opts:pet_selected[p]=min(opts,key=lambda r:(abs((day(r['date'])-t).days),day(r['date'])>t,r['_qc_rank'],r['source_row']))
 pet_lines={int(r['source_row']) for r in pet_selected.values()}
 petvalues={line:value(row.get('CENTILOIDS')) for line,row in records('I11') if line in pet_lines}
 assert len(pet_selected)==335 and all(np.isfinite(v) for v in petvalues.values())
 out=[]
 for r in index.to_dict('records'):
  p=int(r['RID']);t=targets[p];dr=demo[p];age=sex=education=e4=np.nan
  if cov.loc[p,'AGE_AVAILABLE']:
   years={d['year'] for d in dr if d['year']};months={d['month'] for d in dr if d['month']}
   assert len(years)==len(months)==1
   age=t.year-next(iter(years))+(t.month-next(iter(months)))/12
  if cov.loc[p,'SEX_AVAILABLE']:
   s={d['sex'] for d in dr if d['sex'] in ['1','2']};assert len(s)==1;sex=float(next(iter(s))=='2')
  if cov.loc[p,'EDUCATION_365_AVAILABLE']:
   opts=[d for d in dr if d['date'] and np.isfinite(d['education']) and 0<=d['education']<=30 and abs((d['date']-t).days)<=365]
   chosen=min(opts,key=lambda d:(abs((d['date']-t).days),d['date']>t,d['source_row']));education=chosen['education']
  if cov.loc[p,'APOE4_AVAILABLE']:assert len(apo[p])==1;e4=float(next(iter(apo[p])))
  pet=pet_selected.get(p)
  out.append({**r,'age_month_precision':age,'sex_code2':sex,'education':education,'apoe4_carrier':e4,
   'ptau':pvalues[p_rows[p]] if p in p_rows else np.nan,'ptau_batch':flags.loc[p,'PTAU_ANCHOR_BATCH_W0'] if p in p_rows else '',
   'amyloid_CL90':petvalues[int(pet['source_row'])] if pet else np.nan,'PET_tracer':pet['tracer'] if pet else '',
   'AT_operational':flags.loc[p,'PTAU_ANCHOR_AT_W0']})
 out=pd.DataFrame(out)
 out['base_primary']=np.isfinite(out[['ptau','age_month_precision','sex_code2']]).all(axis=1)
 assert int(out.base_primary.sum())==560
 for raw,z in [('ptau','P_RAW'),('ptau','P_LOG'),('amyloid_CL90','P_AMY')]:
  v=np.log2(out[raw]) if z=='P_LOG' else out[raw]
  base=out.base_primary if z!='P_AMY' else np.isfinite(out[['amyloid_CL90','age_month_precision','sex_code2']]).all(axis=1)
  mean=float(v[base].mean());sd=float(v[base].std(ddof=0));out[z]=(v-mean)/sd
  csvout(RESULTS/(z+'_SCALING_AGGREGATE.csv'),[{'variable':z,'n_scaling':int(base.sum()),'transform':'log2' if z=='P_LOG' else 'identity','parameters_location':'restricted only','ddof':0}])
  save_json(RESTRICTED/(z+'_SCALING_RESTRICTED.json'),{'mean':mean,'SD':sd,'n':int(base.sum())})
 csvout(RESTRICTED/'CLINICAL_DESIGN_RESTRICTED.csv',out)
 log('Clinical design prepared: 560 main clinical complete; no protein values opened')

if __name__=='__main__':verify();clinical()
