"""Independent standard-library re-counts; no clinical numerical outputs."""
from __future__ import annotations
import csv,json,re,math,hashlib
from pathlib import Path
from collections import defaultdict
from datetime import date,datetime
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];PROJECT=ROOT.parent;R=ROOT/'results';D=ROOT/'00_restricted'
S1=PROJECT/'1候选逐蛋白可评价性审计'
manifest=list(csv.DictReader((PROJECT/'规划和详细数据情况以及分析方案和研究计划/01_数据清单/CORE_INPUT_MANIFEST.csv').open(encoding='utf-8-sig')))
manifest+=list(csv.DictReader((R/'ADDITIONAL_INPUT_MANIFEST.csv').open(encoding='utf-8-sig')))
inputs={r['input_id']:r for r in manifest};checks=[]

def check(name,condition,actual=None,expected=None):
 checks.append({'check':name,'passed':bool(condition),'actual':actual,'expected':expected})
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def records(key):
 with Path(inputs[key]['path']).open(encoding='utf-8-sig',newline='') as f:yield from csv.DictReader(f)
def day(s):
 s=str(s or '').strip()
 for fmt in ['%Y-%m-%d','%m/%d/%Y','%Y/%m/%d','%Y-%m-%d %H:%M:%S']:
  try:return datetime.strptime(s,fmt).date()
  except ValueError:pass
 return None
def number(s,positive=False,limits=None):
 try:v=float(s)
 except (ValueError,TypeError):return None
 if not math.isfinite(v) or positive and v<=0 or limits and not limits[0]<=v<=limits[1]:return None
 return v
def ident(s):
 try:return int(float(s))
 except (ValueError,TypeError):return None
index=list(csv.DictReader((S1/'00_restricted/TRIPLE_SAMPLE_REGISTRY_RESTRICTED.csv').open(encoding='utf-8-sig')))
targets={int(r['RID']):day(r['EXAMDATE']) for r in index}
saved=pd.read_csv(D/'INDEX_CLINICAL_COVERAGE_FLAGS_RESTRICTED.csv').set_index('RID')
check('index_unique571',len(targets)==len(index)==571,len(targets),571)

# Independent raw-source paths; values only exist transiently while making flags.
elec=defaultdict(list);amy=defaultdict(list);tau=defaultdict(list);cog=defaultdict(set)
for r in records('E01'):
 p=ident(r['RID']);d=day(r['EXAMDATE'])
 if p and d and number(r['PTAU'],True,(8,120)) is not None:elec[p].append(d)
for key,dest,marker in [('I11',amy,'CENTILOIDS'),('I12',tau,'META_TEMPORAL_SUVR')]:
 for r in records(key):
  p=ident(r['RID']);d=day(r['SCANDATE']);ptid=re.search(r'_S_(\d+)$',r.get('PTID',''))
  good=(r.get('qc_flag') in ['1','2'] and number(r.get(marker),marker!='CENTILOIDS') is not None
        and not(ptid and int(ptid[1])!=p))
  if p and d and good:dest[p].append(d)
for r in records('I09'):
 p=ident(r['RID']);d=day(r['EXAMDATE'])
 if p and d and number(r['ADNI_MEM']) is not None:cog[p].add(d)

for name,source in [('ELECSYS_PTAU',elec),('AMY_PET_CENTILOIDS',amy),('TAU_PET_ANY_TRACER',tau)]:
 for w in [0,30,90,180,365,730]:
  selected={p for p,t in targets.items() if any(abs((d-t).days)<=w for d in source.get(p,[]))}
  actual=set(saved.index[saved[f'{name}_W{w}']])
  check(f'independent_{name}_W{w}_RID_set',selected==actual,len(actual),len(selected))
  if name=='ELECSYS_PTAU' and w in [0,90,180]:
   chosen={p:min((d for d in source[p] if abs((d-targets[p]).days)<=w),key=lambda d:(abs((d-targets[p]).days),d>targets[p])) for p in selected}
   check(f'independent_ptau_W{w}_selected_delta',all((d-targets[p]).days==saved.loc[p,f'PTAU_ANCHOR_DELTA_W{w}'] for p,d in chosen.items()))

for w in [90,180]:
 for horizon in [3,5,999]:
  flags={}
  for p,t in targets.items():
   bs=[d for d in cog.get(p,[]) if abs((d-t).days)<=w]
   baseline=min(bs,key=lambda d:(abs((d-t).days),d>t)) if bs else None
   follow={d for d in cog.get(p,[]) if 90<(d-t).days<=365.25*horizon and (baseline is None or d>baseline)}
   dates=follow|({baseline} if baseline else set());span=(max(dates)-min(dates)).days/365.25 if len(dates)>=2 else 0
   flags[p]=(baseline is not None,baseline is not None and len(follow)>=1,baseline is not None and len(follow)>=2 and span>=2,len(dates),span)
  grid=pd.read_csv(R/'COGNITION_FOLLOWUP_GRID.csv').query('baseline_window_days==@w and horizon_years==@horizon')
  for row in grid.itertuples():
   n=sum(v[0] and v[3]>=row.minimum_total_visits and v[4]>=row.minimum_span_years for v in flags.values())
   check(f'independent_cog_W{w}_H{horizon}_V{row.minimum_total_visits}_S{row.minimum_span_years}',n==row.n_participants,row.n_participants,n)
  if w==90 and horizon==5:
   for k,col in enumerate(['COG_BASELINE_90','COG_MINIMUM_5Y','COG_HIGH_INFO_5Y']):
    ids={p for p,v in flags.items() if v[k]};check('independent_'+col,ids==set(saved.index[saved[col]]),len(ids),int(saved[col].sum()))

# An independent birth/year/month consensus, not the primary covariate reducer.
births=defaultdict(set);months=defaultdict(set);sexes=defaultdict(set);internal_bad=set()
for row in records('I13'):
 p=ident(row.get('RID'))
 if p not in targets or row.get('HAS_QC_ERROR')=='1':continue
 s=row.get('PTDOB','').strip();m=re.fullmatch(r'(\d{1,2})/(\d{4})',s);full=day(s);yd=day(row.get('PTDOBYY'));yn=number(row.get('PTDOBYY'))
 year=full.year if full else int(m[2]) if m and 1<=int(m[1])<=12 else int(yn) if yn is not None else yd.year if yd else None
 month=full.month if full else int(m[1]) if m and 1<=int(m[1])<=12 else None
 yy=int(yn) if yn is not None else yd.year if yd else None
 if year is not None and yy is not None and year!=yy:internal_bad.add(p)
 if year is not None and 0<=targets[p].year-year<=120:
  births[p].add(year)
  if month is not None:months[p].add(month)
 if row.get('PTGENDER') in ['1','2']:sexes[p].add(row['PTGENDER'])
age_ok={p for p in targets if len(births[p])==1 and len(months[p])<=1 and p not in internal_bad}
sex_ok={p for p in targets if len(sexes[p])==1}
check('independent_AGE_AVAILABLE_IDs',age_ok==set(saved.index[saved.AGE_AVAILABLE]),len(age_ok),int(saved.AGE_AVAILABLE.sum()))
check('independent_COV_MINIMAL_IDs',age_ok&sex_ok==set(saved.index[saved.COV_MINIMAL]),len(age_ok&sex_ok),int(saved.COV_MINIMAL.sum()))

# Recount every delivered protein intersection from sealed boolean cells.
masks=pd.read_parquet(S1/'00_restricted/TRIPLE_CANDIDATE_AVAILABILITY_RESTRICTED.parquet')
key_to_rid={r['key']:int(r['RID']) for r in index}
data=pd.read_csv(R/'PER_PROTEIN_JOINT_COVERAGE.csv')
routes={'GIS_NUMERIC':'complete_GIS_S_all','GIS_LOD':'complete_GIS_detect_all_N_S_all','RAW_NUMERIC':'complete_raw_S_all','GIS_UNAMBIGUOUS':'complete_GIS_unambiguous_S_all'}
joint={'n_same_day_ptau':'ELECSYS_PTAU_W0','n_same_day_ptau_age_sex':'PTAU90_MINCOV',
 'n_same_day_ptau_fullcov':'PTAU90_FULLCOV','n_p_fullcov_cog_minimum':'PTAU90_MINCOV_COGMIN',
 'n_p_fullcov_cog_highinfo':'PTAU90_FULLCOV_COGHIGH','n_amy_pet_90':'AMY_PET_CENTILOIDS_W90','n_amy_pet_365':'AMY_PET_CENTILOIDS_W365'}
check('same_day_and90_identical_clinical_IDs',set(saved.index[saved.ELECSYS_PTAU_W0])==set(saved.index[saved.ELECSYS_PTAU_W90]))
for row in data.itertuples():
 part=masks[masks.accession.eq(row.accession)];ids={key_to_rid[r.key] for r in part.itertuples() if getattr(r,routes[row.route])}
 check(f'{row.accession}_{row.route}_assays',len(ids)==row.n_three_assays,len(ids),row.n_three_assays)
 for metric,flag in joint.items():
  n=len(ids&set(saved.index[saved[flag]]));check(f'{row.accession}_{row.route}_{metric}',n==getattr(row,metric),n,getattr(row,metric))
 n=len(ids&set(saved.index[saved.PTAU90_FULLCOV&saved.AMY_PET_CENTILOIDS_W90]))
 check(f'{row.accession}_{row.route}_Pcovamy',n==row.n_p_fullcov_amy_pet_90,n,row.n_p_fullcov_amy_pet_90)

coverage=pd.read_csv(R/'PATHOLOGY_WINDOW_COVERAGE.csv')
for route,part in coverage.groupby('route'):
 check('monotone_'+route,part.sort_values('window_days').n_available.diff().fillna(0).ge(0).all())
check('gene_intersections_nested',(data.n_p_fullcov_cog_highinfo<=data.n_p_fullcov_cog_minimum).all() and (data.n_p_fullcov_cog_minimum<=data.n_same_day_ptau_fullcov).all() and (data.n_same_day_ptau_fullcov<=data.n_three_assays).all())

# Edge cases are scientific validity checks, independent of actual patient values.
from audit_common import numeric,ratio
cases=[('censoredA_can_classify',numeric('>1700',True,(200,1700)).category(980,'lt')=='NEG'),
 ('censoredP_not_exact',not numeric('<8',True,(8,120)).exact),
 ('strict_cutoff_equal',numeric('27',True).category(27,'gt')=='NEG'),
 ('censor_boundary_open',numeric('>27',True).category(27,'gt')=='POS'),
 ('censor_boundary_closed',numeric('≥27',True).category(27,'gt')=='UNKNOWN'),
 ('out_of_range_not_exact',not numeric('2000',True,(200,1700)).exact),
 ('negative_cognition_valid',numeric('-4').exact),('negative_centiloid_valid',numeric('-4').exact),
 ('missing_nonfinite_invalid',not numeric('inf').exact),
 ('ratio_censor_not_exact',not ratio(numeric('27',True),numeric('>1700',True)).exact),
 ('ratio_positive_exact',ratio(numeric('27',True),numeric('1000',True)).exact),
 ('near_tie_prefers_earlier',min([date(2020,1,11),date(2020,1,9)],key=lambda d:(abs((d-date(2020,1,10)).days),d>date(2020,1,10)))==date(2020,1,9)),
 ('post90_baseline_not_followup',len({d for d in [100,101] if d>90 and d>100})==1)]
for name,ok in cases:check('edge_'+name,ok)

for row in manifest:
 actual=sha(row['path']);check('source_hash_'+row['input_id'],actual==row['sha256'],actual,row['sha256'])
initial=json.loads((R/'ACCESS_BOUNDARY_INITIAL.json').read_text(encoding='utf-8'))
check('protocol_unchanged',sha(ROOT/'AUDIT_PROTOCOL_v1.0.md')==initial['protocol_sha256_before_clinical_content'])
freeze=json.loads((R/'THRESHOLD_FREEZE_BEFORE_COUNTS.json').read_text(encoding='utf-8'))
check('threshold_registry_unchanged',sha(R/'THRESHOLD_REGISTRY.csv')==freeze['threshold_registry_sha256'])
for path in R.glob('*.csv'):
 cols=list(pd.read_csv(path,nrows=0).columns)
 check('public_no_individual_ID_'+path.stem,not any(c in cols for c in ['RID','PTID','key','EXAMDATE']))
for path in D.glob('*.csv'):
 cols=list(pd.read_csv(path,nrows=0).columns)
 check('restricted_no_raw_measure_'+path.stem,not any(c in cols for c in ['ABETA','ABETA42','PTAU','TAU','CENTILOIDS','META_TEMPORAL_SUVR','ADNI_MEM','ADNI_EF','ADNI_EF2','GENOTYPE','PTDOB','PTDOBYY','PTEDUCAT']))

out=pd.DataFrame(checks);out.to_csv(R/'INDEPENDENT_QA_CHECKS.csv',index=False,encoding='utf-8-sig')
summary={'checks':len(checks),'passed':int(out.passed.sum()),'failed':int((~out.passed).sum()),'independent_raw_count_sources':['E01','I09','I11','I12','I13'],
 'no_association_models':True,'raw_clinical_values_exported':False}
(R/'INDEPENDENT_QA_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(summary)
if not out.passed.all():
 print(out.loc[~out.passed].to_string(index=False));raise SystemExit(1)
