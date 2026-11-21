from __future__ import annotations
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from collections import defaultdict,Counter
import importlib.metadata,platform,sys
from audit_common import *
import clinical_reduce as reduce
from covariate_reduce import covariates,load_diagnosis

def match_records(index,frame,valid_col,window):
 byperson=defaultdict(list)
 if len(frame):
  for row in frame.loc[frame[valid_col].fillna(False)].to_dict('records'):byperson[int(row['RID'])].append(row)
 result={}
 for row in index.itertuples():
  target=parse_date(row.EXAMDATE);options=[]
  for candidate in byperson.get(int(row.RID),[]):
   delta=(candidate['date']-target).days
   if abs(delta)<=window:
    key=(abs(delta),delta>0,candidate.get('_qc_rank',0),candidate.get('batch_priority',0),candidate.get('source_row',0))
    options.append((key,candidate,delta))
  if options:
   _,selected,delta=min(options,key=lambda x:x[0]);result[int(row.RID)]={**selected,'delta_days':delta}
 return result

def cognition_coverage(index,cognition):
 valid=cognition.loc[cognition.ADNI_MEM_EXACT].copy();by=defaultdict(list)
 for row in valid.to_dict('records'):by[int(row['RID'])].append(row)
 grids=[];personrows=[]
 for baseline_window in [90,180]:
  baseline=match_records(index,valid,'ADNI_MEM_EXACT',baseline_window)
  for horizon in [3,5,999]:
   internal=[]
   for person in index.itertuples():
    ridv=int(person.RID);t0=parse_date(person.EXAMDATE);base=baseline.get(ridv)
    follow=sorted({r['date'] for r in by.get(ridv,[]) if (r['date']-t0).days>90 and (not base or r['date']>base['date']) and (r['date']-t0).days<=horizon*365.25})
    dates=([base['date']] if base else [])+follow
    span=(max(dates)-min(dates)).days/365.25 if len(dates)>=2 else 0.
    rec={'RID':ridv,'baseline_window_days':baseline_window,'horizon_years':horizon,
       'baseline_available':base is not None,'baseline_delta_days':base['delta_days'] if base else None,
       'n_followup_dates':len(follow),'n_total_dates':len(dates),'span_years':span,
       'first_followup_days':(follow[0]-t0).days if follow else None,'last_followup_days':(follow[-1]-t0).days if follow else None,
       'minimum_longitudinal_available':base is not None and len(follow)>=1,
       'high_information_available':base is not None and len(follow)>=2 and span>=2}
    internal.append(rec);personrows.append(rec)
   for visits in [2,3,4]:
    for span in [1,2,3]:
     n=sum(r['baseline_available'] and r['n_total_dates']>=visits and r['span_years']>=span for r in internal)
     grids.append({'baseline_window_days':baseline_window,'horizon_years':horizon,'minimum_total_visits':visits,
         'minimum_span_years':span,'n_participants':n,'denominator':len(index),'is_frozen_clinical_eligibility':False})
 write_csv(RESULTS/'COGNITION_FOLLOWUP_GRID.csv',grids)
 data=pd.DataFrame(personrows);write_csv(RESTRICTED/'COGNITION_TIMING_FLAGS_RESTRICTED.csv',data)
 primary=data.loc[data.baseline_window_days.eq(90)&data.horizon_years.eq(5)].set_index('RID')
 return primary

def run():
 log('START coverage audit; no protein-phenotype association and no raw clinical score outputs')
 initial=json.loads((RESULTS/'ACCESS_BOUNDARY_INITIAL.json').read_text(encoding='utf-8'))
 assert sha(ROOT/'AUDIT_PROTOCOL_v1.0.md')==initial['protocol_sha256_before_clinical_content']
 for row in json.loads((RESULTS/'INPUT_HASH_CHECKS.json').read_text(encoding='utf-8')):assert row['passed']
 thresholds=register_thresholds()
 save_json(RESULTS/'SOFTWARE_ENVIRONMENT.json',{'python':sys.version,'platform':platform.platform(),
    'packages':{p:importlib.metadata.version(p) for p in ['pandas','numpy','scipy','pyarrow','openpyxl','matplotlib','pypdf','Pillow']}})
 registry_path=STAGE1/'00_restricted/TRIPLE_SAMPLE_REGISTRY_RESTRICTED.csv'
 mask_path=STAGE1/'00_restricted/TRIPLE_CANDIDATE_AVAILABILITY_RESTRICTED.parquet'
 s1manifest=pd.read_csv(STAGE1/'STAGE1_INPUT_OUTPUT_SHA256.csv')
 for path in [registry_path,mask_path]:
  expected=s1manifest.loc[s1manifest.path.eq(path.relative_to(STAGE1).as_posix())].iloc[0]
  assert sha(path)==expected.sha256
 index=pd.read_csv(registry_path,dtype={'RID':int});assert len(index)==571 and index.RID.is_unique
 frames={'CSF_ELECSYS':reduce.load_csf('E01'),'CSF_ALZBIO3':reduce.load_csf('I10'),
        'AMYLOID_PET':reduce.load_pet('I11'),'TAU_PET':reduce.load_pet('I12')}
 cognition=reduce.load_cognition();dx=load_diagnosis();cov=covariates(index).set_index('RID')
 timing=cognition_coverage(index,cognition)
 for source,frame in frames.items():write_csv(RESTRICTED/(source+'_FLAGS_RESTRICTED.csv'),frame)
 write_csv(RESTRICTED/'COGNITION_DATE_FLAGS_RESTRICTED.csv',cognition)
 write_csv(RESTRICTED/'DIAGNOSIS_DATE_FLAGS_RESTRICTED.csv',dx)
 write_csv(RESTRICTED/'COVARIATE_FLAGS_RESTRICTED.csv',cov.reset_index())
 out=index[['RID','EXAMDATE','key','NULISA_PlateID','SOMASCAN_PlateId','TMT_plex']].copy().set_index('RID')
 for c in cov.columns:out[c]=cov[c]
 out['COG_BASELINE_90']=timing.baseline_available
 out['COG_MINIMUM_5Y']=timing.minimum_longitudinal_available
 out['COG_HIGH_INFO_5Y']=timing.high_information_available
 routes=[]
 elec=frames['CSF_ELECSYS'];master=frames['CSF_ALZBIO3'];amy=frames['AMYLOID_PET'];tau=frames['TAU_PET']
 def route(name,frame,column):routes.append((name,frame,column))
 route('ELECSYS_PTAU',elec,'PTAU_EXACT');route('ELECSYS_ABETA42',elec,'ABETA42_EXACT')
 route('ELECSYS_A4240_RATIO',elec,'A4240_EXACT');route('ELECSYS_PTAU_ABETA42_RATIO',elec,'PTAU_ABETA42_EXACT')
 route('ELECSYS_A_CLASS',elec,'A_VALID');route('ELECSYS_T27_CLASS',elec,'T_VALID');route('ELECSYS_AT_CLASS',elec,'AT_VALID')
 route('ELECSYS_BATCH9_PTAU',elec.loc[elec.batch.eq('UPENNBIOMK9')],'PTAU_EXACT')
 route('ALZBIO3_PTAU_MEDIAN',master.loc[master.batch.eq('MEDIAN')],'PTAU_EXACT')
 route('ALZBIO3_PTAU_ORIGINAL',master.loc[master.batch.ne('MEDIAN')],'PTAU_EXACT')
 route('ALZBIO3_AT_MEDIAN',master.loc[master.batch.eq('MEDIAN')],'AT_VALID')
 route('AMY_PET_CENTILOIDS',amy,'CENTILOIDS_EXACT');route('AMY_PET_A_RELEASE',amy,'A_VALID')
 route('AMY_PET_CENTILOIDS_QC2',amy.loc[amy.qc_full_pass],'CENTILOIDS_EXACT')
 route('TAU_PET_ANY_TRACER',tau,'META_TEMPORAL_SUVR_EXACT')
 for tracer in ['FTP','MK6240','PI2620']:route('TAU_PET_'+tracer,tau.loc[tau.tracer.eq(tracer)],'META_TEMPORAL_SUVR_EXACT')
 route('DX_CATEGORY',dx,'DX_VALID')
 windowrows=[];categoryrows=[];selected={}
 for name,frame,col in routes:
  for window in WINDOWS:
   matches=match_records(index,frame,col,window);selected[(name,window)]=matches
   flag=f'{name}_W{window}';out[flag]=out.index.isin(matches)
   deltas=[r['delta_days'] for r in matches.values()]
   minimal=sum(bool(cov.loc[p,'COV_MINIMAL']) for p in matches)
   extended=sum(bool(cov.loc[p,'COV_COGNITIVE']) for p in matches)
   windowrows.append({'route':name,'window_days':window,'n_available':len(matches),'n_with_age_sex':minimal,
      'n_with_age_sex_education_apoe':extended,'n_same_day':sum(d==0 for d in deltas),
      'n_before_index':sum(d<0 for d in deltas),'n_after_index':sum(d>0 for d in deltas),
      'median_absolute_days':float(np.median(np.abs(deltas))) if deltas else None,'denominator':len(index)})
   if name in ['ELECSYS_AT_CLASS','ALZBIO3_AT_MEDIAN','ELECSYS_PTAU','AMY_PET_A_RELEASE','DX_CATEGORY']:
    catcol='DX_CAT' if name=='DX_CATEGORY' else 'A_CAT' if name=='AMY_PET_A_RELEASE' else 'AT_CAT'
    counter=Counter(r.get(catcol,'UNKNOWN') for r in matches.values())
    for category,n in counter.items():categoryrows.append({'route':name,'window_days':window,'category':category,'n_participants':n})
 write_csv(RESULTS/'PATHOLOGY_WINDOW_COVERAGE.csv',windowrows);write_csv(RESULTS/'OPERATIONAL_CATEGORY_COUNTS.csv',categoryrows)
 # Same-record operational AT labels associated with the selected exact continuous PTAU anchor.
 for window in [0,90,180]:
  matches=selected[('ELECSYS_PTAU',window)]
  out[f'PTAU_ANCHOR_AT_W{window}']=[matches.get(int(p),{}).get('AT_CAT','UNKNOWN') for p in out.index]
  out[f'PTAU_ANCHOR_DELTA_W{window}']=[matches.get(int(p),{}).get('delta_days',np.nan) for p in out.index]
  out[f'PTAU_ANCHOR_BATCH_W{window}']=[matches.get(int(p),{}).get('batch','') for p in out.index]
 out['COV_MINIMAL']=out.COV_MINIMAL.astype(bool);out['COV_COGNITIVE']=out.COV_COGNITIVE.astype(bool)
 out['PTAU90_MINCOV']=out.ELECSYS_PTAU_W90&out.COV_MINIMAL
 out['PTAU90_FULLCOV']=out.ELECSYS_PTAU_W90&out.COV_COGNITIVE
 out['PTAU90_MINCOV_COGMIN']=out.PTAU90_MINCOV&out.COG_MINIMUM_5Y
 out['PTAU90_FULLCOV_COGHIGH']=out.PTAU90_FULLCOV&out.COG_HIGH_INFO_5Y
 write_csv(RESTRICTED/'INDEX_CLINICAL_COVERAGE_FLAGS_RESTRICTED.csv',out.reset_index())
 joins=[]
 for name in ['COG_BASELINE_90','COG_MINIMUM_5Y','COG_HIGH_INFO_5Y','COV_MINIMAL','COV_COGNITIVE','PTAU90_MINCOV','PTAU90_FULLCOV','PTAU90_MINCOV_COGMIN','PTAU90_FULLCOV_COGHIGH']:
  joins.append({'flag':name,'n_participants':int(out[name].sum()),'denominator':571})
 write_csv(RESULTS/'JOINT_COHORT_COVERAGE.csv',joins)
 # Source timing is structural information; no clinical numerical distributions.
 native=[]
 for source,frame in {**frames,'COGNITION':cognition,'DIAGNOSIS':dx}.items():
  native.append({'source':source,'reduced_rows':len(frame),'unique_RID':frame.RID.nunique(),
      'earliest_date':str(frame.date.min()),'latest_date':str(frame.date.max())})
 write_csv(RESULTS/'SOURCE_TIMING_STRUCTURE.csv',native)
 write_csv(RESULTS/'SOURCE_ROW_STRUCTURE.csv',reduce.STRUCTURE);write_csv(RESULTS/'DUPLICATE_AND_CONFLICT_AUDIT.csv',reduce.DUPLICATES)
 write_csv(RESULTS/'CLINICAL_VALUE_STATUS_COUNTS.csv',reduce.VALUE_STATUS);save_json(RESULTS/'FIELD_ACCESS_LOG.json',reduce.FIELD_ACCESS)
 save_json(RESULTS/'COVERAGE_RUN_SUMMARY.json',{'status':'COVERAGE_REDUCTION_COMPLETE','index_n':571,
    'clinical_inputs_reduced_to_flags':['I09','I10','I11','I12','I13','I14','I15','E01'],
    'raw_continuous_clinical_values_saved':False,'protein_quantitative_values_joined':False,'association_models_run':0,
    'threshold_registry_sha256':sha(RESULTS/'THRESHOLD_REGISTRY.csv'),
    'source_frames':{s:len(f) for s,f in frames.items()},'joint_counts':joins})
 log('Coverage flags and count tables completed; no clinical effect models')
 print(pd.DataFrame(windowrows).query("window_days in [0,90,365]")[['route','window_days','n_available','n_with_age_sex']].to_string(index=False))
 print(pd.DataFrame(joins).to_string(index=False))

if __name__=='__main__':run()
