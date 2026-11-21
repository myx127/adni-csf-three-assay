from collections import defaultdict,Counter
from audit_common import *
from run_coverage import match_records

index=pd.read_csv(STAGE1/'00_restricted/TRIPLE_SAMPLE_REGISTRY_RESTRICTED.csv')
cog=pd.read_csv(RESTRICTED/'COGNITION_DATE_FLAGS_RESTRICTED.csv');cog['date']=cog.date.map(parse_date)
out=[]
for domain in ['ADNI_MEM','ADNI_EF','ADNI_EF2','ADNI_LAN','ADNI_VS']:
 flag=domain+'_EXACT';valid=cog.loc[cog[flag]];by=defaultdict(set)
 for row in valid.itertuples():by[int(row.RID)].add(row.date)
 for w in [90,180]:
  baseline=match_records(index,valid,flag,w)
  for horizon in [3,5,999]:
   minimum=high=0
   for row in index.itertuples():
    p=int(row.RID);t=parse_date(row.EXAMDATE);base=baseline.get(p)
    follow={d for d in by[p] if 90<(d-t).days<=365.25*horizon and (not base or d>base['date'])}
    dates=follow|({base['date']} if base else set())
    span=(max(dates)-min(dates)).days/365.25 if len(dates)>=2 else 0
    minimum+=bool(base and follow);high+=bool(base and len(follow)>=2 and span>=2)
   out.append({'domain':domain,'baseline_window_days':w,'horizon_years':horizon,'n_baseline':len(baseline),
       'n_baseline_plus_one_followup':minimum,'n_baseline_plus_two_followups_span2y':high,'denominator':571})
write_csv(RESULTS/'COGNITIVE_DOMAIN_COVERAGE.csv',out)
flags=pd.read_csv(RESTRICTED/'INDEX_CLINICAL_COVERAGE_FLAGS_RESTRICTED.csv')
timing=pd.read_csv(RESTRICTED/'COGNITION_TIMING_FLAGS_RESTRICTED.csv').query('baseline_window_days==90 and horizon_years==5')
structure=[]
for col in ['n_followup_dates','n_total_dates','span_years','first_followup_days','last_followup_days']:
 values=timing[col].dropna();structure.append({'timing_metric':col,'min':float(values.min()),'median':float(values.median()),'max':float(values.max()),'n':len(values)})
write_csv(RESULTS/'COGNITIVE_TIMING_STRUCTURE.csv',structure)
rows_out=[]
for name,mask in [('TRIPLE571',pd.Series(True,index=flags.index)),('PTAU0_MINCOV',flags.PTAU90_MINCOV),
                 ('PTAU0_MINCOV_COGMIN',flags.PTAU90_MINCOV_COGMIN),('PTAU0_FULLCOV_COGHIGH',flags.PTAU90_FULLCOV_COGHIGH)]:
 for cat,n in flags.loc[mask,'PTAU_ANCHOR_AT_W0'].value_counts(dropna=False).items():
  rows_out.append({'cohort':name,'operational_AT_category':str(cat),'n':int(n)})
write_csv(RESULTS/'COHORT_OPERATIONAL_AT_COUNTS.csv',rows_out)
amy=pd.read_csv(RESTRICTED/'AMYLOID_PET_FLAGS_RESTRICTED.csv');amy['date']=amy.date.map(parse_date)
chosen=match_records(index,amy,'CENTILOIDS_EXACT',90)
comparison=Counter((r.get('A_CAT','UNKNOWN'),r.get('CL18_CAT','UNKNOWN'),bool(r.get('CL18_EQUAL',False))) for r in chosen.values())
write_csv(RESULTS/'AMY_RELEASE_VS_CL18_DEFINITIONS.csv',[{'released_category':a,'CL18_category':c,'exactly18CL':e,'n':n,
 'interpretation':'Definition comparison only; not protein association'} for (a,c,e),n in comparison.items()])
print(pd.DataFrame(out).query('baseline_window_days==90 and horizon_years==5').to_string(index=False))
