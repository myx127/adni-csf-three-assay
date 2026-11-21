from config import *
from collections import defaultdict
import pickle

def trajectory(rows,index_date,horizon):
 valid=rows.loc[rows.value.notna()].copy()
 if valid.empty:return None
 valid['day']=(valid.date-index_date).dt.days
 candidates=valid.loc[valid.day.abs().le(90)].copy()
 if candidates.empty:return None
 candidates['after']=candidates.day.gt(0);candidates['absday']=candidates.day.abs()
 base=candidates.sort_values(['absday','after','source_row']).iloc[0]
 follow=valid.loc[valid.day.gt(90)&valid.date.gt(base.date)&valid.day.le(horizon*365.25)]
 part=pd.concat([base.to_frame().T,follow]).sort_values('date').drop_duplicates('date')
 t=part.day.to_numpy(float)/365.25;y=part.value.to_numpy(float)
 if len(t)<2:return {'baseline':float(base.value),'base_day':int(base.day),'n_visits':len(t),'eligible':False,'t':t,'y':y}
 tc=t-t.mean();qc=t*t-np.mean(t*t);w=float(tc@tc)
 if w<=1e-10:return None
 h=float(tc@qc);q=float(qc@qc);sty=float(tc@y);sqy=float(qc@y)
 ar=np.exp(-.6*abs(t[:,None]-t[None,:]));var_ar=float(tc@ar@tc/w)
 return {'baseline':float(base.value),'base_day':int(base.day),'n_visits':len(t),'eligible':True,'n_followups':len(t)-1,'span':float(t.max()-t.min()),'last_year':float(t.max()),'w':w,'h':h,'q':q,'sty':sty,'sqy':sqy,'slope':sty/w,'var_AR_pseudo':var_ar,'highinfo':len(t)>=3 and t.max()-t.min()>=2,'t':t,'y':y}

def run():
 rec=RESULTS/'ACCESS_AND_PROTOCOL_RECORD.json';access=json.loads(rec.read_text(encoding='utf-8'))
 assert access['protocol_sha256']==sha(ROOT/'STAGE10_PROTOCOL_v1.0.md')
 first=RESULTS/'FIRST_NEW_CLINICAL_VALUE_ACCESS.json'
 if not first.exists():save_json(first,{'time':now(),'protocol_sha256':sha(ROOT/'STAGE10_PROTOCOL_v1.0.md'),'purpose':'Assemble fixed cognitive/PET values and design geometry; no candidate clinical coefficient or P value fitted'})
 clinical=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv');assert len(clinical)==571 and clinical.RID.is_unique
 source=pd.read_csv(INPUTS['I09']['path'],usecols=['RID','EXAMDATE','ADNI_MEM','ADNI_EF'],dtype=str)
 source['source_row']=np.arange(len(source))+2;source.RID=pd.to_numeric(source.RID,errors='coerce')
 source['date']=pd.to_datetime(source.EXAMDATE,errors='coerce',format='mixed')
 source=source.loc[source.RID.isin(clinical.RID)&source.date.notna()].copy()
 for marker in ['ADNI_MEM','ADNI_EF']:source[marker]=pd.to_numeric(source[marker],errors='coerce')
 reduced=[];duplicates=[]
 for (rid,day),part in source.groupby(['RID','date'],sort=True):
  row={'RID':int(rid),'date':day,'source_row':int(part.source_row.min())}
  for marker in ['ADNI_MEM','ADNI_EF']:
   values=part[marker].dropna().unique();row[marker]=float(values[0]) if len(values)==1 else np.nan
   if len(part)>1:duplicates.append({'RID':rid,'date':str(day.date()),'marker':marker,'n_rows':len(part),'conflict':len(values)>1})
  reduced.append(row)
 cog=pd.DataFrame(reduced);csvout(PRIVATE/'COGNITIVE_VISITS_SELECTED_RESTRICTED.csv',cog)
 trajectories={};person=[];schedule=[]
 for idx,r in clinical.iterrows():
  part=cog.loc[cog.RID.eq(r.RID)];date=pd.Timestamp(r.EXAMDATE)
  for marker in ['ADNI_MEM','ADNI_EF']:
   for horizon in ([3,5] if marker=='ADNI_MEM' else [5]):
    data=trajectory(part.rename(columns={marker:'value'})[['date','value','source_row']],date,horizon)
    if data:
     ident=(idx,marker,horizon);trajectories[ident]=data
     row={'index':idx,'RID':r.RID,'marker':marker,'horizon':horizon,**{k:v for k,v in data.items() if k not in ['t','y']}}
     person.append(row)
     for t,y in zip(data['t'],data['y']):schedule.append({'index':idx,'RID':r.RID,'marker':marker,'horizon':horizon,'time_years':t,'score':y})
 csvout(PRIVATE/'COGNITIVE_SUFFICIENT_STATISTICS_RESTRICTED.csv',person)
 csvout(PRIVATE/'COGNITIVE_SCHEDULE_RESTRICTED.csv',schedule)
 # Independently recover source PET values using Stage4 metadata selection and compare Stage5 CL values.
 petraw=pd.read_csv(INPUTS['I11']['path'],usecols=['RID','SCANDATE','CENTILOIDS','AMYLOID_STATUS','TRACER'],dtype=str)
 flags=pd.read_csv(S4/'00_restricted/AMYLOID_PET_FLAGS_RESTRICTED.csv')
 flags['date']=pd.to_datetime(flags.date);petrows=[];petchecks=[]
 for idx,r in clinical.iterrows():
  candidates=flags.loc[flags.RID.eq(r.RID)&flags.CENTILOIDS_EXACT.astype(bool)].copy()
  candidates['delta']=(candidates.date-pd.Timestamp(r.EXAMDATE)).dt.days
  candidates=candidates.loc[candidates.delta.abs().le(90)]
  if candidates.empty:continue
  candidates['absdelta']=candidates.delta.abs();candidates['after']=candidates.delta.gt(0)
  sort=['absdelta','after']+[c for c in ['_qc_rank','source_row'] if c in candidates]
  chosen=candidates.sort_values(sort).iloc[0];raw=petraw.iloc[int(chosen.source_row)-2]
  cl=float(raw.CENTILOIDS);category=pd.to_numeric(raw.AMYLOID_STATUS,errors='coerce')
  petrows.append({'index':idx,'RID':r.RID,'CL':cl,'A_release':int(category) if category in [0,1] else np.nan,'tracer':str(raw.TRACER),'delta_days':int(chosen.delta),'source_row':int(chosen.source_row),'qc':str(chosen.qc)})
  petchecks.append({'index_number':idx,'same_stage5_CL':bool(np.isclose(cl,r.amyloid_CL90,rtol=1e-12,atol=1e-12)),'same_tracer':str(raw.TRACER)==str(r.PET_tracer)})
 pet=pd.DataFrame(petrows).set_index('index');assert all(v['same_stage5_CL'] and v['same_tracer'] for v in petchecks)
 csvout(PRIVATE/'AMYLOID_PET_VALUES_RESTRICTED.csv',pet.reset_index());csvout(RESULTS/'PET_SOURCE_REPRODUCTION_CHECKS.csv',petchecks)
 save_json(RESULTS/'CLINICAL_ASSEMBLY_SUMMARY.json',{'analysis_date':'2026-10-29','main_index':571,'cognitive_unique_RID_dates':len(cog),'PET90_n':len(pet),'PET90_A_positive':int(pet.A_release.eq(1).sum()),'PET90_A_negative':int(pet.A_release.eq(0).sum()),'new_associations_fitted':0,'cognitive_version':'UWNPSYCHSUM fixed2023 composites; source file2026','outcome_fields':['ADNI_MEM','ADNI_EF'],'duplicate_date_domain_groups':len(duplicates),'duplicate_conflicts':sum(r['conflict'] for r in duplicates)})
 with (PRIVATE/'ASSEMBLED_CLINICAL_RESTRICTED.pkl').open('wb') as f:pickle.dump({'clinical':clinical,'trajectories':trajectories,'PET':pet},f,pickle.HIGHEST_PROTOCOL)
 coverage=[]
 for marker in ['ADNI_MEM','ADNI_EF']:
  for horizon in ([3,5] if marker=='ADNI_MEM' else [5]):
   ds=[v for (i,m,h),v in trajectories.items() if m==marker and h==horizon]
   eligible=[d for d in ds if d.get('eligible')]
   coverage.append({'marker':marker,'horizon_years':horizon,'baseline':len(ds),'baseline_plus_followup':len(eligible),'highinfo':sum(d['highinfo'] for d in eligible),'median_visits':np.median([d['n_visits'] for d in eligible]),'median_span_years':np.median([d['span'] for d in eligible])})
 csvout(RESULTS/'CLINICAL_COVERAGE_RECOUNT.csv',coverage)
 log('Fixed clinical values/schedules assembled;335 PET source values reproduce Stage5;no candidate clinical coefficient fitted')
 print(pd.DataFrame(coverage).to_string(index=False))
if __name__=='__main__':
 assert_unsealed();run()
