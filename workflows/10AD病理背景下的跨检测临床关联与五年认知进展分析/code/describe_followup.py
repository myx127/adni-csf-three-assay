from config import *
import pickle
def run():
 data=pickle.load((PRIVATE/'ASSEMBLED_CLINICAL_RESTRICTED.pkl').open('rb'));templates=pickle.load((PRIVATE/'ANALYSIS_TEMPLATES_RESTRICTED.pkl').open('rb'));c=data['clinical'];traj=data['trajectories'];pet=data['PET']
 summary=[];yearly=[];selection=[];change=[]
 for acc,gene in CANDIDATES.items():
  t=templates['MEM5Y__'+acc];ids=t['ids'];recs=[traj[(i,'ADNI_MEM',5)] for i in ids];w=np.array([r['w'] for r in recs]);order=np.argsort(w)[::-1]
  summary.append({'gene':gene,'n_people':len(ids),'n_observed_visits':sum(r['n_visits'] for r in recs),'median_visits':np.median([r['n_visits'] for r in recs]),'median_span_years':np.median([r['span'] for r in recs]),'n_highinfo':sum(r['highinfo'] for r in recs),'effective_time_weight_people':w.sum()**2/(w@w),'top10percent_time_weight_fraction':w[order[:max(1,int(np.ceil(.1*len(w))))]].sum()/w.sum(),'maximum_person_weight_fraction':w.max()/w.sum(),'n_last_visit_atleast4point5y':sum(r['last_year']>=4.5 for r in recs),'not_all_followed_full5years':True})
  for year in [1,2,3,4,5]:
   yearly.append({'gene':gene,'year':year,'n_with_last_visit_ge_year_minus_half':sum(r['last_year']>=year-.5 for r in recs),'denominator_main':len(recs)})
   differences=[]
   for r in recs:
    ts=r['t'];ys=r['y'];eligible=np.flatnonzero((ts>90/365.25)&(ts>=year-.5)&(ts<year+.5)&(ts<=5))
    if len(eligible):
     j=eligible[np.argmin(abs(ts[eligible]-year))];differences.append(float(ys[j]-r['baseline']))
   if differences:change.append({'gene':gene,'year_bin_center':year,'n_pairs':len(differences),'mean_observed_change':np.mean(differences),'median_observed_change':np.median(differences),'q25_observed_change':np.quantile(differences,.25),'q75_observed_change':np.quantile(differences,.75),'description':'Observed changing followup sets; no model-based or causal interpretation'})
  for group in ['ALL','HIGHINFO','LOWINFO','PET_AVAILABLE','PET_UNAVAILABLE']:
   sub=[(i,r) for i,r in zip(ids,recs) if group=='ALL' or group=='HIGHINFO' and r['highinfo'] or group=='LOWINFO' and not r['highinfo'] or group=='PET_AVAILABLE' and i in pet.index or group=='PET_UNAVAILABLE' and i not in pet.index]
   ii=[i for i,r in sub];rs=[r for i,r in sub]
   if not ii:continue
   frame=c.iloc[ii]
   selection.append({'gene':gene,'group':group,'n':len(ii),'mean_age':frame.age_month_precision.mean(),'female_fraction':frame.sex_code2.eq(1).mean(),'mean_education':frame.education.mean(),'APOE4_carrier_fraction':frame.apoe4_carrier.mean(),'mean_log2_ptau_globalSD':frame.P_LOG.mean(),'mean_baseline_MEM':np.mean([r['baseline'] for r in rs]),'mean_observed_MEM_slope_descriptive':np.mean([r['slope'] for r in rs]),'median_last_followup_years':np.median([r['last_year'] for r in rs]),'selection_audit_not_adjusted_association':True})
 csvout(RESULTS/'FOLLOWUP_TIME_INFORMATION_SUMMARY.csv',summary);csvout(RESULTS/'FOLLOWUP_REACH_COUNTS.csv',yearly);csvout(RESULTS/'FOLLOWUP_AND_PET_SELECTION_PROFILE.csv',selection);csvout(RESULTS/'OBSERVED_ANNUAL_MEMORY_CHANGE_DESCRIPTIVE.csv',change)
 log('Followup information weights and PET/long-followup selection profiles audited without candidate target regressions')
 print(pd.DataFrame(summary).to_string(index=False))
if __name__=='__main__':
 assert_unsealed();run()
