from config import *
from scipy.linalg import lstsq
import pickle
def run():
 full=pd.read_csv(RESULTS/'CALIBRATION_AND_POWER_FULL.csv');dec=pd.read_csv(RESULTS/'CALIBRATION_QUALIFICATION_DECISIONS.csv').set_index('design_id')
 null=full.loc[full.true_null].copy();rows=[]
 for keys,g in null.groupby(['design_id','gene','branch','formal','test','alpha_name'],sort=False):
  rr=g.loc[g.Wilson95_high.idxmax()]
  rows.append(dict(zip(['design_id','gene','branch','formal','test','alpha_name'],keys))|{'alpha':rr.alpha,'worst_scenario':rr.scenario,'worst_pattern':rr.pattern,'max_false_rate':rr.rate,'max_Wilson95_upper':rr.Wilson95_high,'qualified':bool(dec.loc[keys[0],'qualified'])})
 csvout(RESULTS/'NULL_CALIBRATION_WORST_SUMMARY.csv',rows)
 power=full.loc[full.kind.eq('POWER')].copy();rows=[]
 for keys,g in power.groupby(['design_id','gene','branch','pattern','test','alpha_name'],sort=False):
  rr=g.loc[g.rate.idxmin()]
  rows.append(dict(zip(['design_id','gene','branch','pattern','test','alpha_name'],keys))|{'worst_DGP_rate':rr.rate,'Wilson95_low':rr.Wilson95_low,'Wilson95_high':rr.Wilson95_high,'worst_DGP':rr.scenario,'qualified':bool(dec.loc[keys[0],'qualified']),'sufficient_threshold_not_actual_BY_rank_power':keys[-1]!='NOMINAL'})
 pp=pd.DataFrame(rows);csvout(RESULTS/'CONDITIONAL_POWER_WORST_DGP.csv',pp)
 balanced=pp.loc[pp.branch.eq('MEM5Y')&pp.test.eq('SAME')&pp.pattern.str.startswith('POWER_BALANCED_')&pp.alpha_name.isin(['NOMINAL','BY68_SUFFICIENT'])].copy()
 balanced['working_effect_abs']=balanced.pattern.str.extract(r'POWER_BALANCED_([0-9.]+)_')[0].astype(float)
 rows=[]
 for keys,g in balanced.groupby(['gene','working_effect_abs','alpha_name'],sort=False):
  rr=g.loc[g.worst_DGP_rate.idxmin()];rows.append({'gene':keys[0],'focus14':keys[0] in FOCUS14,'working_effect_abs':keys[1],'alpha_name':keys[2],'worst_DGP_and_sign_rate':rr.worst_DGP_rate,'Wilson95_low':rr.Wilson95_low,'qualified':rr.qualified,'not_a_clinical_MCID':True,'sufficient_not_actual_BY_rank_power':keys[2]!='NOMINAL'})
 csvout(RESULTS/'PRIMARY_WORKING_EFFECT_POWER.csv',rows)
 t=pickle.load((PRIVATE/'ANALYSIS_TEMPLATES_RESTRICTED.pkl').open('rb'));data=pickle.load((PRIVATE/'ASSEMBLED_CLINICAL_RESTRICTED.pkl').open('rb'));clinical=data['clinical'];follow=[];profiles=[]
 for acc,gene in FAMILY34.items():
  for branch,marker in [('MEM5Y','ADNI_MEM'),('EF5Y','ADNI_EF')]:
   tt=t[branch+'__'+acc];records=[data['trajectories'][(int(i),marker,5)] for i in tt['ids']];w=np.array([r['w'] for r in records]);sp=np.array([r['span'] for r in records]);last=np.array([r['last_year'] for r in records]);vis=np.array([r['n_visits'] for r in records])
   follow.append({'gene':gene,'accession':acc,'focus14':gene in FOCUS14,'branch':branch,'n_people':len(records),'visits_total':int(vis.sum()),'median_visits':np.median(vis),'median_span_years':np.median(sp),'span_q25':np.quantile(sp,.25),'span_q75':np.quantile(sp,.75),'minimum_span_years':sp.min(),'maximum_span_years':sp.max(),'high_information_people':sum(r['highinfo'] for r in records),'last_visit_ge4p5years':int((last>=4.5).sum()),'last_visit_ge4p5years_fraction':float((last>=4.5).mean()),'time_weight_effective_people_diagnostic':float(w.sum()**2/(w@w)),'max_person_time_weight_fraction':float(w.max()/w.sum()),'top10pct_time_weight_fraction':float(np.sort(w)[-int(np.ceil(.1*len(w))):].sum()/w.sum())})
  ids=t['MEM5Y__'+acc]['ids'];f=clinical.iloc[ids]
  profiles.append({'gene':gene,'focus14':gene in FOCUS14,'n_people':len(f),'age_mean':f.age_month_precision.mean(),'female_fraction':f.sex_code2.eq(1).mean(),'education_mean':f.education.mean(),'APOE4_carrier_fraction':f.apoe4_carrier.eq(1).mean(),'pTau_global_scaled_mean':f.P_LOG.mean()})
 csvout(RESULTS/'FOLLOWUP_INFORMATION_FAMILY34.csv',follow);csvout(RESULTS/'MAIN_CASE_PROFILE_FAMILY34.csv',profiles)
 # Exact omitted-variable algebra for the fixed14 no-pTau sensitivity; no new tests.
 def solve(X,y):return lstsq(X,y,cond=max(X.shape)*np.finfo(float).eps,lapack_driver='gelsy')[0]
 rows=[]
 for acc,gene in FAMILY34.items():
  if gene not in FOCUS14:continue
  main=t['MEM5Y__'+acc];nop=t['MEM_NO_PTAU__'+acc];assert np.array_equal(main['ids'],nop['ids'])
  w=np.array([data['trajectories'][(int(i),'ADNI_MEM',5)]['w'] for i in main['ids']]);P=clinical.iloc[main['ids']].P_LOG.to_numpy()*np.sqrt(w)
  for j,a in enumerate(['N','S','T']):
   X0=np.c_[nop['ZZ'],nop['targetx'][:,j]];X1=np.c_[X0,P];f0=solve(X0,main['y']);f1=solve(X1,main['y']);pf=solve(X0,P);sigma=main['sigma'][j];change=(f0[-1]-f1[-2])*sigma;product=f1[-1]*pf[-1]*sigma
   rows.append({'gene':gene,'assay':a,'standardized_without_pTau':f0[-1]*sigma,'standardized_with_pTau':f1[-2]*sigma,'change_without_minus_with':change,'delta_times_pi':product,'identity_abs_error':abs(change-product),'descriptive_not_causal_or_mediation':True})
 csvout(RESULTS/'FOCUS14_PTAU_CONDITIONING_DECOMPOSITION.csv',rows);assert max(r['identity_abs_error'] for r in rows)<1e-8
 log('Calibration/power summaries, follow-up information and fixed14 conditioning description completed.')
if __name__=='__main__':assert_unsealed();run()

