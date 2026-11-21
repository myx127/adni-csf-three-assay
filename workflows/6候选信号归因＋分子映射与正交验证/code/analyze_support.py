from config import *
from diagnostics import summarize,parameters
from model import bh
from designs import make
import pickle

def correlations(Y,d,formal=True):
 Q=d['Q'];E=Y-Q@(Q.T@Y);c=E.T@E/d['df'];R=c/np.sqrt(np.outer(np.diag(c),np.diag(c)))
 return {'corr_residual_N_S':R[0,1],'corr_residual_N_targetMS':R[0,2],'corr_residual_S_targetMS':R[1,2]}

def run():
 plan=json.loads((RESULTS/'TARGETED_VERIFICATION_PLAN_FROZEN.json').read_text(encoding='utf-8'))
 first=json.loads((RESULTS/'FIRST_TARGETED_QUANTITATIVE_ACCESS.json').read_text(encoding='utf-8'));assert first['plan_sha256_before_association']==sha(RESULTS/'TARGETED_VERIFICATION_PLAN_FROZEN.json')
 nsc=np.load(RESTRICTED/'SUPPORT_N_S_CACHE_RESTRICTED.npz');t=np.load(RESTRICTED/'TARGETED_LOG_RATIOS_RESTRICTED.npz')
 frames=pickle.load((RESTRICTED/'ORTHOGONAL_SUPPORT_CASE_FRAMES_RESTRICTED.pkl').open('rb'));ds=pickle.load((RESTRICTED/'ORTHOGONAL_SUPPORT_DESIGNS_RESTRICTED.pkl').open('rb'))
 screen=pd.read_csv(RESULTS/'ORTHOGONAL_SUPPORT_DESIGN_SCREEN.csv');cal=pd.read_csv(RESULTS/'ORTHOGONAL_SUPPORT_CALIBRATION_DECISIONS.csv').set_index('design_id')
 assert np.array_equal(nsc['keys'],t['keys']);rows=[];releases=[]
 for r in screen.itertuples():
  frame=frames[r.design_id]['frame'];mask=frames[r.design_id]['mask'];j=list(nsc['accessions']).index(r.accession);k=list(t['feature_keys']).index(r.feature_key)
  Yall=np.c_[nsc['N'][:,j],nsc['S'][:,j],t['Y'][:,k]];assert np.isfinite(Yall[mask]).all()
  base={'design_id':r.design_id,'cohort':r.cohort,'accession':r.accession,'gene':r.gene,'peptide':r.peptide,'feature_key':r.feature_key,'n':r.n,
    'three_discovery_direct_validation':False,'nuisance_conditioning_changed_from_S5':True,'same_parent_ADNI_not_external_cohort':True}
  if r.design_id in ds:
   d=ds[r.design_id];passed=bool(cal.loc[r.design_id,'passed']);stats,_,_,_,_=summarize(Yall[d['indices']],d,0)
   stats={('targetMS_'+k[2:] if k.startswith('T_') else k.replace('corr_residual_NT','corr_residual_N_targetMS').replace('corr_residual_ST','corr_residual_S_targetMS').replace('c_NT','c_NtargetMS').replace('c_ST','c_StargetMS')):v for k,v in stats.items()}
   base.update({**stats,**correlations(Yall[d['indices']],d),'calibration_passed':passed,'analysis_status':'EXPLORATORY_CALIBRATED_MEASUREMENT_SUPPORT' if passed else 'CALIBRATION_FAILED_NO_INFERENTIAL_CLAIM'})
   if not passed:base['p_value_for_family']=1.
   else:base['p_value_for_family']=stats['p_value']
  else:
   d=make(frame,mask,tech=('NULISA_PlateID','SOMASCAN_PlateId','TARGETED_Plate'))
   base.update({'calibration_passed':False,'analysis_status':'DESIGN_INSUFFICIENT_DESCRIPTIVE_ONLY','p_value':np.nan,'p_value_for_family':1.})
   if d is not None and d['df']>0:base.update(correlations(Yall[d['indices']],d))
  base['N_detection_fraction']=float(nsc['N_detect'][mask,j].mean()) if mask.any() else np.nan
  rows.append(base)
  for source in ['V01','V02']:
   sub=mask&frame.TARGETED_Plate.str.startswith(source+'_').to_numpy();dd=make(frame,sub,tech=('NULISA_PlateID','SOMASCAN_PlateId','TARGETED_Plate'))
   if dd is not None and dd['df']>10:
    yy=Yall[dd['indices']];coef=dd['a']@yy;releases.append({'cohort':r.cohort,'accession':r.accession,'gene':r.gene,'peptide':r.peptide,'source':source,'n':dd['n'],'residual_df':dd['df'],
      'N_P_coefficient_descriptive':coef[0],'S_P_coefficient_descriptive':coef[1],'targetMS_P_coefficient_descriptive':coef[2],**correlations(yy,dd),
      'not_separately_calibrated_no_discovery_classification':True})
 out=pd.DataFrame(rows);out['q_BY5_exploratory']=np.nan
 for cohort,part in out.groupby('cohort'):
  assert len(part)==5;out.loc[part.index,'q_BY5_exploratory']=bh(part.p_value_for_family.to_numpy(),True)
 out['relative_constraint_deviation_BY']=out.calibration_passed&out.q_BY5_exploratory.le(.05)
 csvout(RESULTS/'ORTHOGONAL_MEASUREMENT_SUPPORT_RESULTS.csv',out);csvout(RESULTS/'TARGETED_PER_RELEASE_DESCRIPTIVE_RESULTS.csv',releases)
 summary=out.groupby('cohort').agg(n_peptide_hypotheses=('feature_key','size'),n_calibrated=('calibration_passed','sum'),n_BY_relative_deviations=('relative_constraint_deviation_BY','sum'),min_n=('n','min'),max_n=('n','max')).reset_index()
 csvout(RESULTS/'ORTHOGONAL_SUPPORT_COHORT_SUMMARY.csv',summary)
 save_json(RESULTS/'EMORY_PANEL_DIRECT_VALIDATION_STATUS.json',{'scope':'Emory 2022/2023 releases V01 and V02 only; not later legacy extension','NPTX1':'NO_DIRECT_TARGETED_ASSAY_COVERAGE','PDGFRB':'NO_DIRECT_TARGETED_ASSAY_COVERAGE','KLK6':'NO_DIRECT_TARGETED_ASSAY_COVERAGE',
  'direct_three_discovery_validation_completed':False,'reason':'None of three exact accessions in either registered targeted-MS source',
  'other_measured_proteins':'CHI3L1/MDH1/NPTX2/SOD1 measurement support only','new_RID_strict_relative_constraint_inference':'Not performed because residual df below40',
  'external_cohort_validation_completed':False})
 log('Available orthogonal measurement support analyzed; three primary discoveries explicitly remain unvalidated by these new panels')
 print(out[['cohort','gene','peptide','n','analysis_status','p_value','q_BY5_exploratory','corr_residual_N_targetMS','corr_residual_S_targetMS']].to_string(index=False))

if __name__=='__main__':run()
