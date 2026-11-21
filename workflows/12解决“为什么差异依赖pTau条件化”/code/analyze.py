from config import *
from engine import stats
import pickle
from scipy.stats import norm
def run():
 cal=pd.read_csv(RESULTS/'CALIBRATION_QUALIFICATION_DECISIONS.csv').set_index('test_id');summary=json.loads((RESULTS/'CALIBRATION_RUN_SUMMARY.json').read_text(encoding='utf-8'))
 assert summary['cells']==210
 save_json(RESULTS/'FIRST_NEW_EXTENSION_TARGET_ACCESS.json',{'time':now(),'protocol_sha256':sha(ROOT/'STAGE12_PROTOCOL_v1.0.md'),'inference_boundary_sha256':sha(ROOT/'EXECUTION_AND_INFERENCE_BOUNDARY_v1.1.md'),'calibration_summary_sha256':sha(RESULTS/'CALIBRATION_RUN_SUMMARY.json'),'prior_stage11_results_known':True,'calibration_completed_before_new_extension_target_results':True,'selection_after_stage11_not_jointly_calibrated':True})
 T=pickle.load((PRIVATE/'TEST_TEMPLATES_RESTRICTED.pkl').open('rb'));audit=pd.read_csv(RESULTS/'ACTUAL_DESIGN_AUDIT.csv');rows=[];effects=[];contrasts=[];influence=[]
 for rr in audit.to_dict('records'):
  key=rr['test_id'];t=T.get(key);row=rr.copy()
  if t is None:row.update({'p_computed':np.nan,'p_for_family':1.,'calibration_qualified':False,'math_valid':False});rows.append(row);continue
  v=stats(t['y'],t);qual=bool(cal.loc[key,'qualified']);valid=bool(v['valid'][0]);row.update({'p_computed':v['p'][0],'p_for_family':float(v['p'][0]) if qual and valid else 1.,'calibration_qualified':qual,'math_valid':valid,'wald_Q':v['Q'][0]})
  for j,label in enumerate(['S_MINUS_N','T_MINUS_N']):
   h=v['h'][0,j];se=np.sqrt(v['V'][0,j,j]);row[label]=h;row[label+'_SE']=se
   contrasts.append({'test_id':key,'gene':t['gene'],'model':t['model'],'kind':t['kind'],'contrast':label,'estimate':h,'SE_full_joint':se,'point_CI95_low':h-norm.ppf(.975)*se,'point_CI95_high':h+norm.ppf(.975)*se,'calibration_qualified':qual,'not_a_simultaneous_interval':True})
  if t['kind']=='MODEL_DIFF':
   for j,a in enumerate(['N','S','T']):
    gamma=v['gamma'][0,j];se=np.sqrt(v['scaledV'][0,j,j])
    effects.append({'test_id':key,'gene':t['gene'],'model':t['model'],'assay':a,'n_people':t['n_people'],'beta':v['b'][0,j],'fixed_protein_sigma':t['sigma'][j],'gamma':gamma,'SE_full_joint':se,'point_CI95_low':gamma-norm.ppf(.975)*se,'point_CI95_high':gamma+norm.ppf(.975)*se,'calibration_qualified':qual,'units':'ADNI_MEM/year/original_fixed_protein_residual_SD'})
  hscore=np.einsum('ap,gp->ga',t['L'],v['scaled_scores'][0],optimize=True);start=0
  for ax,G in enumerate(t['full_G']):
   u=hscore[start:start+G];energy=(u*u/np.maximum(np.diag(v['V'][0]),1e-100)).sum(axis=1)
   influence.append({'test_id':key,'gene':t['gene'],'model':t['model'],'technical_axis':['NULISA','SOMA','TMT'][ax],'G_scale_reference':G,'max_cluster_energy_fraction':float(energy.max()/energy.sum()) if energy.sum()>0 else np.nan,'effective_cluster_energy_count_diagnostic':float(energy.sum()**2/(energy@energy)) if energy@energy>0 else np.nan});start+=G
  rows.append(row)
 out=pd.DataFrame(rows);out['q_BY_family']=np.nan
 for fam,d in out.groupby('family'):out.loc[d.index,'q_BY_family']=by(d.p_for_family)
 out['q_BY_combined42']=by(out.p_for_family);out['supported_family']=out.q_BY_family.le(.05)&out.calibration_qualified;out['supported_combined42']=out.q_BY_combined42.le(.05)&out.calibration_qualified
 csvout(RESULTS/'ALL42_DIRECTIONAL_EXPLANATION_RESULTS.csv',out);csvout(RESULTS/'PRIMARY9_CONDITIONING_CHANGE_RESULTS.csv',out.loc[out.family.eq('PRIMARY9')]);csvout(RESULTS/'SECONDARY33_RESULTS.csv',out.loc[out.family.eq('SECONDARY33')]);csvout(RESULTS/'MODEL27_CLINICAL_DIFFERENCE_RESULTS.csv',out.loc[out.kind.eq('MODEL_DIFF')]);csvout(RESULTS/'ALL15_PAIRED_CONDITIONING_CHANGE_RESULTS.csv',out.loc[out.kind.eq('CONDITIONING_CHANGE')])
 csvout(RESULTS/'PLATFORM_EFFECTS_ALL27.csv',effects);csvout(RESULTS/'PAIRED_CONTRAST_ESTIMATES_ALL42.csv',contrasts);csvout(RESULTS/'TECHNICAL_CLUSTER_INFLUENCE.csv',influence)
 save_json(RESULTS/'REAL_ANALYSIS_SUMMARY.json',{'time':now(),'supported_primary9':out.loc[out.family.eq('PRIMARY9')&out.supported_family,'test_id'].tolist(),'supported_secondary33':out.loc[out.family.eq('SECONDARY33')&out.supported_family,'test_id'].tolist(),'supported_combined42':out.loc[out.supported_combined42,'test_id'].tolist(),'tests':len(out),'n_primary':9,'n_secondary':33,'pTau_causal_effect_identified':False,'true_protein_reliability_identified':False,'molecular_form_identified':False})
 log('42 fixed target results calculated after calibration; same cases and paired cross-model covariance preserved.')
 print(out[['gene','model','kind','p_computed','q_BY_family','q_BY_combined42','supported_family','calibration_qualified']].to_string(index=False))
if __name__=='__main__':assert_unsealed();run()

