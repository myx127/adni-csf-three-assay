from config import *
from diagnostics import parameters,summarize
from engine import joint_from_z
from model import bh,holm
from freeze_engine import verify
from scipy.special import ndtr
import pickle

def shared(b,se):return float(joint_from_z((np.asarray(b)/np.asarray(se))[None,:])[2][0])

def run():
 verify();runmeta=json.loads((RESULTS/'SIMULATION_RUN.json').read_text(encoding='utf-8'));assert runmeta['independent_noise_repetitions']==21000000
 assert json.loads((RESULTS/'ENGINE_QA_SUMMARY.json').read_text(encoding='utf-8'))['failed']==0
 decisions=pd.read_csv(RESULTS/'CALIBRATION_DECISIONS.csv').set_index('design_id')
 first=RESULTS/'FIRST_NEW_SHARED_STATISTIC_ACCESS.json'
 if not first.exists():save_json(first,{'date':datetime.now().astimezone().isoformat(),'after_complete_frozen_calibration':True,'prior_coefficients_and_results_already_known':True,'engine_freeze_sha256':sha(RESULTS/'ENGINE_FROZEN.json'),'calibration_decisions_sha256':sha(RESULTS/'CALIBRATION_DECISIONS.csv')})
 templates=pickle.load((PRIVATE/'SIMULATION_TEMPLATES_RESTRICTED.pkl').open('rb'));cache=np.load(S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz')
 meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
 previous=pd.read_csv(S6/'results/COVARIATE_FOUR_MODEL_RESULTS.csv').set_index(['accession','model'])
 eq=pd.read_csv(S7/'results/CONSENSUS_INFORMATION_BOTTLENECKS.csv');eq=eq.loc[eq.margin.eq(.5)].set_index('accession')
 source_main=pd.read_csv(S5/'results/PRIMARY_FROZEN_FAMILY_RESULTS.csv').set_index('accession')
 rows=[];platforms=[];moments={};repro=[]
 for key,t in templates.items():
  d=t['d'];acc=t['accession'];j=list(cache['accessions']).index(acc);sj=list(cache['soma_features']).index(meta.loc[acc,'soma_features'])
  Y=np.c_[cache['N'][:,j],cache['S'][:,sj],cache['T_GIS'][:,j]][d['indices']];assert np.isfinite(Y).all()
  stats,theta,IF,V,E=summarize(Y,d,0);b=theta[:3];se=np.sqrt(V.diagonal()[:3]);z=b/se;pp,pm,p,sgn=joint_from_z(z[None,:]);pp,pm,p,sgn=float(pp[0]),float(pm[0]),float(p[0]),int(sgn[0])
  valid=bool(np.isfinite(b).all() and np.isfinite(se).all() and (se>0).all());qualified=bool(decisions.loc[key,'qualified'] and valid)
  selected=ndtr(-sgn*z) if sgn else np.ones(3);weak=['N','S','T'][int(np.argmax(selected))];old=previous.loc[(acc,t['model'])]
  item={'design_id':key,'accession':acc,'gene':t['gene'],'model':t['model'],'in_family34':acc in FAMILY,'n':d['n'],'rank_X':d['rank'],'residual_df':d['df'],
   'N_G':len(d['onehots'][0]),'S_G':len(d['onehots'][1]),'T_G':len(d['onehots'][2]),'calibration_qualified':bool(decisions.loc[key,'qualified']),'math_valid':valid,
   'all_three_observed_signs_agree':bool((b>0).all() or (b<0).all()),'selected_direction':sgn,'weakest_assay':weak,
   'p_plus_computed':pp,'p_minus_computed':pm,'p_shared_computed':p,'p_shared_for_BY':p if qualified else 1.,
   'p_plus_for_directional_BY':pp if qualified else 1.,'p_minus_for_directional_BY':pm if qualified else 1.,
   'nominal_shared_eligible':qualified and p<=.05,'min_directional_scaled_effect_descriptive':float(np.min(sgn*b/np.sqrt(theta[3:6]))),
   'prior_H_p':old.p_value,'prior_H_q_BY_by_model':old.q_BY_exploratory,
   'prior_M0_equivalence_eligible':bool(eq.loc[acc,'equivalence_eligible_Stage6']),'prior_M0_eq_q_at_margin05':eq.loc[acc,'q_BY_Stage6'],
   'prior_M0_factor_structure_compatible':bool(eq.loc[acc,'point_structure_compatible'])}
  for k,letter in enumerate(['N','S','T']):
   item.update({letter+'_b':b[k],letter+'_SE':se[k],letter+'_z':z[k],letter+'_p_plus':float(ndtr(-z[k])),letter+'_p_minus':float(ndtr(z[k]))})
   gradient=np.zeros(9);gradient[k]=1/np.sqrt(theta[3+k]);gradient[3+k]=-.5*b[k]/theta[3+k]**1.5
   scaled=b[k]/np.sqrt(theta[3+k]);scaled_se=np.sqrt(gradient@V@gradient)
   platforms.append({'design_id':key,'accession':acc,'gene':t['gene'],'model':t['model'],'in_family34':acc in FAMILY,'assay':letter,'n':d['n'],'b_per_global_P_SD':b[k],
    'b_SE_SUM3':se[k],'b_CI95_low':b[k]-1.959964*se[k],'b_CI95_high':b[k]+1.959964*se[k],
    'residual_SD':np.sqrt(theta[3+k]),'scaled_b':scaled,'scaled_SE_full_joint_uncertainty':scaled_se,
    'scaled_CI95_low':scaled-1.959964*scaled_se,'scaled_CI95_high':scaled+1.959964*scaled_se})
  err=max(abs(item[letter+'_b']-old[letter+'_b_per_global_P_SD']) for letter in ['N','S','T'])
  errse=max(abs(item[letter+'_SE']-old[letter+'_b_SE']) for letter in ['N','S','T'])
  repro.append({'design_id':key,'coefficient_max_abs_error':err,'SE_max_abs_error':errse,'passed':err<1e-8 and errse<1e-8 and abs(stats['p_value']-old.p_value)<1e-8})
  if t['model']=='M0' and acc in FAMILY:
   assert abs(old.p_value-source_main.loc[acc,'p_value'])<1e-10
  rows.append(item);moments[key]={'theta':theta,'V':V,'IF':IF}
 result=pd.DataFrame(rows)
 for col in ['q_BY_shared34','q_BH_shared34','p_Holm_shared34','q_BY_direction_plus68','q_BY_direction_minus68','q_BY_combined68_shared','q_BY_combined68_H']:result[col]=np.nan
 for model in MODELS:
  ids=result.index[result.model.eq(model)&result.in_family34];assert len(ids)==34
  p=result.loc[ids,'p_shared_for_BY'].to_numpy()
  result.loc[ids,'q_BY_shared34']=bh(p,True);result.loc[ids,'q_BH_shared34']=bh(p);result.loc[ids,'p_Holm_shared34']=holm(p)
  directional=bh(np.r_[result.loc[ids,'p_plus_for_directional_BY'],result.loc[ids,'p_minus_for_directional_BY']],True)
  result.loc[ids,'q_BY_direction_plus68']=directional[:34];result.loc[ids,'q_BY_direction_minus68']=directional[34:]
  combined=bh(np.r_[p,result.loc[ids,'prior_H_p'].to_numpy()],True)
  result.loc[ids,'q_BY_combined68_shared']=combined[:34];result.loc[ids,'q_BY_combined68_H']=combined[34:]
 result['shared_BY34_supported']=result.in_family34&result.q_BY_shared34.le(.05)
 result['prior_H_supported']=result.in_family34&result.prior_H_q_BY_by_model.le(.05)
 result['shared_combined68_supported']=result.in_family34&result.q_BY_combined68_shared.le(.05)
 result['directional68_supported']=result.in_family34&np.where(result.selected_direction.eq(1),result.q_BY_direction_plus68.le(.05),result.q_BY_direction_minus68.le(.05))
 result['evidence_state']=np.select([result.shared_BY34_supported&result.prior_H_supported,result.shared_BY34_supported&~result.prior_H_supported,~result.shared_BY34_supported&result.prior_H_supported],
  ['SHARED_ASSOCIATION_AND_RELATIVE_DEVIATION','SHARED_ASSOCIATION_ONLY_DEVIATION_UNRESOLVED','RELATIVE_DEVIATION_SHARED_ASSOCIATION_UNRESOLVED'],default='UNRESOLVED')
 csvout(RESULTS/'ALL35_THREE_MODEL_SHARED_RESULTS.csv',result);csvout(RESULTS/'PRIMARY_FAMILY34_SHARED_RESULTS.csv',result.loc[result.in_family34&result.model.eq('M0')])
 csvout(RESULTS/'PLATFORM_EFFECTS_ALL_MODELS.csv',platforms);csvout(RESULTS/'PRIOR_PARAMETER_REPRODUCTION.csv',repro);assert all(r['passed'] for r in repro)
 csvout(RESULTS/'MAPT_CONTROL_ONLY.csv',result.loc[~result.in_family34])
 counts=result.loc[result.in_family34].groupby('model').agg(n_family=('accession','size'),n_calibration_qualified=('calibration_qualified','sum'),n_nominal_shared=('nominal_shared_eligible','sum'),
  n_BY34_shared=('shared_BY34_supported','sum'),n_directional68=('directional68_supported','sum'),n_combined68_shared=('shared_combined68_supported','sum'),n_original_H=('prior_H_supported','sum')).reset_index()
 csvout(RESULTS/'MODEL_DISCOVERY_COUNTS.csv',counts)
 csvout(RESULTS/'EVIDENCE_STATE_COUNTS.csv',result.loc[result.in_family34].groupby(['model','evidence_state']).size().rename('n').reset_index())
 comparison=[];by_model={name:result.loc[result.model.eq(name)&result.in_family34].set_index('accession').loc[FAMILY] for name in MODELS}
 for model in ['APOE','BOTH']:
  base=by_model['M0'];new=by_model[model];p0=base.p_shared_for_BY.to_numpy();p1=new.p_shared_for_BY.to_numpy();q0=bh(p0,True);q1=bh(p1,True)
  for index,acc in enumerate(FAMILY):
   r0=base.loc[acc];r1=new.loc[acc];own=p0.copy();own[index]=p1[index];other=p1.copy();other[index]=p0[index]
   q10=bh(own,True)[index];q01=bh(other,True)[index]
   own_contribution=.5*((q10-q0[index])+(q1[index]-q01));other_contribution=.5*((q01-q0[index])+(q1[index]-q10))
   b0=np.array([r0[x+'_b'] for x in ['N','S','T']]);b1=np.array([r1[x+'_b'] for x in ['N','S','T']]);s0=np.array([r0[x+'_SE'] for x in ['N','S','T']]);s1=np.array([r1[x+'_SE'] for x in ['N','S','T']])
   p00=shared(b0,s0);p10=shared(b1,s0);p01=shared(b0,s1);p11=shared(b1,s1)
   comparison.append({'accession':acc,'gene':r0.gene,'comparison':'M0_to_'+model,'n_identical_cases':r0.n,'M0_p':r0.p_shared_computed,'new_p':r1.p_shared_computed,'M0_q':q0[index],'new_q':q1[index],
    'M0_weakest_assay':r0.weakest_assay,'new_weakest_assay':r1.weakest_assay,'q_own_P_contribution':own_contribution,'q_other33_P_contribution':other_contribution,
    'q_identity_error':abs(q1[index]-q0[index]-own_contribution-other_contribution),
    'computed_P_point_contribution':.5*((p10-p00)+(p11-p01)),'computed_P_SE_contribution':.5*((p01-p00)+(p11-p10)),
    'P_identity_error':abs(p11-p00-.5*((p10-p00)+(p11-p01))-.5*((p01-p00)+(p11-p10))),
    'not_causal_mediation':True,'calibration_qualified_both':bool(r0.calibration_qualified and r1.calibration_qualified)})
 csvout(RESULTS/'APOE_SENSITIVITY_COMPARISONS.csv',comparison)
 with (PRIVATE/'REAL_MOMENTS_RESTRICTED.pkl').open('wb') as f:pickle.dump(moments,f,pickle.HIGHEST_PROTOCOL)
 save_json(RESULTS/'REAL_ANALYSIS_SUMMARY.json',{'family':34,'models':MODELS,'MAPT_separate':True,'new_shared_results':counts.to_dict('records'),'prior_19_structure_incompatible_not_excluded':True,'all105_prior_coefficients_SE_reproduced':True,'not_independent_replication':True})
 log('New shared tests, fixed sensitivity, directional68 and combined68 corrections completed after calibration')
 print(counts.to_string(index=False));print(result.loc[result.model.eq('M0')&result.shared_BY34_supported,['gene','p_shared_computed','q_BY_shared34','q_BY_combined68_shared','weakest_assay','evidence_state']].to_string(index=False))
if __name__=='__main__':run()
