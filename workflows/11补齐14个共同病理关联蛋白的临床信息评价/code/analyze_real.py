from config import *
from engine import joint_stats,D
from scipy.stats import norm
import pickle

def run():
 calibration=pd.read_csv(RESULTS/'CALIBRATION_QUALIFICATION_DECISIONS.csv').set_index('design_id')
 summary=json.loads((RESULTS/'CALIBRATION_RUN_SUMMARY.json').read_text(encoding='utf-8'))
 assert summary['cells']==720
 first=RESULTS/'FIRST_REAL_TARGET_ASSOCIATION_ACCESS.json'
 if not first.exists():save_json(first,{'time':now(),'protocol_sha256':sha(ROOT/'STAGE11_PROTOCOL_v1.0.md'),'calibration_decisions_sha256':sha(RESULTS/'CALIBRATION_QUALIFICATION_DECISIONS.csv'),'calibration_completed_before_target_coefficients':True,'prior_pathology_results_already_known':True})
 templates=pickle.load((PRIVATE/'ANALYSIS_TEMPLATES_RESTRICTED.pkl').open('rb'))
 design_results=[];effects=[];contrasts=[];influence=[];allstats={}
 crit=norm.ppf(.975);crit3=norm.ppf(1-.05/6)
 for key,t in templates.items():
  stats=joint_stats(t['y'],t);allstats[key]=stats
  b=stats['b'][0];gamma=stats['scaled'][0];V=stats['Vraw'][0];Vs=stats['Vscaled'][0]
  valid=bool(stats['valid'][0]);qualified=bool(calibration.loc[key,'qualified']);same_sign=bool((b>0).all() or (b<0).all());direction='POSITIVE' if (b>0).all() else 'NEGATIVE' if (b<0).all() else 'MIXED'
  row={'design_id':key,'gene':t['gene'],'accession':t['accession'],'branch':t['branch'],'formal':t['formal'],'n_people':t['n_people'],'information_df':t['summary']['information_df'],'calibration_qualified':qualified,'math_valid':valid,'all_three_observed_signs_agree':same_sign,'direction':direction,'p_same_computed':stats['p_same'][0],'p_diff_computed':stats['p_diff'][0],'p_any_computed':stats['p_any'][0],'Qdiff':stats['Qdiff'][0],'Qany':stats['Qany'][0],'p_plus':stats['p_plus'][0],'p_minus':stats['p_minus'][0]}
  sign=1 if direction=='POSITIVE' else -1 if direction=='NEGATIVE' else 1
  row['weakest_assay']=['N','S','T'][int(np.argmin(sign*stats['z'][0]))]
  for j,assay in enumerate(['N','S','T']):
   se=np.sqrt(V[j,j]);ses=np.sqrt(Vs[j,j])
   effects.append({'design_id':key,'branch':t['branch'],'accession':t['accession'],'gene':t['gene'],'assay':assay,'n_people':t['n_people'],'calibration_qualified':qualified,'b_original_scale':b[j],'b_SE':se,'b_CI95_low':b[j]-crit*se,'b_CI95_high':b[j]+crit*se,'protein_fixed_residual_SD':t['sigma'][j],'standardized_effect':gamma[j],'standardized_SE_full_joint':ses,'standardized_CI95_low':gamma[j]-crit*ses,'standardized_CI95_high':gamma[j]+crit*ses,'units':'CL per protein residualSD' if t['branch']=='AMY_CL' else 'MEM score per protein residualSD' if t['branch']=='MEM_BASE' else ('cognitive score/year per protein residualSD'+(' per50CL modifier' if t['branch']=='MEM_AMY_INTERACTION' else ' APOE-carrier difference' if t['branch']=='MEM_APOE_INTERACTION' else ''))})
   row[assay+'_standardized_effect']=gamma[j];row[assay+'_SE_full_joint']=ses
  for label,c in [('S_MINUS_N',np.array([-1.,1.,0.])),('T_MINUS_N',np.array([-1.,0.,1.])),('T_MINUS_S',np.array([0.,-1.,1.]))]:
   delta=float(c@gamma);se=float(np.sqrt(c@Vs@c))
   contrasts.append({'design_id':key,'branch':t['branch'],'gene':t['gene'],'contrast':label,'difference_standardized_effect':delta,'SE_full_joint':se,'CI95_pointwise_low':delta-crit*se,'CI95_pointwise_high':delta+crit*se,'CI95_within_design_3pair_Bonf_low':delta-crit3*se,'CI95_within_design_3pair_Bonf_high':delta+crit3*se,'pairwise_P_not_used_for_discovery':True})
  rawscore=np.einsum('gjn,n->gj',t['K'],t['y'],optimize=True)
  scaledscore=rawscore*t['sigma'][None,:]+b[None,:]*t['Uv']/(2*t['sigma'][None,:])
  start=0
  for axis,label in enumerate(['NULISA','SOMA','TMT']):
   G=len(np.unique(t['groups'][axis]));ss=scaledscore[start:start+G]
   energy=np.sum((ss@D.T)**2/np.maximum(np.diag(stats['Vdiff'][0]),1e-100),axis=1)
   influence.append({'design_id':key,'branch':t['branch'],'gene':t['gene'],'technical_axis':label,'G_full_scale_cohort':G,'effective_joint_difference_clusters':float(energy.sum()**2/(energy@energy)) if energy@energy>0 else np.nan,'max_cluster_joint_difference_energy_fraction':float(energy.max()/energy.sum()) if energy.sum()>0 else np.nan})
   start+=G
  design_results.append(row)
 out=pd.DataFrame(design_results);registry=pd.read_csv(RESULTS/'FIXED_FORMAL_HYPOTHESIS_REGISTRY.csv');formal=[]
 for r in registry.to_dict('records'):
  row=out.loc[out.design_id.eq(r['design_id'])].iloc[0];p=float(row['p_'+r['test'].lower()+'_computed']);pass_=bool(row.calibration_qualified and row.math_valid)
  formal.append({**r,'n_people':row.n_people,'information_df':row.information_df,'calibration_qualified':bool(row.calibration_qualified),'math_valid':bool(row.math_valid),'p_computed':p,'p_for_family':p if pass_ else 1.,'direction':row.direction,'weakest_assay':row.weakest_assay})
 formal=pd.DataFrame(formal);formal['q_BY_within_family']=np.nan
 for family,part in formal.groupby('family'):formal.loc[part.index,'q_BY_within_family']=by(part.p_for_family)
 formal['q_BY_combined204']=by(formal.p_for_family);formal['supported_within_family']=formal.q_BY_within_family.le(.05)&formal.calibration_qualified;formal['supported_combined204']=formal.q_BY_combined204.le(.05)&formal.calibration_qualified
 csvout(RESULTS/'FORMAL204_CLINICAL_RESULTS.csv',formal);csvout(RESULTS/'PRIMARY688_MEMORY_PROGRESSION_RESULTS.csv',formal.loc[formal.family.eq('PRIMARY68')]);csvout(RESULTS/'SECONDARY136_RESULTS.csv',formal.loc[formal.family.eq('SECONDARY136')])
 sens=[]
 for branch in SENSITIVITIES:
  rows=[]
  for r in out.loc[out.branch.eq(branch)].to_dict('records'):
   for test in ['SAME','DIFF']:
    rows.append({'design_id':r['design_id'],'gene':r['gene'],'branch':branch,'test':test,'n_people':r['n_people'],'calibration_qualified':r['calibration_qualified'],'p_computed':r['p_'+test.lower()+'_computed'],'p_for_diagnostic_family':r['p_'+test.lower()+'_computed'] if r['calibration_qualified'] and r['math_valid'] else 1.,'direction':r['direction']})
  part=pd.DataFrame(rows);part['q_BY28_diagnostic']=by(part.p_for_diagnostic_family);part['not_new_confirmatory_family']=True;sens.extend(part.to_dict('records'))
 csvout(RESULTS/'PRESPECIFIED_SENSITIVITY_RESULTS.csv',sens)
 csvout(RESULTS/'ALL144_DESIGN_CLINICAL_RESULTS.csv',out);csvout(RESULTS/'PLATFORM_CLINICAL_EFFECTS_ALL144.csv',effects);csvout(RESULTS/'PAIRWISE_CLINICAL_CONTRASTS_ALL144.csv',contrasts);csvout(RESULTS/'REAL_CLUSTER_INFLUENCE_DIAGNOSTICS.csv',influence)
 frozen=pd.read_csv(RESULTS/'FROZEN_FAMILY34_AND_FOCUS14.csv');integrated=[]
 for prior in frozen.to_dict('records'):
  gene=prior['gene'];row=prior.copy()
  for branch in FORMAL_BRANCHES:
   for test in ['SAME','DIFF']:
    rr=formal.loc[formal.gene.eq(gene)&formal.branch.eq(branch)&formal.test.eq(test)].iloc[0]
    for field in ['p_computed','q_BY_within_family','q_BY_combined204','supported_within_family','supported_combined204','calibration_qualified','direction','n_people']:
     row[branch+'_'+test+'_'+field]=rr[field]
  integrated.append(row)
 integ=pd.DataFrame(integrated);csvout(RESULTS/'INTEGRATED_PATHOLOGY_CLINICAL_FAMILY34.csv',integ)
 csvout(RESULTS/'FOCUS14_INTEGRATED_EVIDENCE.csv',integ.loc[integ.focus14])
 formal['focus14']=formal.gene.isin(FOCUS14)
 csvout(RESULTS/'FORMAL204_CLINICAL_RESULTS.csv',formal)
 csvout(RESULTS/'PRIMARY68_MEMORY_PROGRESSION_RESULTS.csv',formal.loc[formal.family.eq('PRIMARY68')])
 csvout(RESULTS/'SECONDARY136_RESULTS.csv',formal.loc[formal.family.eq('SECONDARY136')])
 csvout(RESULTS/'FOCUS14_PRIMARY28_FROM_FULL_FAMILY.csv',formal.loc[formal.focus14&formal.family.eq('PRIMARY68')])
 csvout(RESULTS/'FOCUS14_SECONDARY56_FROM_FULL_FAMILY.csv',formal.loc[formal.focus14&formal.family.eq('SECONDARY136')])
 ec=pd.DataFrame(effects);cc=pd.DataFrame(contrasts)
 csvout(RESULTS/'FOCUS14_PLATFORM_EFFECTS.csv',ec.loc[ec.gene.isin(FOCUS14)])
 csvout(RESULTS/'FOCUS14_PAIRED_CONTRASTS.csv',cc.loc[cc.gene.isin(FOCUS14)])
 states=[]
 for branch in FORMAL_BRANCHES:
  for focus,genes in [('FOCUS14',FOCUS14),('OTHER20',[g for g in FAMILY34.values() if g not in FOCUS14])]:
   pp=formal.loc[formal.branch.eq(branch)&formal.gene.isin(genes)]
   for test in ['SAME','DIFF']:
    rr=pp.loc[pp.test.eq(test)]
    states.append({'branch':branch,'subset':focus,'test':test,'n_hypotheses':len(rr),'n_qualified':int(rr.calibration_qualified.sum()),'nominal_P_lt05':int((rr.p_computed.le(.05)&rr.calibration_qualified).sum()),'supported_full_family_BY':int(rr.supported_within_family.sum()),'supported_combined204_BY':int(rr.supported_combined204.sum()),'supported_genes_full_family':';'.join(rr.loc[rr.supported_within_family,'gene']),'supported_genes_combined204':';'.join(rr.loc[rr.supported_combined204,'gene'])})
 csvout(RESULTS/'EVIDENCE_STATE_COUNTS.csv',states)
 old=pd.read_csv(S10/'results/ALL43_DESIGN_CLINICAL_RESULTS.csv').set_index('design_id');repro=[]
 for key,t in templates.items():
  if key not in old.index:continue
  r=out.loc[out.design_id.eq(key)].iloc[0];oo=old.loc[key]
  ds=[abs(float(r['p_'+test+'_computed'])-float(oo['p_'+test+'_computed'])) for test in ['same','diff']]
  errs=[abs(float(r[a+'_standardized_effect'])-float(oo[a+'_standardized_effect'])) for a in ['N','S','T']]
  repro.append({'design_id':key,'gene':t['gene'],'branch':t['branch'],'n_same':int(r.n_people)==int(oo.n_people),'max_abs_P_error':max(ds),'max_abs_effect_error':max(errs),'passed':int(r.n_people)==int(oo.n_people) and max(ds)<1e-8 and max(errs)<1e-8,'q_values_not_expected_equal_new_family_sizes':True})
 csvout(RESULTS/'STAGE10_OVERLAPPING_DESIGN_REPRODUCTION.csv',repro)
 assert all(r['passed'] for r in repro)
 save_json(RESULTS/'REAL_ANALYSIS_SUMMARY.json',{'time':now(),'primary_supported_hypotheses':formal.loc[formal.family.eq('PRIMARY68')&formal.supported_within_family,'hypothesis_id'].tolist(),'secondary_supported_hypotheses':formal.loc[formal.family.eq('SECONDARY136')&formal.supported_within_family,'hypothesis_id'].tolist(),'combined204_supported_hypotheses':formal.loc[formal.supported_combined204,'hypothesis_id'].tolist(),'focus14_primary_supported':formal.loc[formal.focus14&formal.family.eq('PRIMARY68')&formal.supported_within_family,'hypothesis_id'].tolist(),'focus14_secondary_supported':formal.loc[formal.focus14&formal.family.eq('SECONDARY136')&formal.supported_within_family,'hypothesis_id'].tolist(),'fixed_family_genes':34,'fixed_focus_genes':14,'formal_hypotheses':204,'primary_hypotheses':68,'secondary_hypotheses':136,'sensitivity_diagnostic_hypotheses':len(sens),'prior_stage10_designs_reproduced':len(repro),'independent_validation':False,'causal_inference':False,'prediction_model_developed':False})
 log('Clinical analysis complete: fixed34/204 tests; focus14 retains full-family q; Stage10 overlapping designs reproduced.')
 print(pd.DataFrame(states).to_string(index=False))
 print('PRIMARY BY RESULTS');print(formal.loc[formal.family.eq('PRIMARY68')&formal.supported_within_family,['gene','test','direction','p_computed','q_BY_within_family','q_BY_combined204']].to_string(index=False))
if __name__=='__main__':assert_unsealed();run()
