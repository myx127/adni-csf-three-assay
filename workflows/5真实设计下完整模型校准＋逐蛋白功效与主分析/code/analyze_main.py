from config import *
from model import make_design,batch_components,bh,holm
from diagnostics import summarize,covariance_profile
import pickle

def actual_Y(cache,meta,accession,variant='MAIN'):
 j=list(cache['accessions']).index(accession);row=meta.loc[accession];features=str(row.soma_features).split('|');assert len(features)==1
 k=list(cache['soma_features']).index(features[0]);return np.c_[cache['N'][:,j],cache['S'][:,k],cache['T_RAW' if variant=='RAW_TMT' else 'T_GIS'][:,j]]

def adjust(frame,family):
 frame=frame.copy();frame['in_frozen_family']=frame.accession.isin(family);frame['q_BY']=np.nan;frame['q_BH']=np.nan;frame['p_Holm']=np.nan
 wanted=frame.in_frozen_family
 assert int(wanted.sum())==len(family)
 p=frame.loc[wanted,'p_value'].to_numpy();p=np.where(np.isfinite(p),p,1.)
 frame.loc[wanted,'q_BY']=bh(p,True);frame.loc[wanted,'q_BH']=bh(p);frame.loc[wanted,'p_Holm']=holm(p)
 frame['BY_discovery']=wanted&frame.q_BY.le(.05);return frame

def run():
 sap=json.loads((RESULTS/'MAIN_ANALYSIS_SAP_FROZEN.json').read_text(encoding='utf-8'));secondary=json.loads((RESULTS/'SECONDARY_SAP_FROZEN.json').read_text(encoding='utf-8'))
 assert sha(ROOT/'STAGE5_PROTOCOL_v1.0.md')==sap['protocol_sha256']
 for name in ['model.py','diagnostics.py','load_proteins.py','analyze_main.py']:
  assert sha(ROOT/'code'/name)==sap['code_sha256'][name],f'Pre-analysis code changed: {name}'
 clinical=pd.read_csv(RESTRICTED/'CLINICAL_DESIGN_RESTRICTED.csv');meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
 cache=np.load(RESTRICTED/'PROTEIN_TRANSFORMED_RESTRICTED.npz');assert np.array_equal(cache['keys'],clinical.key.to_numpy())
 designs=pickle.load((RESTRICTED/'DESIGNS_MAIN_RESTRICTED.pkl').open('rb'));screen=pd.read_csv(RESULTS/'PRE_EFFECT_PROTEIN_SCREEN.csv')
 rows=[];profile_rows=[];residuals={};details={}
 for m in screen.loc[screen.single_annotation].itertuples():
  base={'accession':m.accession,'gene':m.gene,'role':m.role,'n':m.n,'pre_effect_status':m.design_status,'pre_effect_exclusion_reasons':m.exclusion_reasons,
    'calibration_passed':m.accession in sap['calibration_passed_including_control'],'variance_version':sap['version']}
  if m.accession not in designs:rows.append({**base,'analysis_status':'EXCLUDED_PRE_EFFECT','p_value':np.nan});continue
  d=designs[m.accession];Y=actual_Y(cache,meta,m.accession)[d['indices']];assert np.isfinite(Y).all()
  stats,theta,IF,V,E=summarize(Y,d,sap['version']);assert np.allclose(stats['Q'],batch_components(Y,d)['Q0'][0]/[1.,(d['n']-1)/d['df'],2*(d['n']-1)/d['df']][sap['version']],rtol=1e-8)
  base.update(d['summary']);rows.append({**base,**stats,'analysis_status':'CALIBRATED_MAIN' if base['calibration_passed'] else 'EXPLORATORY_CALIBRATION_FAILED'})
  residuals[m.accession]=E;details[m.accession]=(theta,V)
  for p in covariance_profile(theta,V):profile_rows.append({'accession':m.accession,'gene':m.gene,**p})
 main=adjust(pd.DataFrame(rows),sap['primary_family_accessions']);csvout(RESULTS/'PRIMARY_ALL43_RESULTS.csv',main)
 csvout(RESULTS/'PRIMARY_FROZEN_FAMILY_RESULTS.csv',main.loc[main.in_frozen_family].sort_values('p_value'))
 profile=pd.DataFrame(profile_rows);csvout(RESULTS/'COVARIANCE_OFFSET_GRID125.csv',profile)
 ps=[]
 for acc,p in profile.groupby('accession'):
  q=p.p;ps.append({'accession':acc,'p_grid_min':q.min(),'p_grid_max':q.max(),'all125_valid':bool(p.valid.all()),
   'crosses_nominal05':bool((q<=.05).any() and (q>.05).any()),'not_a_continuous_robust_bound':True,'not_empirically_calibrated_error_range':True})
 csvout(RESULTS/'COVARIANCE_OFFSET_SENSITIVITY_SUMMARY.csv',ps)
 variants=pd.read_csv(RESULTS/'VARIANT_CALIBRATION_DECISIONS.csv').fillna('');extra=pickle.load((RESTRICTED/'DESIGNS_VARIANTS_RESTRICTED.pkl').open('rb'))
 variant_rows=[]
 for r in variants.to_dict('records'):
  ident=r['calibration_design_id'];base={k:r[k] for k in ['variant','accession','gene','role','n','calibration_passed','status','reasons']}
  if not ident:variant_rows.append({**base,'analysis_status':'EXCLUDED_VARIANT_DESIGN','p_value':np.nan});continue
  d=designs[ident[6:]] if ident.startswith('MAIN__') else extra[ident]
  Y=actual_Y(cache,meta,r['accession'],r['variant'])[d['indices']];assert np.isfinite(Y).all()
  stats,_,_,_,_=summarize(Y,d,sap['version']);variant_rows.append({**base,**stats,'analysis_status':'CALIBRATED_SUPPORTIVE' if r['calibration_passed'] else 'EXPLORATORY_UNCALIBRATED'})
 variants=pd.DataFrame(variant_rows);output=[]
 for variant,part in variants.groupby('variant'):
  fam=secondary['amyloid_primary_family'] if variant=='AMY' else sap['primary_family_accessions']
  # Invalid sensitivity designs remain in the fixed family with p=1 for adjustment.
  part=adjust(part,fam);output.append(part)
 variants=pd.concat(output,ignore_index=True);csvout(RESULTS/'ALL_FIXED_VARIANT_RESULTS.csv',variants)
 csvout(RESULTS/'AMYLOID_SECONDARY_RESULTS.csv',variants.loc[variants.variant.eq('AMY')].sort_values('p_value'))
 # Clinical values and per-person residuals remain private; no latent F/D scores.
 np.savez_compressed(RESTRICTED/'MAIN_RESIDUALS_RESTRICTED.npz',**residuals)
 joined=main[['accession','gene','p_value','q_BY','BY_discovery']].merge(pd.DataFrame(ps),on='accession',how='left')
 csvout(RESULTS/'MAIN_INTERPRETATION_CHECKS.csv',joined)
 save_json(RESULTS/'MAIN_RUN_SUMMARY.json',{'primary_family_m':sap['m_primary'],'version':sap['version'],
  'primary_BY_discoveries':int(main.BY_discovery.sum()),'primary_BH_discoveries':int(main.loc[main.in_frozen_family,'q_BH'].le(.05).sum()),
  'primary_Holm_discoveries':int(main.loc[main.in_frozen_family,'p_Holm'].le(.05).sum()),
  'amyloid_family_m':secondary['m_amyloid'],'amyloid_BY_discoveries':int(variants.loc[variants.variant.eq('AMY'),'BY_discovery'].sum()),
  'individual_F_D_scores_created':False,'cognitive_slopes_analyzed':False,'molecular_mechanism_identified':False})
 log(f'Actual main analysis finished: BY discoveries {int(main.BY_discovery.sum())}/{sap["m_primary"]}; full excluded results retained')
 print(main.loc[main.in_frozen_family].sort_values('p_value')[['gene','n','Q','p_value','q_BY','positive_cov_CI_all_three']].head(15).to_string(index=False))

if __name__=='__main__':run()
