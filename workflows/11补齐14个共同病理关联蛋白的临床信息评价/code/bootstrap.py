from config import *
import shutil,pickle

def run():
 assert_unsealed()
 original=pd.read_csv(S8/'results/PRIMARY_FAMILY34_SHARED_RESULTS.csv')
 main=original.loc[original.model.eq('M0')].copy()
 assert dict(zip(main.accession,main.gene))==FAMILY34
 assert set(main.loc[main.shared_BY34_supported,'gene'])==set(FOCUS14)
 rows=[]
 for acc,gene in FAMILY34.items():
  part=original.loc[original.accession.eq(acc)]
  r={'accession':acc,'gene':gene,'focus14':gene in FOCUS14}
  for row in part.itertuples():
   r[row.model+'_pathology_same_q_BY34']=row.q_BY_shared34
   r[row.model+'_relative_H_q_BY34']=row.prior_H_q_BY_by_model
  rows.append(r)
 csvout(RESULTS/'FROZEN_FAMILY34_AND_FOCUS14.csv',rows)
 srcs=[S10/'00_restricted/ASSEMBLED_CLINICAL_RESTRICTED.pkl',S10/'code/prepare_templates.py',S10/'code/engine.py',S10/'code/validate_engine.py',S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz',S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv',S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv',S8/'results/PRIMARY_FAMILY34_SHARED_RESULTS.csv',S10/'results/ALL43_DESIGN_CLINICAL_RESULTS.csv',S10/'results/PLATFORM_CLINICAL_EFFECTS_ALL43.csv']
 rows=[{'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)} for p in srcs]
 csvout(RESULTS/'INPUT_SOURCE_MANIFEST.csv',rows)
 previous=[];cache={}
 prior=pd.read_csv(S10/'results/PRIOR_ALL_NINE_STAGES_SEAL_CHECKS_FINAL.csv')[['path','expected_sha256']].to_dict('records')
 prior += [{'path':str(S10/r.relative_path),'expected_sha256':r.sha256} for r in pd.read_csv(S10/'SHA256_MANIFEST.csv').itertuples()]
 for r in prior:
  p=Path(r['path']).resolve();key=str(p)
  if key not in cache: cache[key]=sha(p) if p.is_file() else ''
  previous.append({'stage':p.relative_to(PROJECT).parts[0] if p.is_relative_to(PROJECT) else 'ORIGINAL_SOURCE','path':key,'expected_sha256':r['expected_sha256'],'actual_sha256':cache[key],'passed':cache[key]==r['expected_sha256']})
 assert all(r['passed'] for r in previous)
 csvout(RESULTS/'PRIOR_TEN_STAGE_SEAL_CHECKS_INITIAL.csv',previous)
 for name in ['prepare_templates.py','engine.py','validate_engine.py']:
  shutil.copyfile(S10/'code'/name,EVIDENCE/('stage10_source_'+name+'.txt'))
 shutil.copyfile(S10/'00_restricted/ASSEMBLED_CLINICAL_RESTRICTED.pkl',PRIVATE/'ASSEMBLED_CLINICAL_RESTRICTED.pkl')
 save_json(RESULTS/'ACCESS_AND_PROTOCOL_RECORD.json',{'time':now(),'protocol_sha256':sha(ROOT/'STAGE11_PROTOCOL_v1.0.md'),'scope':'Fixed34 clinical family; focus14; only CSF','prior_stage10_three_candidate_clinical_results_known':True,'prior_pathology_results_known':True,'clinical_values_previously_accessed':True,'new_noncore_clinical_target_coefficients_fitted':False,'formal_hypotheses':204,'primary':68,'secondary':136,'focus14_not_a_reduced_confirmatory_family':True,'all_prior_stages_read_only':True,'prior_sealed_entries':len(previous),'prior_unique_paths':len(cache),'authorization':'Explicit user instruction to complete Stage11 in this directory'})
 log(f'Protocol fixed; family34/focus14 verified; {len(previous)} prior sealed entries unchanged.')

if __name__=='__main__': run()

