from config import *
import pickle,shutil
def run():
 assert_unsealed();sources=[]
 for k in ['I01','I04','I05','I07']:
  r=INPUTS[k];p=Path(r['path']);actual=sha(p);assert actual==r['sha256'];sources.append({'path':str(p),'bytes':p.stat().st_size,'sha256':actual,'role':k})
 for p in [S11/'00_restricted/ANALYSIS_TEMPLATES_RESTRICTED.pkl',S11/'00_restricted/EXPOSURE_SCALING_RESTRICTED.pkl',S11/'00_restricted/ASSEMBLED_CLINICAL_RESTRICTED.pkl',S11/'results/FROZEN_FAMILY34_AND_FOCUS14.csv',S11/'results/PLATFORM_CLINICAL_EFFECTS_ALL144.csv',S11/'results/ALL144_DESIGN_CLINICAL_RESULTS.csv',S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz',S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv']:
  sources.append({'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p),'role':'SEALED_PRIOR'})
 csvout(RESULTS/'INPUT_SOURCE_MANIFEST.csv',sources)
 rows=pd.read_csv(S11/'results/PRIOR_TEN_STAGE_SEAL_CHECKS_FINAL.csv')[['path','expected_sha256']].to_dict('records')
 rows += [{'path':str(S11/r.relative_path),'expected_sha256':r.sha256} for r in pd.read_csv(S11/'SHA256_MANIFEST.csv').itertuples()]
 cache={};prior=[]
 for r in rows:
  p=Path(r['path'])
  if r['path'] not in cache:cache[r['path']]=sha(p) if p.exists() else ''
  prior.append({**r,'actual_sha256':cache[r['path']],'passed':cache[r['path']]==r['expected_sha256']})
 assert all(r['passed'] for r in prior);csvout(RESULTS/'PRIOR_ELEVEN_STAGE_SEAL_CHECKS_INITIAL.csv',prior)
 old=pickle.load((S11/'00_restricted/ANALYSIS_TEMPLATES_RESTRICTED.pkl').open('rb'))
 selected={k:t for k,t in old.items() if t['accession'] in GENES and t['branch'] in ['MEM5Y','MEM_NO_PTAU']}
 with (PRIVATE/'STAGE11_REFERENCE_TEMPLATES_RESTRICTED.pkl').open('wb') as f:pickle.dump(selected,f,pickle.HIGHEST_PROTOCOL)
 del old
 old=pickle.load((S11/'00_restricted/EXPOSURE_SCALING_RESTRICTED.pkl').open('rb'));selected={k:v for k,v in old.items() if k in GENES}
 with (PRIVATE/'FROZEN_SCALES_RESTRICTED.pkl').open('wb') as f:pickle.dump(selected,f,pickle.HIGHEST_PROTOCOL)
 shutil.copyfile(S11/'00_restricted/ASSEMBLED_CLINICAL_RESTRICTED.pkl',PRIVATE/'ASSEMBLED_CLINICAL_RESTRICTED.pkl')
 for n in ['engine.py','simulation_dgp.py','prepare_templates.py']:shutil.copyfile(S11/'code'/n,EVIDENCE/('stage11_source_'+n+'.txt'))
 save_json(RESULTS/'PROTOCOL_ACCESS_RECORD.json',{'time':now(),'protocol_sha256':sha(ROOT/'STAGE12_PROTOCOL_v1.0.md'),'stage11_results_already_known':True,'result_informed_explanation':True,'new_functional_background_or_paired_change_target_results_fitted':False,'fixed_genes':GENES,'formal_hypotheses':42,'prior_sealed_entries':len(prior),'prior_unique_paths':len(cache),'authorization':'Explicit request to complete Stage12; source stages remain read-only.'})
 log(f'Protocol fixed; {len(prior)} prior sealed entries unchanged. Stage11 models/scales copied read-only.')
if __name__=='__main__':run()

