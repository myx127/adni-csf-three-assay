from config import *
def run():
 science=json.loads((RESULTS/'SCIENTIFIC_QA_SUMMARY.json').read_text(encoding='utf-8'));delivery=json.loads((RESULTS/'DELIVERY_QA_SUMMARY.json').read_text(encoding='utf-8'))
 assert science['failed']==0 and delivery['failed']==0
 real=json.loads((RESULTS/'REAL_ANALYSIS_SUMMARY.json').read_text(encoding='utf-8'));cal=json.loads((RESULTS/'CALIBRATION_RUN_SUMMARY.json').read_text(encoding='utf-8'))
 log('All Stage11 analyses, scientific checks, figure review and delivery checks complete. Sealing.')
 save_json(ROOT/'STAGE11_COMPLETION_STATUS.json',{'stage':11,'status':'COMPLETE','created':now(),'scope':'CSF fixed34 family with focus14; clinical association evaluation','fixed_family':34,'fixed_focus_subset':14,'formal_hypotheses':204,'primary_family':68,'secondary_family':136,'scientific_QA':science,'delivery_QA':delivery,'calibration':cal,'results':real,'prior_stage_files_modified':False,'independent_validation_established':False,'clinical_equivalence_established':False,'causal_inference_established':False,'individual_prediction_model_developed':False,'blood_or_FDG_analyzed':False})
 rows=[]
 for p in sorted(ROOT.rglob('*')):
  if not p.is_file() or p.name=='SHA256_MANIFEST.csv' or '__pycache__' in p.parts or p.suffix in ['.pyc','.tmp'] or p.name.startswith('~$'):continue
  rows.append({'relative_path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p),'restricted':p.is_relative_to(PRIVATE),'category':'SIMULATION' if p.is_relative_to(SIM) else 'RESTRICTED' if p.is_relative_to(PRIVATE) else 'DELIVERY'})
 csvout(ROOT/'SHA256_MANIFEST.csv',rows)
 print(json.dumps({'status':'SEALED','files':len(rows),'bytes':sum(r['bytes'] for r in rows),'primary_focus14':real['focus14_primary_supported'],'secondary_focus14':real['focus14_secondary_supported']},ensure_ascii=False))
if __name__=='__main__':assert_unsealed();run()

