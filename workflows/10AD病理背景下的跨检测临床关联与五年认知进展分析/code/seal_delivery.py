from config import *
def run():
 science=json.loads((RESULTS/'SCIENTIFIC_QA_SUMMARY.json').read_text(encoding='utf-8'));delivery=json.loads((RESULTS/'DELIVERY_QA_SUMMARY.json').read_text(encoding='utf-8'))
 assert science['failed']==0 and delivery['failed']==0
 primary=pd.read_csv(RESULTS/'PRIMARY6_MEMORY_PROGRESSION_RESULTS.csv')
 np_same=primary.loc[primary.gene.eq('NPTX1')&primary.test.eq('SAME')].iloc[0]
 log('Stage10 complete:all required analyses and QA finished;no new clinical BY discoveries;three unqualified amyloid interactions explicitly retained. Sealing delivery.')
 save_json(ROOT/'STAGE10_COMPLETION_STATUS.json',{'stage':10,'status':'COMPLETE','date':'2026-10-29','timezone':'Asia/Shanghai','scientific_QA':science,'delivery_QA':delivery,'core_candidates':['NPTX1','KLK6','PDGFRB'],'three_assay_designs':43,'formal_designs':16,'qualified_formal_designs':13,'sensitivity_designs':27,'qualified_all_designs':40,'independent_noise_draws':2410000,'pattern_evaluations_reusing_noise':89170000,'main_memory_n':{'NPTX1':552,'KLK6':540,'PDGFRB':552},'primary_hypotheses':6,'secondary_hypotheses':26,'new_BY_clinical_discoveries':0,'NPTX1_main_same_p':float(np_same.p_computed),'NPTX1_main_same_BY6_q':float(np_same.q_BY_within_family),'clinical_prognostic_validation_established':False,'amyloid_modification_inference_resolved':False,'clinical_equivalence_established':False,'clinical_risk_prediction_model_fitted':False,'independent_external_validation':False,'prior_stage_files_modified':False,'created':now()})
 rows=[]
 for p in sorted(ROOT.rglob('*')):
  if not p.is_file() or p.name=='SHA256_MANIFEST.csv' or '__pycache__' in p.parts or p.suffix in ['.pyc','.tmp'] or p.name.startswith('~$'):continue
  rows.append({'relative_path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p),'restricted':p.is_relative_to(PRIVATE),'category':'SIMULATION' if p.is_relative_to(SIM) else 'RESTRICTED' if p.is_relative_to(PRIVATE) else 'DELIVERY'})
 csvout(ROOT/'SHA256_MANIFEST.csv',rows)
 print(json.dumps({'status':'SEALED','files':len(rows),'bytes':sum(r['bytes'] for r in rows),'clinical_BY_discoveries':0,'unresolved_amyloid_interaction_designs':3},ensure_ascii=False))
if __name__=='__main__':
 assert_unsealed();run()
