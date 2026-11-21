from config import *
from freeze_engine import verify
import verify_delivery
def run():
 verify();science=json.loads((RESULTS/'INDEPENDENT_QA_SUMMARY.json').read_text(encoding='utf-8'));delivery=json.loads((RESULTS/'DELIVERY_QA_SUMMARY.json').read_text(encoding='utf-8'))
 assert science['failed']==delivery['failed']==0
 for name,obj in [('INDEPENDENT_QA_CHECKS.csv',science),('DELIVERY_QA_CHECKS.csv',delivery)]:
  d=pd.read_csv(RESULTS/name);assert len(d)==obj['checks']>0 and d.passed.all()
 main=pd.read_csv(RESULTS/'PRIMARY_FAMILY34_SHARED_RESULTS.csv');shared=main.loc[main.shared_BY34_supported]
 status={'stage':8,'status':'COMPLETE','date':'2026-10-28','post_prior_results_supplement_not_independent_replication':True,'family_size':34,
 'M0_BY_shared':len(shared),'shared_genes':sorted(shared.gene.tolist()),'same_genes_APOE_BOTH':True,
 'M0_both_shared_and_H':['NPTX1','KLK6'],'PDGFRB_shared_q':float(main.loc[main.gene.eq('PDGFRB'),'q_BY_shared34'].iloc[0]),
 'original_BY_equivalence_count_unchanged':0,'independent_noise_repetitions':21000000,'independent_noise_cells':525,
 'pattern_evaluations_not_independent_noise_repetitions':1113000000,'scientific_QA':science,'delivery_QA':delivery,
 'prior_stage_files_modified':False,'external_individual_replication_performed':False}
 save_json(ROOT/'STAGE8_COMPLETION_STATUS.json',status);rows=[]
 for p in sorted(ROOT.rglob('*')):
  if not p.is_file() or verify_delivery.excluded(p):continue
  name=str(p.relative_to(ROOT)).replace('\\','/');rows.append({'relative_path':name,'bytes':p.stat().st_size,'sha256':sha(p),'restricted':name.startswith('00_restricted/'),'category':name.split('/')[0] if '/' in name else 'delivery'})
 csvout(ROOT/'SHA256_MANIFEST.csv',rows);print({'sealed_files':len(rows),'bytes':sum(r['bytes'] for r in rows),'scientific_QA':science['passed'],'delivery_QA':delivery['passed']});assert verify_delivery.verify()==0
if __name__=='__main__':run()
