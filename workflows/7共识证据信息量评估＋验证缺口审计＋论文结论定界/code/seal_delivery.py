"""Seal only after independent scientific and delivery QA pass."""
from config import *
import verify_delivery

def run():
 science=json.loads((RESULTS/'INDEPENDENT_QA_SUMMARY.json').read_text(encoding='utf-8'))
 delivery=json.loads((RESULTS/'DELIVERY_QA_SUMMARY.json').read_text(encoding='utf-8'))
 assert science['failed']==delivery['failed']==0
 scientific_checks=pd.read_csv(RESULTS/'INDEPENDENT_QA_CHECKS.csv');delivery_checks=pd.read_csv(RESULTS/'DELIVERY_QA_CHECKS.csv')
 assert len(scientific_checks)==science['checks']>0 and scientific_checks.passed.all()
 assert len(delivery_checks)==delivery['checks']>0 and delivery_checks.passed.all()
 status={'stage':7,'status':'COMPLETE_REQUESTED_AUDITS','date':'2026-10-27',
  'scope':['consensus_information_assessment','validation_gap_audit'],
  'post_prior_results_analysis_not_prospective_registration':True,
  'scientific_QA':science,'delivery_QA':delivery,
  'original_family_size':34,'equivalence_eligible':11,'nominal_equivalent_at_widest_working_margin':6,'BY_equivalent_all_working_margins':0,
  'new_independent_simulation_repetitions':1440000,'new_independent_simulation_cells':720,
  '2014_technical_rows_resolved_with_disclosed_precision':303,'2014_verified_specimen_dates':0,
  'conditional2014_geometry_not_formal_replication':{'NPTX1':{'n':144,'df':78},'KLK6':{'n':141,'df':75}},
  'fresh_real_candidate_associations_fitted':0,'formal_external_replication_completed':False,
  'next_required_evidence':'2014 LIMS/laboratory reference to specimen/aliquot and collection date; actual-design dependency audit and calibration',
  'prior_stage_files_modified':False,'individual_data_location':'00_restricted; name does not set OS permissions'}
 save_json(ROOT/'STAGE7_COMPLETION_STATUS.json',status)
 rows=[]
 for path in sorted(ROOT.rglob('*')):
  if not path.is_file() or verify_delivery.excluded(path):continue
  name=str(path.relative_to(ROOT)).replace('\\','/')
  rows.append({'relative_path':name,'bytes':path.stat().st_size,'sha256':sha(path),'restricted':name.startswith('00_restricted/'),'category':name.split('/')[0] if '/' in name else 'delivery'})
 csvout(ROOT/'SHA256_MANIFEST.csv',rows)
 print({'manifest_entries':len(rows),'sealed_bytes':sum(r['bytes'] for r in rows),'scientific_QA_passed':science['passed'],'delivery_QA_passed':delivery['passed']})
 assert verify_delivery.verify()==0
if __name__=='__main__':run()
