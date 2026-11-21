from config import *
from verify_delivery import excluded

def run():
 qa=json.loads((RESULTS/'INDEPENDENT_QA_SUMMARY.json').read_text(encoding='utf-8'));delivery=json.loads((RESULTS/'DELIVERY_QA_SUMMARY.json').read_text(encoding='utf-8'))
 assert qa['failed']==delivery['failed']==0
 save_json(RESULTS/'STAGE6_COMPLETION_STATUS.json',{'status':'COMPLETE_AVAILABLE_SOURCE_ANALYSES_WITH_LIMITED_ORTHOGONAL_SUPPORT',
  'date':'2026-10-27','primary_Stage5_discoveries_unchanged':True,'attribution_models':140,
  'direct_discovery_genes_with_legacy_coverage':2,'direct_peptides':5,'new_subject_disjoint_ADNI_RID':95,
  'formal_relative_disagreement_replications':0,'external_cohort_validations':0,'molecular_mechanisms_verified':False,
  'operational_equivalence_family':34,'equivalence_eligible':11,'BY_equivalent_at_each_working_margin':0,
  'simulation_cells':1010,'independent_simulation_replicates':3160000,'analysis_QA_passed':qa['passed'],'delivery_QA_passed':delivery['passed'],
  'protein_peptide_aggregate_workbook_total_sheets':69,'figures_PNG_PDF_pairs':8,
  '2014_source_scope':'Coverage/source audit only; no verified exact-date/technical crosswalk, therefore no association fits',
  'not_completed_and_not_claimed':['Wet-lab experiments','External cohort replication','Antibody/aptamer epitope identification','Original TMT peptide-level localization','Arbitrary-error-covariance abundance decomposition'],
  'individual_data_directory':'00_restricted; local reproduction only; directory name does not configure OS ACL',
  'manifest_excludes':['itself','__pycache__','*.pyc','*.tmp','*.part']})
 rows=[]
 for p in sorted(ROOT.rglob('*')):
  if not p.is_file() or excluded(p):continue
  rel=p.relative_to(ROOT).as_posix();role='individual_local_only' if rel.startswith('00_restricted/') else 'simulation_replicate' if rel.startswith('results/simulation_replicates/') else 'public_reference_source' if rel.startswith('evidence/') else 'analysis_or_delivery'
  rows.append({'path':rel,'bytes':p.stat().st_size,'sha256':sha(p),'role':role})
 csvout(ROOT/'STAGE6_INPUT_OUTPUT_SHA256.csv',rows)
 # Do not call log(): any writing after the hash would invalidate the manifest.
 print(json.dumps({'sealed_files':len(rows),'analysis_QA':qa['passed'],'delivery_QA':delivery['passed'],'status':'COMPLETE_AVAILABLE_SOURCE_ANALYSES_WITH_LIMITED_ORTHOGONAL_SUPPORT'},ensure_ascii=False))

if __name__=='__main__':run()
