from config import *
from prepare_design import verify
verify()
frozen=json.loads((RESULTS/'SIMULATION_CONFIG_FROZEN.json').read_text(encoding='utf-8'))
assert sha(RESTRICTED/'CLINICAL_DESIGN_RESTRICTED.csv')==frozen['clinical_design_sha256']
baseline=RESULTS/'RESULT_REPRODUCIBILITY_BASELINE.json'
if not baseline.exists():save_json(baseline,{name:sha(RESULTS/name) for name in ['PRIMARY_ALL43_RESULTS.csv','PRIMARY_FROZEN_FAMILY_RESULTS.csv','ALL_FIXED_VARIANT_RESULTS.csv','MAIN_RUN_SUMMARY.json','COVARIANCE_OFFSET_GRID125.csv']})
log('All source hashes and frozen clinical design verified; existing SAP and first-access record retained')
