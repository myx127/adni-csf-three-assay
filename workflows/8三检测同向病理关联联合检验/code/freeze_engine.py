from config import *
def record():
 files=['config.py','engine.py','run_simulations.py']
 return {'protocol_sha256':sha(ROOT/'STAGE8_PROTOCOL_v1.0.md'),'files_sha256':{n:sha(ROOT/'code'/n) for n in files},
 'template_sha256':sha(PRIVATE/'SIMULATION_TEMPLATES_RESTRICTED.pkl'),'reps':REPS,'root_seed':ROOT_SEED,'batch':BATCH,'models':MODELS,'scenarios':SCENARIOS,'patterns':PATTERNS,'before_full_calibration_and_new_real_shared_P':True,'prior_results_known':True}
def verify():
 path=RESULTS/'ENGINE_FROZEN.json';expected=record()
 if path.exists():assert json.loads(path.read_text(encoding='utf-8'))==expected
 else:save_json(path,expected)
if __name__=='__main__':verify();print('Engine/config/runner/templates frozen before synthetic verification and full simulation')
