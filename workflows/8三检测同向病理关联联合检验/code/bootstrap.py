from config import *
import pickle

def verify_prior():
 rows=[];known={}
 sources=[(S1,'STAGE1_INPUT_OUTPUT_SHA256.csv','path'),(S4,'STAGE4_INPUT_OUTPUT_SHA256.csv','relative_path'),(S5,'STAGE5_INPUT_OUTPUT_SHA256.csv','path'),(S6,'STAGE6_INPUT_OUTPUT_SHA256.csv','path'),(S7,'SHA256_MANIFEST.csv','relative_path')]
 for stage,name,col in sources:
  for r in pd.read_csv(stage/name).to_dict('records'):
   p=Path(r[col]);p=p if p.is_absolute() else stage/p;key=str(p.resolve())
   if key not in known:known[key]=sha(p)
   rows.append({'stage':stage.name,'path':str(p),'expected_sha256':r['sha256'],'actual_sha256':known[key],'passed':known[key]==r['sha256']})
 return rows

def run():
 protocol=ROOT/'STAGE8_PROTOCOL_v1.0.md';freeze=RESULTS/'STAGE8_PLAN_FROZEN.json'
 if freeze.exists():assert json.loads(freeze.read_text(encoding='utf-8'))['protocol_sha256']==sha(protocol)
 else:save_json(freeze,{'date':datetime.now().astimezone().isoformat(),'protocol_sha256':sha(protocol),'known_prior_results':True,'before_new_shared_P_computed':True,'family':FAMILY,'models':MODELS,'scenarios':SCENARIOS,'reps_per_noise_cell':REPS,'expected_noise_cells':525,'expected_independent_noise_repetitions':21000000,'alphas':ALPHAS,'patterns':PATTERNS,'root_seed':ROOT_SEED})
 rows=verify_prior();csvout(RESULTS/'PRIOR_SEALED_FILE_CHECKS.csv',rows);assert all(r['passed'] for r in rows)
 refs=[]
 for name in ['model.py','diagnostics.py','designs.py']:
  body=(S7/'code'/name).read_bytes();dest=ROOT/'code'/name
  if dest.exists():assert dest.read_bytes()==body
  else:dest.write_bytes(body)
  refs.append({'file':name,'source':str(S7/'code'/name),'sha256':sha(dest)})
 csvout(RESULTS/'REFERENCE_CODE_ORIGIN.csv',refs)
 source_paths=[S5/'00_restricted/DESIGNS_MAIN_RESTRICTED.pkl',S6/'00_restricted/ATTRIBUTION_MODELS_RESTRICTED.pkl',S6/'results/COVARIATE_FOUR_MODEL_RESULTS.csv',S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz',S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv',S5/'results/PRIMARY_FROZEN_FAMILY_RESULTS.csv',S7/'results/CONSENSUS_INFORMATION_BOTTLENECKS.csv']
 csvout(RESULTS/'SOURCE_MANIFEST.csv',[{'source_id':f'S{i:02d}','path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)} for i,p in enumerate(source_paths,1)])
 log(f'Protocol frozen; {len(rows)} prior sealed entries verified; three reference modules copied unchanged')

if __name__=='__main__':run()
