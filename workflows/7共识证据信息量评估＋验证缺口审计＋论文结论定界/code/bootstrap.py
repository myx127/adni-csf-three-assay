from config import *

def run():
 refs=['model.py','diagnostics.py','equivalence.py','designs.py']
 for name in refs:
  dest=ROOT/'code'/name;body=(S6/'code'/name).read_bytes()
  if dest.exists():assert dest.read_bytes()==body
  else:dest.write_bytes(body)
 csvout(RESULTS/'REFERENCE_CODE_ORIGIN.csv',[{'file':n,'source':str(S6/'code'/n),'sha256':sha(S6/'code'/n)} for n in refs])
 rows=[]
 for stage,manifest,col in [(S1,'STAGE1_INPUT_OUTPUT_SHA256.csv','path'),(S4,'STAGE4_INPUT_OUTPUT_SHA256.csv','relative_path'),(S5,'STAGE5_INPUT_OUTPUT_SHA256.csv','path'),(S6,'STAGE6_INPUT_OUTPUT_SHA256.csv','path')]:
  for r in pd.read_csv(stage/manifest).to_dict('records'):
   p=Path(r[col]);p=p if p.is_absolute() else stage/p;actual=sha(p)
   rows.append({'stage':stage.name,'manifest_path':r[col],'sha256':actual,'expected_sha256':r['sha256'],'passed':actual==r['sha256']})
 csvout(RESULTS/'PRIOR_SEALED_FILE_CHECKS.csv',rows);assert all(r['passed'] for r in rows)
 sources=[{'source_id':k,'path':v['path'],'sha256':v['sha256']} for k,v in INPUTS.items()]
 for manifest in ['ADDITIONAL_INPUT_MANIFEST.csv','LEGACY_INPUT_MANIFEST.csv']:
  sources.extend({'source_id':r['input_id'],'path':r['path'],'sha256':r['sha256']} for r in pd.read_csv(S6/'results'/manifest).to_dict('records'))
 for r in sources:r['bytes']=Path(r['path']).stat().st_size;r['passed']=sha(r['path'])==r['sha256']
 csvout(RESULTS/'INPUT_SOURCE_CHECKS.csv',sources);assert all(r['passed'] for r in sources)
 freeze={'date':datetime.now().astimezone().isoformat(),'protocol_sha256':sha(ROOT/'STAGE7_PROTOCOL_v1.0.md'),'already_seen_Stage6_results':True,
  'family34':FAMILY,'margins_unchanged':MARGINS,'k_grid':K_GRID,'gamma_grid':GAMMAS,'scenarios':SCENARIOS,'reps_per_cell':2000,
  'root_seed':SEED_ROOT,'alpha_BY_rank1_sufficient':ALPHA_BY_SUFFICIENT,'variance_version':0,'prior_manifest_entries_verified':len(rows),'no_new_candidate_association_tests':True}
 dest=RESULTS/'STAGE7_PLAN_FROZEN.json'
 if dest.exists():assert json.loads(dest.read_text(encoding='utf-8'))['protocol_sha256']==freeze['protocol_sha256']
 else:save_json(dest,freeze)
 log(f'Stage7 scope frozen; {len(rows)} prior sealed entries and {len(sources)} raw inputs verified unchanged')

if __name__=='__main__':run()
