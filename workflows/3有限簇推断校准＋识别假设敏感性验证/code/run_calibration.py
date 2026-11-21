from __future__ import annotations
import os
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
import argparse,importlib.metadata,json,platform,sys,time,traceback
import numpy as np
import pandas as pd
from cal_config import *
from cal_data import generate
from cal_inference import infer,sensitivity_profile

CORE=['cal_config.py','cal_data.py','cal_inference.py','run_calibration.py','analyze_calibration.py','reference_stage2_models.py','reference_stage2_data.py','sim_config.py']

def bind_inputs():
 sources=[]
 for directory,manifest_name,names in [(STAGE1,'STAGE1_INPUT_OUTPUT_SHA256.csv',['results/AUDIT_SUMMARY.json','results/TRIPLE_BATCH_SAMPLE_SIZES.csv']),
    (STAGE2,'SIMULATION_INPUT_OUTPUT_SHA256.csv',['code/candidate_models.py','code/synthetic_data.py'])]:
  manifest=pd.read_csv(directory/manifest_name)
  for name in names:
   row=manifest.loc[manifest.path.eq(name)].iloc[0];path=directory/name
   sources.append({'path':str(path),'expected_sha256':row.sha256,'observed_sha256':sha(path),'passed':row.sha256==sha(path)})
 csv(RESULTS/'INPUT_SOURCE_HASH_CHECKS.csv',sources)
 if not all(r['passed'] for r in sources):raise RuntimeError('Frozen source hash mismatch')
 assert sha(ROOT/'code/reference_stage2_models.py')==sha(STAGE2/'code/candidate_models.py')
 assert sha(ROOT/'code/reference_stage2_data.py')==sha(STAGE2/'code/synthetic_data.py')
 sizes=pd.read_csv(STAGE1/'results/TRIPLE_BATCH_SAMPLE_SIZES.csv')
 audited={p:d.n_triple_participants.to_list() for p,d in sizes.groupby('platform')}
 assert {p:len(v) for p,v in audited.items()}=={'NULISA':19,'SOMASCAN':36,'TMT_MS':65}
 assert all(sum(v)==571 for v in audited.values())
 save_json(RESULTS/'AUDITED_BATCH_SIZES.json',audited)
 save_json(RESULTS/'ACCESS_AND_PROTOCOL_RECORD.json',{'clinical_files_opened':0,'individual_protein_files_opened':0,'restricted_files_opened':0,
    'data':'SYNTHETIC_ONLY','source_files':sources,'protocol_sha256_before_experiments':sha(ROOT/'CALIBRATION_PROTOCOL_v1.0.md'),
    'core_code_before_experiments':{f:sha(ROOT/'code'/f) for f in CORE},'client_date':'2026-10-27','timezone':'Asia/Shanghai'})
 return audited

def frozen_core_unchanged(record):return all(sha(ROOT/'code'/name)==value for name,value in record.items())

def main(phase):
 if phase=='DEVELOPMENT':
  audited=bind_inputs()
  registry=[]
  for family,items in [('CALIBRATION',calibration_scenarios()),('SENSITIVITY',sensitivity_scenarios())]:
   for r in items:registry.append({'family':family,**r,'settings':json.dumps(r['settings'],ensure_ascii=False,sort_keys=True)})
  csv(RESULTS/'SCENARIO_REGISTRY.csv',registry);csv(RESULTS/'EXPERIMENT_GRID.csv',grid())
  save_json(RESULTS/'SOFTWARE_ENVIRONMENT.json',{'python':sys.version,'platform':platform.platform(),'threads':1,
    'packages':{p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','pyarrow','openpyxl','matplotlib']}})
 else:
  audited=json.loads((RESULTS/'AUDITED_BATCH_SIZES.json').read_text(encoding='utf-8'))
  access=json.loads((RESULTS/'ACCESS_AND_PROTOCOL_RECORD.json').read_text(encoding='utf-8'))
  assert frozen_core_unchanged(access['core_code_before_experiments']),'Core code changed after development'
  selection=json.loads((RESULTS/'DEVELOPMENT_SELECTION_FROZEN.json').read_text(encoding='utf-8'))
  assert frozen_core_unchanged(selection['core_code_sha256']),'Core code changed after selection'
  assert sha(ROOT/'CALIBRATION_PROTOCOL_v1.0.md')==access['protocol_sha256_before_experiments']
 lookup={r['scenario_id']:r for r in (sensitivity_scenarios() if phase=='SENSITIVITY' else calibration_scenarios())}
 target=RESULTS/'replicates';target.mkdir(exist_ok=True)
 start=time.time();total=0
 log('START '+phase)
 for cell in [c for c in grid() if c['phase']==phase]:
  scenario=lookup[cell['scenario_id']];rows=[];profiles=[]
  for rep in range(cell['repetitions']):
   rng=np.random.default_rng(np.random.SeedSequence([int(cell['seed']),rep,0]))
   bootrng=np.random.default_rng(np.random.SeedSequence([int(cell['seed']),rep,1]))
   data=generate(cell['n'],rng,scenario['settings'],audited)
   results,cache=infer(data,cell['cluster_axis'],cell['bootstrap_B'],bootrng,need_cache=phase=='SENSITIVITY')
   counts=np.bincount(data['groups'],minlength=3)
   for method,result in results.items():
    rows.append({'replicate':rep,'method':method,**result,'n_complete':len(data['Y']),
       'smallest_group_n':int(counts.min()),'clusters_active':int(len(np.unique(data['batches'][:,cell['cluster_axis']]))),
       'minimum_clusters':min(len(np.unique(data['batches'][:,j])) for j in range(3))})
   if phase=='SENSITIVITY':profiles.append({'replicate':rep,**sensitivity_profile(data,cache)})
  frame=pd.DataFrame(rows)
  for column in ['method','status']:frame[column]=frame[column].astype('category')
  frame.to_parquet(target/(cell['cell_id']+'.parquet'),index=False)
  if profiles:pd.DataFrame(profiles).to_parquet(target/(cell['cell_id']+'_sensitivity.parquet'),index=False)
  total+=cell['repetitions']
  log(f"Complete {cell['cell_id']} R={cell['repetitions']} B={cell['bootstrap_B']} phase_total={total} elapsed={time.time()-start:.1f}s")
 save_json(RESULTS/(phase+'_RUN_SUMMARY.json'),{'phase':phase,'cells':sum(c['phase']==phase for c in grid()),'independent_replicates':total,
      'method_rows':total*len(METHODS),'elapsed_seconds':time.time()-start,'clinical_files_opened':0})

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--phase',choices=['DEVELOPMENT','CONFIRMATION','N_GRID','SENSITIVITY'],required=True)
 args=parser.parse_args()
 try:main(args.phase)
 except Exception:
  (LOGS/(args.phase+'_FAILURE.txt')).write_text(traceback.format_exc(),encoding='utf-8');raise
