from __future__ import annotations
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
import json,time
import numpy as np
import pandas as pd
from cal_config import *
from cal_data import generate
from cal_inference import infer

def main():
 audited=json.loads((RESULTS/'AUDITED_BATCH_SIZES.json').read_text(encoding='utf-8'))
 cases={r['scenario_id']:r for r in calibration_scenarios()};rows=[];start=time.time()
 for cell in [c for c in grid() if c['phase']=='CONFIRMATION' and c['scenario_id'] in ['C02','C03','C05','A01']]:
  original=pd.read_parquet(RESULTS/'replicates'/(cell['cell_id']+'.parquet')).set_index(['replicate','method'])
  for rep in range(500):
   data=generate(cell['n'],np.random.default_rng(np.random.SeedSequence([int(cell['seed']),rep,0])),cases[cell['scenario_id']]['settings'],audited)
   out,_=infer(data,cell['cluster_axis'],1999,np.random.default_rng(np.random.SeedSequence([int(cell['seed']),rep,1])))
   for method,result in out.items():
    before=original.loc[(rep,method)]
    rows.append({'scenario_id':cell['scenario_id'],'replicate':rep,'method':method,'original_B':499,'replay_B':1999,
       'original_p':before.p_value,'replay_p':result['p_value'],'original_status':str(before.status),'replay_status':result['status'],
       'classification_changed':bool(np.isfinite(before.p_value) and np.isfinite(result['p_value']) and (before.p_value<=.05)!=(result['p_value']<=.05))})
  log(f'BOOTSTRAP_RESOLUTION {cell["scenario_id"]} reused500 elapsed={time.time()-start:.1f}s')
 frame=pd.DataFrame(rows);frame.to_parquet(RESULTS/'BOOTSTRAP_RESOLUTION_REPLAY.parquet',index=False)
 summary=frame.groupby(['scenario_id','method']).agg(replays=('replicate','size'),classification_change_rate=('classification_changed','mean'),
      original_rejection=('original_p',lambda x:x.le(.05).mean()),high_B_rejection=('replay_p',lambda x:x.le(.05).mean())).reset_index()
 csv(RESULTS/'BOOTSTRAP_RESOLUTION_SUMMARY.csv',summary)
 save_json(RESULTS/'BOOTSTRAP_RESOLUTION_RECORD.json',{'original_confirmation_seeds_replayed':2000,'additional_independent_replicates':0,'method_rows':len(frame),'B':1999})

if __name__=='__main__':main()
