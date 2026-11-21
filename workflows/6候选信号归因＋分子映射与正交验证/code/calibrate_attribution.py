from config import *
from simulation_reference import simulate_cell,NULLS
from concurrent.futures import ProcessPoolExecutor,as_completed
import pickle,time

if __name__=='__main__':
 ds=pickle.load((RESTRICTED/'ATTRIBUTION_CALIBRATION_DESIGNS_RESTRICTED.pkl').open('rb'));tasks=[(key,d,sc,4000,'SECONDARY',0.,1) for key,d in ds.items() for sc in NULLS]
 rows=[];start=time.time()
 with ProcessPoolExecutor(max_workers=6) as pool:
  for i,f in enumerate(as_completed([pool.submit(simulate_cell,t) for t in tasks]),1):
   rows.extend(f.result())
   if i%42==0 or i==len(tasks):log(f'Attribution calibration {i}/{len(tasks)} cells elapsed={time.time()-start:.1f}s')
 out=pd.DataFrame(rows).sort_values(['design_id','scenario','version','alpha']);csvout(RESULTS/'ATTRIBUTION_CALIBRATION_SUMMARY.csv',out)
 selected=out.loc[out.version.eq(0)&out.alpha.eq(.05)];decisions=[]
 for key,part in selected.groupby('design_id'):
  decisions.append({'design_id':key,'passed':bool((part.wilson95_high<=.06).all() and (part.n_failed/part.reps<=.005).all()),
   'max_rate':part.rate.max(),'max_Wilson95_upper':part.wilson95_high.max(),'max_failure_rate':(part.n_failed/part.reps).max()})
 csvout(RESULTS/'ATTRIBUTION_CALIBRATION_DECISIONS.csv',decisions)
 save_json(RESULTS/'ATTRIBUTION_CALIBRATION_RUN.json',{'cells':len(tasks),'independent_replicates':len(tasks)*4000,'seed':SEEDS['SECONDARY'],'selected_variance_version_changed':False})
