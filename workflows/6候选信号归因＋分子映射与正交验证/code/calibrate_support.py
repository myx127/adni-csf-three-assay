from config import *
from simulation_reference import simulate_cell,NULLS
from concurrent.futures import ProcessPoolExecutor,as_completed
import pickle

if __name__=='__main__':
 ds=pickle.load((RESTRICTED/'ORTHOGONAL_SUPPORT_DESIGNS_RESTRICTED.pkl').open('rb'));tasks=[(key,d,sc,4000,'CONFIRM',0.,1) for key,d in ds.items() for sc in NULLS]
 rows=[]
 with ProcessPoolExecutor(max_workers=6) as pool:
  for i,f in enumerate(as_completed([pool.submit(simulate_cell,t) for t in tasks]),1):
   rows.extend(f.result())
   if i%12==0:log(f'Orthogonal-support calibration {i}/{len(tasks)}')
 out=pd.DataFrame(rows).sort_values(['design_id','scenario','version','alpha']);csvout(RESULTS/'ORTHOGONAL_SUPPORT_CALIBRATION_SUMMARY.csv',out)
 decisions=[]
 for key,q in out.loc[out.version.eq(0)&out.alpha.eq(.05)].groupby('design_id'):
  decisions.append({'design_id':key,'passed':bool((q.wilson95_high<=.06).all() and (q.n_failed/q.reps<=.005).all()),
   'max_rate':q.rate.max(),'max_Wilson95_upper':q.wilson95_high.max(),'max_failure_rate':(q.n_failed/q.reps).max()})
 csvout(RESULTS/'ORTHOGONAL_SUPPORT_CALIBRATION_DECISIONS.csv',decisions)
 save_json(RESULTS/'ORTHOGONAL_SUPPORT_CALIBRATION_RUN.json',{'cells':len(tasks),'independent_replicates':len(tasks)*4000,'fixed_variance_version':0})
