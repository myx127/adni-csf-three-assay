from config import *
from simulation import simulate_cell,NULLS
from concurrent.futures import ProcessPoolExecutor,as_completed
import pickle,time

if __name__=='__main__':
 ds=pickle.load((RESTRICTED/'DESIGNS_VARIANTS_RESTRICTED.pkl').open('rb'));rows=[];start=time.time()
 tasks=[(ident,d,sc,4000,'SECONDARY',0.,1) for ident,d in ds.items() for sc in NULLS]
 with ProcessPoolExecutor(max_workers=6) as pool:
  futures=[pool.submit(simulate_cell,t) for t in tasks]
  for i,f in enumerate(as_completed(futures),1):
   rows.extend(f.result())
   if i%24==0 or i==len(tasks):log(f'Variant calibration {i}/{len(tasks)} cells elapsed={time.time()-start:.1f}s')
 csvout(RESULTS/'VARIANT_CALIBRATION_SUMMARY.csv',rows)
 save_json(RESULTS/'VARIANT_CALIBRATION_RUN_SUMMARY.json',{'unique_designs':len(ds),'cells':len(tasks),'independent_replicates':4000*len(tasks),'seconds':time.time()-start})
