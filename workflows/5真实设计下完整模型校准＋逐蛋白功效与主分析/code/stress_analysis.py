from config import *
from model import batch_components,factors,wilson
from simulation import setup,generate
from concurrent.futures import ProcessPoolExecutor,as_completed
from scipy.linalg import qr
import pickle

def task_runner(task):
 acc,d,kind,version=task;rng=np.random.default_rng(np.random.SeedSequence([SEEDS['STRESS'],int(hashlib.sha256(acc.encode()).hexdigest()[:8],16),['WEAK','CORR12','NONLINEAR','OMITTED_TECH'].index(kind)]))
 infer_d=d
 if kind=='OMITTED_TECH':
  infer_d=d.copy();z=d['Z'][:,:3];q,r,piv=qr(z,mode='economic',pivoting=True);rank=int((abs(np.diag(r))>1e-10).sum());qz=q[:,:rank]
  p=d['P'];rP=p-qz@(qz.T@p);ss=rP@rP
  infer_d.update({'Qz':qz,'Q':np.c_[qz,rP/np.sqrt(ss)],'a':rP/ss,'df':d['n']-rank-1,'rank':rank+1})
 spec=setup(d,'C05');qs=[]
 for i in range(0,4000,64):
  Y=generate(d,spec,min(64,4000-i),rng,stress=None if kind=='OMITTED_TECH' else kind)
  qs.extend(batch_components(Y,infer_d)['Q0'])
 qs=np.array(qs);p=np.exp(-qs/(2*factors(infer_d)[version]));valid=np.isfinite(p);k=int(np.sum(valid&(p<=.05)));lo,hi=wilson(k,len(p))
 filename=f'STRESS__{acc}__{kind}.parquet';pd.DataFrame({'replicate':np.arange(len(p)),'Q0':qs,'p_value':p}).to_parquet(SIM/filename,index=False)
 return {'accession':acc,'scenario':kind,'mechanistic_direct_offset_zero':True,'observable_H_null_expected':kind=='WEAK',
  'reps':len(p),'n_failed':int((~valid).sum()),'n_rejected05':k,'rejection_rate':k/len(p),'wilson95_low':lo,'wilson95_high':hi,
  'variance_version':version,'interpretation':'Weak-signal null check' if kind=='WEAK' else 'Assumption stress; rejection does not prove direct pathology measurement effect'}

if __name__=='__main__':
 ds=pickle.load((RESTRICTED/'DESIGNS_MAIN_RESTRICTED.pkl').open('rb'));sap=json.loads((RESULTS/'MAIN_ANALYSIS_SAP_FROZEN.json').read_text(encoding='utf-8'))
 tasks=[(acc,d,kind,sap['version']) for acc,d in ds.items() for kind in ['WEAK','CORR12','NONLINEAR','OMITTED_TECH']]
 rows=[]
 with ProcessPoolExecutor(max_workers=6) as pool:
  for i,f in enumerate(as_completed([pool.submit(task_runner,t) for t in tasks]),1):
   rows.append(f.result())
   if i%28==0:log(f'Predeclared assumption stress {i}/{len(tasks)} cells')
 csvout(RESULTS/'IDENTIFICATION_STRESS_RESULTS.csv',rows)
 save_json(RESULTS/'STRESS_RUN_SUMMARY.json',{'cells':len(tasks),'independent_replicates':4000*len(tasks),'scenario_selection_not_based_on_real_effects':True})
