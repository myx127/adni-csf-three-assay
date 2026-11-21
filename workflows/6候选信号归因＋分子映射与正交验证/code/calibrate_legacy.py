from config import *
from legacy_assay import single_batch
from model import wilson
from concurrent.futures import ProcessPoolExecutor,as_completed
import pickle

SCENARIOS=['L01_IID','L02_HETEROSKEDASTIC','L03_T5','L04_PLATE_RANDOM_SLOPES','L05_RUN_RANDOM_SLOPES','L06_BOTH_RANDOM_SLOPES']

def worker(task):
 ident,d,scenario=task;seed=int(np.random.SeedSequence([2026102908,int(hashlib.sha256(ident.encode()).hexdigest()[:8],16),SCENARIOS.index(scenario)]).generate_state(1,dtype=np.uint64)[0]);rng=np.random.default_rng(seed)
 qrows=[]
 for axis in range(2):
  one=d['onehots'][axis];lab=d['groups'][:,axis];q=d['q_age']-(one@d['q_age']/one.sum(axis=1))[lab];qrows.append((lab,len(one),q))
 values=[]
 for start in range(0,4000,128):
  B=min(128,4000-start);y=rng.standard_t(5,size=(B,d['n']))/np.sqrt(5/3) if scenario=='L03_T5' else rng.normal(size=(B,d['n']))
  if scenario=='L02_HETEROSKEDASTIC':y*=np.exp(.3*np.clip(d['P'],-3,3))[None,:]
  axes=[0] if scenario=='L04_PLATE_RANDOM_SLOPES' else [1] if scenario=='L05_RUN_RANDOM_SLOPES' else [0,1] if scenario=='L06_BOTH_RANDOM_SLOPES' else []
  for axis in axes:
   lab,G,q=qrows[axis];y+=rng.normal(scale=.6 if len(axes)==1 else .4,size=(B,G))[:,lab]*q[None,:]
  values.extend(single_batch(y,d)[2])
 p=np.array(values);valid=np.isfinite(p);k=int((valid&(p<=.05)).sum());lo,hi=wilson(k,len(p));file=f'LEGACY_SINGLE__{ident}__{scenario}.parquet'
 pd.DataFrame({'replicate':np.arange(len(p)),'p_positive':p}).to_parquet(SIM/file,index=False)
 return {'design_id':ident,'scenario':scenario,'reps':len(p),'n_failed':int((~valid).sum()),'n_rejected05':k,'rate':k/len(p),
  'Wilson95_low':lo,'Wilson95_high':hi,'seed':str(seed),'file':file,'variance':'SUM2_CR1','tail':'positive, Student t(min cluster G-1)'}

if __name__=='__main__':
 ds=pickle.load((RESTRICTED/'LEGACY_CLINICAL_DESIGNS_RESTRICTED.pkl').open('rb'));tasks=[(k,d,s) for k,d in ds.items() for s in SCENARIOS]
 save_json(RESULTS/'LEGACY_CALIBRATION_PLAN_FROZEN.json',{'before_legacy_candidate_association_fits':True,'scenarios':SCENARIOS,
  'reps_per_cell':4000,'code_sha256':sha(Path(__file__)),'calibration_pass_rule':'All Wilson95 upper<=.06 and nonfinite fraction<=.005',
  'finite_scenarios_do_not_establish_universal_validity':True,'no_switching_variance_estimator_after_calibration':True})
 rows=[]
 with ProcessPoolExecutor(max_workers=6) as pool:
  for i,f in enumerate(as_completed([pool.submit(worker,t) for t in tasks]),1):
   rows.append(f.result())
   if i%12==0:log(f'Legacy single-assay calibration {i}/{len(tasks)}')
 out=pd.DataFrame(rows).sort_values(['design_id','scenario']);csvout(RESULTS/'LEGACY_SINGLE_ASSAY_NULL_CALIBRATION.csv',out)
 decisions=[]
 for k,part in out.groupby('design_id'):
  decisions.append({'design_id':k,'passed':bool((part.Wilson95_high<=.06).all() and (part.n_failed/part.reps<=.005).all()),
    'max_rate':part.rate.max(),'max_Wilson95_upper':part.Wilson95_high.max(),'max_failure_rate':(part.n_failed/part.reps).max()})
 csvout(RESULTS/'LEGACY_SINGLE_ASSAY_CALIBRATION_DECISIONS.csv',decisions)
 save_json(RESULTS/'LEGACY_SINGLE_ASSAY_CALIBRATION_RUN.json',{'cells':len(tasks),'independent_replicates':len(tasks)*4000})
 log('Legacy single-assay finite-design calibration complete')
