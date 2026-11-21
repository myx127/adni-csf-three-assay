from config import *
from four_assay import moments
from model import wilson
from scipy.stats import norm
from concurrent.futures import ProcessPoolExecutor,as_completed
import pickle

def worker(task):
 ident,d,tol,sign=task;w=np.log(tol);seed=int(np.random.SeedSequence([2026102907,int(hashlib.sha256(ident.encode()).hexdigest()[:8],16),int(tol*100),sign+1]).generate_state(1,dtype=np.uint64)[0]);rng=np.random.default_rng(seed)
 mats=[];qrows=[]
 for axis in range(4):
  lab=d['groups'][:,axis];one=d['onehots'][axis];age=d['q_age'];q=age-(one@age/one.sum(axis=1))[lab];qrows.append((lab,len(one),q));mats.append(one.T*q[:,None])
 A=np.column_stack(mats);A-=d['Q']@(d['Q'].T@A);K=.15**2*np.sum(A*A)/d['df']
 R=np.full((4,4),.4);np.fill_diagonal(R,1.);R[0,3]=R[3,0]=(.4+K)*np.exp(-sign*w)-K;chol=np.linalg.cholesky(R);p=[]
 for first in range(0,2000,64):
  B=min(64,2000-first);Y=rng.normal(size=(B,d['n'],4))@chol.T
  for lab,G,q in qrows:Y+=rng.normal(scale=.15,size=(B,G))[:,lab,None]*q[None,:,None]
  _,_,cc,V,_=moments(Y,d)
  for c,v in zip(cc,V):
   if min(c[1],c[2],c[3],c[4])<=0:p.append(1.);continue
   L=np.log(c[1])+np.log(c[4])-np.log(c[2])-np.log(c[3]);g=np.zeros(6);g[1]=1/c[1];g[4]=1/c[4];g[2]=-1/c[2];g[3]=-1/c[3];se=np.sqrt(g@v@g)
   p.append(max(norm.sf((L+w)/se),norm.sf((w-L)/se)))
 p=np.array(p);k=int((p<=.05).sum());lo,hi=wilson(k,len(p));file=f'FOUR_EQ_BOUNDARY__{ident}__f{int(tol*100)}__s{sign}.parquet'
 pd.DataFrame({'replicate':np.arange(len(p)),'p_equiv':p}).to_parquet(SIM/file,index=False)
 return {'design_id':ident,'ratio_tolerance_factor':tol,'boundary_sign':sign,'reps':len(p),'false_equivalence_count':k,'rate':k/len(p),
  'Wilson95_low':lo,'Wilson95_high':hi,'seed':str(seed),'file':file,'scenario':'Gaussian covariance boundary with four common random slope axes .15','before_strong_cov_gate':True}

if __name__=='__main__':
 ds=pickle.load((RESTRICTED/'FOUR_ASSAY_DESIGNS_RESTRICTED.pkl').open('rb'));tasks=[(k,d,t,s) for k,d in ds.items() for t in [1.25,1.5] for s in [-1,1]]
 save_json(RESULTS/'FOUR_DIRECTION_BOUNDARY_NOTE.json',{'already_seen_four_assay_point_results':True,'not_a_new_discovery_definition':True,'working_tolerance_factors':[1.25,1.5],
  'code_sha256':sha(Path(__file__)),'boundary_scenarios_finite_and_not_universal':True})
 rows=[]
 with ProcessPoolExecutor(max_workers=6) as pool:
  for f in as_completed([pool.submit(worker,t) for t in tasks]):rows.append(f.result())
 csvout(RESULTS/'FOUR_DIRECTION_EQ_BOUNDARY_CALIBRATION.csv',pd.DataFrame(rows).sort_values(['design_id','ratio_tolerance_factor','boundary_sign']))
 save_json(RESULTS/'FOUR_DIRECTION_EQ_CALIBRATION_RUN.json',{'cells':len(tasks),'independent_replicates':len(tasks)*2000})
 log('Four-assay working-direction equivalence boundary checks complete; no primary or exploratory discovery thresholds changed')
