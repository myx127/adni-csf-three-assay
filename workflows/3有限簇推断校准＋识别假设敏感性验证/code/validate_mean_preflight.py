from __future__ import annotations
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
import numpy as np
from cal_config import *
from mean_followup import components,bootstrap_full_mean,followup_grid

def main():
 rng=np.random.default_rng(271771);n=120;P=rng.normal(size=n);Y=rng.normal(size=(n,3))+(P+rng.normal(size=n))[:,None]
 batches=np.c_[np.arange(n)%6,np.arange(n)%8,np.arange(n)%10]
 hstar,W=bootstrap_full_mean(Y,P,batches,7,rng,True);error=0.
 for i,w in enumerate(W):
  X=np.c_[np.ones(n),P];b=np.linalg.solve(X.T@(w[:,None]*X),X.T@(w[:,None]*Y))
  e=Y-X@b;cov=(e*w[:,None]).T@e/w.sum();slope=b[1]
  direct=np.array([slope[1]*cov[0,2]-slope[0]*cov[1,2],slope[2]*cov[0,1]-slope[0]*cov[1,2]])
  error=max(error,float(np.max(abs(hstar[i]-direct))))
 h,IF=components(Y,P);weights=np.ones(n)/n
 def weighted_h(w):
  w=w/w.sum();X=np.c_[np.ones(n),P];b=np.linalg.solve(X.T@(w[:,None]*X),X.T@(w[:,None]*Y));e=Y-X@b
  cov=(e*w[:,None]).T@e;slope=b[1]
  return np.array([slope[1]*cov[0,2]-slope[0]*cov[1,2],slope[2]*cov[0,1]-slope[0]*cov[1,2]])
 if_error=0.
 for i in [0,3,19,78]:
  unit=np.zeros(n);unit[i]=1;step=1e-6
  derivative=(weighted_h((1-step)*weights+step*unit)-weighted_h((1+step)*weights-step*unit))/(2*step)
  if_error=max(if_error,float(np.max(abs(derivative-n*IF[i]))))
 seeds={c['seed'] for c in grid()};new=[c['seed'] for c in followup_grid()]
 checks=[{'check':'weighted_OLS_bootstrap_vs_explicit_loop','passed':error<1e-10,'error':error},
   {'check':'mean_IF_vs_empirical_contamination','passed':if_error<1e-6,'error':if_error},
   {'check':'mean_seeds_disjoint_from_primary_and_each_other','passed':not seeds.intersection(new) and len(new)==len(set(new)),'new_cells':len(new)}]
 result={'status':'PASS' if all(c['passed'] for c in checks) else 'FAIL','checks':checks}
 save_json(RESULTS/'MEAN_PREFLIGHT_VALIDATION.json',result);print(result)
 if result['status']!='PASS':raise RuntimeError('Mean followup preflight failed')

if __name__=='__main__':main()
