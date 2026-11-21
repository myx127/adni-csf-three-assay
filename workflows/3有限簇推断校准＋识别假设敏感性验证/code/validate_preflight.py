from __future__ import annotations
import os
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[key]='1'
import numpy as np
from scipy.stats import f
from cal_config import *
from cal_inference import *
from cal_data import generate
from reference_stage2_models import invariance_constraints

def main():
 checks=[]
 def add(name,passed,detail):checks.append({'check':name,'passed':bool(passed),'detail':detail})
 rng=np.random.default_rng(271990)
 h=rng.normal(size=8);V=np.eye(8);G=19
 result=test(h,V,G);expected=f.sf((G-8)/(8*(G-1))*float(h@h),8,G-8)
 add('Hotelling_type_F_mapping',abs(result['p_value']-expected)<1e-14,expected)
 S=rng.normal(size=(19,8));S-=S.mean(axis=0)
 result,draws=wild_t(h,S,23,np.random.default_rng(99),return_draws=True)
 direct=[]
 for signs in draws['signs']:
  stars=signs[:,None]*S;v=stars.sum(axis=0);centered=stars-stars.mean(axis=0)
  covariance=centered.T@centered*19/18
  direct.append(float(v@np.linalg.solve(covariance,v)))
 error=float(np.max(abs(np.array(direct)-draws['bootQ'])))
 add('wild_fast_formula_vs_explicit_studentization',error<1e-8,error)
 count=int(np.sum(np.array(direct)>=draws['Q']))
 add('wild_plus_one_and_denominator',result['p_value']==(1+count)/24,{'tail_count':count,'B':23})
 n=120;groups=np.arange(n)%3;batches=np.c_[np.arange(n)%6,np.arange(n)%8,np.arange(n)%10];Y=rng.normal(size=(n,3))
 moments,weights,W=pigeonhole_moments(Y,groups,batches,5,rng,return_weights=True)
 maximum=0.
 for b in range(5):
  for g in range(3):
   ids=groups==g;y=Y[ids];w=W[b,ids];mean=(w[:,None]*y).sum(axis=0)/w.sum()
   cov=sum(wi*np.outer(yi-mean,yi-mean) for wi,yi in zip(w,y))/w.sum()
   direct=np.r_[mean,cov[0,1],cov[0,2],cov[1,2]]
   maximum=max(maximum,float(np.max(abs(direct-moments[b,g]))))
 add('full_weighted_moments_vs_explicit_group_loops',maximum<1e-12,maximum)
 maxconstraint=0.
 for kappa in [-.1,0.,.1]:
  batch=constraint_batch(moments,weights,kappa)
  for b in range(5):
   m=moments[b].copy();m[:,3]-=kappa
   reference=invariance_constraints(m,weights[b])[0]
   maxconstraint=max(maxconstraint,float(np.max(abs(reference-batch[b]))))
 add('batched_nonlinear_constraints_vs_frozen_reference',maxconstraint<1e-12,maxconstraint)
 cells=grid();seeds=[c['seed'] for c in cells]
 add('all_phase_cell_seeds_unique',len(seeds)==len(set(seeds)),len(seeds))
 add('plan_counts_74000_592000',sum(c['repetitions'] for c in cells)==74000,{'independent_repeats':74000,'method_rows':592000,'cells':len(cells)})
 sd={p:([30]*19+[1]) if p=='unused' else [] for p in []}
 # Independent covariance/mean corruption can be large but is not a reason to hide failures.
 indefinite=np.eye(8);indefinite[0,1]=indefinite[1,0]=2
 add('negative_covariance_is_explicit_failure',test(h,indefinite)['status']=='INDEFINITE_COVARIANCE',test(h,indefinite)['status'])
 result={'status':'PASS' if all(c['passed'] for c in checks) else 'FAIL','checks':checks,'passed':sum(c['passed'] for c in checks),'total':len(checks)}
 save_json(RESULTS/'PREFLIGHT_VALIDATION.json',result)
 print({k:v for k,v in result.items() if k!='checks'})
 if result['status']!='PASS':raise RuntimeError('Preflight validation failed')

if __name__=='__main__':main()
