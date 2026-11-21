from __future__ import annotations
import numpy as np
from scipy.stats import f
from scipy.special import gammaincc
from reference_stage2_models import invariance_constraints,invariance_moments,multiway_covariance
from cal_config import METHODS

def failed(status):return {'p_value':np.nan,'statistic':np.nan,'status':status,'bootstrap_tail_count':-1}

def inverse_metric(covariance):
 covariance=(covariance+covariance.T)/2
 scale=np.sqrt(np.maximum(np.diag(covariance),0))
 if not np.isfinite(covariance).all() or np.any(scale<=1e-15):return None,'SINGULAR_COVARIANCE'
 corr=covariance/np.outer(scale,scale);eigen=np.linalg.eigvalsh(corr)
 if eigen[0]<-1e-8:return None,'INDEFINITE_COVARIANCE'
 if eigen[0]<=1e-9:return None,'SINGULAR_COVARIANCE'
 return np.linalg.inv(corr)/np.outer(scale,scale),'OK'

def test(h,covariance,G=None):
 inv,status=inverse_metric(covariance)
 if inv is None:return failed(status)
 Q=max(0.,float(h@inv@h));r=len(h)
 if G is None:p=float(gammaincc(r/2,Q/2))
 else:
  if G<=r:return failed('INSUFFICIENT_CLUSTERS_FOR_F')
  p=float(f.sf((G-r)/(r*(G-1))*Q,r,G-r))
 return {'p_value':p,'statistic':Q,'status':'OK','bootstrap_tail_count':-1}

def sums_by_cluster(U,labels):
 _,inverse=np.unique(labels,return_inverse=True);G=int(inverse.max())+1
 sums=np.zeros((G,U.shape[1]));np.add.at(sums,inverse,U)
 sums-=sums.mean(axis=0)
 return sums

def cr1_covariance(sums):
 G=len(sums)
 return sums.T@sums*G/(G-1) if G>1 else np.full((sums.shape[1],)*2,np.nan)

def wild_t(h,sums,B,rng,return_draws=False):
 G,r=sums.shape;cov=cr1_covariance(sums);inv,status=inverse_metric(cov)
 if inv is None:return (failed(status),None) if return_draws else failed(status)
 Q=float(h@inv@h);Ainv=inv*G/(G-1)
 signs=(2*rng.integers(0,2,size=(B,G),dtype=np.int8)-1).astype(float)
 v=signs@sums;t=np.einsum('bi,ij,bj->b',v,Ainv,v)
 denominator=1-t/G
 bootQ=np.divide(t*(G-1)/G,denominator,out=np.full(B,np.inf),where=denominator>1e-12)
 count=int(np.sum(bootQ>=Q));result={'p_value':(1+count)/(B+1),'statistic':Q,'status':'OK','bootstrap_tail_count':count}
 if np.sum(denominator<=1e-12)>max(1,.01*B):result=failed('EXCESS_SINGULAR_WILD_DRAWS')
 return (result,{'signs':signs,'bootQ':bootQ,'sums':sums,'Q':Q}) if return_draws else result

def constraint_batch(m,weights,kappa=0.):
 # Full nonlinear re-estimation of all eight observable restrictions.
 m=m.copy();m[:,:,3]-=kappa
 pooled=np.einsum('bg,bgk->bk',weights,m[:,:,3:]);h=np.zeros((len(m),8))
 for g in [1,2]:
  d=m[:,g,:3]-m[:,0,:3];r=2*(g-1)
  h[:,r]=d[:,1]*pooled[:,1]-d[:,0]*pooled[:,2]
  h[:,r+1]=d[:,2]*pooled[:,0]-d[:,0]*pooled[:,2]
  a,b=m[:,g,3:],m[:,0,3:];r=4+2*(g-1)
  h[:,r]=a[:,0]*b[:,1]-a[:,1]*b[:,0]
  h[:,r+1]=a[:,0]*b[:,2]-a[:,2]*b[:,0]
 return h

def pigeonhole_moments(Y,groups,batches,B,rng,return_weights=False):
 W=np.ones((B,len(Y)))
 for j in range(3):
  weights=rng.exponential(size=(B,int(batches[:,j].max())+1));W*=weights[:,batches[:,j]]
 moments=np.zeros((B,3,6));mass=np.zeros((B,3))
 for g in range(3):
  mask=groups==g;y=Y[mask];wg=W[:,mask];den=wg.sum(axis=1);mass[:,g]=den
  features=np.c_[y,y[:,0]*y[:,1],y[:,0]*y[:,2],y[:,1]*y[:,2]]
  avg=(wg@features)/den[:,None];moments[:,g,:3]=avg[:,:3]
  moments[:,g,3]=avg[:,3]-avg[:,0]*avg[:,1]
  moments[:,g,4]=avg[:,4]-avg[:,0]*avg[:,2]
  moments[:,g,5]=avg[:,5]-avg[:,1]*avg[:,2]
 result=(moments,mass/mass.sum(axis=1)[:,None])
 return (*result,W) if return_weights else result

def bootstrap_tail(h,inv,boot_h):
 if inv is None:return failed('INVALID_BOOTSTRAP_METRIC')
 Q=float(h@inv@h);deviation=boot_h-h
 Qstar=np.einsum('bi,ij,bj->b',deviation,inv,deviation)
 if not np.isfinite(Qstar).all():return failed('NONFINITE_BOOTSTRAP')
 count=int(np.sum(Qstar>=Q));B=len(Qstar)
 return {'p_value':(1+count)/(B+1),'statistic':Q,'status':'OK','bootstrap_tail_count':count}

def prepare(data):
 Y,groups,batches=data['Y'],data['groups'],data['batches'];counts=np.bincount(groups,minlength=3)
 if len(Y)<30 or counts.min()<10:return None
 moment,contributions,counts=invariance_moments(Y,groups)
 return {'moment':moment,'contributions':contributions,'counts':counts,'Y':Y,'groups':groups,'batches':batches}

def at_kappa(state,kappa,axis=2):
 m=state['moment'].copy();m[:,3]-=kappa
 h,jac=invariance_constraints(m,state['counts']/len(state['Y']))
 U=state['contributions']@jac.T
 sums=[sums_by_cluster(U,state['batches'][:,j]) for j in range(3)]
 V=[cr1_covariance(s) for s in sums];Vsum=V[0]+V[1]+V[2]
 return h,U,sums,V,Vsum

def infer(data,axis,B,bootstrap_rng,kappa=0.,need_cache=False):
 state=prepare(data)
 if state is None:
  return {m:failed('INSUFFICIENT_COMPLETE_GROUPS') for m in METHODS},None
 h,U,sums,V,Vsum=at_kappa(state,kappa,axis)
 G=len(sums[axis]);minimum_G=min(len(s) for s in sums)
 results={'IID_CHI2':test(h,U.T@U),'CR1_CHI2':test(h,V[axis]),'CR1_HT_F':test(h,V[axis],G),
          'WILD_SCORE_T':wild_t(h,sums[axis],B,bootstrap_rng),'SUM3_CHI2':test(h,Vsum),
          'SUM3_MIN_F':test(h,Vsum,minimum_G),'CGM3_CHI2':test(h,multiway_covariance(U,state['batches']))}
 boot_m,boot_w=pigeonhole_moments(state['Y'],state['groups'],state['batches'],B,bootstrap_rng)
 inv,status=inverse_metric(Vsum)
 results['PIGEONHOLE_CENTERED']=bootstrap_tail(h,inv,constraint_batch(boot_m,boot_w,kappa)) if inv is not None else failed(status)
 cache={'state':state,'boot_m':boot_m,'boot_w':boot_w} if need_cache else None
 return results,cache

def sensitivity_profile(data,cache):
 if cache is None:return {'status':'INSUFFICIENT_COMPLETE_GROUPS'}
 grid=np.linspace(-.4,.4,9)*.48;values={'CR1_HT_F':[],'SUM3_CHI2':[]}
 for kappa in grid:
  h,U,sums,V,Vsum=at_kappa(cache['state'],float(kappa))
  values['CR1_HT_F'].append(test(h,V[2],len(sums[2]))['p_value'])
  values['SUM3_CHI2'].append(test(h,Vsum)['p_value'])
 out={'status':'OK'}
 true=data['true_kappa12']
 for method,p in values.items():
  p=np.array(p);all_valid=np.isfinite(p).all()
  out[method+'_zero_p']=float(p[4]);out[method+'_grid_min_p']=float(np.min(p)) if all_valid else np.nan
  out[method+'_grid_max_p']=float(np.max(p)) if all_valid else np.nan
  out[method+'_grid_all_valid']=bool(all_valid)
  out[method+'_classification_flips']=bool(np.any(p<=.05) and np.any(p>.05)) if all_valid else False
  on_grid=true is not None and np.min(np.abs(grid-true))<1e-10
  out[method+'_true_fixed_kappa_on_grid']=on_grid
  out[method+'_true_fixed_kappa_p']=float(p[np.argmin(abs(grid-true))]) if on_grid else np.nan
 if 'oracle_quadratic_Y' in data:
  modified={**data,'Y':data['oracle_quadratic_Y']};oracle=prepare(modified)
  if oracle is not None:
   h,U,sums,V,Vsum=at_kappa(oracle,0.)
   out['ORACLE_REMOVE_QUADRATIC_CR1_HT_F_p']=test(h,V[2],len(sums[2]))['p_value']
 return out
