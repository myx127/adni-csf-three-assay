from config import *
from power_engine import setup,generate,expanded_parameters,evaluate
from diagnostics import parameters
from equivalence import actual_equivalence
from scipy.linalg import block_diag,lstsq
import pickle

def run():
 templates=pickle.load((PRIVATE/'POWER_TEMPLATES_RESTRICTED.pkl').open('rb'));checks=[];rng=np.random.default_rng(2026103000)
 def check(name,ok,detail=''):checks.append({'check':name,'passed':bool(ok),'detail':str(detail)})
 for acc,t in templates.items():
  d=t['d'];check('template_cov_reconstruction_'+acc,np.allclose(np.outer(t['lambda'],t['lambda'])+np.diag(t['psi']),t['C'],atol=1e-10))
  Y=generate(setup(t,1,'SLOPES3'),2,rng,.55);theta,V=expanded_parameters(Y,d);p,*_=evaluate(theta,V)
  for i in range(2):
   tt,_,vv,_=parameters(Y[i,0],d);check('k1_theta_V_reference_'+acc+'_'+str(i),np.allclose(theta[i],tt,atol=1e-10) and np.allclose(V[i],vv,atol=1e-10))
   for j,m in enumerate(MARGINS):check('vector_IUT_reference_'+acc+'_'+str(i)+'_'+str(m),np.allclose(p[i,j],actual_equivalence(Y[i,0],d,m)['p_joint_equivalence'],atol=1e-9))
 # Explicit stacked, redundant design SVD: independently checks common P and cluster corrections.
 for acc,t in list(templates.items())[:3]:
  d=t['d'];k=2;spec=setup(t,k,'SLOPES3');Y=generate(spec,1,rng,.55);theta,V=expanded_parameters(Y,d);y=Y[0].reshape(k*d['n'],3)
  Z=block_diag(*([d['Z']]*k));P=np.tile(d['P'],k);X=np.c_[Z,P];coef,_,rank,_=lstsq(X,y,cond=max(X.shape)*np.finfo(float).eps);E=y-X@coef
  cp,_,_,_=lstsq(Z,P,cond=max(Z.shape)*np.finfo(float).eps);rp=P-Z@cp;a=rp/(rp@rp);df=len(y)-rank
  products=np.c_[E[:,0]**2,E[:,1]**2,E[:,2]**2,E[:,0]*E[:,1],E[:,0]*E[:,2],E[:,1]*E[:,2]]
  tt=np.r_[coef[-1],products.sum(axis=0)/df];IF=np.c_[a[:,None]*E,(products-products.mean(axis=0))/df];vv=np.zeros((9,9))
  for axis in range(3):
   g=len(d['onehots'][axis]);labels=np.concatenate([d['groups'][:,axis]+block*g for block in range(k)]);s=np.vstack([IF[labels==j].sum(axis=0) for j in np.unique(labels)]);s-=s.mean(axis=0);vv+=s.T@s*len(s)/(len(s)-1)
  check('expanded_explicit_SVD_theta_'+acc,np.allclose(tt,theta[0],atol=1e-9));check('expanded_explicit_SUM3_V_'+acc,np.allclose(vv,V[0],atol=1e-9));check('expanded_common_P_rank_df_'+acc,df==spec['df'])
  # Expected slope trace computed using the explicit stacked projection and slope matrix.
  mats=[one.T*q[:,None] for one,(_,_,q) in zip(d['onehots'],spec['qrows'])];A=block_diag(*([np.column_stack(mats)]*k));ca,_,_,_=lstsq(X,A,cond=max(X.shape)*np.finfo(float).eps)
  expected=np.sqrt(1+.15**2*np.sum((A-X@ca)**2)/df);check('expanded_slope_normalization_'+acc,np.allclose(expected,spec['normalizer'],atol=1e-10))
 out=pd.DataFrame(checks);csvout(RESULTS/'POWER_ENGINE_INDEPENDENT_CHECKS.csv',out);print({'checks':len(out),'passed':int(out.passed.sum())})
 if not out.passed.all():print(out.loc[~out.passed].to_string(index=False));raise SystemExit(1)

if __name__=='__main__':run()
