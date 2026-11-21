import numpy as np
import pandas as pd
from scipy.linalg import svd
from scipy.stats import norm,chi2

D=np.array([[-1.,1.,0.],[-1.,0.,1.]])
def basis(Z):
 u,s,_=svd(Z,full_matrices=False,check_finite=True)
 rank=int((s>s[0]*max(Z.shape)*np.finfo(float).eps).sum()) if len(s) and s[0]>0 else 0
 return u[:,:rank]
def nuisance(frame,ptau='P_LOG',include_ptau=True,tracer=False,base_score=None,modifier=None):
 blocks=[np.ones((len(frame),1))]
 names=['intercept']
 for c in ['age_month_precision','sex_code2','education','apoe4_carrier']+([ptau] if include_ptau else []):
  v=frame[c].to_numpy(float);sd=v.std();blocks.append(((v-v.mean())/(sd if sd else 1))[:,None]);names.append(c)
 for c in ['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex','ptau_batch']+(['PET_tracer'] if tracer else []):
  dum=pd.get_dummies(frame[c].astype(str),drop_first=True,dtype=float);blocks.append(dum.to_numpy());names.extend(c+'='+v for v in dum.columns)
 if base_score is not None:
  v=np.asarray(base_score,float);blocks.append((v-v.mean())[:,None]);names.append('baseline_MEM')
 Z=np.column_stack(blocks)
 if modifier is not None:
  Z=np.c_[Z,np.asarray(modifier)[:,None]*Z];names=names+['modifier*'+v for v in names]
 return Z,names

def joint_stats(y,template):
 y=np.asarray(y,float)
 if y.ndim==1:y=y[:,None]
 b=template['A']@y
 score=np.einsum('gjn,nb->bgj',template['K'],y,optimize=True)
 V=np.einsum('bgj,bgk->bjk',score,score,optimize=True)
 scaled=b.T*template['sigma'][None,:]
 ss=score*template['sigma'][None,None,:]+b.T[:,None,:]*template['Uv'][None,:,:]/(2*template['sigma'][None,None,:])
 Vs=np.einsum('bgj,bgk->bjk',ss,ss,optimize=True)
 se=np.sqrt(np.maximum(np.diagonal(V,axis1=1,axis2=2),0));z=b.T/se
 pplus=norm.sf(z).max(axis=1);pminus=norm.cdf(z).max(axis=1)
 same=np.minimum(1.,2*np.minimum(pplus,pminus))
 diff=scaled@D.T;Vd=np.einsum('ij,bjk,lk->bil',D,Vs,D)
 Qd,pd_,validd=wald_batch(diff,Vd)
 Qa,pa,valida=wald_batch(b.T,V)
 valid=np.isfinite(se).all(axis=1)&(se>0).all(axis=1)&validd
 return {'b':b.T,'Vraw':V,'scaled':scaled,'Vscaled':Vs,'diff':diff,'Vdiff':Vd,'p_same':same,'p_diff':pd_,'p_any':pa,'Qdiff':Qd,'Qany':Qa,'p_plus':pplus,'p_minus':pminus,'valid':valid,'z':z}

def wald_batch(b,V):
 sd=np.sqrt(np.maximum(np.diagonal(V,axis1=1,axis2=2),0))
 with np.errstate(divide='ignore',invalid='ignore'):
  C=V/(sd[:,:,None]*sd[:,None,:])
 valid=np.isfinite(C).all(axis=(1,2))&(sd>0).all(axis=1)
 if valid.any():valid[valid]&=np.linalg.eigvalsh(C[valid]).min(axis=1)>1e-9
 Q=np.full(len(b),np.nan)
 if valid.any():Q[valid]=np.einsum('bi,bi->b',b[valid],np.linalg.solve(V[valid],b[valid]))
 return Q,chi2.sf(np.maximum(Q,0),b.shape[1]),valid

def stats_from_moments(b,score,t):
 # Only the prespecified marginal IUT and 2df paired contrast determine validity.
 rawvar=np.einsum('bgj,bgj->bj',score,score,optimize=True)
 with np.errstate(divide='ignore',invalid='ignore'):
  zz=b/np.sqrt(rawvar)
 same=np.minimum(1.,2*np.minimum(norm.sf(zz).max(axis=1),norm.cdf(zz).max(axis=1)))
 sigma=t['sigma']; c=t['Uv']/(2*sigma[None,:])
 n=score[:,:,0]*sigma[0]+b[:,None,0]*c[None,:,0]
 s=score[:,:,1]*sigma[1]+b[:,None,1]*c[None,:,1]-n
 u=score[:,:,2]*sigma[2]+b[:,None,2]*c[None,:,2]-n
 v1=np.einsum('bg,bg->b',s,s);v2=np.einsum('bg,bg->b',u,u);v12=np.einsum('bg,bg->b',s,u)
 d1=b[:,1]*sigma[1]-b[:,0]*sigma[0];d2=b[:,2]*sigma[2]-b[:,0]*sigma[0]
 with np.errstate(divide='ignore',invalid='ignore'):
  rho=v12/np.sqrt(v1*v2);x=d1/np.sqrt(v1);y=d2/np.sqrt(v2);stat=(x*x+y*y-2*rho*x*y)/(1-rho*rho)
 valid=np.isfinite(zz).all(axis=1)&(rawvar>0).all(axis=1)&np.isfinite(stat)&(v1>0)&(v2>0)&((1-np.abs(rho))>1e-9)
 diff=np.exp(-np.maximum(stat,0)/2)
 return same,diff,None,valid,zz

def wilson(k,n):
 z=norm.ppf(.975);p=k/n;den=1+z*z/n;center=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 return max(0.,center-half),min(1.,center+half)
