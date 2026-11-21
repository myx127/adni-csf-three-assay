from config import *
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
 valid=np.isfinite(se).all(axis=1)&(se>0).all(axis=1)&validd&valida
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
 V=np.einsum('bgj,bgk->bjk',score,score,optimize=True)
 scaled=b*t['sigma'][None,:]
 ss=score*t['sigma'][None,None,:]+b[:,None,:]*t['Uv'][None,:,:]/(2*t['sigma'][None,None,:])
 Vs=np.einsum('bgj,bgk->bjk',ss,ss,optimize=True)
 se=np.sqrt(np.maximum(np.diagonal(V,axis1=1,axis2=2),0));z=b/se
 p_same=np.minimum(1.,2*np.minimum(norm.sf(z).max(axis=1),norm.cdf(z).max(axis=1)))
 bd=scaled@D.T;Vd=np.einsum('ij,bjk,lk->bil',D,Vs,D)
 _,p_diff,vd=wald_batch(bd,Vd);_,p_any,va=wald_batch(b,V)
 valid=vd&va&np.isfinite(z).all(axis=1)
 return p_same,p_diff,p_any,valid,z

def wilson(k,n):
 z=norm.ppf(.975);p=k/n;den=1+z*z/n;center=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 return max(0.,center-half),min(1.,center+half)
