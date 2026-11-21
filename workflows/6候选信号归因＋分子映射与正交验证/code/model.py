from __future__ import annotations
from config import *
from scipy.linalg import qr
from scipy.stats import norm

def make_design(clinical,mask,variant='MAIN'):
 d=clinical.loc[mask].copy();pcol='P_RAW' if variant=='RAW_P' else 'P_AMY' if variant=='AMY' else 'P_LOG'
 cov=['age_month_precision','sex_code2']+(['education','apoe4_carrier'] if variant=='FULLCOV' else [])
 labels=['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex']
 nuisance=[np.ones((len(d),1))];names=['intercept']
 for c in cov:
  x=d[c].to_numpy(float);sd=x.std();nuisance.append(((x-x.mean())/(sd if sd>0 else 1))[:,None]);names.append(c)
 for c in labels+(['PET_tracer'] if variant=='AMY' else ['ptau_batch']):
  dum=pd.get_dummies(d[c].astype(str),prefix=c,drop_first=True,dtype=float)
  if dum.shape[1]:nuisance.append(dum.to_numpy());names.extend(dum.columns)
 z=np.column_stack(nuisance);q,r,piv=qr(z,mode='economic',pivoting=True)
 tol=np.max(np.abs(np.diag(r)))*max(z.shape)*np.finfo(float).eps
 rank=int((np.abs(np.diag(r))>tol).sum());qz=q[:,:rank]
 p=d[pcol].to_numpy(float);resP=p-qz@(qz.T@p);sxx=float(resP@resP)
 if not np.isfinite(sxx) or sxx<=1e-10:return None
 qfull=np.c_[qz,resP/np.sqrt(sxx)];df=len(d)-qfull.shape[1]
 groups=[];onehots=[];summ={}
 partial=(resP*resP)/sxx;lever=np.sum(qfull*qfull,axis=1)
 for col,letter in zip(labels,['N','S','T']):
  lab,_=pd.factorize(d[col].astype(str),sort=True);G=int(lab.max())+1;groups.append(lab)
  onehot=np.eye(G)[lab].T;onehots.append(onehot)
  part=onehot@partial;cnt=onehot.sum(axis=1)
  summ.update({letter+'_G':G,letter+'_effective_size_clusters':float(cnt.sum()**2/(cnt*cnt).sum()),
    letter+'_P_partial_effective_clusters':float(1/(part*part).sum()),letter+'_P_max_cluster_leverage':float(part.max())})
 return {'indices':d.index.to_numpy(),'keys':d.key.to_numpy(str),'P':p,'Z':z,'Qz':qz,'Q':qfull,'a':resP/sxx,'df':df,'rank':qfull.shape[1],
   'n':len(d),'groups':np.array(groups).T,'onehots':onehots,'q_age':(d.age_month_precision.to_numpy()-d.age_month_precision.mean())/d.age_month_precision.std(ddof=0),
   'summary':{'n':len(d),'rank_X':qfull.shape[1],'residual_df':df,'P_residual_SS':sxx,'P_fraction_remaining':sxx/float(((p-p.mean())**2).sum()),
     'max_hat_leverage':float(lever.max()),'max_P_row_leverage':float(partial.max()),'Z_nominal_columns':z.shape[1],**summ}}

def batch_components(Y,design):
 """Y: B x n x 3; c uses common residual df; IF includes fitted b and residual c."""
 if Y.ndim==2:Y=Y[None,:,:]
 B,n,_=Y.shape;Q=design['Q'];df=design['df'];a=design['a']
 flat=Y.transpose(1,0,2).reshape(n,-1);E=(flat-Q@(Q.T@flat)).reshape(n,B,3).transpose(1,0,2)
 b=np.einsum('n,bnj->bj',a,Y)
 products=np.stack([E[:,:,0]*E[:,:,1],E[:,:,0]*E[:,:,2],E[:,:,1]*E[:,:,2]],axis=2)
 means=products.mean(axis=1);c=products.sum(axis=1)/df;c12,c13,c23=c.T
 h=np.c_[b[:,1]*c13-b[:,0]*c23,b[:,2]*c12-b[:,0]*c23]
 be=a[None,:,None]*E;ci=(products-means[:,None,:])/df
 U=np.stack([-c23[:,None]*be[:,:,0]+c13[:,None]*be[:,:,1]+b[:,1,None]*ci[:,:,1]-b[:,0,None]*ci[:,:,2],
  -c23[:,None]*be[:,:,0]+c12[:,None]*be[:,:,2]+b[:,2,None]*ci[:,:,0]-b[:,0,None]*ci[:,:,2]],axis=2)
 V=np.zeros((B,2,2));uflat=U.transpose(1,0,2).reshape(n,-1)
 for one in design['onehots']:
  G=one.shape[0];s=(one@uflat).reshape(G,B,2).transpose(1,0,2);s-=s.mean(axis=1,keepdims=True)
  if G<=1:V[:]=np.nan;continue
  V+=np.einsum('bgi,bgj->bij',s,s)*G/(G-1)
 det=V[:,0,0]*V[:,1,1]-V[:,0,1]**2
 ratio=det/(V[:,0,0]*V[:,1,1]);ok=np.isfinite(V).all(axis=(1,2))&(V[:,0,0]>0)&(V[:,1,1]>0)&(ratio>1e-9)
 statistic=np.full(B,np.nan)
 statistic[ok]=(h[ok,0]**2*V[ok,1,1]+h[ok,1]**2*V[ok,0,0]-2*h[ok,0]*h[ok,1]*V[ok,0,1])/det[ok]
 return {'Q0':np.maximum(statistic,0),'h':h,'V0':V,'b':b,'c_off':c,'E':E,'U':U,'ok':ok}

def factors(design):
 f=(design['n']-1)/design['df'];return np.array([1.,f,2*f])

def wilson(k,n):
 z=norm.ppf(.975);p=k/n;den=1+z*z/n;center=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 return float(max(0,center-half)),float(min(1,center+half))

def bh(p,dependence=False):
 p=np.asarray(p,float);m=len(p);idx=np.argsort(p);adjust=np.minimum.accumulate((p[idx]*m/np.arange(1,m+1)*(np.sum(1/np.arange(1,m+1)) if dependence else 1))[::-1])[::-1]
 out=np.ones(m);out[idx]=np.minimum(1,adjust);return out

def holm(p):
 p=np.asarray(p,float);m=len(p);idx=np.argsort(p);vals=np.maximum.accumulate(p[idx]*(m-np.arange(m)));out=np.ones(m);out[idx]=np.minimum(1,vals);return out
