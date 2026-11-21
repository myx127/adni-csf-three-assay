from config import *
from scipy.linalg import qr

def make(clinical,mask,extra=(),P='P_LOG',tech=('NULISA_PlateID','SOMASCAN_PlateId','TMT_plex'),additional_factors=('ptau_batch',),cluster_columns=None):
 d=clinical.loc[mask].copy();n=len(d)
 if n<5:return None
 cols=['age_month_precision','sex_code2']+list(extra);blocks=[np.ones((n,1))]
 for c in cols:
  x=d[c].to_numpy(float);sd=x.std();blocks.append(((x-x.mean())/(sd if sd>0 else 1))[:,None])
 for c in list(tech)+list(additional_factors):
  if c in d:blocks.append(pd.get_dummies(d[c].astype(str),drop_first=True,dtype=float).to_numpy())
 Z=np.column_stack(blocks);q,r,piv=qr(Z,mode='economic',pivoting=True);tol=max(Z.shape)*np.finfo(float).eps*abs(np.diag(r)).max()
 rank=int((abs(np.diag(r))>tol).sum());Qz=q[:,:rank];p=d[P].to_numpy(float);rP=p-Qz@(Qz.T@p);ss=rP@rP
 if not np.isfinite(ss) or ss<1e-10:return None
 Q=np.c_[Qz,rP/np.sqrt(ss)];onehots=[];labels=[];summary={}
 if cluster_columns is None:cluster_columns=tech
 for c in cluster_columns:
  lab,_=pd.factorize(d[c].astype(str),sort=True);g=int(lab.max())+1;one=np.eye(g)[lab].T;labels.append(lab);onehots.append(one)
  summary[c+'_G']=g
 age=d.age_month_precision.to_numpy();sd=age.std()
 return {'indices':d.index.to_numpy(),'keys':d.key.to_numpy(str),'n':n,'rank':Q.shape[1],'df':n-Q.shape[1],'P':p,'Z':Z,'Qz':Qz,'Q':Q,
  'a':rP/ss,'q_age':(age-age.mean())/(sd if sd>0 else 1),'groups':np.array(labels).T,'onehots':onehots,
  'summary':{'n':n,'rank_X':Q.shape[1],'residual_df':n-Q.shape[1],'P_residual_SS':float(ss),'max_hat_leverage':float((Q*Q).sum(axis=1).max()),**summary}}
