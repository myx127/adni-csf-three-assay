from config import *
from scipy.linalg import svd
from scipy.stats import norm
D=np.array([[-1.,1.,0.],[-1.,0.,1.]])
def basis(Z):
 u,s,_=svd(Z,full_matrices=False,check_finite=True);rank=int((s>s[0]*max(Z.shape)*np.finfo(float).eps).sum()) if len(s) and s[0]>0 else 0
 return u[:,:rank]
def rcs(x,knots):
 x=np.asarray(x,float);k=np.asarray(knots,float);parts=[x]
 for a in k[:-2]:
  h=np.maximum(x-a,0)**3-np.maximum(x-k[-2],0)**3*(k[-1]-a)/(k[-1]-k[-2])+np.maximum(x-k[-1],0)**3*(k[-2]-a)/(k[-1]-k[-2])
  parts.append(h/(k[-1]-k[0])**2)
 return np.column_stack(parts)
def nuisance(f,model,knots):
 blocks=[np.ones((len(f),1))];names=['intercept']
 for col in ['age_month_precision','sex_code2','education','apoe4_carrier']:
  v=f[col].to_numpy(float);blocks.append(((v-v.mean())/(v.std() or 1))[:,None]);names.append(col)
 for col in ['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex','ptau_batch']:
  dum=pd.get_dummies(f[col].astype(str),drop_first=True,dtype=float);blocks.append(dum.to_numpy());names.extend(col+'='+v for v in dum.columns)
 p=np.empty((len(f),0));pnames=[]
 if model.endswith('LOG_P'):p=f[['P_LOG']].to_numpy();pnames=['P_LOG']
 elif model.endswith('RAW_P'):p=f[['P_RAW']].to_numpy();pnames=['P_RAW']
 elif model.endswith('RCS_P'):p=rcs(f.P_LOG.to_numpy(),knots);pnames=['P_LOG','P_RCS1','P_RCS2']
 for j in range(p.shape[1]):
  v=p[:,j];blocks.append(((v-v.mean())/(v.std() or 1))[:,None]);names.append(pnames[j])
 bgcols=['BG_SOMA','BG_TMT'] if model.startswith('BG2') else ['BG_SOMA','BG_TMT','BG_NULISA'] if model.startswith('BG3') else []
 for col in bgcols:
  v=f[col].to_numpy(float);blocks.append(((v-v.mean())/(v.std() or 1))[:,None]);names.append(col)
 return np.column_stack(blocks),names,p,pnames
def wald2(h,V):
 a=V[:,0,0];b=V[:,1,1];c=V[:,0,1]
 with np.errstate(divide='ignore',invalid='ignore'):
  rho=c/np.sqrt(a*b);x=h[:,0]/np.sqrt(a);y=h[:,1]/np.sqrt(b);Q=(x*x+y*y-2*rho*x*y)/(1-rho*rho)
 good=(a>0)&(b>0)&np.isfinite(Q)&((1-np.abs(rho))>1e-9)
 return np.exp(-np.maximum(Q,0)/2),Q,good
def moment_stats(b,U,t,full=False):
 sigma=t['sigma'];L=t['L'];h=(b*sigma[None,:])@L.T
 scaled=U*sigma[None,None,:]+b[:,None,:]*t['Uv'][None,:,:]/(2*sigma[None,None,:])
 Uh=np.einsum('ap,bgp->bga',L,scaled,optimize=True);V=np.einsum('bga,bgc->bac',Uh,Uh,optimize=True);p,Q,valid=wald2(h,V)
 out={'h':h,'V':V,'p':p,'Q':Q,'valid':valid}
 if full:out.update({'b':b,'rawV':np.einsum('bgp,bgq->bpq',U,U,optimize=True),'gamma':b*sigma[None,:],'scaledV':np.einsum('bgp,bgq->bpq',scaled,scaled,optimize=True),'scores':U,'scaled_scores':scaled})
 return out
def stats(y,t):
 y=np.asarray(y,float)
 if y.ndim==1:y=y[:,None]
 return moment_stats((t['A']@y).T,np.einsum('gpn,nb->bgp',t['K'],y,optimize=True),t,full=True)
def wilson(k,n):
 z=norm.ppf(.975);p=k/n;den=1+z*z/n;c=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 return max(0.,c-half),min(1.,c+half)

