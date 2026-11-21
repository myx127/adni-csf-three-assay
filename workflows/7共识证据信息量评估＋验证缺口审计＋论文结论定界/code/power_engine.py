from config import *
from scipy.stats import norm

def setup(template,k,scenario):
 d=template['d'];qrows=[];mats=[];slopes=0. if scenario=='IID' else .15
 for axis,one in enumerate(d['onehots']):
  lab=d['groups'][:,axis];age=d['q_age'];q=age-(one@age/one.sum(axis=1))[lab];qrows.append((lab,len(one),q));mats.append(one.T*q[:,None])
 A=np.column_stack(mats);AZ=A-d['Qz']@(d['Qz'].T@A);rP=d['a']/(d['a']@d['a']);ss=rP@rP;ploss=np.sum((rP@A)**2)/ss
 df=k*(d['n']-d['Qz'].shape[1])-1;projected=(k*np.sum(AZ*AZ)-ploss)/df
 return {'template':template,'d':d,'k':k,'scenario':scenario,'slope_SD':slopes,'qrows':qrows,'normalizer':np.sqrt(1+slopes**2*projected),'df':df,'rP':rP}

def generate(spec,B,rng,gamma,boundary=None):
 d=spec['d'];k=spec['k'];heavy=spec['scenario']=='T5_SLOPES3'
 def draw(shape):return rng.standard_t(5,size=shape)*np.sqrt(3/5) if heavy else rng.normal(size=shape)
 W=draw((B,k,d['n'],4))
 if spec['slope_SD']:
  for lab,G,q in spec['qrows']:W+=draw((B,k,G,4))[:,:,lab,:]*spec['slope_SD']*q[None,None,:,None]
 W/=spec['normalizer'];lam=spec['template']['lambda'];psi=spec['template']['psi']
 Y=lam[None,None,None,:]*(gamma*d['P'][None,None,:,None]+W[:,:,:,0,None])+np.sqrt(psi)[None,None,None,:]*W[:,:,:,1:]
 if boundary:
  margin,target,sign=boundary;Y[:,:,:,target]+=sign*margin*np.sqrt(spec['template']['C'][target,target])*d['P'][None,None,:]
 return Y

def expanded_parameters(Y,d):
 B,k,n,_=Y.shape;flat=Y.transpose(2,0,1,3).reshape(n,-1);R=flat-d['Qz']@(d['Qz'].T@flat)
 b=(d['a']@flat).reshape(B,k,3).mean(axis=1);rP=d['a']/(d['a']@d['a']);E=R.reshape(n,B,k,3).transpose(1,2,0,3)-rP[None,None,:,None]*b[:,None,None,:]
 df=k*(n-d['Qz'].shape[1])-1
 prod=np.stack([E[:,:,:,0]**2,E[:,:,:,1]**2,E[:,:,:,2]**2,E[:,:,:,0]*E[:,:,:,1],E[:,:,:,0]*E[:,:,:,2],E[:,:,:,1]*E[:,:,:,2]],axis=3)
 c=prod.sum(axis=(1,2))/df;IF=np.concatenate([E*d['a'][None,None,:,None]/k,(prod-prod.mean(axis=(1,2),keepdims=True))/df],axis=3)
 flat=IF.transpose(2,0,1,3).reshape(n,-1);V=np.zeros((B,9,9))
 for one in d['onehots']:
  G=len(one);s=(one@flat).reshape(G,B,k,9).transpose(1,2,0,3).reshape(B,k*G,9);s-=s.mean(axis=1,keepdims=True)
  V+=np.einsum('bgi,bgj->bij',s,s)*(k*G)/(k*G-1)
 return np.c_[b,c],V

def evaluate(theta,V):
 B=len(theta);b=theta[:,:3];v=theta[:,3:6];c=theta[:,6:9];c12,c13,c23=c.T
 A=np.c_[b[:,1]*c13-b[:,0]*c23,b[:,2]*c12-b[:,0]*c23];JA=np.zeros((B,2,9))
 JA[:,0,0]=-c23;JA[:,0,1]=c13;JA[:,0,7]=b[:,1];JA[:,0,8]=-b[:,0]
 JA[:,1,0]=-c23;JA[:,1,2]=c12;JA[:,1,6]=b[:,2];JA[:,1,8]=-b[:,0]
 sv=np.sqrt(v);den=np.c_[c13*sv[:,1],c12*sv[:,2]];JB=np.zeros((B,2,9))
 JB[:,0,7]=sv[:,1];JB[:,0,4]=.5*c13/sv[:,1];JB[:,1,6]=sv[:,2];JB[:,1,5]=.5*c12/sv[:,2]
 valid=np.isfinite(theta).all(axis=1)&np.isfinite(V).all(axis=(1,2))&(v>0).all(axis=1)
 positive=(c>0).all(axis=1);commonvar=np.c_[c12*c13/c23,c12*c23/c13,c13*c23/c12]
 struct=positive&(v-commonvar>=-1e-12).all(axis=1)
 covstrong=(c-1.959964*np.sqrt(np.maximum(0,np.diagonal(V,axis1=1,axis2=2)[:,6:9]))>0).all(axis=1)
 vb=np.einsum('bji,bik,bjk->bj',JB,V,JB);bounded=((den>0)&(den*den>norm.isf(.05)**2*vb)).all(axis=1)
 eligible=valid&struct&covstrong&bounded;pvalues=[]
 for margin in MARGINS:
  ps=[]
  for j in range(2):
   for sign in [-1,1]:
    g=margin*JB[:,j,:]+sign*JA[:,j,:];vv=np.einsum('bi,bij,bj->b',g,V,g);distance=margin*den[:,j]+sign*A[:,j]
    ps.append(norm.sf(distance/np.sqrt(np.maximum(vv,1e-100))))
  p=np.max(ps,axis=0);p[(den<=0).any(axis=1)]=1.;p[~valid]=np.nan;pvalues.append(p)
 return np.array(pvalues).T,eligible,valid,struct,covstrong,bounded
