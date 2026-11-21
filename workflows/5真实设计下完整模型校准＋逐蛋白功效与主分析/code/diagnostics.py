from config import *
from model import factors
from itertools import product

def parameters(Y,d):
 Q=d['Q'];E=Y-Q@(Q.T@Y);b=d['a']@Y;df=d['df']
 products=np.column_stack([E[:,0]**2,E[:,1]**2,E[:,2]**2,E[:,0]*E[:,1],E[:,0]*E[:,2],E[:,1]*E[:,2]])
 c=products.sum(axis=0)/df;IF=np.c_[d['a'][:,None]*E,(products-products.mean(axis=0))/df]
 V=np.zeros((9,9))
 for one in d['onehots']:
  s=one@IF;s-=s.mean(axis=0);V+=s.T@s*len(s)/(len(s)-1)
 return np.r_[b,c],IF,V,E

def constraints_jacobian(theta,rho=(0.,0.,0.)):
 b=theta[:3];v=theta[3:6];off=theta[6:9].copy();D=np.zeros((3,6))
 for k,(i,j) in enumerate([(0,1),(0,2),(1,2)]):
  off[k]-=rho[k]*np.sqrt(v[i]*v[j]);D[k,i]=-.5*rho[k]*np.sqrt(v[j]/v[i]);D[k,j]=-.5*rho[k]*np.sqrt(v[i]/v[j]);D[k,3+k]=1
 c12,c13,c23=off;h=np.array([b[1]*c13-b[0]*c23,b[2]*c12-b[0]*c23])
 Jb=np.array([[-c23,c13,0],[-c23,0,c12]]);Jo=np.array([[0,b[1],-b[0]],[b[2],0,-b[0]]])
 return h,np.c_[Jb,Jo@D]

def wald(h,V):
 sd=np.sqrt(np.maximum(np.diag(V),0))
 if not np.isfinite(V).all() or np.min(sd)<=0:return np.nan,np.nan
 corr=V/np.outer(sd,sd)
 if np.linalg.eigvalsh(corr)[0]<=1e-9:return np.nan,np.nan
 q=max(0.,float(h@np.linalg.solve(V,h)));return q,float(np.exp(-q/2))

def summarize(Y,d,version):
 theta,IF,V,E=parameters(Y,d);V*=factors(d)[version];h,J=constraints_jacobian(theta);Q,p=wald(h,J@V@J.T)
 result={'Q':Q,'p_value':p,'H1':h[0],'H2':h[1],'residual_cov_positive':bool((theta[6:9]>0).all())}
 for j,label in enumerate(['N','S','T']):
  b=theta[j];v=theta[3+j];gradient=np.zeros(9);gradient[j]=1/np.sqrt(v);gradient[3+j]=-.5*b/v**1.5
  se=float(np.sqrt(max(0,gradient@V@gradient)));scaled=b/np.sqrt(v)
  result.update({label+'_b_per_global_P_SD':b,label+'_b_SE':float(np.sqrt(V[j,j])),label+'_residual_SD':float(np.sqrt(v)),
    label+'_b_over_residual_SD':scaled,label+'_scaled_b_CI95_low':scaled-1.959964*se,label+'_scaled_b_CI95_high':scaled+1.959964*se})
 cov_low=[]
 for k,(i,j,label) in enumerate([(0,1,'NS'),(0,2,'NT'),(1,2,'ST')]):
  c=theta[6+k];se=float(np.sqrt(V[6+k,6+k]));corr=c/np.sqrt(theta[3+i]*theta[3+j]);low=c-1.959964*se;cov_low.append(low)
  result.update({'c_'+label:c,'corr_residual_'+label:corr,'c_'+label+'_CI95_low':low,'c_'+label+'_CI95_high':c+1.959964*se})
 result['positive_cov_CI_all_three']=bool(min(cov_low)>0)
 for axis,label in enumerate(['N','S','T']):
  s=d['onehots'][axis]@(IF@J.T);energy=np.sum(s*s/np.maximum(np.diag(J@V@J.T),1e-100),axis=1)
  result[label+'_IF_effective_clusters']=float(energy.sum()**2/(energy@energy)) if energy@energy>0 else np.nan
  result[label+'_max_IF_energy_share']=float(energy.max()/energy.sum()) if energy.sum()>0 else np.nan
 return result,theta,IF,V,E

def covariance_profile(theta,V):
 rows=[]
 for rho in product([-.2,-.1,0.,.1,.2],repeat=3):
  h,J=constraints_jacobian(theta,rho);Q,p=wald(h,J@V@J.T)
  rows.append({'rho_NS':rho[0],'rho_NT':rho[1],'rho_ST':rho[2],'Q':Q,'p':p,'valid':bool(np.isfinite(p))})
 return rows
