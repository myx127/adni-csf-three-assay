from config import *
from model import batch_components
from diagnostics import parameters,constraints_jacobian
from scipy.stats import norm

MARGINS=[.2,.35,.5]

def equivalence_batch(Y,d,margin):
 out=batch_components(Y,d);E=out['E'];b=out['b'];df=d['df'];B,n,_=E.shape
 variances=np.sum(E*E,axis=1)/df;sv=np.sqrt(variances)
 off=np.stack([E[:,:,0]*E[:,:,1],E[:,:,0]*E[:,:,2],E[:,:,1]*E[:,:,2]],axis=2)
 ci=(off-off.mean(axis=1,keepdims=True))/df
 diag=(E*E-(E*E).mean(axis=1,keepdims=True))/df
 c12,c13,c23=out['c_off'].T;A=out['h'];den=np.c_[c13*sv[:,1],c12*sv[:,2]]
 UB=np.stack([sv[:,1,None]*ci[:,:,1]+.5*c13[:,None]/sv[:,1,None]*diag[:,:,1],
              sv[:,2,None]*ci[:,:,0]+.5*c12[:,None]/sv[:,2,None]*diag[:,:,2]],axis=2)
 U=np.concatenate([out['U'],UB],axis=2);flat=U.transpose(1,0,2).reshape(n,-1);V=np.zeros((B,4,4))
 for one in d['onehots']:
  G=len(one);s=(one@flat).reshape(G,B,4).transpose(1,0,2);s-=s.mean(axis=1,keepdims=True);V+=np.einsum('bgi,bgj->bij',s,s)*G/(G-1)
 vA=np.stack([V[:,0,0],V[:,1,1]],axis=1);vB=np.stack([V[:,2,2],V[:,3,3]],axis=1);vAB=np.stack([V[:,0,2],V[:,1,3]],axis=1)
 vplus=vA+margin**2*vB+2*margin*vAB;vminus=vA+margin**2*vB-2*margin*vAB
 lower=norm.sf((A+margin*den)/np.sqrt(np.maximum(vplus,1e-100)));upper=norm.sf((margin*den-A)/np.sqrt(np.maximum(vminus,1e-100)))
 p=np.maximum(lower,upper).max(axis=1);p[np.any(den<=0,axis=1)]=1.
 return p

def actual_equivalence(Y,d,margin):
 theta,IF,V,E=parameters(Y,d);A,JA=constraints_jacobian(theta);v=theta[3:6]
 den=np.array([theta[7]*np.sqrt(v[1]),theta[6]*np.sqrt(v[2])]);JB=np.zeros((2,9))
 JB[0,7]=np.sqrt(v[1]);JB[0,4]=.5*theta[7]/np.sqrt(v[1]);JB[1,6]=np.sqrt(v[2]);JB[1,5]=.5*theta[6]/np.sqrt(v[2])
 rows=[];ps=[];z=norm.ppf(.95)
 for j in range(2):
  gl=JA[j]+margin*JB[j];gu=JA[j]-margin*JB[j]
  lower=norm.sf((A[j]+margin*den[j])/np.sqrt(gl@V@gl));upper=norm.sf((margin*den[j]-A[j])/np.sqrt(gu@V@gu));ps.extend([lower,upper])
  va=JA[j]@V@JA[j];vb=JB[j]@V@JB[j];cab=JA[j]@V@JB[j]
  aa=den[j]**2-z*z*vb;bb=-2*(A[j]*den[j]-z*z*cab);cc=A[j]**2-z*z*va;discr=bb*bb-4*aa*cc
  bounded=aa>0 and discr>=0
  roots=sorted([(-bb-np.sqrt(discr))/(2*aa),(-bb+np.sqrt(discr))/(2*aa)]) if bounded else [np.nan,np.nan]
  rows.append({'contrast':'S_relative_to_N' if j==0 else 'T_relative_to_N','delta_point':A[j]/den[j] if den[j] else np.nan,
    'fieller90_bounded':bounded,'fieller90_lower':roots[0],'fieller90_upper':roots[1],'p_lower_boundary':float(lower),'p_upper_boundary':float(upper)})
 return {'p_joint_equivalence':float(max(ps)) if (den>0).all() else 1.,'two_contrasts':rows}
