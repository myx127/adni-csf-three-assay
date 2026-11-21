from __future__ import annotations
import numpy as np
from scipy.special import ndtri,expit
from reference_stage2_data import batch_labels

def generate(n,rng,settings,audited):
 c={'sd':[.6,.8,1.2],'beta':.55,'uvar':[1.,1.,1.],'batch_sd':[0.,0.,0.],**settings}
 batches=batch_labels(n,rng,audited)
 probs=np.array(c.get('group_probs',[1/3]*3));counts=np.floor(n*probs).astype(int)
 counts[np.argsort(-(n*probs-counts))[:n-counts.sum()]]+=1
 groups=np.repeat(np.arange(3),counts);rng.shuffle(groups)
 if 'exposure_cluster_axis' in c:
  labels=batches[:,c['exposure_cluster_axis']];G=int(labels.max())+1
  sizes=np.floor(G*probs).astype(int);sizes[np.argsort(-(G*probs-sizes))[:G-sizes.sum()]]+=1
  group_by_cluster=np.repeat(np.arange(3),sizes);rng.shuffle(group_by_cluster)
  groups=group_by_cluster[labels]
 technical=c.get('technical','none');plexshift=rng.normal(scale=.75,size=batches[:,2].max()+1)
 if technical in ['raw','gis']:
  score=plexshift[batches[:,2]]+rng.normal(scale=.25,size=n)
  groups=np.empty(n,dtype=int);groups[np.argsort(score)]=np.minimum(2,np.arange(n)*3//n)
 P=groups.astype(float)-1.
 T=c['beta']*P+np.sqrt(np.array(c['uvar'])[groups])*rng.normal(size=n)
 z=rng.normal(size=(n,3));rho=c.get('rho',0.)
 rho=np.array(rho)[groups] if isinstance(rho,list) else rho
 z[:,1]=rho*z[:,0]+np.sqrt(1-np.asarray(rho)**2)*z[:,1]
 eps=z*np.array(c['sd']);eps[:,2]*=np.exp(c.get('hetero',0.)*P)
 response=np.repeat(T[:,None],3,axis=1)
 q=c.get('quadratic',0.);response[:,2]+=q*T*T
 if c.get('saturation',False):response[:,2]=1.8*np.tanh(T/1.8)
 Y=response+eps+P[:,None]*np.array(c.get('delta',[0,0,0]))+(P*T)[:,None]*np.array(c.get('load_delta',[0,0,0]))
 for j,sd in enumerate(c['batch_sd']):
  if sd:
   shift=rng.normal(scale=sd,size=batches[:,j].max()+1);Y[:,j]+=shift[batches[:,j]]
 if technical=='raw':Y[:,2]+=plexshift[batches[:,2]]
 if technical=='gis':Y[:,2]-=rng.normal(scale=.6,size=len(plexshift))[batches[:,2]]
 missing=c.get('missing','none')
 if missing=='mcar':Y[rng.uniform(size=Y.shape)<.2]=np.nan
 if missing in ['lod30','lod70','floor30']:
  fraction=.7 if missing=='lod70' else .3
  threshold=np.sqrt(c['beta']**2*2/3+np.mean(c['uvar'])+c['sd'][0]**2)*ndtri(fraction)
  if missing=='floor30':Y[:,0]=np.maximum(Y[:,0],threshold)
  else:Y[Y[:,0]<threshold,0]=np.nan
 if missing=='selective':
  s=c.get('selection_strength',.7)
  Y[rng.uniform(size=n)>expit(1.6-s*Y[:,2]+s/.7*.6*P),2]=np.nan
 good=np.isfinite(Y).all(axis=1)
 # Oracle removes only a known fixed quadratic component using true T; no extra data in deployed paths.
 oracle=Y.copy();oracle[:,2]-=q*T*T
 return {'Y':Y[good],'groups':groups[good],'batches':batches[good],'P':P[good],
         'T':T[good],'oracle_quadratic_Y':oracle[good],'n_original':n,
         'true_kappa12':None if isinstance(c.get('rho',0.),list) else float(c.get('rho',0.)*c['sd'][0]*c['sd'][1])}
