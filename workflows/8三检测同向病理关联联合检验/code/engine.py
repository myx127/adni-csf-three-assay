"""Fixed-X coefficient conjunction engine; no latent-factor decomposition."""
from config import *
from scipy.special import ndtr

def setup(template,scenario):
 d=template['d'];Q=d['Q'];a=d['a'];operators=[a[None,:]];blocks=[];start=1
 for one in d['onehots']:
  weighted=one*a[None,:];K=weighted-(weighted@Q)@Q.T
  operators.append(K);blocks.append((start,start+len(one),len(one)));start+=len(one)
 O=np.vstack(operators);slopes=[];slope_sd=.35 if scenario in ['SLOPE_N','SLOPES3','T5_SLOPES3'] else 0.
 axes=[0] if scenario=='SLOPE_N' else range(3) if slope_sd else []
 for axis in axes:
  one=d['onehots'][axis];lab=d['groups'][:,axis];age=d['q_age'];q=age-(one@age/one.sum(axis=1))[lab]
  slopes.append((lab,len(one),q))
 hetero=np.ones(d['n'])
 if scenario=='HETERO':
  pc=d['P']-d['P'].mean();hetero=np.sqrt(.5+.5*pc*pc/np.mean(pc*pc))
 trace=float(np.sum((1-np.sum(Q*Q,axis=1))*hetero**2))
 for lab,G,q in slopes:
  A=np.eye(G)[lab]*q[:,None];AQ=A-Q@(Q.T@A);trace+=slope_sd**2*np.sum(AQ*AQ)
 scale=np.sqrt(trace/d['df']);assert np.isfinite(scale) and scale>0
 return {'template':template,'d':d,'scenario':scenario,'O':O,'blocks':blocks,'hetero':hetero,'slopes':slopes,'slope_sd':slope_sd,'normalizer':scale}

def generate(spec,B,rng):
 d=spec['d'];heavy=spec['scenario']=='T5_SLOPES3'
 def draw(shape):return rng.standard_t(5,size=shape)*np.sqrt(3/5) if heavy else rng.normal(size=shape)
 W=draw((B,d['n'],3))*spec['hetero'][None,:,None]
 for lab,G,q in spec['slopes']:W+=draw((B,G,3))[:,lab,:]*q[None,:,None]*spec['slope_sd']
 return (W/spec['normalizer'])@spec['template']['L'].T

def fit_noise(Y,spec):
 B,n,_=Y.shape;projected=(spec['O']@Y.transpose(1,0,2).reshape(n,-1)).reshape(-1,B,3).transpose(1,0,2)
 coef=projected[:,0,:];variance=np.zeros((B,3))
 for first,last,G in spec['blocks']:
  scores=projected[:,first:last,:];scores=scores-scores.mean(axis=1,keepdims=True)
  variance+=np.sum(scores*scores,axis=1)*G/(G-1)
 sd=spec['template']['sigma'];se=np.sqrt(np.maximum(variance,0))
 valid=np.isfinite(coef).all(axis=1)&np.isfinite(se).all(axis=1)&(se>0).all(axis=1)
 return coef/sd,se/sd,valid

def joint_from_z(z):
 positive=np.max(ndtr(-z),axis=-1);negative=np.max(ndtr(z),axis=-1)
 shared=np.minimum(1.,2*np.minimum(positive,negative))
 direction=np.where(positive<negative,1,np.where(negative<positive,-1,0))
 return positive,negative,shared,direction

def evaluate(noise,se,effect,valid=None):
 z=(noise+np.asarray(effect)[None,:])/se
 plus,minus,p,direction=joint_from_z(z)
 if valid is not None:p=np.where(valid,p,np.nan)
 return plus,minus,p,direction
