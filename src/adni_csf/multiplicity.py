import numpy as np

def by(p):
 p=np.asarray(p,float); m=len(p); order=np.argsort(p); rank=np.arange(1,m+1)
 q=np.minimum.accumulate((p[order]*m*np.sum(1/rank)/rank)[::-1])[::-1]
 out=np.ones(m); out[order]=np.minimum(1,q); return out
