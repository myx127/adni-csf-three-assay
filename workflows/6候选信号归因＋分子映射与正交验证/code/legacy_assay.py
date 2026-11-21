"""Separate single-assay pTau association; this does not test relative constraints."""
from config import *
from scipy.stats import t

def single_batch(Y,d):
 Y=np.asarray(Y,float)
 if Y.ndim==1:Y=Y[None,:]
 E=Y-(Y@d['Q'])@d['Q'].T;b=Y@d['a'];u=E*d['a'][None,:];V=np.zeros(len(Y))
 for one in d['onehots']:
  sums=u@one.T;sums-=sums.mean(axis=1,keepdims=True)
  V+=(sums*sums).sum(axis=1)*len(one)/(len(one)-1)
 se=np.sqrt(V);df_t=min(len(o) for o in d['onehots'])-1
 p=t.sf(b/se,df_t);p[~np.isfinite(p)]=np.nan
 return b,se,p,E,df_t

def single_summary(y,d):
 b,se,p,E,df=single_batch(y,d);sd=np.sqrt((E[0]@E[0])/d['df']);z=t.ppf(.975,df)
 return {'b_P':b[0],'se_SUM2_CR1':se[0],'t_df_min_Gminus1':df,'p_positive_direction':p[0],
  'b_P_95low':b[0]-z*se[0],'b_P_95high':b[0]+z*se[0],
  'b_P_per_residual_SD':b[0]/sd,'CI_per_residual_SD_conditional_on_SD_low':(b[0]-z*se[0])/sd,
  'CI_per_residual_SD_conditional_on_SD_high':(b[0]+z*se[0])/sd,'residual_SD':sd}
