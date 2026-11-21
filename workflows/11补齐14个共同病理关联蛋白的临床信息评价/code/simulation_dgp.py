from config import *
def patterns(t):
 unit=t['target_unit'];rows=[]
 def add(name,v,nullsame=False,nulldiff=False,nullany=False,kind='NULL'):
  rows.append({'name':name,'effect':np.array(v,float),'null_same':nullsame,'null_diff':nulldiff,'null_any':nullany,'kind':kind})
 add('ZERO',[0,0,0],True,True,True)
 for mag in [.35,.70]:
  for sign in [-1,1]:add(f'COMMON_{mag}_{sign}',[sign*mag*unit]*3,False,True)
 for axis in range(3):
  for sign in [-1,1]:
   for typ,other,one in [('ONE_ZERO',.7,0.),('TWO_ZERO',0.,.7),('MIXED',.35,-.02)]:
    v=np.full(3,other);v[axis]=one;add(f'{typ}_{axis}_{sign}',v*sign*unit,True,False)
 if t['branch']=='AMY_CL':mags=[5.,10.,20.]
 elif t['branch']=='MEM_BASE':mags=[.1,.25,.5]
 else:mags=[.02,.05,.10]
 for mag in mags:
  for sign in [-1,1]:add(f'POWER_BALANCED_{mag}_{sign}',[sign*mag]*3,False,True,False,'POWER')
 for axis in range(3):
  for sign in [-1,1]:
   v=np.full(3,mags[1]);v[axis]=mags[0];add(f'POWER_WEAK_{axis}_{sign}',v*sign,False,False,False,'POWER')
 for sign in [-1,1]:add(f'POWER_DIFFERENTIAL_{sign}',np.array([0.,mags[1],mags[2]])*sign,True,False,False,'POWER')
 return rows

def innovations(rng,t,scenario,B):
 n=t['n_people'];M=t['A'].shape[1];idx=t['case_local'];scale=t['noise_scale']
 z=rng.standard_t(5,size=(M,B))/np.sqrt(5/3) if scenario=='T5_TECHNICAL_SLOPES' else rng.normal(size=(M,B))
 if scenario=='AR1_AND_PERSON_SLOPE':
  if t['is_quadratic']:
   v1=t['ar_sd'][::2]**2;v2=t['ar_sd'][1::2]**2;cross=t['ar_cross'];a=np.sqrt(np.maximum(v1,1e-12));b=cross/a;c=np.sqrt(np.maximum(v2-b*b,0));first=z[::2].copy();second=z[1::2].copy();z[::2]=a[:,None]*first;z[1::2]=b[:,None]*first+c[:,None]*second
  else:z*=t['ar_sd'][:,None]
  z+=.35*t['slope_load'][:,None]*rng.normal(size=(n,B))[idx]
 else:z*=t['visit_sd'][:,None]
 if scenario=='HETERO_P':z*=np.exp(.3*np.clip(t['P'][idx],-3,3))[:,None]
 if scenario in ['TECHNICAL_AGE_SLOPES','T5_TECHNICAL_SLOPES']:
  age=t['age'];age=(age-age.mean())/(age.std() if age.std()>0 else 1)
  for labels in t['case_groups']:
   G=int(labels.max())+1;mean=np.bincount(labels,weights=age)/np.bincount(labels);q=age-mean[labels]
   z+=.35*t['slope_load'][:,None]*q[idx,None]*rng.normal(size=(G,B))[labels[idx]]
 return z*scale
