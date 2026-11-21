from __future__ import annotations
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
import argparse,json,time
import numpy as np
import pandas as pd
from cal_config import *
from cal_data import generate
from cal_inference import test,wild_t,sums_by_cluster,cr1_covariance,inverse_metric,bootstrap_tail
from reference_stage2_models import multiway_covariance
from analyze_calibration import wilson

def followup_grid():
 result=[]
 for phase,root,reps,B in [('MEAN_DEV',2026102704,1000,199),('MEAN_CONFIRM',2026102705,2000,499)]:
  for s in calibration_scenarios():
   seed=int(np.random.SeedSequence([root,s['index'],571]).generate_state(1,dtype=np.uint64)[0])
   result.append({'cell_id':f'{phase}_{s["scenario_id"]}_n571','phase':phase,'scenario_id':s['scenario_id'],
                  'scenario_index':s['index'],'n':571,'repetitions':reps,'bootstrap_B':B,'seed':str(seed),'cluster_axis':s['cluster_axis']})
 return result

def components(Y,P,kappa=0.):
 n=len(Y);gc=P-P.mean();sxx=gc@gc
 b=gc@Y/sxx;e=Y-Y.mean(axis=0)-gc[:,None]*b
 products=np.c_[e[:,0]*e[:,1],e[:,0]*e[:,2],e[:,1]*e[:,2]]
 c12,c13,c23=products.mean(axis=0);c12-=kappa
 h=np.array([b[1]*c13-b[0]*c23,b[2]*c12-b[0]*c23])
 J=np.array([[-c23,c13,0.,0.,b[1],-b[0]],[-c23,0.,c12,b[2],0.,-b[0]]])
 IF=np.c_[gc[:,None]*e/sxx,(products-products.mean(axis=0))/n]@J.T
 return h,IF

def bootstrap_full_mean(Y,P,batches,B,rng,return_weights=False):
 W=np.ones((B,len(Y)))
 for j in range(3):W*=rng.exponential(size=(B,int(batches[:,j].max())+1))[:,batches[:,j]]
 total=W.sum(axis=1);meanP=(W@P)/total;varP=(W@(P*P))/total-meanP**2
 features=np.c_[Y,P[:,None]*Y,Y[:,0]*Y[:,1],Y[:,0]*Y[:,2],Y[:,1]*Y[:,2]]
 avg=(W@features)/total[:,None];mean=avg[:,:3]
 slope=(avg[:,3:6]-meanP[:,None]*mean)/varP[:,None]
 cov=avg[:,6:]-np.c_[mean[:,0]*mean[:,1],mean[:,0]*mean[:,2],mean[:,1]*mean[:,2]]-varP[:,None]*np.c_[slope[:,0]*slope[:,1],slope[:,0]*slope[:,2],slope[:,1]*slope[:,2]]
 h=np.c_[slope[:,1]*cov[:,1]-slope[:,0]*cov[:,2],slope[:,2]*cov[:,0]-slope[:,0]*cov[:,2]]
 return (h,W) if return_weights else h

def mean_infer(data,axis,B,rng):
 from cal_inference import failed
 if len(data['Y'])<30 or np.bincount(data['groups'],minlength=3).min()<10:return {m:failed('INSUFFICIENT_COMPLETE_GROUPS') for m in METHODS}
 Y,P,batches=data['Y'],data['P'],data['batches'];h,IF=components(Y,P)
 S=[sums_by_cluster(IF,batches[:,j]) for j in range(3)];V=[cr1_covariance(s) for s in S];Vsum=sum(V);G=len(S[axis]);Gmin=min(map(len,S))
 results={'IID_CHI2':test(h,IF.T@IF),'CR1_CHI2':test(h,V[axis]),'CR1_HT_F':test(h,V[axis],G),
    'WILD_SCORE_T':wild_t(h,S[axis],B,rng),'SUM3_CHI2':test(h,Vsum),'SUM3_MIN_F':test(h,Vsum,Gmin),
    'CGM3_CHI2':test(h,multiway_covariance(IF,batches))}
 inv,status=inverse_metric(Vsum)
 results['PIGEONHOLE_CENTERED']=bootstrap_tail(h,inv,bootstrap_full_mean(Y,P,batches,B,rng)) if inv is not None else failed(status)
 return results

def profile(data):
 out={};gridk=np.linspace(-.4,.4,9)*.48
 for method in ['CR1_HT_F','SUM3_CHI2']:
  ps=[]
  for kappa in gridk:
   h,IF=components(data['Y'],data['P'],kappa);S=[sums_by_cluster(IF,data['batches'][:,j]) for j in range(3)];V=[cr1_covariance(s) for s in S]
   ps.append(test(h,V[2],len(S[2]))['p_value'] if method=='CR1_HT_F' else test(h,sum(V))['p_value'])
  ps=np.array(ps);valid=np.isfinite(ps).all();true=data['true_kappa12'];on_grid=true is not None and min(abs(gridk-true))<1e-10
  out[method+'_zero_p']=float(ps[4]);out[method+'_grid_max_p']=float(max(ps)) if valid else np.nan
  out[method+'_classification_flips']=bool(np.any(ps<=.05)&np.any(ps>.05)) if valid else False
  out[method+'_true_fixed_kappa_p']=float(ps[np.argmin(abs(gridk-true))]) if on_grid else np.nan
  out[method+'_true_fixed_kappa_on_grid']=on_grid
 h,IF=components(data['oracle_quadratic_Y'],data['P']);S=sums_by_cluster(IF,data['batches'][:,2])
 out['ORACLE_REMOVE_QUADRATIC_CR1_HT_F_p']=test(h,cr1_covariance(S),len(S))['p_value']
 return out

def run(phase):
 audited=json.loads((RESULTS/'AUDITED_BATCH_SIZES.json').read_text(encoding='utf-8'))
 if phase=='MEAN_DEV':
  csv(RESULTS/'MEAN_FOLLOWUP_GRID.csv',followup_grid())
  save_json(RESULTS/'MEAN_FOLLOWUP_FREEZE.json',{'protocol_sha256':sha(ROOT/'LOW_DIMENSIONAL_FOLLOWUP_PROTOCOL_v1.0.md'),
     'core':{f:sha(ROOT/'code'/f) for f in ['mean_followup.py','cal_data.py','cal_inference.py','cal_config.py']},
     'eight_df_development_already_seen':True,'mean_validation_not_yet_generated':True})
 elif phase=='MEAN_CONFIRM':
  frozen=json.loads((RESULTS/'MEAN_SELECTION_FROZEN.json').read_text(encoding='utf-8'))
  assert all(sha(ROOT/'code'/f)==v for f,v in frozen['core'].items())
  assert sha(ROOT/'LOW_DIMENSIONAL_FOLLOWUP_PROTOCOL_v1.0.md')==frozen['protocol_sha256']
 lookup={s['scenario_id']:s for s in (sensitivity_scenarios() if phase=='MEAN_SENS_REPLAY' else calibration_scenarios())}
 cells=[{**c,'original_cell_id':c['cell_id'],'cell_id':c['cell_id'].replace('SENSITIVITY','MEAN_SENS_REPLAY'),'phase':'MEAN_SENS_REPLAY'} for c in grid() if c['phase']=='SENSITIVITY'] if phase=='MEAN_SENS_REPLAY' else [c for c in followup_grid() if c['phase']==phase]
 start=time.time();total=0
 for cell in cells:
  rows=[];profiles=[]
  for rep in range(cell['repetitions']):
   rng=np.random.default_rng(np.random.SeedSequence([int(cell['seed']),rep,0]));brng=np.random.default_rng(np.random.SeedSequence([int(cell['seed']),rep,2]))
   data=generate(cell['n'],rng,lookup[cell['scenario_id']]['settings'],audited);est=mean_infer(data,cell['cluster_axis'],cell['bootstrap_B'],brng)
   for method,result in est.items():rows.append({'replicate':rep,'method':method,**result,'n_complete':len(data['Y'])})
   if phase=='MEAN_SENS_REPLAY':
    valid=len(data['Y'])>=30 and np.bincount(data['groups'],minlength=3).min()>=10
    profiles.append({'replicate':rep,'status':'OK' if valid else 'INSUFFICIENT_COMPLETE_GROUPS',**(profile(data) if valid else {})})
  pd.DataFrame(rows).to_parquet(RESULTS/'replicates'/(cell['cell_id']+'.parquet'),index=False)
  if profiles:pd.DataFrame(profiles).to_parquet(RESULTS/'replicates'/(cell['cell_id']+'_sensitivity.parquet'),index=False)
  total+=cell['repetitions'];log(f'MEAN {cell["cell_id"]} total={total} elapsed={time.time()-start:.1f}s')
 save_json(RESULTS/(phase+'_RUN_SUMMARY.json'),{'phase':phase,'repeats':total,'independent_new_repeats':0 if phase=='MEAN_SENS_REPLAY' else total,'method_rows':total*8})

def summarize(select=False):
 rows=[]
 for c in followup_grid():
  path=RESULTS/'replicates'/(c['cell_id']+'.parquet')
  if not path.exists():continue
  for method,q in pd.read_parquet(path).groupby('method'):
   ok=q.status.eq('OK');k=int((ok&q.p_value.le(.05)).sum());lo,hi=wilson(k,len(q))
   rows.append({**c,'method':method,'n_replicates':len(q),'n_failed':int((~ok).sum()),'n_rejected':k,'rejection_rate':k/len(q),
      'wilson95_lower':lo,'wilson95_upper':hi,'failure_rate':float((~ok).mean())})
 out=pd.DataFrame(rows);csv(RESULTS/'MEAN_FOLLOWUP_SUMMARY.csv',out)
 if select:
  dev=out.loc[out.phase.eq('MEAN_DEV')];decisions=[];scores=[]
  for scope,cases in {'ONE_WAY_SCOPE':['C01','C02','C03','C04','C05','C08','C09','C10'],'CROSSED_SCOPE':[f'C{i:02d}' for i in range(1,11)]}.items():
   eligible=[]
   for method in NEW_METHODS:
    q=dev.loc[dev.scenario_id.isin(cases)&dev.method.eq(method)];passed=len(q)==len(cases) and q.wilson95_upper.le(.08).all() and q.failure_rate.le(.01).all()
    power=dev.loc[dev.scenario_id.isin(['A01','A02'])&dev.method.eq(method),'rejection_rate'].mean()
    scores.append({'scope':scope,'method':method,'passed':bool(passed),'maximum_null_rate':float(q.rejection_rate.max()),'power_A01_A02':float(power)})
    if passed:eligible.append((power,-NEW_METHODS.index(method),method))
   decisions.append({'scope':scope,'selected_method':max(eligible)[2] if eligible else None})
  frozen=json.loads((RESULTS/'MEAN_FOLLOWUP_FREEZE.json').read_text(encoding='utf-8'))
  save_json(RESULTS/'MEAN_SELECTION_FROZEN.json',{**frozen,'decisions':decisions,'confirmation_generated_before_selection':False})
  csv(RESULTS/'MEAN_SELECTION_SCORES.csv',scores);print(pd.DataFrame(scores).to_string(index=False));print(decisions)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--phase',choices=['MEAN_DEV','MEAN_CONFIRM','MEAN_SENS_REPLAY']);p.add_argument('--select',action='store_true');args=p.parse_args()
 if args.phase:run(args.phase)
 else:summarize(args.select)
