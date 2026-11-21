from config import *
from engine import stats_from_moments,wilson
from simulation_dgp import innovations,patterns
from scipy.stats import norm
from concurrent.futures import ProcessPoolExecutor,as_completed
import pickle,time

def signature():
 return hashlib.sha256(''.join(sha(ROOT/'code'/n) for n in ['config.py','engine.py','simulation_dgp.py','run_calibration.py']).encode()+sha(PRIVATE/'ANALYSIS_TEMPLATES_RESTRICTED.pkl').encode()).hexdigest()

def worker(task):
 t,scenario,sig=task;key=t['design_id'];reps=FORMAL_REPS if t['formal'] else SENS_REPS
 seed=int(np.random.SeedSequence([ROOT_SEED,int(hashlib.sha256(key.encode()).hexdigest()[:8],16),SCENARIOS.index(scenario)]).generate_state(1,dtype=np.uint64)[0])
 file=SIM/(key+'__'+scenario+'.parquet');meta=SIM/(key+'__'+scenario+'.json')
 if file.exists() and meta.exists():
  old=json.loads(meta.read_text(encoding='utf-8'))
  if old.get('seed')==str(seed) and old.get('reps')==reps and old.get('run_signature')==sig and old.get('file_sha256')==sha(file):return old
 start_time=time.perf_counter();rng=np.random.default_rng(seed);pp=patterns(t);As=t['sigma'][:,None]*t['A'];Gram=As@As.T;means=[]
 for row in pp:
  mu=As.T@np.linalg.solve(Gram,row['effect']);means.append((t['A']@mu,np.einsum('gjn,n->gj',t['K'],mu,optimize=True)))
 arrays={f'{kind}_{p["name"]}':np.empty(reps) for p in pp for kind in ['SAME','DIFF']}
 failures=np.zeros(reps,bool)
 component_counts={p['name']:np.zeros(3,dtype=np.int64) for p in pp if p['kind']=='NULL' and np.any(p['effect']==0)}
 for start in range(0,reps,BATCH):
  B=min(BATCH,reps-start);Y=innovations(rng,t,scenario,B);b0=(t['A']@Y).T;U0=np.einsum('gjn,nb->bgj',t['K'],Y,optimize=True)
  for p,(mu_b,mu_u) in zip(pp,means):
   b=b0+mu_b;U=U0+mu_u[None,:,:]
   same,diff,unused,valid,z=stats_from_moments(b,U,t);failures[start:start+B]|=~valid
   arrays['SAME_'+p['name']][start:start+B]=np.where(valid,same,1.)
   arrays['DIFF_'+p['name']][start:start+B]=np.where(valid,diff,1.)
   if p['name'] in component_counts:component_counts[p['name']]+=(norm.sf(z)<=.025).sum(axis=0)
 pd.DataFrame({'replicate':np.arange(reps),'numerical_failure':failures,**arrays}).to_parquet(file,index=False,compression='zstd')
 selected_alphas={k:v for k,v in ALPHAS.items() if k!='BY28_DIAGNOSTIC_SUFFICIENT'} if t['formal'] else {k:v for k,v in ALPHAS.items() if k in ['NOMINAL','BY28_DIAGNOSTIC_SUFFICIENT']}
 rows=[];comp=[]
 for p in pp:
  for test,null in [('SAME',p['null_same']),('DIFF',p['null_diff'])]:
   for name,alpha in selected_alphas.items():
    k=int((arrays[test+'_'+p['name']]<=alpha).sum());lo,hi=wilson(k,reps)
    rows.append({'design_id':key,'gene':t['gene'],'focus14':t['gene'] in FOCUS14,'branch':t['branch'],'formal':t['formal'],'scenario':scenario,'pattern':p['name'],'kind':p['kind'],'test':test,'true_null':null,'alpha_name':name,'alpha':alpha,'reps':reps,'rejected':k,'rate':k/reps,'Wilson95_low':lo,'Wilson95_high':hi})
  if p['name'] in component_counts:
   for j in np.flatnonzero(p['effect']==0):
    k=int(component_counts[p['name']][j]);lo,hi=wilson(k,reps)
    comp.append({'design_id':key,'scenario':scenario,'pattern':p['name'],'assay':['N','S','T'][j],'reps':reps,'rejected025':k,'rate':k/reps,'Wilson95_high':hi})
 output={'design_id':key,'scenario':scenario,'formal':t['formal'],'reps':reps,'seed':str(seed),'n_failed':int(failures.sum()),'file':file.name,'file_sha256':sha(file),'run_signature':sig,'seconds':time.perf_counter()-start_time,'summary':rows,'component':comp}
 save_json(meta,output);return output

def decide(templates,summary,comp,inventory):
 decisions=[]
 for key,t in templates.items():
  part=summary.loc[summary.design_id.eq(key)]
  if part.empty:
   decisions.append({'design_id':key,'gene':t['gene'],'branch':t['branch'],'formal':t['formal'],'qualified':False,'reason':'DESIGN_GATE_FAILED'});continue
  null=part.loc[part.true_null];nom=null.loc[null.alpha_name.eq('NOMINAL')];tail=null.loc[~null.alpha_name.eq('NOMINAL')]
  cc=comp.loc[comp.design_id.eq(key)];inv=inventory.loc[inventory.design_id.eq(key)]
  failure=(inv.n_failed/inv.reps).max();np_=nom.Wilson95_high.le(.06).all();tp=tail.Wilson95_high.le(tail.alpha).all();mp=cc.Wilson95_high.le(.035).all()
  decisions.append({'design_id':key,'gene':t['gene'],'branch':t['branch'],'formal':t['formal'],'qualified':bool(t['eligible'] and np_ and tp and mp and failure<=.005),'nominal_pass':bool(np_),'tail_pass':bool(tp),'marginal_pass':bool(mp),'max_nominal_false_rate':nom.rate.max(),'max_nominal_Wilson95_upper':nom.Wilson95_high.max(),'max_tail_Wilson_over_alpha':(tail.Wilson95_high/tail.alpha).max(),'max_component025_Wilson95_upper':cc.Wilson95_high.max(),'max_failure_fraction':failure,'reason':'FIXED_FINITE_SCENARIO_RULE'})
 return pd.DataFrame(decisions)

def run():
 assert json.loads((RESULTS/'ENGINE_VALIDATION_SUMMARY.json').read_text(encoding='utf-8'))['failed']==0
 templates=pickle.load((PRIVATE/'ANALYSIS_TEMPLATES_RESTRICTED.pkl').open('rb'));targets=[t for t in templates.values() if t['eligible']]
 for a,n in [(min(ALPHAS[k] for k in ['BY68_SUFFICIENT','BY136_SUFFICIENT','BY204_SUFFICIENT']),FORMAL_REPS),(ALPHAS['BY28_DIAGNOSTIC_SUFFICIENT'],SENS_REPS)]:assert wilson(0,n)[1]<a
 sig=signature();grid=[]
 for t in targets:
  for p in patterns(t):grid.append({'design_id':t['design_id'],'pattern':p['name'],'effect_N':p['effect'][0],'effect_S':p['effect'][1],'effect_T':p['effect'][2],'null_same':p['null_same'],'null_diff':p['null_diff'],'kind':p['kind']})
 csvout(RESULTS/'SIMULATION_MEAN_PATTERN_GRID.csv',grid)
 save_json(RESULTS/'CALIBRATION_FREEZE_BEFORE_RUN.json',{'time':now(),'formal_reps':FORMAL_REPS,'sensitivity_reps':SENS_REPS,'scenarios':SCENARIOS,'alphas':ALPHAS,'designs':len(targets),'patterns_per_design':37,'run_signature':sig,'code_sha256':{n:sha(ROOT/'code'/n) for n in ['config.py','engine.py','simulation_dgp.py','run_calibration.py']},'new_noncore_clinical_target_coefficients_fitted':False,'interpretation':'Conditional fixed-design finite calibration; not a universal guarantee or missingness solution'})
 tasks=[(t,scenario,sig) for t in targets for scenario in SCENARIOS];outputs=[]
 with ProcessPoolExecutor(max_workers=WORKERS) as pool:
  futures=[pool.submit(worker,t) for t in tasks]
  for i,f in enumerate(as_completed(futures),1):
   o=f.result();outputs.append(o);log(f'Calibration {i}/{len(tasks)} {o["design_id"]} {o["scenario"]}: {o["reps"]} noise draws, {o["seconds"]:.1f}s')
 outputs.sort(key=lambda x:(x['design_id'],x['scenario']))
 full=pd.DataFrame([r for o in outputs for r in o['summary']]);comp=pd.DataFrame([r for o in outputs for r in o['component']])
 inv=pd.DataFrame([{k:v for k,v in o.items() if k not in ['summary','component']} for o in outputs])
 csvout(RESULTS/'CALIBRATION_AND_POWER_FULL.csv',full);csvout(RESULTS/'COMPONENT_NULL_CALIBRATION.csv',comp);csvout(RESULTS/'SIMULATION_FILE_INVENTORY.csv',inv)
 decisions=decide(templates,full,comp,inv);csvout(RESULTS/'CALIBRATION_QUALIFICATION_DECISIONS.csv',decisions)
 save_json(RESULTS/'CALIBRATION_RUN_SUMMARY.json',{'time':now(),'designs':len(targets),'cells':len(tasks),'independent_noise_draws':int(inv.reps.sum()),'pattern_evaluations_reusing_noise':int(inv.reps.sum()*37),'formal_designs':sum(t['formal'] for t in targets),'sensitivity_designs':sum(not t['formal'] for t in targets),'qualified':int(decisions.qualified.sum()),'run_signature':sig})
 log(f'Calibration complete: {int(decisions.qualified.sum())}/{len(templates)} designs qualified.')
 print(decisions[['gene','branch','qualified','max_nominal_false_rate','max_tail_Wilson_over_alpha']].to_string(index=False))
if __name__=='__main__':assert_unsealed();run()

