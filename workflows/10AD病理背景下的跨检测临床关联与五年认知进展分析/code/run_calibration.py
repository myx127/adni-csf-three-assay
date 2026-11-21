from config import *
from engine import stats_from_moments,wilson
from scipy.stats import norm
from concurrent.futures import ProcessPoolExecutor,as_completed
import pickle,time

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

def worker(task):
 t,scenario=task;key=t['design_id'];reps=FORMAL_REPS if t['formal'] else SENS_REPS
 seed=int(np.random.SeedSequence([ROOT_SEED,int(hashlib.sha256(key.encode()).hexdigest()[:8],16),SCENARIOS.index(scenario)]).generate_state(1,dtype=np.uint64)[0])
 file=SIM/(key+'__'+scenario+'.parquet');meta=SIM/(key+'__'+scenario+'.json')
 if file.exists() and meta.exists():
  old=json.loads(meta.read_text(encoding='utf-8'))
  if old.get('seed')==str(seed) and old.get('reps')==reps and old.get('file_sha256')==sha(file):return old
 rng=np.random.default_rng(seed);pp=patterns(t);As=t['sigma'][:,None]*t['A'];Gram=As@As.T;means=[]
 for row in pp:
  mu=As.T@np.linalg.solve(Gram,row['effect']);b=t['A']@mu;score=np.einsum('gjn,n->gj',t['K'],mu,optimize=True);means.append((b,score))
 arrays={f'{kind}_{p["name"]}':np.empty(reps) for p in pp for kind in ['SAME','DIFF','ANY']}
 failures=np.zeros(reps,bool);components={p['name']:np.zeros((reps,3),np.float64) for p in pp if p['kind']=='NULL' and np.any(p['effect']==0)}
 for start in range(0,reps,BATCH):
  B=min(BATCH,reps-start);Y=innovations(rng,t,scenario,B);b0=(t['A']@Y).T;U0=np.einsum('gjn,nb->bgj',t['K'],Y,optimize=True)
  for p,(mu_b,mu_u) in zip(pp,means):
   b=b0+mu_b;U=U0+mu_u[None,:,:]
   same,diff,any_,valid,z=stats_from_moments(b,U,t);failures[start:start+B]|=~valid
   for name,vals in [('SAME',same),('DIFF',diff),('ANY',any_)]:arrays[name+'_'+p['name']][start:start+B]=np.where(valid,vals,1.)
   if p['name'] in components:components[p['name']][start:start+B]=norm.sf(z)
 frame=pd.DataFrame({'replicate':np.arange(reps),'numerical_failure':failures,**arrays});frame.to_parquet(file,index=False,compression='zstd')
 rows=[]
 for p in pp:
  for test,null in [('SAME',p['null_same']),('DIFF',p['null_diff']),('ANY',p['null_any'])]:
   for alpha_name,alpha in ALPHAS.items():
    if not t['formal'] and alpha_name=='BY32_SUFFICIENT':continue
    vals=arrays[test+'_'+p['name']];k=int((vals<=alpha).sum());lo,hi=wilson(k,reps)
    rows.append({'design_id':key,'gene':t['gene'],'branch':t['branch'],'formal':t['formal'],'scenario':scenario,'pattern':p['name'],'kind':p['kind'],'test':test,'true_null':null,'alpha_name':alpha_name,'alpha':alpha,'reps':reps,'rejected':k,'rate':k/reps,'Wilson95_low':lo,'Wilson95_high':hi})
 comp=[]
 for p in pp:
  if p['name'] not in components:continue
  for j in np.flatnonzero(p['effect']==0):
   values=components[p['name']][:,j];k=int((values<=.025).sum());lo,hi=wilson(k,reps)
   comp.append({'design_id':key,'scenario':scenario,'pattern':p['name'],'assay':['N','S','T'][j],'reps':reps,'rejected025':k,'rate':k/reps,'Wilson95_high':hi})
 output={'design_id':key,'scenario':scenario,'formal':t['formal'],'reps':reps,'seed':str(seed),'n_failed':int(failures.sum()),'file':file.name,'file_sha256':sha(file),'summary':rows,'component':comp}
 save_json(meta,output);return output

def run():
 assert json.loads((RESULTS/'ENGINE_VALIDATION_SUMMARY.json').read_text(encoding='utf-8'))['failed']==0
 templates=pickle.load((PRIVATE/'ANALYSIS_TEMPLATES_RESTRICTED.pkl').open('rb'));targets=[t for t in templates.values() if t['eligible']]
 grid=[]
 for t in targets:
  for p in patterns(t):grid.append({'design_id':t['design_id'],'pattern':p['name'],'effect_N':p['effect'][0],'effect_S':p['effect'][1],'effect_T':p['effect'][2],'null_same':p['null_same'],'null_diff':p['null_diff'],'null_any':p['null_any'],'kind':p['kind']})
 csvout(RESULTS/'SIMULATION_MEAN_PATTERN_GRID.csv',grid)
 save_json(RESULTS/'CALIBRATION_FREEZE_BEFORE_RUN.json',{'time':now(),'formal_reps':FORMAL_REPS,'sensitivity_reps':SENS_REPS,'scenarios':SCENARIOS,'designs':len(targets),'formal_designs':sum(t['formal'] for t in targets),'patterns_per_design':37,'code_sha256':{p.name:sha(p) for p in [Path(__file__),ROOT/'code/engine.py',ROOT/'code/config.py']},'real_clinical_target_coefficients_fitted':False,'fixed_observed_exposures_and_schedules':True,'interpretation':'Conditional finite scenarios; standardization uncertainty propagated analytically, not a universal random-exposure calibration'})
 tasks=[(t,scenario) for t in targets for scenario in SCENARIOS];allout=[]
 with ProcessPoolExecutor(max_workers=WORKERS) as pool:
  fs=[pool.submit(worker,task) for task in tasks]
  for i,f in enumerate(as_completed(fs),1):
   out=f.result();allout.append(out)
   log(f'Calibration {i}/{len(tasks)}: {out["design_id"]} {out["scenario"]};{out["reps"]} draws')
 summary=pd.DataFrame([r for o in allout for r in o['summary']]);comp=pd.DataFrame([r for o in allout for r in o['component']])
 csvout(RESULTS/'CALIBRATION_AND_POWER_FULL.csv',summary);csvout(RESULTS/'COMPONENT_NULL_CALIBRATION.csv',comp)
 inventory=[{k:v for k,v in o.items() if k not in ['summary','component']} for o in allout];csvout(RESULTS/'SIMULATION_FILE_INVENTORY.csv',inventory)
 decisions=[]
 for key,t in templates.items():
  rows=summary.loc[summary.design_id.eq(key)];cc=comp.loc[comp.design_id.eq(key)];inv=pd.DataFrame(inventory);ii=inv.loc[inv.design_id.eq(key)]
  if rows.empty:
   decisions.append({'design_id':key,'qualified':False,'reason':'DESIGN_GATE_FAILED'});continue
  relevant=['SAME','DIFF'] if 'INTERACTION' not in t['branch'] else ['ANY','DIFF']
  null=rows.loc[rows.true_null&rows.test.isin(relevant)]
  nominal=null.loc[null.alpha_name.eq('NOMINAL')]
  tail=null.loc[~null.alpha_name.eq('NOMINAL')]
  failures=(ii.n_failed/ii.reps).max();nompass=nominal.Wilson95_high.le(.06).all();tailpass=tail.Wilson95_high.le(tail.alpha).all();margpass=cc.Wilson95_high.le(.035).all() if 'SAME' in relevant else True
  decisions.append({'design_id':key,'gene':t['gene'],'branch':t['branch'],'formal':t['formal'],'qualified':bool(t['eligible'] and nompass and tailpass and margpass and failures<=.005),'nominal_pass':bool(nompass),'tail_pass':bool(tailpass),'marginal_pass':bool(margpass),'max_nominal_false_rate':nominal.rate.max(),'max_nominal_Wilson95_upper':nominal.Wilson95_high.max(),'max_tail_Wilson_over_alpha':(tail.Wilson95_high/tail.alpha).max(),'max_failure_fraction':failures,'reason':'FIXED_FINITE_SCENARIO_RULE'})
 csvout(RESULTS/'CALIBRATION_QUALIFICATION_DECISIONS.csv',decisions)
 save_json(RESULTS/'CALIBRATION_RUN_SUMMARY.json',{'time':now(),'designs':len(targets),'cells':len(tasks),'independent_noise_draws':sum(o['reps'] for o in allout),'pattern_evaluations_reuse_noise':sum(o['reps']*37 for o in allout),'qualified':sum(r['qualified'] for r in decisions),'formal_designs':sum(t['formal'] for t in targets),'sensitivity_designs':sum(not t['formal'] for t in targets)})
 log(f'Calibration completed:{sum(r["qualified"] for r in decisions)}/{len(decisions)} qualified')
 print(pd.DataFrame(decisions)[['gene','branch','qualified','max_nominal_false_rate','max_nominal_Wilson95_upper','max_tail_Wilson_over_alpha']].to_string(index=False))
if __name__=='__main__':
 assert_unsealed();run()
