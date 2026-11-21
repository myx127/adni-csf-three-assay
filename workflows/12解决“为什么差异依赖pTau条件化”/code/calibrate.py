from config import *
from engine import moment_stats,wilson
from concurrent.futures import ProcessPoolExecutor,as_completed
import pickle,time
def innovations(rng,t,scenario,B):
 n=t['n_people'];z=rng.standard_t(5,size=(n,B))/np.sqrt(5/3) if scenario=='T5_TECHNICAL_SLOPES' else rng.normal(size=(n,B))
 if scenario=='AR1_AND_PERSON_SLOPE':z=z*t['ar_sd'][:,None]+.35*t['slope_load'][:,None]*rng.normal(size=(n,B))
 if scenario=='HETERO_P':z*=np.exp(.3*np.clip(t['P'],-3,3))[:,None]
 if scenario in ['TECHNICAL_AGE_SLOPES','T5_TECHNICAL_SLOPES']:
  age=(t['age']-np.mean(t['age']))/(np.std(t['age']) or 1)
  for labels in t['groups']:
   G=int(labels.max())+1;count=np.bincount(labels,minlength=G);mean=np.divide(np.bincount(labels,weights=age,minlength=G),count,out=np.zeros(G),where=count>0);q=age-mean[labels]
   z+=.35*t['slope_load'][:,None]*q[:,None]*rng.normal(size=(G,B))[labels]
 return z*t['noise_scale']
def patterns(t):
 A=t['Ah'];Gram=A@A.T;out=[{'name':'ZERO','kind':'NULL','mu':np.zeros(t['n_people']),'target':np.zeros(2)}]
 for j in range(6):
  v=t['prototypes'][:,j].copy();v-=t['Qcommon']@(t['Qcommon'].T@v);v-=A.T@np.linalg.solve(Gram,A@v);sd=np.sqrt(np.mean(v*v))
  if sd>1e-10:v*=.7*t['noise_scale']/sd
  else:v[:]=0
  assert abs(A@v).max()<1e-8
  out.append({'name':'NULL_PROTOTYPE_'+str(j),'kind':'NULL','mu':v,'target':np.zeros(2)})
 for mag in [.02,.05,.10]:
  for j,d in enumerate([np.array([1.,0]),np.array([0.,1]),np.array([1.,1])/np.sqrt(2),np.array([1.,-1])/np.sqrt(2)]):
   for sign in [-1,1]:
    target=d*mag*sign;mu=A.T@np.linalg.solve(Gram,target);out.append({'name':f'POWER_{mag}_{j}_{sign}','kind':'POWER','mu':mu,'target':target})
 assert len(out)==31
 return out
def worker(task):
 t,scenario,sig=task;seed=int(np.random.SeedSequence([ROOT_SEED,int(hashlib.sha256(t['test_id'].encode()).hexdigest()[:8],16),SCENARIOS.index(scenario)]).generate_state(1,dtype=np.uint64)[0]);f=SIM/(t['test_id']+'__'+scenario+'.parquet');jf=f.with_suffix('.json')
 if f.exists() and jf.exists():
  old=json.loads(jf.read_text(encoding='utf-8'))
  if old.get('run_signature')==sig and old.get('seed')==str(seed) and old.get('file_sha256')==sha(f):return old
 clock=time.perf_counter();rng=np.random.default_rng(seed);ps=patterns(t);means=[(t['A']@p['mu'],np.einsum('gpn,n->gp',t['K'],p['mu'],optimize=True)) for p in ps]
 arrays={p['name']:np.empty(REPS) for p in ps};bad=np.zeros(REPS,bool)
 for start in range(0,REPS,BATCH):
  B=min(BATCH,REPS-start);Y=innovations(rng,t,scenario,B);b0=(t['A']@Y).T;U0=np.einsum('gpn,nb->bgp',t['K'],Y,optimize=True)
  for pp,(bm,um) in zip(ps,means):
   v=moment_stats(b0+bm[None,:],U0+um[None,:,:],t);bad[start:start+B]|=~v['valid'];arrays[pp['name']][start:start+B]=np.where(v['valid'],v['p'],1.)
 pd.DataFrame({'replicate':np.arange(REPS),'numerical_failure':bad,**arrays}).to_parquet(f,index=False,compression='zstd');summary=[]
 for p in ps:
  for name,a in ALPHAS.items():
   k=int((arrays[p['name']]<=a).sum());lo,hi=wilson(k,REPS)
   summary.append({'test_id':t['test_id'],'gene':t['gene'],'model':t['model'],'kind':t['kind'],'scenario':scenario,'pattern':p['name'],'mean_kind':p['kind'],'alpha_name':name,'alpha':a,'reps':REPS,'rejected':k,'rate':k/REPS,'Wilson95_low':lo,'Wilson95_high':hi})
 out={'test_id':t['test_id'],'scenario':scenario,'seed':str(seed),'reps':REPS,'n_failed':int(bad.sum()),'file':f.name,'file_sha256':sha(f),'run_signature':sig,'seconds':time.perf_counter()-clock,'summary':summary};save_json(jf,out);return out
def qualification(audit,full,inv):
 out=[]
 for r in audit.itertuples():
  d=full.loc[full.test_id.eq(r.test_id)&full.mean_kind.eq('NULL')];ii=inv.loc[inv.test_id.eq(r.test_id)]
  if d.empty:out.append({'test_id':r.test_id,'gene':r.gene,'model':r.model,'qualified':False,'reason':'DESIGN_GATE_FAILED'});continue
  nom=d.loc[d.alpha_name.eq('NOMINAL')];tail=d.loc[~d.alpha_name.eq('NOMINAL')];np_=nom.Wilson95_high.le(.06).all();tp=tail.Wilson95_high.le(tail.alpha).all();fail=(ii.n_failed/ii.reps).max()
  out.append({'test_id':r.test_id,'gene':r.gene,'model':r.model,'kind':r.kind,'family':r.family,'qualified':bool(r.eligible and np_ and tp and fail<=.005),'nominal_pass':bool(np_),'tail_pass':bool(tp),'max_null_rate':nom.rate.max(),'max_null_Wilson95_upper':nom.Wilson95_high.max(),'max_tail_Wilson_over_alpha':(tail.Wilson95_high/tail.alpha).max(),'max_failure_fraction':fail,'reason':'FIXED_FINITE_SCENARIO_RULE'})
 return pd.DataFrame(out)
def run():
 assert json.loads((RESULTS/'ENGINE_VALIDATION_SUMMARY.json').read_text(encoding='utf-8'))['failed']==0
 T=pickle.load((PRIVATE/'TEST_TEMPLATES_RESTRICTED.pkl').open('rb'));audit=pd.read_csv(RESULTS/'ACTUAL_DESIGN_AUDIT.csv');targets=[t for t in T.values() if t['eligible']]
 signature=hashlib.sha256((''.join(sha(ROOT/'code'/n) for n in ['config.py','engine.py','calibrate.py'])+sha(PRIVATE/'TEST_TEMPLATES_RESTRICTED.pkl')+sha(ROOT/'STAGE12_PROTOCOL_v1.0.md')).encode()).hexdigest()
 save_json(RESULTS/'CALIBRATION_FREEZE_BEFORE_RUN.json',{'time':now(),'run_signature':signature,'reps':REPS,'patterns':31,'scenarios':SCENARIOS,'alphas':ALPHAS,'new_real_extension_target_coefficients_fitted':False,'code_sha256':{n:sha(ROOT/'code'/n) for n in ['config.py','engine.py','calibrate.py']}})
 grid=[]
 for t in targets:
  for p in patterns(t):grid.append({'test_id':t['test_id'],'pattern':p['name'],'mean_kind':p['kind'],'target_S_minus_N':p['target'][0],'target_T_minus_N':p['target'][1],'mean_h_max_error':float(abs(t['Ah']@p['mu']-p['target']).max())})
 csvout(RESULTS/'CALIBRATION_MEAN_PATTERN_GRID.csv',grid);outputs=[];tasks=[(t,s,signature) for t in targets for s in SCENARIOS]
 with ProcessPoolExecutor(max_workers=WORKERS) as pool:
  futures=[pool.submit(worker,x) for x in tasks]
  for i,f in enumerate(as_completed(futures),1):
   r=f.result();outputs.append(r);log(f'Calibration {i}/{len(tasks)} {r["test_id"]} {r["scenario"]}: {r["seconds"]:.1f}s')
 outputs.sort(key=lambda x:(x['test_id'],x['scenario']));full=pd.DataFrame([r for o in outputs for r in o['summary']]);inv=pd.DataFrame([{k:v for k,v in o.items() if k!='summary'} for o in outputs]);dec=qualification(audit,full,inv)
 csvout(RESULTS/'CALIBRATION_AND_POWER_FULL.csv',full);csvout(RESULTS/'SIMULATION_FILE_INVENTORY.csv',inv);csvout(RESULTS/'CALIBRATION_QUALIFICATION_DECISIONS.csv',dec)
 save_json(RESULTS/'CALIBRATION_RUN_SUMMARY.json',{'time':now(),'tests':len(audit),'eligible_tests':len(targets),'cells':len(tasks),'independent_noise_draws':int(inv.reps.sum()),'pattern_evaluations_reusing_noise':int(inv.reps.sum()*31),'qualified':int(dec.qualified.sum()),'run_signature':signature})
 log(f'Calibration finished {int(dec.qualified.sum())}/{len(audit)} tests qualified.')
 print(dec.to_string(index=False))
if __name__=='__main__':assert_unsealed();run()

