from __future__ import annotations
from config import *
from model import batch_components,factors,wilson
from concurrent.futures import ProcessPoolExecutor,as_completed
import pickle,time,argparse

NULLS=['C01','C02','C03','C04','C05','C06']
SCENARIO_SPEC={
 'C01':{'label':'IID equal precision','sigma':[.8,.8,.8],'slopes':0.,'axes':[]},
 'C02':{'label':'IID unequal precision','sigma':[.6,.8,1.2],'slopes':0.,'axes':[]},
 'C03':{'label':'Independent heteroskedastic errors','sigma':[.6,.8,1.2],'slopes':0.,'axes':[],'hetero':.35},
 'C04':{'label':'19-plate residual random slopes','sigma':[.6,.8,1.2],'slopes':.6,'axes':[0]},
 'C05':{'label':'Three-way residual random slopes','sigma':[.6,.8,1.2],'slopes':.35,'axes':[0,1,2]},
 'C06':{'label':'t5 errors and three-way slopes','sigma':[.6,.8,1.2],'slopes':.35,'axes':[0,1,2],'t5':True}}

def setup(design,scenario):
 spec=SCENARIO_SPEC.get(scenario,SCENARIO_SPEC['C05']).copy();A=[];q=[]
 for axis in spec['axes']:
  lab=design['groups'][:,axis];one=design['onehots'][axis];v=design['q_age'];center=v-(one@v/one.sum(axis=1))[lab]
  q.append((axis,center));A.append(one.T*center[:,None])
 projected=0.
 if A:
  mat=np.column_stack(A);mat-=design['Q']@(design['Q'].T@mat);projected=float(np.sum(mat*mat)/design['df'])
 spec['q']=q;spec['scale']=np.sqrt(1+np.array(spec['sigma'])**2+2*spec['slopes']**2*projected)
 nuisanceP=design['Qz']@(design['Qz'].T@design['P'])
 spec['mu']=nuisanceP[:,None]*np.array([.1,.9,-.7])+design['q_age'][:,None]*np.array([-.08,.02,.14])
 return spec

def generate(design,spec,B,rng,delta=0.,target=1,stress=None):
 n=design['n'];heavy=spec.get('t5',False)
 def draw(shape):return rng.standard_t(5,size=shape)*np.sqrt(3/5) if heavy else rng.normal(size=shape)
 U=draw((B,n));e=draw((B,n,3))*np.array(spec['sigma'])
 if spec.get('hetero'):
  mult=np.exp(np.clip(spec['hetero']*design['P'],-3,3));mult/=np.sqrt(np.mean(mult*mult));e*=mult[None,:,None]
 if stress=='CORR12':
  base=e[:,:,0]/spec['sigma'][0];e[:,:,1]=spec['sigma'][1]*(.35*base+np.sqrt(1-.35**2)*e[:,:,1]/spec['sigma'][1])
 for axis,center in spec['q']:
  lab=design['groups'][:,axis];G=design['onehots'][axis].shape[0];coeff=draw((B,G,4))*spec['slopes']
  U+=coeff[:,lab,0]*center[None,:]
  e+=coeff[:,lab,1:]*center[None,:,None]
 T=.55*design['P'][None,:]+U
 Y=(T[:,:,None]+e)/spec['scale'][None,None,:]+spec['mu'][None,:,:]
 if stress=='NONLINEAR':Y[:,:,2]+=.25*T*T/spec['scale'][2]
 if stress=='WEAK':Y=(.1*T[:,:,None]+e)/spec['scale'][None,None,:]+spec['mu'][None,:,:]
 if delta:Y[:,:,target]+=delta*design['P'][None,:]
 return Y

def simulate_cell(task):
 ident,design,scenario,reps,phase,delta,target=task
 seed=int(np.random.SeedSequence([SEEDS[phase],int(hashlib.sha256(ident.encode()).hexdigest()[:8],16),NULLS.index(scenario) if scenario in NULLS else 99,int(round(delta*100)),target]).generate_state(1,dtype=np.uint64)[0])
 rng=np.random.default_rng(seed);spec=setup(design,scenario);stats=[];start=time.time()
 for first in range(0,reps,64):
  B=min(64,reps-first);Y=generate(design,spec,B,rng,delta,target)
  stats.extend(batch_components(Y,design)['Q0'].tolist())
 name=f'{phase}__{ident}__{scenario}__d{int(round(delta*100)):03d}__M{target+1}'
 pd.DataFrame({'replicate':np.arange(reps),'Q0':stats}).to_parquet(SIM/(name+'.parquet'),index=False)
 rows=[];q=np.array(stats)
 for version,mult in enumerate(factors(design)):
  p=np.exp(-q/(2*mult));valid=np.isfinite(p)
  for alpha in [.05,.005,.001]:
   k=int(np.sum(valid&(p<=alpha)));lo,hi=wilson(k,reps)
   rows.append({'phase':phase,'design_id':ident,'scenario':scenario,'delta':delta,'target':target+1,'version':version,'alpha':alpha,'reps':reps,
    'n_failed':int((~valid).sum()),'n_rejected':k,'rate':k/reps,'wilson95_low':lo,'wilson95_high':hi,'seed':str(seed),'seconds':time.time()-start,'file':name+'.parquet'})
 return rows

def main(phase,workers=6,only=None):
 designs=pickle.load((RESTRICTED/'DESIGNS_MAIN_RESTRICTED.pkl').open('rb'))
 if only:designs={k:v for k,v in designs.items() if k==only}
 tasks=[]
 if phase in ['DEV','CONFIRM']:
  reps=1000 if phase=='DEV' else 4000
  for ident,d in designs.items():
   for scenario in NULLS:tasks.append((ident,d,scenario,reps,phase,0.,1))
 elif phase=='POWER':
  selected=json.loads((RESULTS/'MAIN_ANALYSIS_SAP_FROZEN.json').read_text(encoding='utf-8'))
  designs={k:v for k,v in designs.items() if k in selected['calibration_passed_including_control']}
  for ident,d in designs.items():
   for target in [1,2]:
    for delta in [.2,.35,.5,.75,1.,1.5]:tasks.append((ident,d,'C05',2000,phase,delta,target))
 engine=RESULTS/'ENGINE_FREEZE_BEFORE_SIMULATION.json'
 hashes={p.name:sha(p) for p in [ROOT/'code/model.py',ROOT/'code/simulation.py',ROOT/'code/config.py']}
 if engine.exists():assert json.loads(engine.read_text(encoding='utf-8'))['code']==hashes,'Frozen engine changed'
 else:save_json(engine,{'code':hashes,'spec':SCENARIO_SPEC,'batch_size':64,'before_first_synthetic_Y':True,'real_protein_Y_read':False})
 rows=[];t0=time.time()
 with ProcessPoolExecutor(max_workers=workers) as executor:
  futures={executor.submit(simulate_cell,t):t[0:3] for t in tasks}
  for num,future in enumerate(as_completed(futures),1):
   result=future.result();rows.extend(result)
   if num%12==0 or num==len(tasks):log(f'{phase} cells {num}/{len(tasks)} elapsed={time.time()-t0:.1f}s')
 csvout(RESULTS/(phase+'_SIMULATION_SUMMARY'+('_'+only if only else '')+'.csv'),rows)
 save_json(RESULTS/(phase+'_RUN_SUMMARY'+('_'+only if only else '')+'.json'),{'phase':phase,'cells':len(tasks),'independent_replicates':sum(t[3] for t in tasks),'seconds':time.time()-t0,'workers':workers,'only':only})

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--phase',choices=['DEV','CONFIRM','POWER'],required=True);p.add_argument('--workers',type=int,default=6);p.add_argument('--only');args=p.parse_args();main(args.phase,args.workers,args.only)
