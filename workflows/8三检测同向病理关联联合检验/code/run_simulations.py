from config import *
from engine import setup,generate,fit_noise,evaluate
from model import wilson
from scipy.special import ndtr
from concurrent.futures import ProcessPoolExecutor,as_completed
from freeze_engine import verify
import pickle,time
TEMPLATES=None
def initialize():
 global TEMPLATES
 TEMPLATES=pickle.load((PRIVATE/'SIMULATION_TEMPLATES_RESTRICTED.pkl').open('rb'))
def worker(task):
 key,scenario=task;t=TEMPLATES[key];label=key+'__'+scenario
 seed=int(np.random.SeedSequence([ROOT_SEED,int(hashlib.sha256(label.encode()).hexdigest()[:16],16)]).generate_state(1,dtype=np.uint64)[0])
 file='SIM__'+label+'.parquet';checkpoint=SIM/('SIM__'+label+'.json')
 if checkpoint.exists():
  obj=json.loads(checkpoint.read_text(encoding='utf-8'))
  assert obj['seed']==str(seed) and obj['reps']==REPS and sha(SIM/file)==obj['parquet_sha256']
  return obj
 rng=np.random.default_rng(seed);spec=setup(t,scenario);noise=[];stderr=[];valid=[]
 for start in range(0,REPS,BATCH):
  b,se,ok=fit_noise(generate(spec,min(BATCH,REPS-start),rng),spec);noise.append(b);stderr.append(se);valid.append(ok)
 noise=np.concatenate(noise);stderr=np.concatenate(stderr);valid=np.concatenate(valid)
 frame=pd.DataFrame({'replicate':np.arange(REPS),'valid':valid})
 for j,letter in enumerate(['N','S','T']):frame[letter+'_noise_over_templateSD']=noise[:,j];frame[letter+'_SE_over_templateSD']=stderr[:,j]
 frame.to_parquet(SIM/file,index=False)
 common={'design_id':key,'accession':t['accession'],'gene':t['gene'],'model':t['model'],'scenario':scenario,'in_family34':t['accession'] in FAMILY,'n':t['d']['n'],'reps':REPS,'failure_count':int((~valid).sum()),'seed':str(seed),'file':file}
 rows=[]
 for pattern in PATTERNS:
  _,_,p,direction=evaluate(noise,stderr,pattern['effect'],valid)
  for threshold,alpha in ALPHAS.items():
   success=valid&np.isfinite(p)&(p<=alpha);count=int(success.sum());lo,hi=wilson(count,REPS)
   rows.append({**common,**{k:v for k,v in pattern.items() if k!='effect'},'bN_true_scaled':pattern['effect'][0],'bS_true_scaled':pattern['effect'][1],'bT_true_scaled':pattern['effect'][2],
    'threshold':threshold,'alpha':alpha,'success_count':count,'rate':count/REPS,'Wilson95_low':lo,'Wilson95_high':hi,
    'positive_rejection_count':int((success&(direction==1)).sum()),'negative_rejection_count':int((success&(direction==-1)).sum()),
    'wrong_direction_count':int((success&(direction!=pattern['sign'])).sum()) if pattern['kind']=='POWER' else np.nan})
 components=[]
 for j,letter in enumerate(['N','S','T']):
  z=noise[:,j]/stderr[:,j]
  for sign in [-1,1]:
   p=ndtr(-sign*z);count=int((valid&(p<=.025)).sum());lo,hi=wilson(count,REPS)
   components.append({**common,'assay':letter,'direction':sign,'one_sided_alpha':.025,'reject_count':count,'rate':count/REPS,'Wilson95_low':lo,'Wilson95_high':hi})
 obj={'seed':str(seed),'reps':REPS,'file':file,'parquet_sha256':sha(SIM/file),'rows':rows,'components':components,'normalizer':spec['normalizer']}
 save_json(checkpoint,obj);return obj
if __name__=='__main__':
 verify();initialize();tasks=[(key,scenario) for key in TEMPLATES for scenario in SCENARIOS];assert len(tasks)==525
 started=time.time();rows=[];components=[];inventory=[]
 with ProcessPoolExecutor(max_workers=WORKERS,initializer=initialize) as pool:
  for i,future in enumerate(as_completed([pool.submit(worker,t) for t in tasks]),1):
   obj=future.result();rows.extend(obj['rows']);components.extend(obj['components']);inventory.append({'file':obj['file'],'seed':obj['seed'],'reps':obj['reps'],'sha256':obj['parquet_sha256']})
   if i%15==0 or i==len(tasks):log(f'Noise cells {i}/{len(tasks)}; elapsed {time.time()-started:.1f}s; patterns share each independent noise repeat')
 frame=pd.DataFrame(rows).sort_values(['design_id','scenario','kind','pattern','threshold'])
 csvout(RESULTS/'NULL_CALIBRATION_SUMMARY.csv',frame.loc[frame.kind.eq('NULL')]);csvout(RESULTS/'POWER_SCENARIO_SUMMARY.csv',frame.loc[frame.kind.eq('POWER')])
 csvout(RESULTS/'COMPONENT_NULL_CALIBRATION.csv',pd.DataFrame(components).sort_values(['design_id','scenario','assay','direction']))
 csvout(RESULTS/'SIMULATION_FILE_INVENTORY.csv',pd.DataFrame(inventory).sort_values('file'))
 save_json(RESULTS/'SIMULATION_RUN.json',{'independent_noise_cells':len(tasks),'independent_noise_repetitions':len(tasks)*REPS,'null_patterns':25,'power_patterns':28,'pattern_evaluations_not_independent_repetitions':len(tasks)*REPS*len(PATTERNS),'elapsed_seconds':time.time()-started})
