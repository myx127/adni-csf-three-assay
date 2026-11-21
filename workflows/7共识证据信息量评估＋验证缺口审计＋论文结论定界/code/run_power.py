from config import *
from power_engine import setup,generate,expanded_parameters,evaluate
from model import wilson
from concurrent.futures import ProcessPoolExecutor,as_completed
import pickle,time,argparse

def worker(task):
 kind,acc,template,k,scenario,gamma,boundary=task
 label=f'{kind}__{acc}__K{k:02d}__{scenario}__g{int(gamma*100):03d}'
 if boundary:label+=f'__d{int(boundary[0]*100):02d}__M{boundary[1]+1}__s{boundary[2]}'
 seed=int(np.random.SeedSequence([SEED_ROOT,int(hashlib.sha256(label.encode()).hexdigest()[:16],16)]).generate_state(1,dtype=np.uint64)[0]);rng=np.random.default_rng(seed);spec=setup(template,k,scenario);records=[]
 for start in range(0,2000,32):
  B=min(32,2000-start);theta,V=expanded_parameters(generate(spec,B,rng,gamma,boundary),template['d']);p,gate,valid,struct,covstrong,bounded=evaluate(theta,V)
  records.append(pd.DataFrame({'replicate':np.arange(start,start+B),'eligible':gate,'valid':valid,'structure':struct,'covariance_information':covstrong,'bounded_Fieller90':bounded,
   **{'p_equiv_d'+str(int(m*100)):p[:,j] for j,m in enumerate(MARGINS)}}))
 output=pd.concat(records,ignore_index=True);file=label+'.parquet';output.to_parquet(SIM/file,index=False);rows=[]
 for j,m in enumerate(MARGINS):
  if boundary and m!=boundary[0]:continue
  p=output['p_equiv_d'+str(int(m*100))].to_numpy();gate=output.eligible.to_numpy(bool)
  for threshold,alpha in [('NOMINAL',.05),('BY_RANK1_SUFFICIENT',ALPHA_BY_SUFFICIENT)]:
   raw=int((np.isfinite(p)&(p<=alpha)).sum());gated=int((np.isfinite(p)&gate&(p<=alpha)).sum());lo,hi=wilson(gated,len(p))
   rows.append({'kind':kind,'accession':acc,'gene':template['gene'],'in_family34':acc in FAMILY,'k_independent_panels':k,'template_n':template['d']['n'],
    'future_n_scenario':template['d']['n']*k,'scenario':scenario,'gamma_common_P':gamma,'margin':m,'threshold':threshold,'alpha':alpha,'reps':len(p),
    'raw_success_count':raw,'gate_success_count':gated,'rate':gated/len(p),'Wilson95_low':lo,'Wilson95_high':hi,'gate_pass_rate':gate.mean(),'failure_count':int((~output.valid).sum()),
    'boundary_target':boundary[1]+1 if boundary else np.nan,'boundary_sign':boundary[2] if boundary else np.nan,'seed':str(seed),'file':file,
    'not_actual_BY_ranking_power':threshold!='NOMINAL'})
 return rows

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--kind',choices=['POWER','BOUNDARY'],required=True);opt=parser.parse_args();templates=pickle.load((PRIVATE/'POWER_TEMPLATES_RESTRICTED.pkl').open('rb'))
 freeze=RESULTS/'POWER_ENGINE_FROZEN.json';obj={'before_any_Stage7_power_or_boundary_simulation':True,'already_known_Stage6_results':True,
  'engine_sha256':sha(ROOT/'code/power_engine.py'),'runner_sha256':sha(Path(__file__)),'template_sha256':sha(PRIVATE/'POWER_TEMPLATES_RESTRICTED.pkl'),
  'execution_supplement_sha256':sha(ROOT/'POWER_EXECUTION_DETAIL_v1.1.md'),'templates':list(templates),'reps':2000,'k':K_GRID,'gamma':GAMMAS,'scenarios':SCENARIOS,'batch_size':32,'root_seed':SEED_ROOT}
 if freeze.exists():assert json.loads(freeze.read_text(encoding='utf-8'))==obj
 else:save_json(freeze,obj)
 tasks=[]
 if opt.kind=='POWER':tasks=[('POWER',a,t,k,s,g,None) for a,t in templates.items() for k in K_GRID for s in SCENARIOS for g in GAMMAS]
 else:tasks=[('BOUNDARY',a,t,k,'SLOPES3',.55,(m,target,sign)) for a,t in templates.items() for k in [1,8] for m in MARGINS for target in [1,2] for sign in [-1,1]]
 rows=[];start=time.time()
 with ProcessPoolExecutor(max_workers=6) as pool:
  for i,f in enumerate(as_completed([pool.submit(worker,t) for t in tasks]),1):
   rows.extend(f.result())
   if i%12==0 or i==len(tasks):log(f'{opt.kind}: {i}/{len(tasks)} cells; elapsed={time.time()-start:.1f}s')
 csvout(RESULTS/(opt.kind+'_SIMULATION_SUMMARY.csv'),pd.DataFrame(rows).sort_values(['accession','k_independent_panels','scenario','gamma_common_P','margin','threshold','boundary_target','boundary_sign']))
 save_json(RESULTS/(opt.kind+'_RUN.json'),{'cells':len(tasks),'independent_repetitions':len(tasks)*2000,'elapsed_seconds':time.time()-start})
