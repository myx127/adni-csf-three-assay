from config import *
from equivalence import equivalence_batch,MARGINS
from simulation_reference import setup,generate
from model import wilson
from concurrent.futures import ProcessPoolExecutor,as_completed
import pickle

def worker(task):
 acc,d,margin,target,sign=task;seed=int(np.random.SeedSequence([SEEDS['CONFIRM'],int(hashlib.sha256(acc.encode()).hexdigest()[:8],16),int(margin*100),target,sign+1]).generate_state(1,dtype=np.uint64)[0])
 rng=np.random.default_rng(seed);spec=setup(d,'C05');ps=[]
 for i in range(0,2000,64):
  Y=generate(d,spec,min(64,2000-i),rng,margin*sign,target);ps.extend(equivalence_batch(Y,d,margin))
 p=np.array(ps);k=int((p<=.05).sum());lo,hi=wilson(k,len(p));name=f'EQUIV_BOUNDARY__{acc}__d{int(margin*100)}__M{target+1}__s{sign}.parquet'
 pd.DataFrame({'replicate':np.arange(len(p)),'p_equiv':p}).to_parquet(SIM/name,index=False)
 return {'accession':acc,'margin':margin,'boundary_target':target+1,'boundary_sign':sign,'reps':len(p),'false_equivalence_count':k,
  'rate':k/len(p),'Wilson95_low':lo,'Wilson95_high':hi,'seed':str(seed),'scenario':'C05 three-way slopes','file':name,'qualified_gate_not_applied_upper_bound_check':True}

if __name__=='__main__':
 ds=pickle.load((S5/'00_restricted/DESIGNS_MAIN_RESTRICTED.pkl').open('rb'));tasks=[(acc,d,m,t,s) for acc,d in ds.items() for m in MARGINS for t in [1,2] for s in [-1,1]]
 save_json(RESULTS/'EQUIVALENCE_PLAN_FREEZE.json',{'margins':MARGINS,'reference_margin':.2,'post_S5_exploratory_working_margins_not_biological_tolerance':True,
  'before_stage6_equivalence_results_computed':True,'boundary_simulation_reps':2000,'cells':len(tasks),'variance_version':0,
  'code_sha256':{p.name:sha(p) for p in [ROOT/'code/equivalence.py',Path(__file__)]}})
 rows=[]
 with ProcessPoolExecutor(max_workers=6) as pool:
  for i,f in enumerate(as_completed([pool.submit(worker,t) for t in tasks]),1):
   rows.append(f.result())
   if i%42==0:log(f'Equivalence boundary calibration {i}/{len(tasks)}')
 csvout(RESULTS/'EQUIVALENCE_BOUNDARY_CALIBRATION.csv',pd.DataFrame(rows).sort_values(['accession','margin','boundary_target','boundary_sign']))
 save_json(RESULTS/'EQUIVALENCE_CALIBRATION_RUN.json',{'independent_replicates':len(tasks)*2000,'cells':len(tasks),'definition':'IUT of four one-sided boundary tests; tests small relative response, not true abundance consensus'})
