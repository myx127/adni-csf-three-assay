from config import *
from model import batch_components
from simulation import setup,generate
import pickle
main=pickle.load((RESTRICTED/'DESIGNS_MAIN_RESTRICTED.pkl').open('rb'));extra=pickle.load((RESTRICTED/'DESIGNS_VARIANTS_RESTRICTED.pkl').open('rb'))
rows=[]
for table in ['DEV_SIMULATION_SUMMARY.csv','CONFIRM_SIMULATION_SUMMARY.csv','VARIANT_CALIBRATION_SUMMARY.csv','POWER_SIMULATION_SUMMARY.csv']:
 data=pd.read_csv(RESULTS/table,dtype={'seed':str});data=data.loc[data.version.eq(0)&data.alpha.eq(.05)].sort_values('design_id')
 selected=data.groupby(['target','delta'] if table.startswith('POWER') else ['scenario'],sort=True).head(1)
 for r in selected.itertuples():
  d=main[r.design_id] if r.design_id in main else extra[r.design_id];rng=np.random.default_rng(int(r.seed));spec=setup(d,r.scenario);values=[]
  for first in range(0,r.reps,64):
   y=generate(d,spec,min(64,r.reps-first),rng,r.delta,r.target-1);values.extend(batch_components(y,d)['Q0'])
  saved=pd.read_parquet(SIM/r.file).Q0.to_numpy();ok=np.allclose(values,saved,rtol=1e-11,atol=1e-11,equal_nan=True)
  rows.append({'summary_table':table,'file':r.file,'reps_replayed':r.reps,'max_abs_difference':float(np.nanmax(abs(np.array(values)-saved))),'passed':bool(ok)})
csvout(RESULTS/'SEEDED_REPLAY_VALIDATION.csv',rows)
save_json(RESULTS/'SEEDED_REPLAY_SUMMARY.json',{'cells':len(rows),'replayed_replicates':sum(r['reps_replayed'] for r in rows),'passed':sum(r['passed'] for r in rows),
 'replays_counted_as_new_independent_simulations':False,'saved_simulation_files_not_rewritten':True})
assert all(r['passed'] for r in rows)
log(f'Read-only seeded replay {len(rows)}/{len(rows)} cells passed; not counted as independent new simulations')
