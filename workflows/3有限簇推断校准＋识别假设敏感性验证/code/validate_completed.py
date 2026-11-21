from __future__ import annotations
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
import json
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from cal_config import *
from cal_data import generate
from cal_inference import infer
from mean_followup import followup_grid,mean_infer

def main():
 checks=[]
 def add(name,passed,detail):checks.append({'check':name,'passed':bool(passed),'detail':detail})
 for name in ['PREFLIGHT_VALIDATION.json','MEAN_PREFLIGHT_VALIDATION.json']:
  record=json.loads((RESULTS/name).read_text(encoding='utf-8'));add(name,record['status']=='PASS',record)
 access=json.loads((RESULTS/'ACCESS_AND_PROTOCOL_RECORD.json').read_text(encoding='utf-8'))
 add('primary_protocol_unchanged',sha(ROOT/'CALIBRATION_PROTOCOL_v1.0.md')==access['protocol_sha256_before_experiments'],access['protocol_sha256_before_experiments'])
 add('primary_core_unchanged',all(sha(ROOT/'code'/n)==v for n,v in access['core_code_before_experiments'].items()),list(access['core_code_before_experiments']))
 mf=json.loads((RESULTS/'MEAN_SELECTION_FROZEN.json').read_text(encoding='utf-8'))
 add('mean_followup_protocol_and_core_unchanged',sha(ROOT/'LOW_DIMENSIONAL_FOLLOWUP_PROTOCOL_v1.0.md')==mf['protocol_sha256'] and all(sha(ROOT/'code'/n)==v for n,v in mf['core'].items()),list(mf['core']))
 sources=pd.read_csv(RESULTS/'INPUT_SOURCE_HASH_CHECKS.csv')
 add('all_source_inputs_unchanged',all(sha(r.path)==r.expected_sha256 for r in sources.itertuples()),len(sources))
 for selection_name,prefix in [('DEVELOPMENT_SELECTION_FROZEN.json','CONFIRMATION_'),('MEAN_SELECTION_FROZEN.json','MEAN_CONFIRM_')]:
  times=[p.stat().st_mtime for p in (RESULTS/'replicates').glob(prefix+'*.parquet') if not p.name.endswith('_sensitivity.parquet')]
  add('selection_precedes_data_'+prefix,len(times)==14 and min(times)>(RESULTS/selection_name).stat().st_mtime,{'confirmation_files':len(times),'selection_file':selection_name})
 revalidation=json.loads((RESULTS/'MEAN_DEVELOPMENT_REVALIDATION_QA.json').read_text(encoding='utf-8'))
 add('mean_interface_fix_did_not_change14_development_outputs',revalidation['all14_unchanged'],revalidation)
 audited=json.loads((RESULTS/'AUDITED_BATCH_SIZES.json').read_text(encoding='utf-8'))
 cells=[{**c,'dimension':8,'new_independent':True} for c in grid()]+[{**c,'dimension':2,'new_independent':True} for c in followup_grid()]
 cells +=[{**c,'cell_id':c['cell_id'].replace('SENSITIVITY','MEAN_SENS_REPLAY'),'dimension':2,'new_independent':False} for c in grid() if c['phase']=='SENSITIVITY']
 cellchecks=[];total=0;failcounts={};independent=0
 for cell in cells:
  frame=pd.read_parquet(RESULTS/'replicates'/(cell['cell_id']+'.parquet'));N=cell['repetitions'];ok=frame.status.eq('OK')
  passed=len(frame)==N*8 and not frame.duplicated(['replicate','method']).any() and set(frame.method.astype(str))==set(METHODS)
  passed=passed and frame.groupby('method',observed=True).size().eq(N).all() and frame.loc[ok,'p_value'].between(0,1).all() and frame.loc[~ok,'p_value'].isna().all()
  bootstrap=frame.loc[ok&frame.method.isin(['WILD_SCORE_T','PIGEONHOLE_CENTERED'])]
  passed=passed and np.allclose(bootstrap.p_value,(1+bootstrap.bootstrap_tail_count)/(cell['bootstrap_B']+1),atol=1e-14)
  cellchecks.append({'cell_id':cell['cell_id'],'passed':bool(passed),'rows':len(frame)})
  total+=len(frame);independent+=N if cell['new_independent'] else 0
  for status,k in frame.loc[~ok,'status'].value_counts().items():
   if k:failcounts[str(status)]=failcounts.get(str(status),0)+int(k)
  if (cell['phase']=='CONFIRMATION' and cell['scenario_id'] in ['C02','C05','A01']) or (cell['phase']=='MEAN_CONFIRM' and cell['scenario_id'] in ['C02','C08','A04']):
   case={s['scenario_id']:s for s in calibration_scenarios()}[cell['scenario_id']]
   rep=7;data=generate(cell['n'],np.random.default_rng(np.random.SeedSequence([int(cell['seed']),rep,0])),case['settings'],audited)
   rng=np.random.default_rng(np.random.SeedSequence([int(cell['seed']),rep,1 if cell['dimension']==8 else 2]))
   out=infer(data,cell['cluster_axis'],cell['bootstrap_B'],rng)[0] if cell['dimension']==8 else mean_infer(data,cell['cluster_axis'],cell['bootstrap_B'],rng)
   stored=frame.loc[frame.replicate.eq(rep)].set_index('method')
   for method,result in out.items():
    before=stored.loc[method];same=str(before.status)==result['status'] and (np.isnan(before.p_value) and np.isnan(result['p_value']) or abs(before.p_value-result['p_value'])<1e-12)
    add(f'seed_replay_{cell["cell_id"]}_{method}',same,'Full generator/inference replay; not independent scientific validation')
 add('all112_primary_followup_and_sensitivity_method_cells',len(cellchecks)==112 and all(c['passed'] for c in cellchecks),cellchecks)
 add('116000_independent_and1120000_method_rows_before_resolution',independent==116000 and total==1120000,{'independent':independent,'method_rows':total})
 envelope_checks=[]
 for prefix in ['SENSITIVITY','MEAN_SENS_REPLAY']:
  for p in (RESULTS/'replicates').glob(prefix+'*_sensitivity.parquet'):
   q=pd.read_parquet(p)
   for method in ['CR1_HT_F','SUM3_CHI2']:
    zero=q[method+'_zero_p'];maximum=q[method+'_grid_max_p'];valid=maximum.notna()&zero.notna()
    truth=q[method+'_true_fixed_kappa_p'];vtruth=maximum.notna()&truth.notna()
    envelope_checks.append((maximum[valid]>=zero[valid]-1e-12).all() and (maximum[vtruth]>=truth[vtruth]-1e-12).all())
 add('all_discrete_nuisance_envelopes_include_zero_and_known_true_grid_points',len(envelope_checks)==96 and all(envelope_checks),len(envelope_checks))
 mean_summary=pd.read_csv(RESULTS/'MEAN_FOLLOWUP_SUMMARY.csv');summary=pd.read_csv(RESULTS/'ALL_EXPERIMENT_SUMMARY.csv')
 summary_ok=[];wilson_error=0.
 for frame in [summary,mean_summary]:
  for row in frame.itertuples():
   raw=pd.read_parquet(RESULTS/'replicates'/(row.cell_id+'.parquet'),columns=['method','status','p_value'])
   raw=raw.loc[raw.method.eq(row.method)];k=int((raw.status.eq('OK')&raw.p_value.le(.05)).sum())
   summary_ok.append(k==row.n_rejected and len(raw)==row.n_replicates and int(raw.status.ne('OK').sum())==row.n_failed)
  # Independently solve the Wilson score equation, rather than call production interval code.
  for row in frame.iloc[::31].itertuples():
   phat=row.n_rejected/row.n_replicates;z=1.959963984540054
   def equation(p):return row.n_replicates*(phat-p)**2-z*z*p*(1-p)
   if phat==0:lower,upper=0.,z*z/(row.n_replicates+z*z)
   elif phat==1:lower,upper=row.n_replicates/(row.n_replicates+z*z),1.
   else:lower,upper=brentq(equation,0,phat),brentq(equation,phat,1)
   wilson_error=max(wilson_error,abs(lower-row.wilson95_lower),abs(upper-row.wilson95_upper))
 add('all704_main_and_mean_summary_counts_match_parquet',len(summary_ok)==704 and all(summary_ok),len(summary_ok))
 add('Wilson_intervals_independent_score_equation_roots',wilson_error<1e-8,wilson_error)
 resolution=pd.read_parquet(RESULTS/'BOOTSTRAP_RESOLUTION_REPLAY.parquet')
 analytic=resolution.loc[~resolution.method.isin(['WILD_SCORE_T','PIGEONHOLE_CENTERED'])]
 same=(analytic.original_status==analytic.replay_status)&(np.isclose(analytic.original_p,analytic.replay_p,atol=1e-12,rtol=1e-12)|analytic.original_p.isna()&analytic.replay_p.isna())
 add('bootstrap_resolution_reuses_data_and_preserves_analytic_tests',len(resolution)==16000 and same.all(),{'rows':len(resolution),'analytic_rows':len(analytic)})
 seeds=[c['seed'] for c in grid()]+[c['seed'] for c in followup_grid()]
 add('all88_independent_cell_seeds_disjoint',len(seeds)==88 and len(set(seeds))==88,len(seeds))
 result={'status':'PASS' if all(c['passed'] for c in checks) else 'FAIL','passed':sum(c['passed'] for c in checks),'total':len(checks),
     'clinical_files_opened':0,'individual_protein_files_opened':0,'expected_method_failures':failcounts,'checks':checks,
     'interpretation':'Implementation, plan and reproducibility QA; does not turn failed calibration or structural assumptions into valid inference.'}
 save_json(RESULTS/'COMPLETED_VALIDATION_QA.json',result);print({k:v for k,v in result.items() if k!='checks'})
 if result['status']!='PASS':raise RuntimeError('Completed validation failed')

if __name__=='__main__':main()
