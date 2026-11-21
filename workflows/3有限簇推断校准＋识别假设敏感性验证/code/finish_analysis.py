from __future__ import annotations
import json,math
import numpy as np
import pandas as pd
from cal_config import *
from mean_followup import followup_grid,summarize
from analyze_calibration import aggregate,wilson

def confirmation_gates(summary,phase,record,dimension):
 out=[]
 scopes={'ONE_WAY_SCOPE':['C01','C02','C03','C04','C05','C08','C09','C10'],'CROSSED_SCOPE':[f'C{i:02d}' for i in range(1,11)]}
 for decision in record['decisions']:
  method=decision['selected_method'];scope=decision['scope'];cells=scopes[scope]
  q=summary.loc[summary.phase.eq(phase)&summary.method.eq(method)&summary.scenario_id.isin(cells)] if method else summary.iloc[:0]
  power=summary.loc[summary.phase.eq(phase)&summary.method.eq(method)&summary.scenario_id.eq('A04')] if method else summary.iloc[:0]
  passed=method is not None and len(q)==len(cells) and q.wilson95_upper.le(.075).all() and q.failure_rate.le(.01).all()
  adequate=bool(len(power)==1 and power.wilson95_lower.iloc[0]>=.80)
  out.append({'dimension_df':dimension,'scope':scope,'development_selected_method':method,'confirmation_null_gate_pass':bool(passed),
      'maximum_null_rate':float(q.rejection_rate.max()) if len(q) else np.nan,
      'maximum_null_upper95':float(q.wilson95_upper.max()) if len(q) else np.nan,
      'maximum_failure':float(q.failure_rate.max()) if len(q) else np.nan,
      'strong_A04_power':float(power.rejection_rate.iloc[0]) if len(power) else np.nan,
      'strong_A04_lower95':float(power.wilson95_lower.iloc[0]) if len(power) else np.nan,
      'strong_power_gate_pass':adequate,'control_and_strong_power_pass':bool(passed and adequate),
      'decision':'CONDITIONAL_CANDIDATE' if passed and adequate else 'CONTROL_ONLY_POWER_INSUFFICIENT' if passed else 'NOT_CONFIRMED'})
 return out

def sensitivity_summary():
 rows=[];profile=[]
 for dimension,prefix in [(8,'SENSITIVITY'),(2,'MEAN_SENS_REPLAY')]:
  for cell in [c for c in grid() if c['phase']=='SENSITIVITY']:
   name=cell['cell_id'].replace('SENSITIVITY',prefix)
   frame=pd.read_parquet(RESULTS/'replicates'/(name+'.parquet'))
   for method,q in frame.groupby('method',observed=True):
    ok=q.status.eq('OK');k=int((ok&q.p_value.le(.05)).sum());lo,hi=wilson(k,len(q))
    rows.append({'dimension_df':dimension,'scenario_id':cell['scenario_id'],'method':str(method),'n_replicates':len(q),
       'n_rejected':k,'n_failed':int((~ok).sum()),'rejection_rate':k/len(q),'wilson95_lower':lo,'wilson95_upper':hi,
       'failure_rate':float((~ok).mean()),'n_complete_median':float(q.n_complete.median())})
   q=pd.read_parquet(RESULTS/'replicates'/(name+'_sensitivity.parquet'));N=len(q)
   for method in ['CR1_HT_F','SUM3_CHI2']:
    base=q[method+'_zero_p'];envelope=q[method+'_grid_max_p'];truth=q[method+'_true_fixed_kappa_p']
    profile.append({'dimension_df':dimension,'scenario_id':cell['scenario_id'],'method':method,'n_replicates':N,
       'zero_assumption_rejection':float(base.le(.05).mean()),'grid_envelope_rejection':float(envelope.le(.05).mean()),
       'classification_flip_rate':float(q[method+'_classification_flips'].fillna(False).mean()),
       'valid_grid_fraction':float(envelope.notna().mean()),'known_true_fixed_kappa_rejection':float(truth.le(.05).mean()) if truth.notna().all() else np.nan,
       'known_true_kappa_available_fraction':float(truth.notna().mean()),
       'oracle_remove_quadratic_rejection':float(q.ORACLE_REMOVE_QUADRATIC_CR1_HT_F_p.le(.05).mean())})
 csv(RESULTS/'SENSITIVITY_METHOD_SUMMARY.csv',rows);csv(RESULTS/'ASSUMPTION_GRID_ENVELOPE_SUMMARY.csv',profile)

def main():
 summary=aggregate();summarize(False)
 mean=pd.read_csv(RESULTS/'MEAN_FOLLOWUP_SUMMARY.csv')
 original=json.loads((RESULTS/'DEVELOPMENT_SELECTION_FROZEN.json').read_text(encoding='utf-8'))
 supplement=json.loads((RESULTS/'MEAN_SELECTION_FROZEN.json').read_text(encoding='utf-8'))
 gates=confirmation_gates(summary,'CONFIRMATION',original,8)+confirmation_gates(mean,'MEAN_CONFIRM',supplement,2)
 csv(RESULTS/'INDEPENDENT_CONFIRMATION_GATES.csv',gates)
 tables=[]
 for dimension,frame,phase in [(8,summary,'CONFIRMATION'),(2,mean,'MEAN_CONFIRM')]:
  q=frame.loc[frame.phase.eq(phase)].copy();q['dimension_df']=dimension;tables.append(q)
 confirmation=pd.concat(tables,ignore_index=True);csv(RESULTS/'CONFIRMATION_METHOD_COMPARISON.csv',confirmation)
 sensitivity_summary()
 save_json(RESULTS/'ANALYSIS_SUMMARY.json',{'primary_independent_replicates':74000,'mean_followup_independent_replicates':42000,
       'total_independent_replicates':116000,'mean_sensitivity_replays_not_new':24000,'resolution_replays_not_new':2000,
       'primary_method_rows':592000,'followup_method_rows':336000,'mean_sensitivity_method_rows':192000,'resolution_method_rows':16000,
       'all_method_rows_including_replays':1136000,'clinical_files_opened':0,'formal_SAP_frozen':False,'independent_confirmation':gates})
 print(pd.DataFrame(gates).to_string(index=False))

if __name__=='__main__':main()
