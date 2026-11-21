from __future__ import annotations
import argparse,json,math
import numpy as np
import pandas as pd
from datetime import datetime
from cal_config import *

def wilson(k,n):
 if not n:return np.nan,np.nan
 p=k/n;z=1.959963984540054;den=1+z*z/n
 mid=(p+z*z/(2*n))/den;delta=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 return mid-delta,mid+delta

def aggregate(phase=None):
 rows=[];states=[]
 lookup={r['scenario_id']:r for r in calibration_scenarios()}
 sens={r['scenario_id']:r for r in sensitivity_scenarios()}
 for cell in grid():
  if phase is not None and cell['phase']!=phase:continue
  path=RESULTS/'replicates'/(cell['cell_id']+'.parquet')
  if not path.exists():continue
  frame=pd.read_parquet(path);case=(sens if cell['phase']=='SENSITIVITY' else lookup)[cell['scenario_id']]
  for method,q in frame.groupby('method',observed=True):
   ok=q.status.eq('OK');reject=ok&q.p_value.le(.05);k=int(reject.sum());n=len(q);lo,hi=wilson(k,n)
   rows.append({**cell,'method':str(method),'label_zh':case['label_zh'],'compatible_null':case['compatible_null'],
      'n_replicates':n,'n_rejected':k,'n_failed':int((~ok).sum()),'rejection_rate':k/n,
      'wilson95_lower':lo,'wilson95_upper':hi,'failure_rate':float((~ok).mean()),
      'rejection_successful_only':k/int(ok.sum()) if ok.any() else np.nan,
      'mc_se':math.sqrt(k/n*(1-k/n)/n),'n_complete_median':float(q.n_complete.median()),
      'smallest_group_n_median':float(q.smallest_group_n.median())})
   for status,count in q.status.value_counts().items():
    if count:states.append({'cell_id':cell['cell_id'],'method':str(method),'status':str(status),'count':int(count)})
 summary=pd.DataFrame(rows)
 csv(RESULTS/('DEVELOPMENT_SUMMARY.csv' if phase=='DEVELOPMENT' else 'ALL_EXPERIMENT_SUMMARY.csv'),summary)
 csv(RESULTS/('DEVELOPMENT_STATUS_COUNTS.csv' if phase=='DEVELOPMENT' else 'METHOD_STATUS_COUNTS.csv'),states)
 return summary

def select_development():
 summary=aggregate('DEVELOPMENT');decisions=[];scores=[]
 scopes={'ONE_WAY_SCOPE':['C01','C02','C03','C04','C05','C08','C09','C10'],
         'CROSSED_SCOPE':[f'C{i:02d}' for i in range(1,11)]}
 for scope,cases in scopes.items():
  eligible=[]
  for method in NEW_METHODS:
   q=summary.loc[summary.scenario_id.isin(cases)&summary.method.eq(method)]
   passed=len(q)==len(cases) and q.wilson95_upper.le(.080).all() and q.failure_rate.le(.01).all()
   power=summary.loc[summary.scenario_id.isin(['A01','A02'])&summary.method.eq(method),'rejection_rate'].mean()
   scores.append({'scope':scope,'method':method,'passed_all_development_nulls':bool(passed),'maximum_null_rate':float(q.rejection_rate.max()),
       'maximum_null_upper95':float(q.wilson95_upper.max()),'maximum_failure':float(q.failure_rate.max()),'mean_A01_A02_power':float(power)})
   if passed:eligible.append((float(power),-NEW_METHODS.index(method),method))
  chosen=max(eligible)[2] if eligible else None
  decisions.append({'scope':scope,'selected_method':chosen,'eligible_methods':[r[2] for r in sorted(eligible,reverse=True)],
                    'selection_rule':'All predeclared nulls upper Wilson <=.080 and failure<=.01; then A01/A02 average power'})
 csv(RESULTS/'DEVELOPMENT_SELECTION_SCORES.csv',scores)
 access=json.loads((RESULTS/'ACCESS_AND_PROTOCOL_RECORD.json').read_text(encoding='utf-8'))
 save_json(RESULTS/'DEVELOPMENT_SELECTION_FROZEN.json',{'frozen_at_local':datetime.now().astimezone().isoformat(),
    'decisions':decisions,'core_code_sha256':access['core_code_before_experiments'],
    'protocol_sha256':access['protocol_sha256_before_experiments'],'confirmation_data_generated_before_selection':False})
 print(pd.DataFrame(scores).to_string(index=False));print(json.dumps(decisions,ensure_ascii=False))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--select-development',action='store_true');args=p.parse_args()
 if args.select_development:select_development()
 else:print(aggregate().shape)
