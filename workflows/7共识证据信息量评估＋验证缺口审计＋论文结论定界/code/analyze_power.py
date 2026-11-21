from config import *
from model import wilson
import pickle

def run():
 power=pd.read_csv(RESULTS/'POWER_SIMULATION_SUMMARY.csv');boundary=pd.read_csv(RESULTS/'BOUNDARY_SIMULATION_SUMMARY.csv');checks=[]
 for r in boundary.itertuples():
  lo,hi=wilson(r.raw_success_count,r.reps);checks.append({'accession':r.accession,'gene':r.gene,'k':r.k_independent_panels,'margin':r.margin,'target':r.boundary_target,'sign':r.boundary_sign,
   'threshold':r.threshold,'alpha':r.alpha,'reps':r.reps,'raw_false_equiv_count':r.raw_success_count,'raw_rate':r.raw_success_count/r.reps,'raw_Wilson95_low':lo,'raw_Wilson95_high':hi,
   'gated_false_equiv_count':r.gate_success_count,'gated_rate':r.rate,'gated_Wilson95_high':r.Wilson95_high,'failure_count':r.failure_count,'file':r.file})
 checks=pd.DataFrame(checks);csvout(RESULTS/'REAL_TEMPLATE_BOUNDARY_CHECKS.csv',checks);decisions=[]
 for (acc,k,m),part in checks.loc[checks.threshold.eq('NOMINAL')].groupby(['accession','k','margin']):
  decisions.append({'accession':acc,'k':k,'margin':m,'n_boundary_cells':len(part),'raw_reference_pass':bool(len(part)==4 and part.raw_Wilson95_high.le(.06).all() and (part.failure_count/part.reps).le(.005).all()),
   'max_raw_false_equiv_rate':part.raw_rate.max(),'max_raw_Wilson95_upper':part.raw_Wilson95_high.max(),'max_gated_false_equiv_rate':part.gated_rate.max()})
 decisions=pd.DataFrame(decisions);csvout(RESULTS/'REAL_TEMPLATE_BOUNDARY_DECISIONS.csv',decisions)
 worst=power.groupby(['accession','gene','in_family34','k_independent_panels','template_n','future_n_scenario','gamma_common_P','margin','threshold']).agg(
  worst_scenario_rate=('rate','min'),best_scenario_rate=('rate','max'),worst_Wilson95_lower=('Wilson95_low','min'),worst_gate_pass_rate=('gate_pass_rate','min')).reset_index()
 ok=decisions.groupby(['accession','margin']).raw_reference_pass.all().to_dict();worst['both_endpoint_finite_boundary_checks_passed']=[ok[(r.accession,r.margin)] for r in worst.itertuples()]
 worst['k2_k4_boundary_not_separately_simulated']=worst.k_independent_panels.isin([2,4]);csvout(RESULTS/'POWER_WORST_SCENARIO_INFORMATION.csv',worst)
 requirements=[]
 for (acc,gene,infam,gamma,m,threshold),part in worst.groupby(['accession','gene','in_family34','gamma_common_P','margin','threshold']):
  for target in [.8,.9]:
   passed=part.loc[part.worst_Wilson95_lower.ge(target)&part.both_endpoint_finite_boundary_checks_passed].sort_values('k_independent_panels');r=passed.iloc[0] if len(passed) else None
   requirements.append({'accession':acc,'gene':gene,'in_family34':infam,'gamma_common_P':gamma,'margin':m,'threshold':threshold,'target_probability':target,
    'first_tested_k_with_all3_scenario_Wilson_lower_above_target':r.k_independent_panels if r is not None else np.nan,
    'future_n_scenario_at_that_k':r.future_n_scenario if r is not None else np.nan,'maximum_tested_k':8,
    'status':'REACHED_IN_FIXED_TEMPLATE_GRID' if r is not None else 'NOT_REACHED_THROUGH8_PANELS_OR_ENDPOINT_CALIBRATION_UNPASSED',
    'not_unconditional_sample_size_recommendation':True})
 csvout(RESULTS/'POWER80_90_INDEPENDENT_PANEL_REQUIREMENTS.csv',requirements)
 counts=[]
 for (gamma,m,k,threshold),part in worst.loc[worst.in_family34].groupby(['gamma_common_P','margin','k_independent_panels','threshold']):
  for target in [.8,.9]:counts.append({'gamma_common_P':gamma,'margin':m,'k':k,'threshold':threshold,'target':target,'n_templates':len(part),
   'n_Wilson_lower_meets_target':int(part.worst_Wilson95_lower.ge(target).sum()),'median_worst_scenario_rate':part.worst_scenario_rate.median(),
   'minimum_worst_scenario_rate':part.worst_scenario_rate.min(),'maximum_worst_scenario_rate':part.worst_scenario_rate.max(),
   'boundary_checked_only_SLOPES3_gamma_point55_at_k1_and8':True})
 csvout(RESULTS/'POWER_TARGET_ATTAINMENT_COUNTS.csv',counts)
 # Inspect precision of the previous finite boundary checks at the smaller sufficient threshold.
 tails=[]
 for r in pd.read_csv(S6/'results/EQUIVALENCE_BOUNDARY_CALIBRATION.csv').itertuples():
  p=pd.read_parquet(S6/'results/simulation_replicates'/r.file).p_equiv.to_numpy()
  for label,alpha in [('NOMINAL',.05),('BY_RANK1_SUFFICIENT',ALPHA_BY_SUFFICIENT)]:
   k=int((np.isfinite(p)&(p<=alpha)).sum());lo,hi=wilson(k,len(p))
   tails.append({'accession':r.accession,'margin':r.margin,'threshold':label,'alpha':alpha,'reps':len(p),'false_equiv_count':k,'rate':k/len(p),'Wilson95_low':lo,'Wilson95_high':hi,
    'Wilson_upper_below_alpha':hi<=alpha,'zero_event_upper_if_zero_events':wilson(0,len(p))[1],'source_file':r.file})
 csvout(RESULTS/'STAGE6_BOUNDARY_EXTREME_TAIL_PRECISION.csv',tails)
 summary={'templates':power.accession.nunique(),'family34_templates':int(power[['accession','in_family34']].drop_duplicates().in_family34.sum()),
  'power_independent_cells':power.file.nunique(),'power_independent_repetitions':int(power[['file','reps']].drop_duplicates().reps.sum()),
  'boundary_independent_cells':boundary.file.nunique(),'boundary_independent_repetitions':int(boundary[['file','reps']].drop_duplicates().reps.sum()),
  'nominal_boundary_design_margin_endpoints':len(decisions),'nominal_boundary_passed':int(decisions.raw_reference_pass.sum()),
  'max_raw_nominal_false_equivalence_rate':float(checks.loc[checks.threshold.eq('NOMINAL'),'raw_rate'].max()),
  'max_raw_nominal_Wilson95_upper':float(checks.loc[checks.threshold.eq('NOMINAL'),'raw_Wilson95_high'].max()),
  'alpha_BY_rank1_sufficient':ALPHA_BY_SUFFICIENT,'this_threshold_does_not_equal_actual_BY_power':True,
  'same_cluster_sample_size_not_assumed_to_give_V_over_k':True}
 save_json(RESULTS/'POWER_INFORMATION_SUMMARY.json',summary);log('Power scenarios, endpoint-boundary checks, prospective panel requirements and extreme-tail precision summarized')
 print(summary)
 print(worst.loc[worst.in_family34&worst.gamma_common_P.eq(.55)&worst.margin.eq(.5)].groupby(['k_independent_panels','threshold']).worst_scenario_rate.describe()[['min','50%','max']].to_string())

if __name__=='__main__':run()
