from config import *
def run():
 null=pd.read_csv(RESULTS/'NULL_CALIBRATION_SUMMARY.csv');component=pd.read_csv(RESULTS/'COMPONENT_NULL_CALIBRATION.csv');decisions=[]
 for key,part in null.groupby('design_id'):
  nominal=part.loc[part.threshold.eq('NOMINAL')];tail=part.loc[~part.threshold.eq('NOMINAL')];cp=component.loc[component.design_id.eq(key)]
  nominal_pass=bool(len(nominal)==125 and nominal.Wilson95_high.le(.06).all())
  marginal_pass=bool(len(cp)==30 and cp.Wilson95_high.le(.03).all())
  tailpass=bool(len(tail)==250 and tail.Wilson95_high.le(tail.alpha).all())
  failure_pass=bool((part.failure_count/part.reps).le(.005).all())
  first=part.iloc[0]
  decisions.append({'design_id':key,'accession':first.accession,'gene':first.gene,'model':first.model,'in_family34':first.in_family34,'nominal_pass':nominal_pass,
   'marginal_one_sided_pass':marginal_pass,'two_tail_thresholds_finite_grid_pass':tailpass,'failure_pass':failure_pass,'qualified':nominal_pass and marginal_pass and tailpass and failure_pass,
   'max_nominal_false_rate':nominal.rate.max(),'max_nominal_Wilson95_high':nominal.Wilson95_high.max(),'max_component_Wilson95_high':cp.Wilson95_high.max(),
   'max_tail_Wilson_to_alpha':(tail.Wilson95_high/tail.alpha).max(),'n_null_evaluations_not_independent_cells':len(part)})
 decisions=pd.DataFrame(decisions);csvout(RESULTS/'CALIBRATION_DECISIONS.csv',decisions)
 summary=decisions.groupby('model').agg(designs=('qualified','size'),qualified=('qualified','sum'),max_nominal_rate=('max_nominal_false_rate','max'),max_nominal_Wilson_upper=('max_nominal_Wilson95_high','max'),max_tail_upper_to_alpha=('max_tail_Wilson_to_alpha','max')).reset_index()
 csvout(RESULTS/'CALIBRATION_MODEL_SUMMARY.csv',summary)
 power=pd.read_csv(RESULTS/'POWER_SCENARIO_SUMMARY.csv')
 worst=power.groupby(['design_id','accession','gene','model','in_family34','pattern','weak_axis','magnitude','sign','threshold']).agg(worst_DGP_rate=('rate','min'),worst_Wilson95_lower=('Wilson95_low','min'),best_DGP_rate=('rate','max'),max_wrong_direction_count=('wrong_direction_count','max')).reset_index()
 worst=worst.merge(decisions[['design_id','qualified']],on='design_id');csvout(RESULTS/'POWER_WORST_DGP_SUMMARY.csv',worst)
 balanced=worst.loc[worst.pattern.str.startswith('BALANCED')&worst.in_family34]
 csvout(RESULTS/'BALANCED_POWER_FAMILY_SUMMARY.csv',balanced.groupby(['model','magnitude','sign','threshold']).agg(n_templates=('gene','size'),median_worst_rate=('worst_DGP_rate','median'),minimum_worst_rate=('worst_DGP_rate','min'),maximum_worst_rate=('worst_DGP_rate','max'),n_Wilson_lower_atleast80=('worst_Wilson95_lower',lambda x:int(x.ge(.8).sum()))).reset_index())
 save_json(RESULTS/'CALIBRATION_SCOPE.json',{'strict_tail_gate':True,'tail_thresholds':ALPHAS,'Wilson_intervals_are_cellwise_not_simultaneous':True,'finite_DGP_grid_not_universal_P_superuniformity_proof':True,'family_size_unchanged':34,'models':summary.to_dict('records')})
 log('All finite-grid calibration gates evaluated without method replacement');print(summary.to_string(index=False))
if __name__=='__main__':run()
