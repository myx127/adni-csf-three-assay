from config import *

def freeze():
 version=json.loads((RESULTS/'DEVELOPMENT_SELECTION_FROZEN.json').read_text(encoding='utf-8'))['selected_version']
 if version is None:version=2
 main=pd.read_csv(RESULTS/'CONFIRM_SIMULATION_SUMMARY.csv');extra=pd.read_csv(RESULTS/'VARIANT_CALIBRATION_SUMMARY.csv')
 lookup=pd.read_csv(RESULTS/'VARIANT_DESIGN_SCREEN.csv').fillna('');rows=[];families=[]
 for r in lookup.to_dict('records'):
  ident=r['calibration_design_id'];q=pd.DataFrame();ok=False
  if ident.startswith('MAIN__'):q=main.loc[main.design_id.eq(ident[6:])&main.version.eq(version)&main.alpha.eq(.05)]
  elif ident:q=extra.loc[extra.design_id.eq(ident)&extra.version.eq(version)&extra.alpha.eq(.05)]
  if len(q)==6:ok=bool((q.wilson95_high<=.06).all() and (q.n_failed/q.reps<=.005).all())
  rows.append({**r,'variance_version':version,'calibration_passed':ok,'max_null_rate':float(q.rate.max()) if len(q) else np.nan,
    'max_null_wilson95_upper':float(q.wilson95_high.max()) if len(q) else np.nan,'max_failure_rate':float((q.n_failed/q.reps).max()) if len(q) else np.nan})
  if ok and r['variant']=='AMY' and r['role']=='PRIMARY_DISCOVERY_CANDIDATE':families.append(r['accession'])
 csvout(RESULTS/'VARIANT_CALIBRATION_DECISIONS.csv',rows)
 save_json(RESULTS/'SECONDARY_SAP_FROZEN.json',{'variance_version':version,'amyloid_primary_family':sorted(families),'m_amyloid':len(families),
   'before_real_protein_Y_access':True,'fixed_window_days':90,'primary_adjustment':'BY .05 separate secondary family',
   'raw_P_raw_TMT_LOD_Batch9_fullcov':'supportive sensitivity, no primary family reselection','decisions_sha256':sha(RESULTS/'VARIANT_CALIBRATION_DECISIONS.csv')})
 log(f'Variant calibration frozen before actual protein effects: {sum(r["calibration_passed"] for r in rows)}/{len(rows)} design rows pass; amyloid family {len(families)}')

if __name__=='__main__':freeze()
