from config import *
import argparse

def select():
 data=pd.read_csv(RESULTS/'DEV_SIMULATION_SUMMARY.csv');rows=[];selected=None
 for version in [0,1,2]:
  q=data.loc[data.version.eq(version)&data.alpha.eq(.05)]
  ok=(q.wilson95_high<=.06)&(q.n_failed/q.reps<=.005)
  passed=bool(len(q)==35*6 and ok.all())
  rows.append({'version':version,'all_design_scenario_cells_passed':passed,'cells':len(q),'n_failed_cells':int((~ok).sum()),'max_rate':float(q.rate.max()),'max_wilson95_high':float(q.wilson95_high.max())})
  if passed and selected is None:selected=version
 csvout(RESULTS/'DEVELOPMENT_VARIANCE_SELECTION.csv',rows)
 save_json(RESULTS/'DEVELOPMENT_SELECTION_FROZEN.json',{'selected_version':selected,'selection_order':[0,1,2],
  'dev_summary_sha256':sha(RESULTS/'DEV_SIMULATION_SUMMARY.csv'),'engine_freeze_sha256':sha(RESULTS/'ENGINE_FREEZE_BEFORE_SIMULATION.json'),
  'confirmation_generated_before_selection':False,'real_protein_values_read':False})
 print(pd.DataFrame(rows).to_string(index=False));print('SELECTED',selected)

def freeze_sap():
 selected=json.loads((RESULTS/'DEVELOPMENT_SELECTION_FROZEN.json').read_text(encoding='utf-8'))['selected_version']
 data=pd.read_csv(RESULTS/'CONFIRM_SIMULATION_SUMMARY.csv');screen=pd.read_csv(RESULTS/'PRE_EFFECT_PROTEIN_SCREEN.csv').set_index('accession')
 if selected is None:
  version=2;mode='EXPLORATORY_NO_DEVELOPMENT_VERSION_PASSED'
 else:version=selected;mode='CALIBRATION_CONDITIONALLY_PASSED'
 data=data.loc[data.version.eq(version)&data.alpha.eq(.05)];decisions=[];passed=[];family=[]
 for accession,q in data.groupby('design_id'):
  ok=bool(selected is not None and len(q)==6 and (q.wilson95_high<=.06).all() and (q.n_failed/q.reps<=.005).all())
  role=screen.loc[accession,'role']
  if ok:passed.append(accession)
  if ok and role=='PRIMARY_DISCOVERY_CANDIDATE':family.append(accession)
  decisions.append({'accession':accession,'gene':screen.loc[accession,'gene'],'role':role,'version':version,'passed_independent_confirmation':ok,
   'max_null_rate':float(q.rate.max()),'max_null_wilson95_upper':float(q.wilson95_high.max()),'max_failure_rate':float((q.n_failed/q.reps).max())})
 csvout(RESULTS/'INDEPENDENT_CONFIRMATION_DECISIONS.csv',decisions)
 old=RESULTS/'MAIN_ANALYSIS_SAP_FROZEN.json'
 frozen={'mode':mode,'version':version,'primary_P':'log2 same-day Elecsys PTAU, standardized across clinical560',
  'primary_model':'Common full nuisance design: age, sex, all three technical factors, Elecsys batch',
  'calibration_passed_including_control':sorted(passed),'primary_family_accessions':sorted(family),'m_primary':len(family),
  'primary_adjustment':'Benjamini-Yekutieli .05; assumes valid marginal p-values; finite catalogue calibration only',
  'BH_secondary':True,'Holm_secondary':True,'n_design_qualified':35,'MAPT_outside_primary_family':True,
  'real_protein_Y_read_before_freeze':False,'protocol_sha256':sha(ROOT/'STAGE5_PROTOCOL_v1.0.md'),
  'pre_effect_screen_sha256':sha(RESULTS/'PRE_EFFECT_PROTEIN_SCREEN.csv'),
  'confirmation_summary_sha256':sha(RESULTS/'CONFIRM_SIMULATION_SUMMARY.csv'),
  'code_sha256':{p.name:sha(p) for p in (ROOT/'code').glob('*.py')}}
 if old.exists():assert json.loads(old.read_text(encoding='utf-8'))['primary_family_accessions']==sorted(family),'Do not redefine family after analysis'
 else:save_json(old,frozen)
 log(f'SAP frozen before protein connection: version V{version}, primary family {len(family)}, passed designs {len(passed)}/35')
 print(pd.DataFrame(decisions).to_string(index=False))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('step',choices=['select','freeze']);args=p.parse_args();select() if args.step=='select' else freeze_sap()
