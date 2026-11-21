from __future__ import annotations
import csv as csvlib,json,re
from datetime import datetime
from pathlib import Path
import pandas as pd
from openpyxl import load_workbook
from cal_config import *

def main():
 checks=[]
 def add(name,passed,detail):checks.append({'check':name,'passed':bool(passed),'detail':detail})
 qa=json.loads((RESULTS/'COMPLETED_VALIDATION_QA.json').read_text(encoding='utf-8'))
 add('all_completed_validation_checks',qa['status']=='PASS' and qa['passed']==qa['total'],{'passed':qa['passed'],'total':qa['total']})
 summary=json.loads((RESULTS/'ANALYSIS_SUMMARY.json').read_text(encoding='utf-8'))
 add('independent_replicate_and_row_accounting',summary['total_independent_replicates']==116000 and summary['all_method_rows_including_replays']==1136000,summary)
 gates=pd.read_csv(RESULTS/'INDEPENDENT_CONFIRMATION_GATES.csv')
 g2=gates.loc[gates.dimension_df.eq(2)];g8=gates.loc[gates.dimension_df.eq(8)]
 add('method_memo_matches_confirmation',g2.control_and_strong_power_pass.all() and not g8.strong_power_gate_pass.any(),gates.to_dict('records'))
 sources=pd.read_csv(RESULTS/'INPUT_SOURCE_HASH_CHECKS.csv')
 add('original_bound_sources_still_unchanged',all(sha(r.path)==r.expected_sha256 for r in sources.itertuples()),len(sources))
 access=json.loads((RESULTS/'ACCESS_AND_PROTOCOL_RECORD.json').read_text(encoding='utf-8'))
 add('primary_frozen_code_and_protocol',all(sha(ROOT/'code'/n)==v for n,v in access['core_code_before_experiments'].items()) and sha(ROOT/'CALIBRATION_PROTOCOL_v1.0.md')==access['protocol_sha256_before_experiments'],'Core changed only in separately declared followup, not original frozen engine')
 mean=json.loads((RESULTS/'MEAN_SELECTION_FROZEN.json').read_text(encoding='utf-8'))
 add('supplement_frozen_code_and_protocol',all(sha(ROOT/'code'/n)==v for n,v in mean['core'].items()) and sha(ROOT/'LOW_DIMENSIONAL_FOLLOWUP_PROTOCOL_v1.0.md')==mean['protocol_sha256'],'Checked')
 wb=load_workbook(ROOT/'有限簇校准与识别敏感性总表.xlsx',read_only=True,data_only=False);schema=pd.read_csv(RESULTS/'WORKBOOK_SCHEMA.csv');counts=[]
 for row in schema.itertuples():
  ws=wb[row.sheet];counts.append({'sheet':row.sheet,'rows':ws.max_row-1,'columns':ws.max_column,'passed':ws.max_row==row.rows+1 and ws.max_column==row.columns})
 add('workbook_all17_sheets',len(wb.sheetnames)==17 and all(c['passed'] for c in counts),counts)
 plan=wb['主实验种子计划'];rows=list(plan.values);j=list(rows[0]).index('seed');values=[r[j] for r in rows[1:]]
 add('primary_seeds_exact_text',values==[c['seed'] for c in grid()] and all(isinstance(v,str) for v in values),len(values))
 from mean_followup import followup_grid
 plan=wb['补充种子计划'];rows=list(plan.values);j=list(rows[0]).index('seed');values=[r[j] for r in rows[1:]]
 add('followup_seeds_exact_text',values==[c['seed'] for c in followup_grid()] and all(isinstance(v,str) for v in values),len(values))
 wb.close()
 required=['FINITE_CLUSTER_AND_ASSUMPTION_REPORT_ZH.md','NEXT_METHOD_DECISION_MEMO.md','README_阅读顺序.md','code/run_all.ps1']
 for stem in ['01_INDEPENDENT_NULL_CALIBRATION','02_CONTROL_POWER_TRADEOFF','03_ERROR_CORRELATION_SENSITIVITY','04_BOOTSTRAP_RESOLUTION_REPLAY']:
  required += [f'figures/{stem}.png',f'figures/{stem}.pdf']
 add('reports_figures_and_reproduction_script',all((ROOT/f).is_file() and (ROOT/f).stat().st_size>0 for f in required),required)
 missing=[]
 for path in ROOT.glob('*.md'):
  for link in re.findall(r'\]\(([^)]+)\)',path.read_text(encoding='utf-8')):
   if link.startswith('http') or link in ['results/DELIVERY_QA.json','CALIBRATION_INPUT_OUTPUT_SHA256.csv']:continue
   if not (ROOT/link).exists():missing.append({'file':path.name,'link':link})
 add('all_local_report_links_resolve',not missing,missing)
 add('no_patient_input_or_formal_SAP_claim',access['clinical_files_opened']==0 and access['individual_protein_files_opened']==0 and not summary['formal_SAP_frozen'],'Bound input/code evidence, not an OS file-access forensic audit')
 passed=all(c['passed'] for c in checks)
 result={'status':'COMPLETE_WITH_CONDITIONAL_MEAN_TEST_CANDIDATE' if passed else 'DELIVERY_FAILED',
    'completed_at_Asia_Shanghai':datetime.now().astimezone().isoformat(),'passed':sum(c['passed'] for c in checks),'total':len(checks),
    'independent_synthetic_replicates':116000,'method_rows_including_replays':1136000,'clinical_files_opened':0,
    'individual_protein_files_opened':0,'final_SAP_frozen':False,'confirmed_conditional_candidate':'2df relative mean + SUM3_CHI2',
    'scientific_limits':['Conservative finite-scenario control, not universal exact calibration','Not all loading changes covered','Error correlation and nonlinearity identification assumptions remain','Technical confounding not removed by robust variance'],
    'checks':checks}
 save_json(RESULTS/'DELIVERY_QA.json',result)
 if not passed:raise RuntimeError('Delivery check failed')
 (LOGS/'FINALIZATION.log').write_text(datetime.now().astimezone().isoformat()+f' {len(checks)}/{len(checks)} delivery checks passed; manifest follows.\n',encoding='utf-8')
 target=ROOT/'CALIBRATION_INPUT_OUTPUT_SHA256.csv';records=[]
 for row in sources.itertuples():
  path=Path(row.path);records.append({'type':'FROZEN_SOURCE_INPUT','path':str(path),'bytes':path.stat().st_size,'sha256':sha(path)})
 for path in [STAGE1/'STAGE1_INPUT_OUTPUT_SHA256.csv',STAGE2/'SIMULATION_INPUT_OUTPUT_SHA256.csv']:
  records.append({'type':'SOURCE_HASH_BINDING_MANIFEST','path':str(path),'bytes':path.stat().st_size,'sha256':sha(path)})
 for path in sorted(ROOT.rglob('*')):
  if path.is_file() and path!=target and '__pycache__' not in path.parts:
   records.append({'type':'STAGE3_OUTPUT','path':path.relative_to(ROOT).as_posix(),'bytes':path.stat().st_size,'sha256':sha(path)})
 csv(target,records)
 with target.open(encoding='utf-8-sig',newline='') as f:written=list(csvlib.DictReader(f))
 for row in written:
  path=ROOT/row['path'] if row['type']=='STAGE3_OUTPUT' else Path(row['path'])
  if path.stat().st_size!=int(row['bytes']) or sha(path)!=row['sha256']:raise RuntimeError('Manifest hash verification failed')
 print(json.dumps({'status':result['status'],'delivery_checks':f'{len(checks)}/{len(checks)}','manifest_entries_verified':len(written)},ensure_ascii=False))

if __name__=='__main__':main()
