from __future__ import annotations
import csv,json
from audit_common import *
from openpyxl import load_workbook
from PIL import Image

checks=[]
def check(name,condition,detail=''):
 checks.append({'check':name,'passed':bool(condition),'detail':detail})

qa=json.loads((RESULTS/'INDEPENDENT_QA_SUMMARY.json').read_text(encoding='utf-8'))
check('independent_validation_passed',qa['failed']==0,str(qa))
registry=pd.read_csv(RESULTS/'WORKBOOK_SHEET_COUNTS.csv')
book=load_workbook(ROOT/'真实预处理与病理认知覆盖审计总表.xlsx',read_only=True,data_only=True)
check('workbook_sheet_count',len(book.sheetnames)==len(registry),str(len(book.sheetnames)))
for row in registry.itertuples():
 ws=book[row.sheet];check('sheet_rows_'+row.sheet,ws.max_row==row.rows+1,f'{ws.max_row-1}/{row.rows}')
 check('sheet_columns_'+row.sheet,ws.max_column==row.columns,f'{ws.max_column}/{row.columns}')
 header=[c.value for c in next(ws.iter_rows())]
 check('sheet_no_individual_ID_'+row.sheet,not any(c in header for c in ['RID','PTID','key','EXAMDATE']))
book.close()
report=(ROOT/'STAGE4_AUDIT_REPORT_ZH.md').read_text(encoding='utf-8')
check('report_QA_matches',f"{qa['passed']}/{qa['checks']}" in report)
rules=pd.read_csv(RESULTS/'PREPROCESSING_DECISIONS.csv')
check('31_rules_unique',len(rules)==31 and rules.rule_id.is_unique)
overlap=pd.read_csv(RESULTS/'CLINICAL_CONSTRUCT_OVERLAP_FLAGS.csv')
check('MAPT_construct_overlap_explicit',overlap.overlaps_primary_ptau_parent_protein.sum()==1 and overlap.loc[overlap.overlaps_primary_ptau_parent_protein,'nulisa_target'].iloc[0]=='MAPT')
check('70x4_protein_routes',len(pd.read_csv(RESULTS/'PER_PROTEIN_JOINT_COVERAGE.csv'))==280)
for p in sorted(FIGURES.glob('*.png')):
 with Image.open(p) as im:
  check('figure_valid_'+p.stem,im.width>1000 and im.height>700,f'{im.width}x{im.height}')
 check('figure_PDF_'+p.stem,p.with_suffix('.pdf').stat().st_size>10000)
check('four_figure_pairs',len(list(FIGURES.glob('*.png')))==len(list(FIGURES.glob('*.pdf')))==4)
for path in [ROOT/'README_审计总览.md',ROOT/'PREPROCESSING_DECISION_MEMO_ZH.md',ROOT/'code/run_all.ps1',ROOT/'logs/EXECUTION_NOTES_ZH.md']:
 check('delivery_exists_'+path.name,path.is_file() and path.stat().st_size>100)
# Verify sealed Stage1 lookup files, without opening individual protein quantities.
manifest=pd.read_csv(STAGE1/'STAGE1_INPUT_OUTPUT_SHA256.csv')
for rel in ['00_restricted/TRIPLE_SAMPLE_REGISTRY_RESTRICTED.csv','00_restricted/TRIPLE_CANDIDATE_AVAILABILITY_RESTRICTED.parquet']:
 expected=manifest.loc[manifest.path.eq(rel),'sha256'].iloc[0]
 check('sealed_stage1_'+Path(rel).name,sha(STAGE1/rel)==expected)
summary=json.loads((RESULTS/'COVERAGE_RUN_SUMMARY.json').read_text(encoding='utf-8'))
check('no_association_models',summary['association_models_run']==0)
check('no_raw_clinical_value_output',summary['raw_continuous_clinical_values_saved'] is False)
check('no_protein_quantitative_join',summary['protein_quantitative_values_joined'] is False)
check('restricted_no_temporary_signatures',all(not any(c.startswith('_sig_') for c in pd.read_csv(p,nrows=0).columns) for p in RESTRICTED.glob('*.csv')))
for name in ['PREPROCESSING_RULES_v1.0.json','FIELD_ACCESS_LOG.json','ACCESS_BOUNDARY_INITIAL.json','THRESHOLD_FREEZE_BEFORE_COUNTS.json','EVIDENCE_REGISTRY.csv','NEXT_ACTIONS.csv']:
 check('trace_exists_'+name,(RESULTS/name).exists())
out=pd.DataFrame(checks);write_csv(RESULTS/'STAGE4_COMPLETION_CHECKLIST.csv',out)
state={'status':'COMPLETE_WITH_OPEN_METHOD_DECISIONS' if out.passed.all() else 'QA_FAILED',
 'completion_checks':len(out),'completion_passed':int(out.passed.sum()),'completion_failed':int((~out.passed).sum()),
 'independent_checks':qa['checks'],'independent_passed':qa['passed'],'workbook_sheets':len(registry),'figure_pairs':4,
 'rules':31,'raw_continuous_clinical_values_saved':False,'association_models_run':0,
 'formal_SAP_frozen':False,'remaining_issues':'Actual-design model calibration/power/multiplicity; age conflicts8; source version3447/3907; cognition representation'}
save_json(RESULTS/'AUDIT_COMPLETION_QA.json',state)
if not out.passed.all():
 print(out.loc[~out.passed].to_string(index=False));raise SystemExit(1)
log(f"Final delivery QA {int(out.passed.sum())}/{len(out)} passed; independent QA {qa['passed']}/{qa['checks']}; no clinical effects estimated")
files=[]
excluded={'STAGE4_INPUT_OUTPUT_SHA256.csv'}
for p in sorted(ROOT.rglob('*')):
 if not p.is_file() or '__pycache__' in p.parts or p.name in excluded:continue
 files.append({'relative_path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p),
  'delivery_class':'RESTRICTED_FLAGS_NO_RAW_VALUES' if '00_restricted' in p.parts else 'CODE' if 'code' in p.parts else 'AGGREGATE_OR_DOCUMENT'})
for row in INPUTS.values():files.append({'relative_path':row['path'],'bytes':Path(row['path']).stat().st_size,'sha256':sha(row['path']),'delivery_class':'INPUT_UNMODIFIED'})
write_csv(ROOT/'STAGE4_INPUT_OUTPUT_SHA256.csv',files)
print(json.dumps(state,ensure_ascii=False))
