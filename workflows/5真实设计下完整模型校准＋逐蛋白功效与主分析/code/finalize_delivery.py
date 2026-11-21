from config import *
from openpyxl import load_workbook
from PIL import Image
checks=[]
def check(name,ok,detail=''):checks.append({'check':name,'passed':bool(ok),'detail':str(detail)})
qa=json.loads((RESULTS/'INDEPENDENT_QA_SUMMARY.json').read_text(encoding='utf-8'));check('independent_QA_all_passed',qa['failed']==0,qa)
check('all5930000_recounted',qa['recounted_independent_replicates']==5930000 and qa['simulation_files_recounted']==1850)
replay=json.loads((RESULTS/'SEEDED_REPLAY_SUMMARY.json').read_text(encoding='utf-8'));check('seeded_replay_passed',replay['passed']==replay['cells'],replay)
base=json.loads((RESULTS/'RESULT_REPRODUCIBILITY_BASELINE.json').read_text(encoding='utf-8'))
for name,digest in base.items():check('actual_result_reproducibility_'+name,sha(RESULTS/name)==digest)
registry=pd.read_csv(RESULTS/'WORKBOOK_SHEET_COUNTS.csv');book=load_workbook(ROOT/'真实设计模型校准_逐蛋白功效与主分析总表.xlsx',read_only=True,data_only=True)
check('workbook_sheet_count',len(book.sheetnames)==len(registry))
for row in registry.itertuples():
 sheet=book[row.sheet];check('sheet_dimensions_'+row.sheet,sheet.max_row==row.rows+1 and sheet.max_column==row.columns)
 cols=[c.value for c in next(sheet.iter_rows())];check('sheet_no_person_ID_'+row.sheet,not any(c in cols for c in ['RID','PTID','key','EXAMDATE','GENOTYPE','PTDOB']))
book.close()
for p in sorted(FIGURES.glob('*.png')):
 with Image.open(p) as im:check('figure_dimensions_'+p.stem,im.width>1000 and im.height>650)
 check('figure_PDF_'+p.stem,p.with_suffix('.pdf').exists() and p.with_suffix('.pdf').stat().st_size>10000)
check('six_figure_pairs',len(list(FIGURES.glob('*.png')))==len(list(FIGURES.glob('*.pdf')))==6)
report=(ROOT/'STAGE5_MAIN_REPORT_ZH.md').read_text(encoding='utf-8');check('report_QA_current',f"{qa['passed']}/{qa['checks']}" in report)
for text in ['教育/APOE','Heywood','不是独立稳健性证据','5930000','不能解释为已经识别','认知斜率']:check('report_limit_'+text,text in report)
sap=json.loads((RESULTS/'MAIN_ANALYSIS_SAP_FROZEN.json').read_text(encoding='utf-8'));check('frozen_family34',sap['m_primary']==34 and sap['MAPT_outside_primary_family'])
main=pd.read_csv(RESULTS/'PRIMARY_FROZEN_FAMILY_RESULTS.csv');check('main_three_reported',set(main.loc[main.BY_discovery,'gene'])=={'PDGFRB','NPTX1','KLK6'})
check('excluded_candidates_retained',len(pd.read_csv(RESULTS/'PRIMARY_ALL43_RESULTS.csv'))==43 and len(pd.read_csv(RESULTS/'PRE_EFFECT_PROTEIN_SCREEN.csv'))==70)
variants=pd.read_csv(RESULTS/'ALL_FIXED_VARIANT_RESULTS.csv');check('extended_cov_no_BY',not variants.loc[variants.variant.eq('FULLCOV'),'BY_discovery'].any())
check('amyloid_no_BY',not variants.loc[variants.variant.eq('AMY'),'BY_discovery'].any())
check('no_individual_latent_scores',json.loads((RESULTS/'MAIN_RUN_SUMMARY.json').read_text(encoding='utf-8'))['individual_F_D_scores_created'] is False)
for name in ['STAGE5_PROTOCOL_v1.0.md','METHOD_IMPLEMENTATION_ZH.md','METHODS_RESULTS_WORDING_EN.md','README_研究结果入口.md','code/run_all.ps1','logs/EXECUTION_NOTES_ZH.md']:
 check('delivery_exists_'+name,(ROOT/name).exists())
out=pd.DataFrame(checks);csvout(RESULTS/'DELIVERY_COMPLETION_CHECKS.csv',out)
status={'status':'COMPLETE_CONDITIONAL_MAIN_ANALYSIS' if out.passed.all() else 'QA_FAILED','completion_checks':len(out),'completion_passed':int(out.passed.sum()),
 'completion_failed':int((~out.passed).sum()),'independent_QA':qa,'seeded_replay':replay,'workbook_sheets':len(registry),'figure_pairs':6,
 'independent_synthetic_replicates':5930000,'primary_family34_BY_discoveries':3,'expanded_covariate_BY_discoveries':0,'amyloid_BY_discoveries':0,
 'observational_mechanism_identified':False,'cognitive_slopes_analyzed':False,'old_SAP_and_stages_overwritten':False}
save_json(RESULTS/'AUDIT_COMPLETION_QA.json',status)
if not out.passed.all():print(out.loc[~out.passed].to_string(index=False));raise SystemExit(1)
log(f'Final delivery QA {len(out)}/{len(out)}; independent QA {qa["passed"]}/{qa["checks"]}; 5930000 independent simulations; conditional main analysis sealed')
files=[]
for p in sorted(ROOT.rglob('*')):
 if not p.is_file() or '__pycache__' in p.parts or p.name=='STAGE5_INPUT_OUTPUT_SHA256.csv':continue
 files.append({'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p),'class':'RESTRICTED_INDIVIDUAL_VALUES' if '00_restricted' in p.parts else 'SYNTHETIC_Q_ONLY' if 'simulation_replicates' in p.parts else 'CODE' if 'code' in p.parts else 'AGGREGATE_OR_DOCUMENT'})
for row in INPUTS.values():files.append({'path':row['path'],'bytes':Path(row['path']).stat().st_size,'sha256':sha(row['path']),'class':'INPUT_UNMODIFIED'})
csvout(ROOT/'STAGE5_INPUT_OUTPUT_SHA256.csv',files);print(json.dumps(clean(status),ensure_ascii=False))
