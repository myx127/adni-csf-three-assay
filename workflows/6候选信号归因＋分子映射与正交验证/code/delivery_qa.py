from config import *
from openpyxl import load_workbook
from PIL import Image
from pypdf import PdfReader
import ast

def run():
 checks=[]
 def check(name,ok,detail=''):checks.append({'check':name,'passed':bool(ok),'detail':str(detail)})
 qa=json.loads((RESULTS/'INDEPENDENT_QA_SUMMARY.json').read_text(encoding='utf-8'));check('independent_QA_all_pass',qa['failed']==0 and qa['checks']>=570)
 build=json.loads((RESULTS/'WORKBOOK_BUILD.json').read_text(encoding='utf-8'));book=load_workbook(ROOT/build['filename'],read_only=True,data_only=True)
 index=pd.read_csv(RESULTS/'WORKBOOK_SHEET_INDEX.csv');check('workbook_sheet_count',len(book.sheetnames)==build['total_sheets']==len(index)+2)
 check('workbook_has_readme_and_index',book.sheetnames[:2]==['README','INDEX'])
 allvalues=[]
 for r in index.itertuples():
  f=pd.read_csv(ROOT/r.source_csv);ws=book[r.sheet];rows=ws.iter_rows(values_only=True);header=next(rows)
  check('workbook_dimensions_'+r.sheet,ws.max_row==len(f)+1 and ws.max_column==len(f.columns))
  check('workbook_headers_'+r.sheet,list(header)==list(f.columns));ok=True
  for expected,actual in zip(f.itertuples(index=False,name=None),rows):
   for e,a in zip(expected,actual):
    if pd.isna(e):ok&=a is None
    elif isinstance(e,(bool,np.bool_)):ok&=bool(a)==bool(e)
    elif isinstance(e,(int,float,np.integer,np.floating)):
     try:ok&=bool(np.isclose(float(e),float(a),rtol=1e-12,atol=1e-12))
     except (ValueError,TypeError):ok=False
    else:ok&=str(e)==str(a)
  check('workbook_values_'+r.sheet,ok);allvalues.append(ok)
 check('workbook_all67_source_tables_match',len(allvalues)==67 and all(allvalues));book.close()
 figs=pd.read_csv(RESULTS/'FIGURE_INDEX.csv');check('eight_scientific_figures',len(figs)==8)
 for r in figs.itertuples():
  with Image.open(ROOT/r.PNG) as img:check('figure_pixels_'+r.figure,img.width>=2500 and img.height>=1000,(img.width,img.height))
  pdf=PdfReader(ROOT/r.PDF);check('figure_vector_PDF_'+r.figure,len(pdf.pages)==1)
 for p in (ROOT/'code').glob('*.py'):
  try:ast.parse(p.read_text(encoding='utf-8'));valid=True
  except SyntaxError:valid=False
  check('python_syntax_'+p.name,valid)
 files=['STAGE6_MAIN_REPORT_ZH.md','DECISION_MEMO_ZH.md','README_研究结果入口.md','FIELD_GUIDE_ZH.md','METHODS_RESULTS_WORDING_EN.md',
  'EXECUTION_SCOPE_FINAL_v1.5.md','REFERENCES_AND_METHOD_SOURCES.md','REPRODUCIBILITY_ZH.md','NUMERICAL_QA_NOTE.md']
 for name in files:
  content=(ROOT/name).read_text(encoding='utf-8');check('document_nonempty_'+name,len(content)>100 and '\ufffd' not in content)
  # Check actual relative markdown targets, not links used as code in examples.
  import re
  links=re.findall(r'\]\(([^)]+)\)',content)
  for link in links:
   if '://' not in link and not link.startswith('#'):check('document_relative_link_'+name+'_'+link,(ROOT/link.split('#')[0]).exists())
 report=(ROOT/'STAGE6_MAIN_REPORT_ZH.md').read_text(encoding='utf-8');check('report_QA_count_matches',f"{qa['passed']}/{qa['checks']}" in report)
 for phrase in ['3,160,000','169,651','95人','不是外部队列复现','未通过预先冻结的40门槛']:
  check('report_required_evidence_'+phrase,phrase in report)
 status=json.loads((RESULTS/'DIRECT_DISCOVERY_VALIDATION_STATUS.json').read_text(encoding='utf-8'))
 check('final_status_not_Emory_initial',status['direct_discovery_accessions_covered_count']==2 and status['relative_disagreement_replications_completed']==0)
 check('original_Emory_scope_snapshot_retained',(RESULTS/'EMORY_PANEL_DIRECT_VALIDATION_STATUS.json').exists())
 for p in RESULTS.glob('*.csv'):
  cols=pd.read_csv(p,nrows=0).columns;check('no_public_individual_columns_'+p.stem,not set(cols)&{'RID','PTID','key','EXAMDATE','SEQUENCE_NO','REPLICATE1','REPLICATE2','GENOTYPE','PTDOB'})
 out=pd.DataFrame(checks);csvout(RESULTS/'DELIVERY_QA_CHECKS.csv',out)
 summary={'checks':len(out),'passed':int(out.passed.sum()),'failed':int((~out.passed).sum()),'workbook_result_sheets':len(index),
  'scientific_figures':len(figs),'PNG_and_PDF':2*len(figs),'independent_analysis_QA':qa['checks']}
 save_json(RESULTS/'DELIVERY_QA_SUMMARY.json',summary);print(summary)
 if not out.passed.all():print(out.loc[~out.passed].to_string(index=False));raise SystemExit(1)

if __name__=='__main__':run()
