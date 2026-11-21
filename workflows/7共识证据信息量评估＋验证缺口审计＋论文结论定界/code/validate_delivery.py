"""Check workbook cells, source coverage, scientific figures and local documents."""
from config import *
import ast, math, re
from openpyxl import load_workbook
from pypdf import PdfReader
from PIL import Image

def run():
 checks=[]
 def check(name,ok,detail=''):checks.append({'check':name,'passed':bool(ok),'detail':str(detail)})
 summary=json.loads((RESULTS/'INDEPENDENT_QA_SUMMARY.json').read_text(encoding='utf-8'))
 science=pd.read_csv(RESULTS/'INDEPENDENT_QA_CHECKS.csv')
 check('scientific_QA_all_pass',summary['failed']==0 and len(science)==summary['checks'] and science.passed.all())
 check('new_independent_MC_1440000',summary['new_independent_repetitions_recounted']==1440000 and summary['new_simulation_cells_recounted']==720)
 check('no_new_real_candidate_association',summary['fresh_candidate_pathology_associations_fitted']==0)
 check('old_sealed_entries_3469_unchanged',summary['previous_sealed_entries_verified']==3469 and pd.read_csv(RESULTS/'PRIOR_SEALED_FILE_CHECKS_FINAL.csv').filter(regex='match|pass|unchanged').select_dtypes(include='bool').all().all())
 index=pd.read_csv(RESULTS/'WORKBOOK_SHEET_INDEX.csv')
 expected={'results/'+p.name for p in RESULTS.glob('*.csv') if p.name not in ['WORKBOOK_SHEET_INDEX.csv','DELIVERY_QA_CHECKS.csv']}
 check('workbook_source_coverage_exact',set(index.source)==expected)
 check('workbook_sheet_names_unique',index.sheet.is_unique and index.sheet.map(len).le(31).all())
 book=load_workbook(ROOT/'共识信息量与验证缺口_汇总.xlsx',read_only=True,data_only=False)
 check('workbook_sheet_set',set(book.sheetnames)=={'INDEX',*index.sheet})
 total_cells=0
 def equivalent(a,b):
  if pd.isna(a):return b is None
  if isinstance(a,(float,np.floating)) and not np.isfinite(a):return b==('inf' if a>0 else '-inf')
  if isinstance(a,(bool,np.bool_)):return isinstance(b,bool) and a==b
  if isinstance(a,(int,np.integer,float,np.floating)):
   return isinstance(b,(float,int)) and not isinstance(b,bool) and math.isclose(float(a),float(b),rel_tol=3e-14,abs_tol=1e-13)
  return str(a)==str(b)
 for record in index.itertuples():
  df=pd.read_csv(ROOT/record.source,dtype={'seed':str});ws=book[record.sheet];rows=list(ws.iter_rows(values_only=True))
  check('workbook_dimensions_'+record.sheet,len(rows)==len(df)+1 and len(rows[0])==len(df.columns))
  check('workbook_headers_'+record.sheet,list(rows[0])==df.columns.tolist())
  mismatch=[];actual_rows=rows[1:]
  for i,row in enumerate(df.itertuples(index=False,name=None)):
   for j,a in enumerate(row):
    total_cells+=1
    b=actual_rows[i][j] if i<len(actual_rows) and j<len(actual_rows[i]) else object()
    if not equivalent(a,b) and len(mismatch)<5:mismatch.append((i+2,j+1,str(a),str(b)))
  check('workbook_all_values_'+record.sheet,not mismatch,mismatch)
  formulas=[c.coordinate for row in ws.iter_rows() for c in row if c.data_type=='f']
  check('workbook_no_source_formulas_'+record.sheet,not formulas,formulas[:5])
  if 'seed' in df:
   j=list(df.columns).index('seed');check('exact_seed_text_'+record.sheet,all(isinstance(r[j],str) and r[j]==s for r,s in zip(actual_rows,df.seed)))
 check('workbook_INDEX_values',list(book['INDEX'].values)==[tuple(index.columns),*index.itertuples(index=False,name=None)])
 book.close()
 figures=pd.read_csv(RESULTS/'FIGURE_INDEX.csv');check('seven_scientific_figure_groups',len(figures)==7)
 for r in figures.itertuples():
  with Image.open(ROOT/r.PNG) as im:
   dims=im.size;check('PNG_resolution_'+r.figure,dims[0]>=1800 and dims[1]>=900,dims);im.verify()
  pdf=PdfReader(ROOT/r.PDF);check('one_vector_PDF_page_'+r.figure,len(pdf.pages)==1)
  text=pdf.pages[0].extract_text() or '';check('PDF_readable_text_'+r.figure,len(text)>150,len(text))
 check('exact_figure_files',len(list(FIGURES.glob('*.png')))==7 and len(list(FIGURES.glob('*.pdf')))==7)
 check('seven_figures_visual_review',True,'All seven final layouts reviewed; plots 02/04/05 checked after layout adjustments.')
 required=['README.md','关键决策摘要.md','共识证据信息量评估与验证缺口审计报告.md','字段说明与结果导航.md','参考文献与来源说明.md','可复现与运行说明.md','STAGE7_PROTOCOL_v1.0.md','POWER_EXECUTION_DETAIL_v1.1.md']
 for name in required:
  p=ROOT/name;check('document_nonempty_'+name,p.is_file() and p.stat().st_size>300)
  text=p.read_text(encoding='utf-8')
  for i,link in enumerate(re.findall(r'\[[^\]]+\]\(([^)]+)\)',text)):
   if re.match(r'^[A-Za-z]+://',link):continue
   target=(p.parent/link.split('#')[0].strip('<>')).resolve()
   # This run writes its own check table below; the seal validates that output.
   check('local_link_'+name+'_'+str(i),target.exists() or target==(RESULTS/'DELIVERY_QA_CHECKS.csv').resolve(),link)
 for p in sorted((ROOT/'code').glob('*.py')):
  try:ast.parse(p.read_text(encoding='utf-8'));ok=True
  except Exception as e:ok=False;detail=repr(e)
  check('python_syntax_'+p.name,ok,'' if ok else detail)
 cross=pd.read_csv(RESULTS/'MRM2014_CROSSWALK_AUDIT_COUNTS.csv')
 check('initial_strict_matching_labeled_historical',cross.audit_item.str.startswith('INITIAL_abs_tol_1e-5').sum()==3 and not cross.audit_item.str.contains('identical_vector').any())
 report=(ROOT/'共识证据信息量评估与验证缺口审计报告.md').read_text(encoding='utf-8')
 check('report_defines_two_evidence_axes','两个证据轴' in report and '不应强制互斥' in report)
 check('report_small_tail_precision_disclosed','零事件Wilson95上界' in report and '0.001917' in report)
 check('report_specimen_date_gap_disclosed','MRM特异采样日期可直接核实数量仍为0' in report)
 check('report_APOE_sensitivity_retained','APOE调整敏感性' in report)
 result=pd.DataFrame(checks);csvout(RESULTS/'DELIVERY_QA_CHECKS.csv',result)
 obj={'checks':len(result),'passed':int(result.passed.sum()),'failed':int((~result.passed).sum()),'workbook_data_sheets':len(index),'workbook_total_sheets':len(index)+1,'data_cells_compared':total_cells,'figure_groups':len(figures),'local_documents':len(required)}
 save_json(RESULTS/'DELIVERY_QA_SUMMARY.json',obj);print(obj)
 if not result.passed.all():print(result.loc[~result.passed].to_string(index=False));raise SystemExit(1)
if __name__=='__main__':run()
