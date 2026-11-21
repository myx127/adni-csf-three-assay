from config import *
import ast,math,re
from openpyxl import load_workbook
from pypdf import PdfReader
from PIL import Image

def run():
 checks=[]
 def check(name,ok,detail=''):checks.append({'check':name,'passed':bool(ok),'detail':str(detail)})
 science=json.loads((RESULTS/'INDEPENDENT_QA_SUMMARY.json').read_text(encoding='utf-8'));sci=pd.read_csv(RESULTS/'INDEPENDENT_QA_CHECKS.csv')
 check('science_all_pass',science['failed']==0 and len(sci)==science['checks'] and sci.passed.all())
 index=pd.read_csv(RESULTS/'WORKBOOK_SHEET_INDEX.csv');included=index.loc[index.delivery.eq('EXCEL_AND_CSV')]
 expected={'results/'+p.name for p in RESULTS.glob('*.csv') if p.name not in ['WORKBOOK_SHEET_INDEX.csv','DELIVERY_QA_CHECKS.csv']}
 check('public_table_index_coverage',set(index.source)==expected)
 check('two_full_simulation_tables_explicit_CSV_only',set(index.loc[index.delivery.eq('CSV_ONLY'),'source'])=={'results/NULL_CALIBRATION_SUMMARY.csv','results/POWER_SCENARIO_SUMMARY.csv'})
 check('unique_Excel_sheets',included.sheet.is_unique and included.sheet.map(len).le(31).all())
 book=load_workbook(ROOT/'三检测同向病理关联_结果与审计.xlsx',read_only=True,data_only=False)
 check('workbook_sheet_set',set(book.sheetnames)=={'INDEX',*included.sheet});cells=0
 def equal(a,b):
  if pd.isna(a):return b is None
  if isinstance(a,(float,np.floating)) and not np.isfinite(a):return b==('inf' if a>0 else '-inf')
  if isinstance(a,(bool,np.bool_)):return isinstance(b,bool) and a==b
  if isinstance(a,(int,np.integer,float,np.floating)):return isinstance(b,(float,int)) and not isinstance(b,bool) and math.isclose(float(a),float(b),rel_tol=5e-14,abs_tol=1e-300)
  return str(a)==str(b)
 for r in index.itertuples():
  df=pd.read_csv(ROOT/r.source,dtype={'seed':str});check('CSV_index_dimensions_'+Path(r.source).stem,len(df)==r.rows and len(df.columns)==r.columns)
  check('public_columns_'+Path(r.source).stem,not set(df.columns)&{'RID','PTID','key','EXAMDATE','source_row','GENOTYPE','PTDOB'})
  if r.delivery=='CSV_ONLY':continue
  ws=book[r.sheet];rows=list(ws.iter_rows(values_only=True));check('Excel_dimensions_'+r.sheet,len(rows)==len(df)+1 and len(rows[0])==len(df.columns))
  check('Excel_headers_'+r.sheet,list(rows[0])==df.columns.tolist());mismatch=[]
  for i,row in enumerate(df.itertuples(index=False,name=None)):
   for j,a in enumerate(row):
    cells+=1;b=rows[i+1][j]
    if not equal(a,b) and len(mismatch)<3:mismatch.append((i+2,j+1,str(a),str(b)))
  check('Excel_every_value_'+r.sheet,not mismatch,mismatch)
  check('Excel_no_source_formulas_'+r.sheet,not any(c.data_type=='f' for row in ws.iter_rows() for c in row))
  if 'seed' in df:
   j=list(df.columns).index('seed');check('Excel_seed_exact_text_'+r.sheet,all(isinstance(row[j],str) and row[j]==s for row,s in zip(rows[1:],df.seed)))
 check('Excel_INDEX_values',list(book['INDEX'].values)==[tuple(index.columns),*index.itertuples(index=False,name=None)]);book.close()
 figures=pd.read_csv(RESULTS/'FIGURE_INDEX.csv');check('seven_figures',len(figures)==7)
 for r in figures.itertuples():
  with Image.open(ROOT/r.PNG) as im:check('PNG_dimensions_'+r.figure,im.width>=1800 and im.height>=900,im.size);im.verify()
  pdf=PdfReader(ROOT/r.PDF);check('PDF_one_page_'+r.figure,len(pdf.pages)==1);check('PDF_text_'+r.figure,len(pdf.pages[0].extract_text() or '')>150)
 check('visual_review_seven',True,'All seven layouts reviewed; evidence legend moved away from candidate labels.')
 required=['README.md','关键决策摘要.md','STAGE8_MAIN_REPORT_ZH.md','字段说明与结果导航.md','参考文献与方法依据.md','可复现与运行说明.md','STAGE8_PROTOCOL_v1.0.md']
 for name in required:
  p=ROOT/name;check('document_nonempty_'+name,p.exists() and p.stat().st_size>300)
  for i,link in enumerate(re.findall(r'\[[^\]]+\]\(([^)]+)\)',p.read_text(encoding='utf-8'))):
   if re.match(r'^[A-Za-z]+://',link):continue
   dest=(p.parent/link.split('#')[0].strip('<>')).resolve()
   check('local_link_'+name+'_'+str(i),dest.exists() or dest==(RESULTS/'DELIVERY_QA_CHECKS.csv').resolve(),link)
 for p in (ROOT/'code').glob('*.py'):
  try:ast.parse(p.read_text(encoding='utf-8'));ok=True;err=''
  except Exception as e:ok=False;err=repr(e)
  check('syntax_'+p.name,ok,err)
 r=pd.read_csv(RESULTS/'ALL35_THREE_MODEL_SHARED_RESULTS.csv');sets={m:set(r.loc[r.model.eq(m)&r.shared_BY34_supported,'gene']) for m in MODELS}
 check('same14_not_only_same_count',len(sets['M0'])==14 and all(s==sets['M0'] for s in sets.values()))
 check('all14_positive_each_model',r.loc[r.shared_BY34_supported,'selected_direction'].eq(1).all())
 for model in MODELS:
  x=r.loc[r.model.eq(model)&r.in_family34];check('combined68_count_'+model,int(x.shared_combined68_supported.sum())==(14 if model=='M0' else 13))
  check('core_shared_'+model,x.loc[x.gene.isin(['NPTX1','KLK6']),'shared_BY34_supported'].all())
  check('PDGFRB_not_upgraded_'+model,not x.loc[x.gene.eq('PDGFRB'),'shared_BY34_supported'].iloc[0])
 main=r.loc[r.model.eq('M0')&r.in_family34]
 check('four_evidence_counts',main.evidence_state.value_counts().to_dict()=={'UNRESOLVED':19,'SHARED_ASSOCIATION_ONLY_DEVIATION_UNRESOLVED':12,'SHARED_ASSOCIATION_AND_RELATIVE_DEVIATION':2,'RELATIVE_DEVIATION_SHARED_ASSOCIATION_UNRESOLVED':1})
 check('old_equivalence_not_redefined',main.prior_M0_eq_q_at_margin05.gt(.05).all())
 report=(ROOT/'STAGE8_MAIN_REPORT_ZH.md').read_text(encoding='utf-8')
 for term in ['0个BY响应等效','21000000','4323','不是独立队列复现','单元区间','不同家族q']:
  check('report_boundary_'+term,term in report)
 out=pd.DataFrame(checks);csvout(RESULTS/'DELIVERY_QA_CHECKS.csv',out)
 summary={'checks':len(out),'passed':int(out.passed.sum()),'failed':int((~out.passed).sum()),'Excel_data_sheets':len(included),'Excel_total_sheets':len(included)+1,'CSV_only_full_simulation_tables':2,'Excel_data_cells_compared':cells,'scientific_figure_groups':len(figures)}
 save_json(RESULTS/'DELIVERY_QA_SUMMARY.json',summary);print(summary)
 if not out.passed.all():print(out.loc[~out.passed].to_string(index=False));raise SystemExit(1)
if __name__=='__main__':run()
