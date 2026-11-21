from config import *
import openpyxl
from PIL import Image
from pypdf import PdfReader
import argparse

def run(allow_unreviewed=False):
 checks=[]
 def add(name,ok,detail=''):checks.append({'check':name,'passed':bool(ok),'detail':str(detail)})
 idx=pd.read_csv(RESULTS/'WORKBOOK_TABLE_INDEX.csv',keep_default_na=False);build=json.loads((RESULTS/'DELIVERY_BUILD_RECORD.json').read_text(encoding='utf-8'));book=ROOT/build['workbook'];wb=openpyxl.load_workbook(book,read_only=True,data_only=True);ncell=0
 def same(a,b):
  if b=='':return a is None or a==''
  if isinstance(b,(bool,np.bool_)):return a==bool(b)
  if isinstance(b,(int,float,np.integer,np.floating)):
   return isinstance(a,(int,float)) and np.isclose(a,b,rtol=1e-12,atol=1e-14,equal_nan=True)
  return str(a)==str(b)
 for t in idx.to_dict('records'):
  p=RESULTS/t['file'];add('TABLE_SHA/'+t['file'],sha(p)==t['sha256']);d=pd.read_csv(p,keep_default_na=False);add('TABLE_SHAPE/'+t['file'],d.shape==(t['rows'],t['columns']))
  if not t['sheet']:continue
  ws=wb[t['sheet']];add('EXCEL_SHAPE/'+t['sheet'],ws.max_row==len(d)+1 and ws.max_column==len(d.columns));actual=ws.iter_rows(values_only=True);add('EXCEL_HEADER/'+t['sheet'],next(actual)==tuple(d.columns));bad=[]
  for j,(a,b) in enumerate(zip(actual,d.itertuples(index=False,name=None)),2):
   ncell+=len(b)
   if len(a)!=len(b) or not all(same(x,y) for x,y in zip(a,b)):bad.append(j)
  add('EXCEL_ALL_VALUES/'+t['sheet'],len(bad)==0,str(bad[:5]))
 add('EXCEL_SHEET_COUNT',len(wb.sheetnames)==build['sheets_including_index']);wb.close()
 figs=pd.read_csv(RESULTS/'FIGURE_INVENTORY.csv');add('FIGURES16',len(figs)==16)
 for r in figs.itertuples():
  p=FIGURES/r.file;add('FIG_SHA/'+r.file,sha(p)==r.sha256)
  if p.suffix=='.png':
   with Image.open(p) as im:add('PNG_SIZE/'+r.file,im.width>=1600 and im.height>=650,f'{im.width}x{im.height}');im.verify()
  else:
   doc=PdfReader(str(p));text=' '.join(x.extract_text() or '' for x in doc.pages);add('PDF_ONE_PAGE/'+r.file,len(doc.pages)==1);add('PDF_VECTOR_TEXT/'+r.file,len(text)>100,len(text))
 reviewfile=RESULTS/'FIGURE_REVIEW_RECORD.csv';reviewed=False
 if reviewfile.exists():
  review=pd.read_csv(reviewfile);reviewed=len(review)==8 and review.passed.all() and all(sha(FIGURES/r.file)==r.png_sha256 for r in review.itertuples())
 add('FIGURE_VISUAL_REVIEW',reviewed or allow_unreviewed,'Visually reviewed PNGs required for sealing; fresh reproductions may defer this review.')
 required=['README.md','00_先读_关键结论与决策.md','01_完整分析报告.md','02_原理与论文表述.md','03_Methods_EN.md','04_参考文献与证据用途.md','05_表格与字段说明.md','STAGE12_PROTOCOL_v1.0.md','EXECUTION_AND_INFERENCE_BOUNDARY_v1.1.md','run_all.ps1']
 for f in required:
  text=(ROOT/f).read_text(encoding='utf-8');add('DOC_NONEMPTY/'+f,len(text)>100);add('DOC_UNICODE/'+f,'\ufffd' not in text)
 full=(ROOT/'01_完整分析报告.md').read_text(encoding='utf-8');add('REPORT_SELECTION_BOUNDARY','不能宣称选择后' in full);add('REPORT_NOT_H','不是第5阶段' in full);add('REPORT_INFORMATION_LOSS','信息可进一步减少' in full);add('REPORT_MECHANISM_BOUNDARY','没有因果机制' in (ROOT/'README.md').read_text(encoding='utf-8'));add('REPORT_SCOPE_LABELS','q33' in full and 'q68' in full)
 qa=json.loads((RESULTS/'SCIENTIFIC_QA_SUMMARY.json').read_text(encoding='utf-8'));add('SCIENTIFIC_QA',qa['checks']==qa['passed'] and qa['failed']==0);d=pd.read_csv(RESULTS/'DECISION_SUMMARY.csv');add('THREE_CANDIDATES',set(d.gene)==set(GENES.values()));add('ALL_PRIMARY_FORM_CHANGES',d.three_P_forms_change_supported_scope.all());add('NO_CAUSAL_OR_SELECTIVE_CLAIM',not d.post_selection_confirmatory_FDR_claimed.any())
 for p in (ROOT/'code').glob('*.py'):
  try:compile(p.read_text(encoding='utf-8'),str(p),'exec');ok=True
  except Exception:ok=False
  add('CODE_COMPILES/'+p.name,ok)
 csvout(RESULTS/'DELIVERY_QA_CHECKS.csv',checks);out={'time':now(),'checks':len(checks),'passed':sum(x['passed'] for x in checks),'failed':sum(not x['passed'] for x in checks),'excel_cells_compared':ncell,'tables':len(idx),'figure_groups':8,'visual_review_completed':reviewed,'allow_unreviewed_copy':allow_unreviewed,'workbook_sha256':sha(book)};save_json(RESULTS/'DELIVERY_QA_SUMMARY.json',out)
 assert out['failed']==0,[x for x in checks if not x['passed']];log(f'Delivery QA {out["passed"]}/{out["checks"]}; {ncell:,} Excel cells match CSV.')
if __name__=='__main__':
 assert_unsealed();ap=argparse.ArgumentParser();ap.add_argument('--allow-unreviewed',action='store_true');run(ap.parse_args().allow_unreviewed)
