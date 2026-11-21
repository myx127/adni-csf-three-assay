from config import *
from openpyxl import load_workbook
from PIL import Image
import ast,math,re
def run():
 checks=[];cells=0
 def add(name,ok,detail=''):checks.append({'check':name,'passed':bool(ok),'detail':str(detail)})
 index=pd.read_csv(RESULTS/'WORKBOOK_SHEET_INDEX.csv');wb=load_workbook(ROOT/'AD病理背景与五年认知进展_结果与审计.xlsx',read_only=True,data_only=False)
 add('WORKBOOK_INDEX',wb.sheetnames[0]=='INDEX')
 for r in index.itertuples():
  p=RESULTS/r.source_CSV;d=pd.read_csv(p,dtype={c:str for c in pd.read_csv(p,nrows=0).columns if 'seed' in c.lower()})
  add('CSV_DIMENSIONS/'+p.name,len(d)==r.rows and len(d.columns)==r.columns)
  if r.storage=='CSV_ONLY':continue
  ws=wb[r.sheet];rows=ws.iter_rows(values_only=True);header=next(rows);add('HEADER/'+p.name,list(header)==list(d));add('DIMENSIONS/'+p.name,ws.max_row==len(d)+1 and ws.max_column==len(d.columns))
  ok=True;bad=''
  for number,(actual,expected) in enumerate(zip(rows,d.itertuples(index=False,name=None)),2):
   for j,(a,e) in enumerate(zip(actual,expected),1):
    cells+=1
    if isinstance(e,np.generic):e=e.item()
    if e is None or isinstance(e,float) and not math.isfinite(e):one=a is None
    elif isinstance(e,bool):one=isinstance(a,bool) and a==e
    elif isinstance(e,(int,float)):
     one=isinstance(a,(int,float)) and not isinstance(a,bool) and math.isclose(a,e,rel_tol=1e-12,abs_tol=1e-15)
     if e!=0 and abs(e)<1e-6:one=one and a!=0
    else:one=str(a)==str(e)
    if not one:ok=False;bad=f'row{number}/col{j}: {a!r} != {e!r}';break
   if not ok:break
  add('ALL_EXCEL_VALUES/'+p.name,ok,bad)
 add('SHEET_COUNT',len(wb.sheetnames)==1+int(index.storage.eq('EXCEL_AND_CSV').sum()));wb.close()
 allowed={'WORKBOOK_SHEET_INDEX.csv','DELIVERY_QA_CHECKS.csv'}
 missing={p.name for p in RESULTS.glob('*.csv')}-set(index.source_CSV)-allowed;add('ALL_CSV_INDEXED',not missing,sorted(missing))
 for p in (ROOT/'code').glob('*.py'):
  try:ast.parse(p.read_text(encoding='utf-8-sig'));ok=True
  except SyntaxError:ok=False
  add('PYTHON_SYNTAX/'+p.name,ok)
 for p in ROOT.glob('*.md'):
  text=p.read_text(encoding='utf-8-sig');add('UTF8/'+p.name,'\ufffd' not in text)
  for target in re.findall(r'\]\(([^)]+)\)',text):
   if not target.startswith(('http:','https:','#')):add('LINK/'+p.name+target,(ROOT/target.strip('<>')).is_file())
 figures=pd.read_csv(RESULTS/'FIGURE_INDEX.csv')
 for r in figures.itertuples():
  with Image.open(ROOT/r.PNG) as im:add('PNG/'+r.figure,im.width>=1000 and im.height>=700)
  add('PDF/'+r.figure,(ROOT/r.PDF).read_bytes().startswith(b'%PDF-'))
 review=pd.read_csv(RESULTS/'FIGURE_REVIEW_RECORD.csv');add('ALL_EIGHT_VISUALLY_REVIEWED',len(review)==8 and review.reviewed.all())
 science=json.loads((RESULTS/'SCIENTIFIC_QA_SUMMARY.json').read_text(encoding='utf-8'));add('SCIENTIFIC_QA_ZERO_FAILURES',science['failed']==0)
 temp=json.loads((RESULTS/'TEMPLATE_FREEZE_RECORD.json').read_text(encoding='utf-8'));add('TEMPLATE_FREEZE_UNCHANGED',temp['templates_sha256']==sha(PRIVATE/'ANALYSIS_TEMPLATES_RESTRICTED.pkl'))
 form=pd.read_csv(RESULTS/'FORMAL32_CLINICAL_RESULTS.csv');add('CLINICAL_NO_BY_FINDINGS',not form.supported_within_family.any() and not form.supported_combined32.any())
 unqual=pd.read_csv(RESULTS/'UNQUALIFIED_INTERACTION_POINTS_ONLY.csv');add('UNQUALIFIED_POINTS_ONLY',len(unqual)==9 and not any('CI' in c or c.lower().startswith('p_') for c in unqual))
 add('SIX_UNQUALIFIED_FORMAL_P1',int((~form.calibration_qualified&form.p_for_family.eq(1)).sum())==6)
 add('NO_REAL_PERSON_FILES_IN_WORKBOOK',not any('RESTRICTED' in f for f in index.source_CSV))
 add('WORKFLOW_FILES',all((ROOT/'code'/n).is_file() for n in ['verify_delivery.py','reproduce_in_copy.py','run_all.ps1','seal_delivery.py']))
 csvout(RESULTS/'DELIVERY_QA_CHECKS.csv',checks)
 result={'time':now(),'checks':len(checks),'passed':sum(r['passed'] for r in checks),'failed':sum(not r['passed'] for r in checks),'Excel_data_cells_compared':cells,'Excel_total_sheets':1+int(index.storage.eq('EXCEL_AND_CSV').sum()),'indexed_tables':len(index),'CSV_ONLY_tables':int(index.storage.eq('CSV_ONLY').sum()),'figure_groups':len(figures)}
 save_json(RESULTS/'DELIVERY_QA_SUMMARY.json',result)
 assert all(r['passed'] for r in checks),[r for r in checks if not r['passed']][:10]
 log(f'Delivery QA:{len(checks)}/{len(checks)} passed;{cells} Excel data cells compared')
if __name__=='__main__':
 assert_unsealed();run()
