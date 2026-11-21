from config import *
import openpyxl
from PIL import Image
def run():
 checks=[];cells=0
 def add(name,ok,detail=''):checks.append({'check':name,'passed':bool(ok),'detail':str(detail)})
 idx=pd.read_csv(RESULTS/'WORKBOOK_SHEET_INDEX.csv');book=ROOT/'14个共同病理关联蛋白_临床信息评价_结果与审计.xlsx'
 w=openpyxl.load_workbook(book,read_only=True,data_only=True)
 for row in idx.itertuples():
  d=pd.read_csv(RESULTS/row.csv);add('CSV_SHAPE/'+row.csv,len(d)==row.rows and len(d.columns)==row.columns)
  if row.sheet=='CSV_ONLY':continue
  add('WORKBOOK_SHEET/'+row.sheet,row.sheet in w.sheetnames);s=w[row.sheet]
  it=s.iter_rows(values_only=True);header=next(it);add('HEADER/'+row.sheet,list(header)==list(d.columns))
  correct=True;maxerr=0.;nrows=0
  for expected,actual in zip(d.itertuples(index=False,name=None),it):
   nrows+=1
   for a,b in zip(expected,actual):
    cells+=1
    if pd.isna(a):ok=b is None or b==''
    elif isinstance(a,(bool,np.bool_)):ok=isinstance(b,bool) and bool(a)==b
    elif isinstance(a,(int,float,np.number)):
     ok=b is not None and isinstance(b,(int,float)) and np.isclose(a,b,rtol=1e-11,atol=1e-13)
     if b is not None and isinstance(b,(int,float)):maxerr=max(maxerr,abs(float(a)-float(b)))
    else:ok=str(a)==str(b)
    if not ok:correct=False
  add('WORKBOOK_ALL_CELLS/'+row.sheet,correct and nrows==len(d),f'{nrows} rows; max numeric error={maxerr:.3g}')
 w.close()
 figures=pd.read_csv(RESULTS/'FIGURE_INDEX.csv');review=pd.read_csv(RESULTS/'FIGURE_REVIEW_RECORD.csv')
 add('EIGHT_FIGURES',len(figures)==8);add('EIGHT_VISUAL_REVIEWS',len(review)==8 and review.passed.all())
 for r in figures.itertuples():
  for col in ['png','pdf']:
   f=ROOT/getattr(r,col);add('FIGURE_FILE/'+f.name,f.is_file() and f.stat().st_size>5000)
  with Image.open(ROOT/r.png) as im:add('PNG_RESOLUTION/'+r.basename,im.width>=1000 and im.height>=700,str(im.size))
 for name in ['README.md','STAGE11_PROTOCOL_v1.0.md','STAGE11_MAIN_REPORT_ZH.md','关键决策摘要.md','模型原理与结果解释.md','参考文献与来源说明.md','可复现与运行说明.md','字段说明与结果导航.md']:
  add('DOCUMENT/'+name,(ROOT/name).exists() and (ROOT/name).stat().st_size>200)
 f=pd.read_csv(RESULTS/'FORMAL204_CLINICAL_RESULTS.csv');add('ALL34_RETAINED',f.gene.nunique()==34);add('ALL14_RETAINED',f.loc[f.focus14].gene.nunique()==14)
 add('ALL204_RETAINED',len(f)==204);add('NO_MAPT',not f.gene.eq('MAPT').any())
 add('ALL_SCIENCE_PASSED',json.loads((RESULTS/'SCIENTIFIC_QA_SUMMARY.json').read_text(encoding='utf-8'))['failed']==0)
 out={'time':now(),'checks':len(checks),'passed':sum(r['passed'] for r in checks),'failed':sum(not r['passed'] for r in checks),'Excel_cells_compared':cells,'workbook_sheets':len(idx.loc[idx.sheet.ne('CSV_ONLY')])+1,'indexed_tables':len(idx),'csv_only_tables':idx.loc[idx.sheet.eq('CSV_ONLY'),'csv'].tolist(),'figure_groups':8}
 csvout(RESULTS/'DELIVERY_QA_CHECKS.csv',checks);save_json(RESULTS/'DELIVERY_QA_SUMMARY.json',out)
 assert out['failed']==0,[r for r in checks if not r['passed']]
 log(f'Delivery QA passed {len(checks)}/{len(checks)}; {cells} Excel cells compared.')
if __name__=='__main__':assert_unsealed();run()

