from config import *
from openpyxl import load_workbook
from PIL import Image
import ast,re,math

def run():
    checks=[];compared=0
    def add(name,ok,detail=''):checks.append({'check':name,'passed':bool(ok),'detail':str(detail)})
    index=pd.read_csv(RESULTS/'WORKBOOK_SHEET_INDEX.csv')
    book=load_workbook(ROOT/'核心候选证据整合与验证_结果与审计.xlsx',read_only=True,data_only=False)
    add('WORKBOOK_INDEX',book.sheetnames[0]=='INDEX')
    for entry in index.itertuples():
        p=RESULTS/entry.source_CSV
        d=pd.read_csv(p,dtype={c:str for c in pd.read_csv(p,nrows=0).columns if 'seed' in c.lower()})
        add('CSV_EXISTS/'+p.name,p.is_file())
        add('CSV_DIMENSIONS/'+p.name,len(d)==entry.rows and len(d.columns)==entry.columns)
        if entry.storage=='CSV_ONLY':continue
        ws=book[entry.sheet];rows=ws.iter_rows(values_only=True);header=next(rows)
        add('EXCEL_HEADER/'+p.name,list(header)==list(d.columns))
        add('EXCEL_DIMENSIONS/'+p.name,ws.max_row==len(d)+1 and ws.max_column==len(d.columns))
        ok=True;bad=''
        for number,(actual,expected) in enumerate(zip(rows,d.itertuples(index=False,name=None)),2):
            for col,(a,e) in enumerate(zip(actual,expected),1):
                compared+=1
                if isinstance(e,np.generic):e=e.item()
                if e is None or (isinstance(e,float) and not math.isfinite(e)):one=a is None
                elif isinstance(e,bool):one=isinstance(a,bool) and a==e
                elif isinstance(e,(int,float)) and not isinstance(e,bool):
                    one=isinstance(a,(int,float)) and not isinstance(a,bool) and math.isclose(a,e,rel_tol=1e-12,abs_tol=1e-15)
                    if e!=0 and abs(e)<1e-6:one=one and a!=0
                else:one=str(a)==str(e)
                if not one:ok=False;bad=f'row{number}/col{col}: {a!r} != {e!r}';break
            if not ok:break
        add('EXCEL_CELL_VALUES/'+p.name,ok,bad)
    add('EXCEL_SHEET_COUNT',len(book.sheetnames)==1+int(index.storage.eq('EXCEL_AND_CSV').sum()))
    book.close()
    indexed=set(index.source_CSV);allowed={'WORKBOOK_SHEET_INDEX.csv','DELIVERY_QA_CHECKS.csv'}
    missing={p.name for p in RESULTS.glob('*.csv')}-indexed-allowed
    add('ALL_CSV_INDEXED_OR_QA_SELF',not missing,sorted(missing))
    for p in (ROOT/'code').glob('*.py'):
        try:ast.parse(p.read_text(encoding='utf-8-sig'));ok=True
        except SyntaxError as exc:ok=False
        add('PYTHON_SYNTAX/'+p.name,ok)
    for p in ROOT.glob('*.md'):
        text=p.read_text(encoding='utf-8-sig')
        add('UTF8_REPORT/'+p.name,'\ufffd' not in text)
        for target in re.findall(r'\]\(([^)]+)\)',text):
            if target.startswith(('http:','https:','#')):continue
            target=target.strip('<>');add('REPORT_LINK/'+p.name+'/'+target,(ROOT/target).is_file())
    figures=pd.read_csv(RESULTS/'FIGURE_INDEX.csv')
    for r in figures.itertuples():
        with Image.open(ROOT/r.PNG) as im:
            add('FIGURE_PNG/'+r.figure,im.width>=1000 and im.height>=700)
        add('FIGURE_PDF/'+r.figure,(ROOT/r.PDF).read_bytes().startswith(b'%PDF-'))
    review=pd.read_csv(RESULTS/'FIGURE_REVIEW_RECORD.csv')
    add('SEVEN_FIGURES_REVIEWED',len(review)==7 and review.PNG_and_PDF_present.all())
    protocol=json.loads((RESULTS/'EXECUTION_PLAN_RECORD.json').read_text(encoding='utf-8'))
    add('ORIGINAL_STAGE9_PROTOCOL_UNCHANGED',protocol['protocol_sha256']==sha(ROOT/'STAGE9_PROTOCOL_v1.0.md'))
    science=json.loads((RESULTS/'SCIENTIFIC_QA_SUMMARY.json').read_text(encoding='utf-8'))
    add('SCIENTIFIC_QA_PASSED',science['failed']==0)
    core=pd.read_csv(RESULTS/'THREE_CORE_CANDIDATE_EVIDENCE_MATRIX.csv')
    add('NO_UNSUPPORTED_VALIDATION_UPGRADE',not core.strict_relative_pattern_validation_completed.any() and not core.independent_external_cohort_validation_completed.any() and not core.specific_molecular_form_or_mechanism_verified.any())
    add('RESTRICTED_TEMPLATE_NOT_IN_WORKBOOK',not any('REQUEST_TEMPLATE' in str(s) for s in index.source_CSV))
    add('VERIFIER_AND_SAFE_REPRODUCER_EXIST',all((ROOT/'code'/n).exists() for n in ['verify_delivery.py','reproduce_in_copy.py','run_all.ps1']))
    csvout(RESULTS/'DELIVERY_QA_CHECKS.csv',checks)
    save_json(RESULTS/'DELIVERY_QA_SUMMARY.json',{'checks':len(checks),'passed':sum(x['passed'] for x in checks),'failed':sum(not x['passed'] for x in checks),'Excel_total_sheets':1+int(index.storage.eq('EXCEL_AND_CSV').sum()),'Excel_data_cells_compared':compared,'indexed_tables':len(index),'CSV_only_tables':int(index.storage.eq('CSV_ONLY').sum()),'figure_groups':len(figures),'time':now()})
    assert all(x['passed'] for x in checks),[x for x in checks if not x['passed']][:10]
    log(f'Delivery QA:{len(checks)}/{len(checks)} passed;{compared} Excel cells compared')

if __name__=='__main__':
    assert_unsealed()
    run()
