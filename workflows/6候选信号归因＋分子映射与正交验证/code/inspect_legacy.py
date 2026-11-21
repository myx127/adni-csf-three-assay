from config import *
from pypdf import PdfReader
import zipfile,io
from openpyxl import load_workbook

base=Path(INPUTS['I10']['path']).parent
names=[('L01','CSFMRM_15May2026.csv'),('L02','ADNI_CSFMRM_CONSOLIDATED_DATA_10Apr2014.xlsx'),('L03','ADNI_Bimarkers_Consortium_MRM_Data_Primer04Apr2014.pdf'),
 ('L04','csfproteomics_20190403.zip'),('L05','BiomarkersConsortiumCSF_ADNI1_02Aug2013.zip'),('L06','BiomarkersConsortiumCSF_Proteomics_02Aug2013.zip')]
rows=[{'input_id':k,'path':str(base/n),'bytes':(base/n).stat().st_size,'sha256':sha(base/n)} for k,n in names]
registry=RESULTS/'LEGACY_INPUT_MANIFEST.csv'
if registry.exists():
 prev=pd.read_csv(registry).set_index('input_id');assert all(prev.loc[r['input_id'],'sha256']==r['sha256'] for r in rows)
else:csvout(registry,rows)
headers=[];members=[]
for key,name in names:
 p=base/name
 if p.suffix=='.csv':
  with p.open(encoding='utf-8-sig',newline='') as f:h=next(csv.reader(f))
  headers.append({'source':key,'member':name,'columns':h})
 elif p.suffix=='.pdf':
  t='\n'.join('[PDF page '+str(i+1)+']\n'+(page.extract_text() or '') for i,page in enumerate(PdfReader(p).pages));(EVIDENCE/(key+'_METHODS_EXTRACTED.txt')).write_text(t,encoding='utf-8')
 elif p.suffix=='.xlsx':
  book=load_workbook(p,read_only=True,data_only=True)
  for sheet in book:
   first=next(sheet.iter_rows(values_only=True));headers.append({'source':key,'member':sheet.title,'columns':[str(v) for v in first]})
  book.close()
 elif p.suffix=='.zip':
  with zipfile.ZipFile(p) as z:
   for member in z.namelist():
    members.append({'source':key,'member':member,'bytes':z.getinfo(member).file_size})
    if member.startswith('__MACOSX/') or Path(member).name.startswith('._'):continue
    if member.lower().endswith('.csv'):
     with z.open(member) as f:
      text=io.TextIOWrapper(f,encoding='utf-8-sig',errors='replace',newline='');h=next(csv.reader(text))
     headers.append({'source':key,'member':member,'columns':h})
    if member.lower().endswith('.pdf'):
     t='\n'.join('[PDF page '+str(i+1)+']\n'+(page.extract_text() or '') for i,page in enumerate(PdfReader(io.BytesIO(z.read(member))).pages))
     (EVIDENCE/(key+'_'+Path(member).stem+'_METHODS.txt')).write_text(t,encoding='utf-8')
    if 'dict' in member.lower() and member.lower().endswith('.csv'):
     with z.open(member) as f:raw=f.read()
     (EVIDENCE/(key+'_'+Path(member).name)).write_bytes(raw)
save_json(RESULTS/'LEGACY_SCHEMAS.json',headers);csvout(RESULTS/'LEGACY_ARCHIVE_MEMBERS.csv',members)
for h in headers:
 print(h['source'],h['member'],'ncols',len(h['columns']),'leading',h['columns'][:9],
  'candidate_headers',[c for c in h['columns'] if any(s in c.upper() for s in ['NPTX1','PDGFRB','KLK6','Q15818','P09619','Q92876'])])
log('Six older CSF sources registered before metadata inspection; no legacy candidate effects viewed')
