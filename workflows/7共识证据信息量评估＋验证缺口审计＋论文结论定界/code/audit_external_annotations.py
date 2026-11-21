from config import *
from openpyxl import load_workbook
import urllib.request,re

def run():
 url='https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13195-022-01113-5/MediaObjects/13195_2022_1113_MOESM2_ESM.xlsx'
 path=PRIVATE/'PUBLIC_EMORY2022_SUPPLEMENT2_ORIGINAL.xlsx';manifest=RESULTS/'EXTERNAL_ANNOTATION_SOURCE_MANIFEST.json'
 if not path.exists():
  request=urllib.request.Request(url,headers={'User-Agent':'Research-metadata-audit/1.0'})
  with urllib.request.urlopen(request,timeout=45) as response:body=response.read()
  assert body[:2]==b'PK';path.write_bytes(body)
  save_json(manifest,{'source_id':'X01','url':url,'publisher_article':'https://doi.org/10.1186/s13195-022-01113-5','path':str(path.relative_to(ROOT)),
   'sha256':sha(path),'bytes':len(body),'retrieved':datetime.now().astimezone().isoformat(),'permitted_use_this_stage':'Published gene/platform overlap annotation only; no individual candidate associations or correlation effects analyzed'})
 else:assert sha(path)==json.loads(manifest.read_text(encoding='utf-8'))['sha256']
 book=load_workbook(path,read_only=True,data_only=True);print('Public supplementary sheet names:',book.sheetnames)
 schemas=[];hits=[];coverage=[]
 for ws in book:
  number=re.search(r'^(?:Table[ _-]*)?(\d+)[.]?',ws.title,re.I)
  if not number or int(number[1]) not in [4,5]:continue
  data=list(ws.iter_rows(values_only=True));schemas.append({'sheet':ws.title,'rows':len(data),'columns':ws.max_column,'annotation_only':True})
  print('Overlap sheet',ws.title,'header rows:',data[3:8])
  for i,row in enumerate(data):
   for j,v in enumerate(row):
    if isinstance(v,str):
     for acc,gene in CANDIDATES.items():
      if re.search(r'(?<![A-Z0-9])'+gene+r'(?![A-Z0-9])',v):
       hits.append({'accession':acc,'gene':gene,'sheet':ws.title,'row_1based':i+1,'column_1based':j+1,'annotation_value':v,
        'column_header_firstrow':str(data[0][j]),'not_proof_all_assays_or_unfiltered_detection':True,'not_a_replication_result':True})
       print('Target annotation row:',gene,row)
       if j==0 and v==gene:
        assert data[4][:4]==('Gene Symbol','SomaScan.CSF','Olink.CSF','MS.CSF.Depleted' if int(number[1])==4 else 'MS.CSF.Undepleted')
        coverage.append({'accession_from_current_study':acc,'gene':gene,'sheet':ws.title,'SomaScan_filtered_gene_present':bool(row[1]),'Olink_filtered_gene_present':bool(row[2]),
          'MS_filtered_gene_present':bool(row[3]),'NULISA_in_this_external_study':False,'reported_cohort_n':36,'three_assay_subject_upper_bound':35,
          'per_gene_complete_case_n_not_verified':True,'MS_peptide_identity_not_audited':True,'no_external_association_reanalysis':True})
 book.close();csvout(RESULTS/'EXTERNAL_EMORY2022_ANNOTATION_SHEETS.csv',schemas)
 csvout(RESULTS/'EXTERNAL_EMORY2022_CANDIDATE_ANNOTATION_HITS.csv',pd.DataFrame(hits,columns=['accession','gene','sheet','row_1based','column_1based','annotation_value','column_header_firstrow','not_proof_all_assays_or_unfiltered_detection','not_a_replication_result']))
 csvout(RESULTS/'EXTERNAL_EMORY2022_FILTERED_GENE_PLATFORM_COVERAGE.csv',coverage)
 log('External published Emory annotation checked only; no individual data analysis or candidate effect selection')

if __name__=='__main__':run()
