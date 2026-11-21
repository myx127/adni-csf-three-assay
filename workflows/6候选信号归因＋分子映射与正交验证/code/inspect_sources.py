from config import *
from collections import Counter
from pypdf import PdfReader
import logging
logging.getLogger('pypdf').setLevel(logging.ERROR)

def run():
 base=Path(INPUTS['I10']['path']).parent
 extras=[('V01',base/'EMORY_PEPTIDERATIOS_15May2026.csv'),('V02',base/'EMORY_PEPTIDERATIOS_SET2_15May2026.csv'),
  ('VM01',base/'ADNI_MethodsReport_EmoryRevised_20221019.pdf'),('VM02',base/'ADNI_Targeted_MS_Methods_EmoryJuly2023release.pdf')]
 checks=[]
 for k,row in INPUTS.items():
  digest=sha(row['path']);checks.append({'id':k,'path':row['path'],'sha256':digest,'expected_sha256':row['sha256'],'passed':digest==row['sha256']})
 m=pd.read_csv(S5/'STAGE5_INPUT_OUTPUT_SHA256.csv')
 for rel in ['results/MAIN_ANALYSIS_SAP_FROZEN.json','results/SECONDARY_SAP_FROZEN.json','results/PRIMARY_ALL43_RESULTS.csv',
  '00_restricted/CLINICAL_DESIGN_RESTRICTED.csv','00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz','00_restricted/DESIGNS_MAIN_RESTRICTED.pkl','code/model.py','code/diagnostics.py','code/simulation.py']:
  expected=m.loc[m.path.eq(rel),'sha256'].iloc[0];digest=sha(S5/rel)
  checks.append({'id':'S5/'+rel,'path':str(S5/rel),'sha256':digest,'expected_sha256':expected,'passed':digest==expected})
 csvout(RESULTS/'SOURCE_IDENTITY_CHECKS.csv',checks);assert all(c['passed'] for c in checks)
 rows=[{'input_id':k,'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)} for k,p in extras]
 manifest=RESULTS/'ADDITIONAL_INPUT_MANIFEST.csv'
 if manifest.exists():
  previous=pd.read_csv(manifest).set_index('input_id');assert all(previous.loc[r['input_id'],'sha256']==r['sha256'] for r in rows)
 else:csvout(manifest,rows)
 boundary=RESULTS/'ACCESS_BOUNDARY_INITIAL.json'
 if not boundary.exists():save_json(boundary,{'protocol_sha256':sha(ROOT/'STAGE6_PROTOCOL_v1.0.md'),'new_sources_registered_before_row_scan':True,'targeted_MS_associations_not_read':True,'Stage5_discoveries_already_known':True})
 else:assert json.loads(boundary.read_text(encoding='utf-8'))['protocol_sha256']==sha(ROOT/'STAGE6_PROTOCOL_v1.0.md')
 metadata=[];schema=[]
 for k,p in extras:
  if p.suffix.lower()=='.pdf':
   reader=PdfReader(p);text='\n'.join(f'[PDF page {i+1}]\n'+(page.extract_text() or '') for i,page in enumerate(reader.pages))
   (EVIDENCE/(k+'_METHODS_EXTRACTED.txt')).write_text(text,encoding='utf-8');continue
  with p.open(encoding='utf-8-sig',newline='') as f:
   reader=csv.DictReader(f);schema.append({'input_id':k,'columns':reader.fieldnames});counts=Counter()
   fields=['TYPE','TESTNAME','UNITS'] if k=='V01' else ['PEPTIDE','PEPTIDE_MOD_SEQUENCE','PROTEIN']
   n=0
   for row in reader:
    n+=1;counts[tuple(row.get(c,'').strip() for c in fields)]+=1
   for key,count in counts.items():metadata.append({'source':k,**dict(zip(fields,key)),'n_rows':count})
   print(k,'rows',n,'annotation_combinations',len(counts),'columns',len(reader.fieldnames))
 save_json(RESULTS/'TARGETED_FILE_SCHEMAS.json',schema);csvout(RESULTS/'TARGETED_ANNOTATION_INVENTORY.csv',metadata)
 dictionary=pd.read_csv(INPUTS['I03']['path'],dtype=str,keep_default_na=False)
 d=dictionary.loc[dictionary.TBLNAME.str.upper().str.contains('EMORY.*PEPTIDE',regex=True)]
 csvout(RESULTS/'TARGETED_DICTIONARY_FIELDS.csv',d)
 for name in ['model.py','diagnostics.py']:
  target=ROOT/'code'/name
  if not target.exists():target.write_bytes((S5/'code'/name).read_bytes())
  assert sha(target)==sha(S5/'code'/name)
 ref=ROOT/'code/simulation_reference.py'
 if not ref.exists():ref.write_bytes((S5/'code/simulation.py').read_bytes())
 assert sha(ref)==sha(S5/'code/simulation.py')
 save_json(RESULTS/'REFERENCE_CODE_ORIGIN.json',{'reference_stage':'Stage5','new_config_redirects_all_writes_to_Stage6':True,
  'source_sha256':{n:sha(S5/'code'/n) for n in ['model.py','diagnostics.py','simulation.py']}})
 log('New targeted sources registered and annotation-only inventory complete; no targeted quantitative associations examined')

if __name__=='__main__':run()
