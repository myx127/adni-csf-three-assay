from pathlib import Path
import csv,json,hashlib,re,logging
import pandas as pd
from pypdf import PdfReader
from inspect_inputs import ROOT,PROJECT,RESULTS,MANIFEST

def main():
 inputs={r['input_id']:r for r in csv.DictReader(MANIFEST.open(encoding='utf-8-sig'))}
 source=Path(inputs['I10']['path']).parent.parent
 extras=[('E01',source/'CSF/UPENNBIOMK_ROCHE_ELECSYS_15May2026.csv'),
 ('M01',source/'CSF/UPENNBIOMK_MASTER_METHODS_20231109.txt'),
 ('M02',source/'CSF/UPENNBIOMK_ROCHE_ELECSYS_METHODS_20231109.txt'),
 ('M03',source/'PET_UC Berkeley_Analysis/ADNI_UCBerkeley_AmyloidPET_Methods_v3_2025-06-30.pdf'),
 ('M04',source/'PET_UC Berkeley_Analysis/ADNI_UCBerkeley_TauPET_Methods_2025-06-30.pdf'),
 ('M05',source/'Neuropsychological/ADNI_Methods_UWNPSYCHSUM_20231018.pdf')]
 rows=[];notes=[]
 logging.getLogger('pypdf').setLevel(logging.ERROR)
 for key,path in extras:
  if not path.exists():raise FileNotFoundError(path)
  rows.append({'input_id':key,'path':str(path),'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
 registry=RESULTS/'ADDITIONAL_INPUT_MANIFEST.csv'
 if registry.exists():
  prior={r['input_id']:r for r in csv.DictReader(registry.open(encoding='utf-8-sig'))}
  assert all(r['sha256']==prior[r['input_id']]['sha256'] and int(prior[r['input_id']]['bytes'])==r['bytes'] for r in rows),'New source changed; do not silently rebind'
 else:pd.DataFrame(rows).to_csv(registry,index=False,encoding='utf-8-sig')
 # Register or verify all additions before accessing their headers/content.
 for key,path in extras:
  if key=='E01':
   with path.open(encoding='utf-8-sig',newline='') as f:header=next(csv.reader(f))
   (RESULTS/'ELECSYS_SCHEMA.json').write_text(json.dumps(header,indent=2),encoding='utf-8');print('E01',header)
   continue
  if path.suffix=='.txt':
   text=path.read_text(encoding='utf-8-sig',errors='replace')
  else:
   reader=PdfReader(path);text='\n'.join(f'[PDF page {i+1}]\n'+(p.extract_text() or '') for i,p in enumerate(reader.pages))
  (RESULTS/(key+'_METHODS_EXTRACTED.txt')).write_text(text,encoding='utf-8')
  lines=text.splitlines();terms=['cutoff','cut-off','threshold','QC_FLAG','QC flag','qc_flag','centiloid','META_TEMPORAL','rescal','repeat','median','first','last update','PHC','2023','-4','1.3','23','192']
  selected=set()
  for i,line in enumerate(lines):
   if any(term.lower() in line.lower() for term in terms):selected.update(range(max(0,i-1),min(len(lines),i+3)))
  excerpt='\n'.join(lines[i] for i in sorted(selected))
  notes.append('# '+key+' '+path.name+'\n\n'+excerpt)
  print(key,'extracted_characters',len(text),'evidence_lines',len(selected))
 (RESULTS/'METHODS_EVIDENCE_NOTES.md').write_text('\n\n'.join(notes),encoding='utf-8')
 dictionary=pd.read_csv(inputs['I03']['path'],dtype=str,keep_default_na=False)
 tables=['UPENNBIOMK_MASTER','UPENNBIOMK_ROCHE_ELECSYS','UCBERKELEY_AMY_6MM','UCBERKELEY_TAU_6MM','UWNPSYCHSUM','PTDEMOG','APOERES','DXSUM']
 fields=['ABETA','PTAU','TAU','ABETA42','PTAU181','BATCH','KIT','STDS','DRAWDTE','RUNDATE','AMYLOID_STATUS','qc_flag','QC_FLAG','TRACER','META_TEMPORAL_SUVR','CENTILOIDS','PTDOB','PTDOBYY','PTGENDER','PTEDUCAT','GENOTYPE','APUSABLE','DIAGNOSIS','HAS_QC_ERROR','ADNI_MEM']
 print('Dictionary columns',dictionary.columns.to_list())
 subset=dictionary.loc[dictionary.TBLNAME.str.upper().isin(tables)]
 subset.to_csv(RESULTS/'RELEVANT_DICTIONARY_FIELDS.csv',index=False,encoding='utf-8-sig')
 print('Dictionary tables',subset.TBLNAME.value_counts().to_dict())

if __name__=='__main__':main()
