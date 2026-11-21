from config import *
from openpyxl import load_workbook
import tarfile,zipfile,re,importlib.util

def run():
 legacy=pd.read_csv(S6/'results/LEGACY_INPUT_MANIFEST.csv').set_index('input_id');base=Path(legacy.loc['L01','path']).parents[1]
 wanted=[('A01',base/'ADNIMERGE/ADNIMERGE2.tar.gz'),('A02',base/'Biospecimen_Results.zip')]
 sources=[];members=[]
 for ident,path in wanted:
  if not path.exists():continue
  sources.append({'source_id':ident,'path':str(path),'bytes':path.stat().st_size,'sha256':sha(path)})
 csvout(RESULTS/'ADDITIONAL_METADATA_INPUTS.csv',sources)
 schema=[]
 book=load_workbook(legacy.loc['L02','path'],read_only=True,data_only=True)
 for ws in book:
  it=ws.iter_rows(values_only=True);header=next(it);schema.append({'sheet':ws.title,'rows':ws.max_row,'columns':ws.max_column,'leading_header_types':[type(x).__name__ for x in header[:7]],
   'has_date_named_column':any(re.search('DATE|COLLECT',str(x),re.I) for x in header)})
  if ws.title.replace(' ','')=='SampleInfo':
   info=pd.DataFrame(list(it),columns=header);info.to_csv(PRIVATE/'MRM2014_SAMPLEINFO_RESTRICTED.csv',index=False,encoding='utf-8-sig')
   csvout(RESULTS/'MRM2014_SAMPLEINFO_FIELDS.csv',[{'field':str(c),'nonmissing':int(info[c].notna().sum()),'n_unique':int(info[c].nunique()),
    'is_date_named':bool(re.search('DATE|COLLECT',str(c),re.I))} for c in info])
  if ws.title=='CSV-Export':
   export=pd.DataFrame(list(it),columns=header);export.to_csv(PRIVATE/'MRM2014_CSV_EXPORT_RESTRICTED.csv',index=False,encoding='utf-8-sig')
 book.close();save_json(RESULTS/'MRM2014_WORKBOOK_STRUCTURE.json',schema)
 for source,path in wanted:
  if not path.exists():continue
  if path.suffix=='.zip':
   with zipfile.ZipFile(path) as z:
    members.extend({'source':source,'member':m,'bytes':z.getinfo(m).file_size} for m in z.namelist())
  else:
   with tarfile.open(path,'r:gz') as t:members.extend({'source':source,'member':m.name,'bytes':m.size} for m in t.getmembers())
 csvout(RESULTS/'METADATA_ARCHIVE_MEMBERS.csv',members)
 # Read only relevant fields of existing source dictionaries and CSF metadata.
 fields=[];extra=[]
 for name in ['UPENNBIOMK_MASTER_15May2026.csv','UPENNBIOMK_ROCHE_ELECSYS_15May2026.csv','UPENNBIOMKADNIDIAN2017_15May2026.csv']:
  path=Path(legacy.loc['L01','path']).parent/name
  if path.exists():
   with path.open(encoding='utf-8-sig',newline='') as f:h=next(csv.reader(f))
   fields.append({'source':name,'columns':h,'collection_or_specimen_fields':[c for c in h if re.search('LIMS|SAMP|SPEC|EXAMDATE|VISCODE|RID',c,re.I)]})
   if name not in [Path(x['path']).name for x in INPUTS.values()]:extra.append({'source_id':'META_'+name,'path':str(path),'bytes':path.stat().st_size,'sha256':sha(path)})
 if extra:csvout(RESULTS/'ADDITIONAL_METADATA_INPUTS.csv',sources+extra)
 save_json(RESULTS/'CSF_METADATA_FIELD_SCHEMAS.json',fields)
 print(pd.DataFrame(schema)[['sheet','rows','columns','has_date_named_column']].to_string(index=False))
 print('SampleInfo fields:',list(info.columns));print('CSV-Export leading:',list(export.columns)[:9]);print('pyreadr available:',bool(importlib.util.find_spec('pyreadr')))
 print('Relevant package members:',[r['member'] for r in members if re.search('MRM|BIOSPEC|biomk|SAMPLE|specimen|adnimerge2\.rda|csf',r['member'],re.I)][:60])
 log('2014 workbook and additional local metadata registered/audited without fitting candidate pathology effects')

if __name__=='__main__':run()
