from config import *
import urllib.request,time

def download(acc):
 url=f'https://rest.uniprot.org/uniprotkb/{acc}.json';path=EVIDENCE/(acc+'_UNIPROT.json')
 if not path.exists():
  request=urllib.request.Request(url,headers={'User-Agent':'ADNI-measurement-audit/1.0','Accept':'application/json'})
  for attempt in range(3):
   try:
    with urllib.request.urlopen(request,timeout=30) as r:body=r.read()
    obj=json.loads(body);assert obj['primaryAccession']==acc
    path.write_bytes(body);break
   except Exception:
    if attempt==2:raise
    time.sleep(1)
 return {'accession':acc,'url':url,'retrieved_local':datetime.now().astimezone().isoformat(),'path':path.relative_to(ROOT).as_posix(),'sha256':sha(path),'bytes':path.stat().st_size}

if __name__=='__main__':
 inventory=pd.read_csv(RESULTS/'TARGETED_PROTEIN_COVERAGE_ALL70.csv');wanted=set(CANDIDATES)|set(inventory.loc[inventory.targeted_directly_covered,'accession'])
 rows=[]
 for acc in sorted(wanted):
  rows.append(download(acc));log('Official UniProt sequence/annotation cached: '+acc)
 csvout(RESULTS/'UNIPROT_SOURCE_MANIFEST.csv',rows)
