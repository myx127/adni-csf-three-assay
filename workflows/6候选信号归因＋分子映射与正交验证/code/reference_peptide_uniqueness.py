"""Sequence specificity within the official human reference proteome, including isoforms.
Not proof of analytical specificity, proteolytic form, modification or epitope identity.
"""
from config import *
import urllib.request,urllib.parse,re,gzip

def run():
 path=EVIDENCE/'HUMAN_UP000005640_CANONICAL_AND_ISOFORMS.fasta.gz'
 params={'query':'proteome:UP000005640','format':'fasta','includeIsoform':'true','compressed':'true'}
 url='https://rest.uniprot.org/uniprotkb/stream?'+urllib.parse.urlencode(params)
 metadata=RESULTS/'HUMAN_REFERENCE_PROTEOME_MANIFEST.json'
 if not path.exists():
  req=urllib.request.Request(url,headers={'User-Agent':'ADNI-sequence-specificity-audit/1.0'})
  with urllib.request.urlopen(req,timeout=60) as r:
   headers=dict(r.headers);parts=[]
   while True:
    block=r.read(1024*1024)
    if not block:break
    parts.append(block)
  body=b''.join(parts);assert body[:2]==b'\x1f\x8b';path.write_bytes(body)
  save_json(metadata,{'url':url,'retrieved_local':datetime.now().astimezone().isoformat(),'path':path.relative_to(ROOT).as_posix(),
    'sha256':sha(path),'bytes':len(body),'response_headers':headers,
    'scope':'UniProt human reference proteome UP000005640 canonical + isoform FASTA; not all species, contaminants, or all TrEMBL entries outside this reference proteome'})
 else:assert sha(path)==json.loads(metadata.read_text(encoding='utf-8'))['sha256']
 modern=pd.read_csv(RESULTS/'TARGETED_FEATURE_VALIDITY_COVERAGE.csv')[['accession','peptide']].drop_duplicates()
 legacy=pd.read_csv(RESULTS/'LEGACY_CANDIDATE_TRANSITION_REGISTRY.csv')[['accession','peptide']].drop_duplicates()
 canonical={p.name.split('_')[0]:json.loads(p.read_text(encoding='utf-8')) for p in EVIDENCE.glob('*_UNIPROT.json')}
 targets=pd.concat([modern.assign(panel='EMORY_2022_2023'),legacy.assign(panel='CAPRION_2019')]).drop_duplicates()
 targets=targets.loc[targets.accession.isin(canonical)].sort_values(['accession','peptide','panel'])
 entries=[];header=None;seq=[]
 with gzip.open(path,'rt',encoding='utf-8') as f:
  for line in f:
   if line.startswith('>'):
    if header is not None:entries.append((header,''.join(seq)))
    header=line[1:].strip();seq=[]
   else:seq.append(line.strip())
 if header is not None:entries.append((header,''.join(seq)))
 rows=[];hits=[]
 for row in targets.itertuples():
  gene=canonical[row.accession]['genes'][0]['geneName']['value'];seq=canonical[row.accession]['sequence']['value'];starts=[m.start()+1 for m in re.finditer('(?='+re.escape(row.peptide)+')',seq)]
  matched={'EXACT':[],'I_L_EQUIVALENT':[]}
  for h,s in entries:
   acc=h.split('|')[1];m=re.search(r'\bGN=(\S+)',h);g=m[1] if m else '';base=acc.split('-')[0]
   for mode,ss,pp in [('EXACT',s,row.peptide),('I_L_EQUIVALENT',s.replace('I','L'),row.peptide.replace('I','L'))]:
    if pp in ss:
     matched[mode].append((acc,g));hits.append({'target_accession':row.accession,'target_gene':gene,'peptide':row.peptide,'match_mode':mode,
      'matched_accession':acc,'matched_gene':g,'same_annotated_gene':g==gene,'same_base_accession':base==row.accession})
  summary={'accession':row.accession,'gene':gene,'peptide':row.peptide,'panel':row.panel,'canonical_start_1based':';'.join(map(str,starts)),
   'canonical_end_1based':';'.join(str(x+len(row.peptide)-1) for x in starts),'reference_proteome_sequences_screened':len(entries),
   'source_reference_proteome':'UP000005640 canonical + isoforms','not_proof_of_epitope_or_mature_species':True}
  for mode,found in matched.items():
   other=sorted({g or '[NO_GENE_ANNOTATION]:'+a for a,g in found if g!=gene})
   summary[mode+'_sequence_matches']=len(found);summary[mode+'_other_genes_or_unannotated_entries']=';'.join(other)
   summary[mode+'_gene_specific_within_reference']=bool(found) and not other
   summary[mode+'_unique_isoform']=len(found)==1
  rows.append(summary)
 csvout(RESULTS/'HUMAN_REFERENCE_PEPTIDE_SPECIFICITY.csv',rows);csvout(RESULTS/'HUMAN_REFERENCE_PEPTIDE_ALL_MATCHES.csv',hits)
 log(f'Official reference-proteome peptide scan complete: {len(entries)} sequences, {len(rows)} peptide targets; sequence specificity does not establish assay specificity')
 print(pd.DataFrame(rows)[['gene','peptide','EXACT_sequence_matches','I_L_EQUIVALENT_sequence_matches','I_L_EQUIVALENT_gene_specific_within_reference','I_L_EQUIVALENT_unique_isoform']].to_string(index=False))

if __name__=='__main__':run()
