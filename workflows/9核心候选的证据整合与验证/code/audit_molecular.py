from config import *
import gzip,re,requests
from concurrent.futures import ThreadPoolExecutor

def download(acc):
    url=f'https://rest.uniprot.org/uniprotkb/{acc}.json';out=EVIDENCE/f'{acc}_UNIPROT_REFRESH.json'
    record={'accession':acc,'url':url,'retrieved_time':now()}
    if out.exists() and '--refresh' not in sys.argv:
        obj=json.loads(out.read_text(encoding='utf-8'))
        snapshot=EVIDENCE/'UNIPROT_FETCH_SNAPSHOT.json'
        original=next((r for r in json.loads(snapshot.read_text(encoding='utf-8')) if r['accession']==acc),{}) if snapshot.exists() else {}
        record.update({'status':'FROZEN_STAGE9_OFFICIAL_REST_SNAPSHOT','path':str(out),'bytes':out.stat().st_size,'sha256':sha(out),'retrieved_time':original.get('retrieved_time'),'loaded_time':now(),'release':original.get('release')})
        return acc,obj,record
    try:
        response=requests.get(url,timeout=20);response.raise_for_status();obj=response.json()
        assert obj['primaryAccession']==acc;out.write_bytes(response.content)
        record.update({'status':'REFRESHED_OFFICIAL_REST','path':str(out),'bytes':out.stat().st_size,'sha256':sha(out),'release':response.headers.get('X-UniProt-Release'),'release_date':response.headers.get('X-UniProt-Release-Date')})
    except Exception as e:
        out=S6/'evidence'/f'{acc}_UNIPROT.json';obj=json.loads(out.read_text(encoding='utf-8'))
        record.update({'status':'FROZEN_STAGE6_OFFICIAL_CACHE_USED','path':str(out),'bytes':out.stat().st_size,'sha256':sha(out),'error':str(e)})
    return acc,obj,record

def run():
    downloaded=list(ThreadPoolExecutor(max_workers=3).map(download,CANDIDATES))
    entries={a:o for a,o,_ in downloaded};csvout(RESULTS/'UNIPROT_REFRESH_SOURCE_MANIFEST.csv',[r for _,_,r in downloaded])
    audit=[];regions=[];peptide=[]
    for acc,obj in entries.items():
        old=json.loads((S6/'evidence'/f'{acc}_UNIPROT.json').read_text(encoding='utf-8'))
        seq=obj['sequence']['value'];audit.append({'accession':acc,'gene':CANDIDATES[acc],'length_aa':len(seq),'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),'same_sequence_as_stage6':seq==old['sequence']['value'],**obj.get('entryAudit',{})})
        for f in obj.get('features',[]):
            if f['type'] in ['Signal','Propeptide','Chain','Topological domain','Transmembrane','Domain','Active site','Alternative sequence']:
                loc=f['location'];regions.append({'accession':acc,'gene':CANDIDATES[acc],'feature_type':f['type'],'start':loc['start'].get('value'),'end':loc['end'].get('value'),'description':f.get('description',''),'annotation_not_assay_binding_site':True})
        for pep in PEPTIDES.get(CANDIDATES[acc],[]):
            positions=[m.start()+1 for m in re.finditer('(?='+pep+')',seq)]
            assert positions
            features=[]
            for f in regions:
                if f['accession']==acc and isinstance(f['start'],int) and isinstance(f['end'],int) and any(f['start']<=p and p+len(pep)-1<=f['end'] for p in positions):features.append(f['feature_type']+':'+f['description'])
            peptide.append({'accession':acc,'gene':CANDIDATES[acc],'peptide':pep,'canonical_positions_1based':'|'.join(str(p) for p in positions),'canonical_end_1based':'|'.join(str(p+len(pep)-1) for p in positions),'overlapping_annotated_regions':'|'.join(features),'activation_or_cleavage_junction_verified':False,'intact_protein_or_activity_verified':False})
    csvout(RESULTS/'UNIPROT_SEQUENCE_VERSION_AUDIT.csv',audit);csvout(RESULTS/'CORE_UNIPROT_REGIONS.csv',regions)
    # Independently stream every sequence in the frozen human canonical+isoform reference.
    reference=S6/'evidence/HUMAN_UP000005640_CANONICAL_AND_ISOFORMS.fasta.gz'
    wanted=sum(PEPTIDES.values(),[]);matches=[];n=0
    def consume(header,seq):
        nonlocal n
        if not header:return
        n+=1;parts=header.split('|');acc=parts[1] if len(parts)>2 else header.split()[0];gm=re.search(r'\bGN=([^\s]+)',header);gene=gm.group(1) if gm else ''
        for mode,s in [('EXACT',seq),('IL_EQUIVALENT',seq.replace('I','L'))]:
            for pep in wanted:
                target=pep if mode=='EXACT' else pep.replace('I','L')
                positions=[m.start()+1 for m in re.finditer('(?='+target+')',s)]
                if positions:matches.append({'peptide':pep,'match_mode':mode,'matched_accession':acc,'matched_gene':gene,'positions_1based':'|'.join(map(str,positions))})
    with gzip.open(reference,'rt',encoding='utf-8') as f:
        header='';parts=[]
        for line in f:
            if line.startswith('>'):
                consume(header,''.join(parts));header=line[1:].strip();parts=[]
            else:parts.append(line.strip())
        consume(header,''.join(parts))
    csvout(RESULTS/'INDEPENDENT_HUMAN_REFERENCE_ALL_MATCHES.csv',matches)
    mm=pd.DataFrame(matches);summ=[]
    previous=pd.read_csv(S6/'results/HUMAN_REFERENCE_PEPTIDE_SPECIFICITY.csv')
    checks=[]
    for row in peptide:
        pep=row['peptide'];expected_gene=row['gene']
        for mode in ['EXACT','IL_EQUIVALENT']:
            part=mm.loc[mm.peptide.eq(pep)&mm.match_mode.eq(mode)]
            row[mode+'_reference_sequence_matches']=len(part)
            row[mode+'_other_or_unannotated_genes']='|'.join(sorted(set(part.matched_gene)-{expected_gene}))
            row[mode+'_gene_specific_in_frozen_reference']=bool(len(part)>0 and part.matched_gene.eq(expected_gene).all())
            old=previous.loc[previous.peptide.eq(pep)].iloc[0];col='EXACT_sequence_matches' if mode=='EXACT' else 'I_L_EQUIVALENT_sequence_matches'
            checks.append({'check':pep+'_'+mode+'_reference_matches','actual':len(part),'expected':old[col],'passed':len(part)==old[col]})
        row['reference_sequences_screened']=n;row['reference_is_stage6_frozen_version_not_claimed_current_release']=True
    csvout(RESULTS/'CORE_PEPTIDE_SEQUENCE_REGION_SPECIFICITY.csv',peptide)
    csvout(RESULTS/'INDEPENDENT_MOLECULAR_CHECKS.csv',checks)
    assert all(x['passed'] for x in checks)
    save_json(RESULTS/'FROZEN_REFERENCE_SCAN_RECORD.json',{'path':str(reference),'sha256':sha(reference),'sequences':n,'all_five_peptides_both_modes_scanned':True,'timestamp':now(),'scope':'UniProt human UP000005640 frozen Stage6 canonical+isoform reference; specificity limited to this set'})
    log(f'Molecular audit:3 official entries;{n} frozen human sequences scanned independently;five peptides x2 modes reproduced')
    print(pd.DataFrame(audit).to_string(index=False))
    print(pd.DataFrame(peptide)[['gene','peptide','canonical_positions_1based','EXACT_reference_sequence_matches','IL_EQUIVALENT_reference_sequence_matches']].to_string(index=False))

if __name__=='__main__':
    assert_unsealed()
    run()
