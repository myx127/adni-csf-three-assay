from config import *
import re,zipfile

def positions(sequence,peptide):
 return [m.start()+1 for m in re.finditer('(?='+re.escape(peptide)+')',sequence)] if peptide else []

def run():
 meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv');annot=pd.read_csv(INPUTS['I05']['path'],dtype=str,keep_default_na=False)
 selected=[]
 for acc,gene in CANDIDATES.items():
  m=meta.loc[meta.accession.eq(acc)].iloc[0];fields=str(m.soma_features).split('|');matched=annot.loc[annot.Analytes.isin(fields)]
  for row in matched.to_dict('records'):selected.append({'candidate_accession':acc,'candidate_gene':gene,**row})
 csvout(RESULTS/'CANDIDATE_SOMA_ORIGINAL_ANNOTATION.csv',selected)
 uniprot={p.name.split('_')[0]:json.loads(p.read_text(encoding='utf-8')) for p in EVIDENCE.glob('*_UNIPROT.json')}
 features=pd.read_csv(RESULTS/'TARGETED_FEATURE_VALIDITY_COVERAGE.csv').drop_duplicates(['accession','feature_key'])
 mappings=[];regions=[];entryrows=[]
 for acc,obj in uniprot.items():
  seq=obj['sequence']['value'];audit=obj.get('entryAudit',{});names=obj.get('genes',[{}]);gene=names[0].get('geneName',{}).get('value','')
  entryrows.append({'accession':acc,'gene':gene,'length_aa':len(seq),'entry_version':audit.get('entryVersion'),'sequence_version':audit.get('sequenceVersion'),
   'sequence_update_date':audit.get('lastSequenceUpdateDate'),'annotation_update_date':audit.get('lastAnnotationUpdateDate'),
   'UniProt_entry':f'https://www.uniprot.org/uniprotkb/{acc}/entry','canonical_sequence_sha256':hashlib.sha256(seq.encode()).hexdigest()})
  for f in obj.get('features',[]):
   if f['type'] in ['Signal','Propeptide','Chain','Topological domain','Transmembrane','Domain','Region','Active site','Alternative sequence']:
    loc=f['location'];regions.append({'accession':acc,'gene':gene,'feature_type':f['type'],'start':loc['start'].get('value'),
       'end':loc['end'].get('value'),'description':f.get('description',''),'is_annotation_not_assay_binding_evidence':True})
 for f in features.to_dict('records'):
  acc=f['accession']
  if acc not in uniprot:continue
  pep=f['peptide'];seq=uniprot[acc]['sequence']['value'];exact=positions(seq,pep);IL=positions(seq.replace('I','L'),pep.replace('I','L'))
  other=[]
  for otheracc,obj in uniprot.items():
   if otheracc!=acc and positions(obj['sequence']['value'],pep):other.append(otheracc)
  mappings.append({'accession':acc,'peptide':pep,'feature_key':f['feature_key'],'canonical_exact_n_matches':len(exact),'canonical_exact_start_1based':';'.join(map(str,exact)),
    'canonical_exact_end_1based':';'.join(str(p+len(pep)-1) for p in exact),'I_L_equivalent_n_matches':len(IL),
    'other_downloaded_entries_exact_match':';'.join(other),'whole_human_proteome_uniqueness_tested':False,
    'sequence_match_does_not_prove_exclusive_isoform':True})
 csvout(RESULTS/'UNIPROT_ENTRY_AUDIT.csv',entryrows);csvout(RESULTS/'UNIPROT_REGIONS_AND_PROCESSING.csv',regions);csvout(RESULTS/'TARGETED_PEPTIDE_CANONICAL_MAPPING.csv',mappings)
 matrix=[]
 legacy=pd.read_csv(RESULTS/'LEGACY_CANDIDATE_TRANSITION_REGISTRY.csv') if (RESULTS/'LEGACY_CANDIDATE_TRANSITION_REGISTRY.csv').exists() else pd.DataFrame(columns=['accession','peptide'])
 for acc,gene in CANDIDATES.items():
  m=meta.loc[meta.accession.eq(acc)].iloc[0];obj=uniprot[acc];ft=[f['type']+': '+f.get('description','') for f in obj.get('features',[]) if f['type'] in ['Signal','Propeptide','Transmembrane','Topological domain']]
  matrix.append({'accession':acc,'gene':gene,'nulisa_label':m.nulisa_target,'nulisa_raw_accession':m.nulisa_raw_accession_fields,
   'soma_features':m.soma_features,'soma_raw_accession':m.soma_raw_accession_fields,'TMT_aggregate_column':m.tmt_columns,
   'Emory_2022_2023_direct_peptides':0,'Caprion_2019_direct_peptides':int(legacy.loc[legacy.accession.eq(acc),'peptide'].nunique()),'antibody_epitope_available':False,'aptamer_binding_region_available':False,
   'original_TMT_peptide_level_evidence_available':False,'specific_mature_fragment_or_modification_verified':False,
   'official_sequence_length':obj['sequence']['length'],'annotated_processing_or_topology':'; '.join(ft),
   'molecular_interpretation':'Protein accession matches are annotation-level; actual bound/digested molecular species not identified',
   'direct_new_assay_validation':'SEE_FINAL_DIRECT_DISCOVERY_VALIDATION_STATUS; no formal relative-pattern replication completed'})
 csvout(RESULTS/'THREE_DISCOVERY_MOLECULAR_EVIDENCE_MATRIX.csv',matrix)
 with zipfile.ZipFile(INPUTS['I07']['path']) as z:
  csvout(RESULTS/'TMT_ARCHIVE_MEMBER_INVENTORY.csv',[{'member':n,'bytes':z.getinfo(n).file_size,'is_protein_aggregate_matrix':Path(n).name=='EMORY_CSF_TMT_MS.csv'} for n in z.namelist()])
 log('Official sequence mapping and three-candidate molecular evidence matrix complete; assay epitopes/peptide evidence not invented')

if __name__=='__main__':run()
