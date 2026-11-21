from config import *
from collections import Counter,defaultdict
from datetime import datetime
import re

def dt(x):
 s=str(x or '').strip()
 for fmt in ['%Y-%m-%d','%m/%d/%Y','%Y/%m/%d','%Y-%m-%d %H:%M:%S']:
  try:return datetime.strptime(s,fmt).date()
  except ValueError:pass
 return None
def rid(x):
 try:
  v=float(x);return int(v) if np.isfinite(v) and v>0 and v.is_integer() else None
 except (ValueError,TypeError):return None
def numeric(x):
 try:v=float(x);return v if np.isfinite(v) else np.nan
 except (ValueError,TypeError):return np.nan
def parse_feature(row,key):
 if key=='V01':
  name=row['TESTNAME'].strip();parts=name.split('|');peptide=parts[0];protein='|'.join(parts[1:]);modified=peptide
 else:peptide=row['PEPTIDE'].strip();modified=row['PEPTIDE_MOD_SEQUENCE'].strip();protein=row['PROTEIN'].strip()
 match=re.search(r'(?:sp|tr)\|([^|]+)\|',protein);acc=match[1] if match else ''
 seq=re.sub(r'\[[^]]+\]','',peptide).replace('.','').strip()
 # Only fixed carbamidomethyl C is treated as a common peptide identity across releases.
 mod_class='FIXED_C_OR_UNMODIFIED' if re.sub(r'C\[\+57(?:\.0+)?\]','C',modified)==seq else modified
 return {'accession':acc,'protein_label':protein,'peptide':seq,'modified_sequence':modified,'mod_class':mod_class,
  'feature_key':acc+'|'+seq+'|'+mod_class}
def plate_from_id(s):
 text=str(s or '')
 for pat in [r'(?:^|_)B(\d+)[A-H]\d+',r'(?:^|_)P0*(\d+)[A-H]\d+',r'(?:^|_)SmpP0*(\d+)',r'(?:^|_)QCsP0*(\d+)']:
  m=re.search(pat,text,re.I)
  if m:return 'P'+str(int(m[1])).zfill(2)
 return ''

def run():
 extra=pd.read_csv(RESULTS/'ADDITIONAL_INPUT_MANIFEST.csv').set_index('input_id');inventory=[];private=[];structures=[];formats=[]
 clinical=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv');oldkeys=set(clinical.key);oldrids=set(clinical.RID)
 for source in ['V01','V02']:
  byfeature=defaultdict(list);count=Counter();formatcounts=Counter()
  with Path(extra.loc[source,'path']).open(encoding='utf-8-sig',newline='') as f:
   for line,row in enumerate(csv.DictReader(f),2):
    count['rows']+=1;feature=parse_feature(row,source);person=rid(row.get('RID'));day=dt(row.get('EXAMDATE'))
    typev=row.get('TYPE','');valid_person=person is not None and person!=999999 and day is not None and (source!='V01' or typev=='CSF')
    if valid_person:
     count['clinical_rows']+=1;k=str(person)+'|'+day.isoformat()
    else:count['QC_or_invalid_ID_date_rows']+=1;k=''
    # Values are reduced to positive/finite flags only, not associated with clinical P here.
    value=numeric(row.get('TESTVALUE' if source=='V01' else 'TOTAL_AREA_RATIO'));valid=np.isfinite(value) and value>0
    reps=[row.get('REPLICATE1','').strip(),row.get('REPLICATE2','').strip()]
    for i,rep in enumerate(reps,1):
     if rep:count[f'replicateID{i}_nonblank']+=1;count[f'replicateID{i}_numeric']+=int(np.isfinite(numeric(rep)))
    rawplate=row.get('PLATEPOSITION','') if source=='V01' else reps[0]
    plate=plate_from_id(rawplate);plate2=plate_from_id(reps[1]);plate_conflict=bool(plate and plate2 and plate!=plate2)
    if not plate:plate=plate2
    count['plate_label_available_rows']+=bool(plate);count['inhouse_ID_plate_conflict_rows']+=plate_conflict
    formatcounts['plate_field_shape_'+re.sub(r'\d','d',re.sub(r'[A-Za-z]','a',rawplate[:40]))]+=1
    out={'source':source,'source_row':line,'RID':person,'date':str(day or ''),'key':k,'valid_participant':valid_person,
      **feature,'positive_finite':valid,'in_stage5_key':k in oldkeys,'in_stage5_RID':person in oldrids,
      'plate':plate,'plate2':plate2,'plate_conflict':plate_conflict,'run_date':str(dt(row.get('RUNDATE')) or ''),
      'replicateID1_present':bool(reps[0]),'replicateID2_present':bool(reps[1]),'VISCODE2':row.get('VISCODE2',''),
      'PTID_RID_match':not re.search(r'_S_(\d+)$',row.get('PTID','')) or int(re.search(r'_S_(\d+)$',row['PTID'])[1])==person}
    private.append(out);byfeature[feature['feature_key']].append(out)
  for feat,rows in byfeature.items():
   first=rows[0];sample=[r for r in rows if r['valid_participant']];valid=[r for r in sample if r['positive_finite']]
   inventory.append({'source':source,**{k:first[k] for k in ['accession','protein_label','peptide','modified_sequence','mod_class','feature_key']},
    'rows_all':len(rows),'n_sample_rows':len(sample),'n_positive_rows':len(valid),'unique_RID':len({r['RID'] for r in sample}),
    'unique_sample_keys':len({r['key'] for r in sample}),'n_stage5_same_day_keys':len({r['key'] for r in valid if r['in_stage5_key']}),
    'n_new_RID':len({r['RID'] for r in valid if not r['in_stage5_RID']}),'n_plates':len({r['plate'] for r in sample if r['plate']})})
  structures.append({'source':source,**count})
  formats.extend({'source':source,'field_shape':k,'n_rows':v} for k,v in formatcounts.items())
 frame=pd.DataFrame(private);csvout(RESTRICTED/'TARGETED_VALIDITY_AND_TIMING_RESTRICTED.csv',frame)
 csvout(RESULTS/'TARGETED_FEATURE_VALIDITY_COVERAGE.csv',inventory);csvout(RESULTS/'TARGETED_SOURCE_STRUCTURE.csv',structures)
 # Only aggregate formatting patterns; remove rare potentially identifying patterns.
 csvout(RESULTS/'TARGETED_PLATE_FORMAT_COUNTS.csv',[r for r in formats if r['n_rows']>=50])
 meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv');s5=pd.read_csv(S5/'results/PRE_EFFECT_PROTEIN_SCREEN.csv').set_index('accession');features=pd.DataFrame(inventory)
 status=[]
 for r in meta.itertuples():
  matched=features.loc[features.accession.eq(r.accession)]
  status.append({'accession':r.accession,'gene':r.nulisa_target,'single_annotation':r.annotation_strict_candidate,
    's5_design_qualified':s5.loc[r.accession,'design_status']=='DESIGN_QUALIFIED','is_three_discoveries':r.accession in CANDIDATES,
    'targeted_directly_covered':len(matched)>0,'targeted_unique_features':matched.feature_key.nunique(),'targeted_unmodified_peptides':matched.peptide.nunique(),
    'targeted_sources':';'.join(sorted(matched.source.unique())),'direct_discovery_validation_status':'NOT_COVERED_NO_SUBSTITUTION' if r.accession in CANDIDATES and matched.empty else 'ANNOTATION_COVERED_NEEDS_DESIGN' if len(matched) else 'NOT_COVERED'})
 csvout(RESULTS/'TARGETED_PROTEIN_COVERAGE_ALL70.csv',status)
 duplicate=frame.loc[frame.valid_participant].groupby(['source','key','feature_key']).size()
 csvout(RESULTS/'TARGETED_DUPLICATE_STRUCTURE.csv',[{'source':source,'duplicate_key_feature_groups':int((part>1).sum()),'max_rows_per_sample_feature':int(part.max())} for source,part in duplicate.groupby(level=0)])
 log('Targeted validity/date/label coverage audited; targeted clinical associations still not examined')
 print(pd.DataFrame(status).query('targeted_directly_covered or is_three_discoveries')[['gene','s5_design_qualified','is_three_discoveries','targeted_directly_covered','targeted_unique_features']].to_string(index=False))

if __name__=='__main__':run()
