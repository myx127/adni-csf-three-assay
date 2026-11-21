from config import *
from targeted_inventory import dt,rid,numeric
from collections import Counter,defaultdict
import zipfile,re
from openpyxl import load_workbook

legacy=pd.read_csv(RESULTS/'LEGACY_INPUT_MANIFEST.csv').set_index('input_id');old=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv')
oldkeys=set(old.key);oldrids=set(old.RID);labels=[];structures=[];private=[]
with zipfile.ZipFile(legacy.loc['L04','path']) as z:
 member='csfproteomics/adni_csfproteomics.csv'
 with z.open(member) as f:header=next(csv.reader(__import__('io').TextIOWrapper(f,encoding='utf-8-sig')))
 candidates=[c for c in header if re.match(r'(NPTX1|KLK6)_HUMAN\.',c)]
 metadata=[c for c in header if '_HUMAN.' not in c and c not in ['CMGA','FABPH','NPTX2','SCG2','VGF']]
 with z.open(member) as f:frame=pd.read_csv(f,usecols=list(dict.fromkeys(metadata+candidates)),low_memory=False)
 for col in candidates:
  parts=col.split('.');gene=parts[0].split('_')[0];pep=parts[1];acc='Q15818' if gene=='NPTX1' else 'Q92876'
  v=pd.to_numeric(frame[col],errors='coerce');labels.append({'source':'L04','accession':acc,'gene':gene,'peptide':pep,'transition_field':col,'released_scale':'normalized log2 per dictionary; primer also says natural-log; use directly',
    'n_finite':int(np.isfinite(v).sum()),'literal_zero_count':int(v.eq(0).sum()),'literal_minus4_count':int(v.eq(-4).sum()),'both_are_not_automatically_missing_without_documented_rule':True})
 counts=Counter();cats=defaultdict(Counter)
 for line,row in enumerate(frame.to_dict('records'),2):
  r=rid(row.get('RID'));date=dt(row.get('EXAMDATE'));valid=r is not None and r!=999999 and date is not None
  k=str(r)+'|'+date.isoformat() if valid else ''
  for field in ['INSTRUMENT','RUN.NUMBER','SUBRUN','MS.PLATE.INDEX','PINK']:cats[field][str(row.get(field,''))]+=1
  out={field:row.get(field,'') for field in metadata};out.update({'source_row':line,'RID':r,'date':str(date or ''),'key':k,'valid_participant':valid,
    'in_stage5_key':k in oldkeys,'in_stage5_RID':r in oldrids})
  for col in candidates:out[col+'_finite']=np.isfinite(numeric(row.get(col)))
  private.append(out);counts['rows']+=1;counts['valid_participant_rows']+=valid
 for field,c in cats.items():
  for value,n in c.items():structures.append({'source':'L04','field':field,'category':value,'n_rows':n})
csvout(RESULTS/'LEGACY_CANDIDATE_TRANSITION_REGISTRY.csv',labels);csvout(RESULTS/'LEGACY_MRM_TECHNICAL_METADATA_COUNTS.csv',structures)
csvout(RESTRICTED/'LEGACY_L04_VALIDITY_TIMING_RESTRICTED.csv',private)
data=pd.DataFrame(private);summary=[{'source':'L04','rows':len(data),'valid_participant_rows':int(data.valid_participant.sum()),'unique_RID':data.loc[data.valid_participant,'RID'].nunique(),
 'unique_RID_date':data.loc[data.valid_participant,'key'].nunique(),'n_stage5_same_day_keys':data.loc[data.in_stage5_key,'key'].nunique(),
 'n_new_RID':data.loc[data.valid_participant&~data.in_stage5_RID,'RID'].nunique(),'duplicate_RID_date_groups':int((data.loc[data.valid_participant].groupby('key').size()>1).sum())}]
with Path(legacy.loc['L01','path']).open(encoding='utf-8-sig',newline='') as f:
 reader=csv.DictReader(f);h=reader.fieldnames;rr=list(reader)
 summary.append({'source':'L01','rows':len(rr),'unique_RID':len({r['RID'] for r in rr}),'date_in_release':False,'direct_candidate_fields':';'.join(c for c in h if c.startswith('KLK6_') or c.startswith('NPTX1_'))})
book=load_workbook(legacy.loc['L02','path'],read_only=True,data_only=True)
trans=list(book['Transitions'].iter_rows(values_only=True));header=trans[0];tt=[dict(zip(header,r)) for r in trans[1:]]
csvout(RESULTS/'LEGACY_2014_CANDIDATE_TRANSITION_ANNOTATIONS.csv',[r for r in tt if 'NPTX1' in str(r.get('Protein')) or 'KLK6' in str(r.get('Protein'))])
book.close()
with zipfile.ZipFile(legacy.loc['L05','path']) as z:
 member=next(n for n in z.namelist() if 'QC Multiplex' in n and n.endswith('.csv'))
 with z.open(member) as f:h=next(csv.reader(__import__('io').TextIOWrapper(f,encoding='utf-8-sig')))
 csvout(RESULTS/'LEGACY_MULTIPLEX_PDGF_LABEL_CHECK.csv',[{'source':'L05','published_analyte_label':c,'is_PDGF_ligand_not_PDGFRB_receptor':'Growth Factor' in c or 'PDGF-BB' in c,'direct_PDGFRB_confirmed':False} for c in h if 'PDGF' in c.upper() or 'PLATELET' in c.upper()])
csvout(RESULTS/'LEGACY_SOURCE_SAMPLE_STRUCTURE.csv',summary)
log('Legacy candidate peptide/transition and sample metadata audit complete; no legacy pathology association fitted')
print(pd.DataFrame(summary).to_string(index=False));print(pd.DataFrame(structures).query("field in ['RUN.NUMBER','MS.PLATE.INDEX','PINK']").to_string(index=False))
