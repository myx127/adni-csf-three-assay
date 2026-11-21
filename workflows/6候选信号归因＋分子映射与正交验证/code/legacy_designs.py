from config import *
from designs import make
from targeted_inventory import dt,rid,numeric
from collections import defaultdict
import pickle,re

def clinical_info(keys):
 old=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv').set_index('key');demo=defaultdict(list);ptau=defaultdict(list)
 for row in csv.DictReader(Path(INPUTS['I13']['path']).open(encoding='utf-8-sig')):
  r=rid(row.get('RID'))
  if not r or row.get('HAS_QC_ERROR')=='1':continue
  m=re.fullmatch(r'(\d{1,2})/(\d{4})',row.get('PTDOB','').strip());full=dt(row.get('PTDOB'));yd=dt(row.get('PTDOBYY'));yn=numeric(row.get('PTDOBYY'))
  y=full.year if full else int(m[2]) if m else int(yn) if np.isfinite(yn) else yd.year if yd else None;mo=full.month if full else int(m[1]) if m else None;yy=int(yn) if np.isfinite(yn) else yd.year if yd else None
  demo[r].append((y,mo,row.get('PTGENDER'),bool(y and yy and y!=yy)))
 for line,row in enumerate(csv.DictReader(Path(INPUTS['E01']['path']).open(encoding='utf-8-sig')),2):
  r=rid(row.get('RID'));day=dt(row.get('EXAMDATE'));p=numeric(row.get('PTAU'));batch=row.get('BATCH','');m=re.search(r'(\d+)$',batch)
  if r and day and 8<=p<=120:ptau[str(r)+'|'+day.isoformat()].append((int(m[1]) if m else 0,line,p,batch))
 scale=json.loads((S5/'00_restricted/P_LOG_SCALING_RESTRICTED.json').read_text(encoding='utf-8'));rows=[]
 for key in keys:
  r=int(key.split('|')[0]);day=dt(key.split('|')[1]);d=demo[r];years={v[0] for v in d if v[0]};months={v[1] for v in d if v[1]};sex={v[2] for v in d if v[2] in ['1','2']}
  age=day.year-next(iter(years))+(day.month-next(iter(months)))/12 if len(years)==len(months)==1 and not any(v[3] for v in d) else np.nan
  p=min(ptau[key]) if ptau[key] else None
  row={'key':key,'RID':r,'date':day.isoformat(),'age_month_precision':age,'sex_code2':float(next(iter(sex))=='2') if len(sex)==1 else np.nan,
   'P_LOG':(np.log2(p[2])-scale['mean'])/scale['SD'] if p else np.nan,'ptau_batch':p[3] if p else '', 'is_stage5_RID':r in set(old.RID)}
  rows.append(row)
 return pd.DataFrame(rows)

def run():
 metadata=pd.read_csv(RESTRICTED/'LEGACY_L04_VALIDITY_TIMING_RESTRICTED.csv',low_memory=False).fillna('');metadata=metadata.loc[metadata.valid_participant].copy()
 # Select reported replicate by metadata, never by target value or clinical effect.
 metadata=metadata.sort_values(['key','RUN.NUMBER','SUBRUN','RUN.INJECTION.NUMBER','source_row']).drop_duplicates('key',keep='first')
 clinical=clinical_info(metadata.key.tolist()).set_index('key');frame=metadata.set_index('key').join(clinical.drop(columns=['RID','date']),how='left').reset_index()
 frame['MRM_PLATE']=frame['MS.PLATE.INDEX'].astype(str);frame['MRM_RUN']=frame['RUN.NUMBER'].astype(str)+'_'+frame['SUBRUN'].astype(str)
 frame['MRM_INSTRUMENT']=frame['INSTRUMENT'].astype(str);frame['PINK_FLAG']=pd.to_numeric(frame.PINK,errors='coerce').eq(1)
 old=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv').set_index('key');oldidx={k:i for i,k in enumerate(old.index)}
 frame['S5_index']=[oldidx.get(k,-1) for k in frame.key]
 for field in ['NULISA_PlateID','SOMASCAN_PlateId']:frame[field]=[str(old.loc[k,field]) if k in old.index else '' for k in frame.key]
 frame=frame.sort_values(['RID','date']).reset_index(drop=True);csvout(RESTRICTED/'LEGACY_L04_CLINICAL_DESIGN_RESTRICTED.csv',frame)
 trans=pd.read_csv(RESULTS/'LEGACY_CANDIDATE_TRANSITION_REGISTRY.csv');peptides=[]
 for (acc,pep),part in trans.groupby(['accession','peptide']):
  assert len(part)==2;peptides.append({'accession':acc,'gene':part.gene.iloc[0],'peptide':pep,'transition_fields':part.transition_field.tolist()})
 raw=np.load(S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz');meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
 rows=[];ds={};clinicalds={};masks={};first_by_rid=frame.drop_duplicates('RID',keep='first').index
 for p in peptides:
  finite=np.ones(len(frame),bool)
  for field in p['transition_fields']:finite&=frame[field+'_finite'].astype(bool).to_numpy()
  idx=frame.S5_index.to_numpy(int);j=list(raw['accessions']).index(p['accession']);sf=list(raw['soma_features']).index(meta.loc[p['accession'],'soma_features'])
  nsvalid=np.array([np.isfinite(raw['N'][i,j]) and np.isfinite(raw['S'][i,sf]) if i>=0 else False for i in idx])
  cc=np.isfinite(frame[['P_LOG','age_month_precision','sex_code2']]).all(axis=1).to_numpy()
  for cohort,mask in [('S5_SAME_DATE',cc&finite&nsvalid&frame.in_stage5_key.to_numpy()),('S5_NO_PINK',cc&finite&nsvalid&frame.in_stage5_key.to_numpy()&~frame.PINK_FLAG.to_numpy())]:
   d=make(frame,mask,tech=('NULISA_PlateID','SOMASCAN_PlateId','MRM_PLATE','MRM_RUN'),additional_factors=('ptau_batch','MRM_INSTRUMENT'))
   key=cohort+'__'+p['accession']+'__'+p['peptide'];reason=[]
   if d is None:reason.append('P_NOT_ESTIMABLE')
   else:
    if d['n']<80:reason.append('N_LT80')
    if d['df']<40:reason.append('DF_LT40')
    if min(len(one) for one in d['onehots'][:3])<6:reason.append('N_S_PLATE_CLUSTERS_LT6')
   rows.append({'design_id':key,'cohort':cohort,**p,'n':int(mask.sum()),'eligible_strict_relative_test':not reason,'reason':';'.join(reason),**(d['summary'] if d else {})})
   masks[key]=mask
   if not reason:ds[key]=d
  # Subject-disjoint clinical peptide–pTau question is separate from relative-pattern replication.
  for cohort,mask in [('NEW_RID_CLINICAL',cc&finite&~frame.is_stage5_RID.to_numpy()&frame.index.isin(first_by_rid)),('NEW_RID_CLINICAL_NO_PINK',cc&finite&~frame.is_stage5_RID.to_numpy()&frame.index.isin(first_by_rid)&~frame.PINK_FLAG.to_numpy())]:
   d=make(frame,mask,tech=('MRM_PLATE','MRM_RUN'),additional_factors=('ptau_batch','MRM_INSTRUMENT'))
   key=cohort+'__'+p['accession']+'__'+p['peptide'];eligible=d is not None and d['n']>=50 and d['df']>=20 and min(len(o) for o in d['onehots'])>=3
   rows.append({'design_id':key,'cohort':cohort,**p,'n':int(mask.sum()),'eligible_strict_relative_test':False,'eligible_single_assay_clinical':eligible,
    'reason':'Different hypothesis: subject-disjoint peptide association, not relative-pattern validation',**(d['summary'] if d else {})})
   masks[key]=mask
   if eligible:clinicalds[key]=d
 csvout(RESULTS/'LEGACY_DIRECT_DESIGN_SCREEN.csv',rows)
 with (RESTRICTED/'LEGACY_RELATIVE_DESIGNS_RESTRICTED.pkl').open('wb') as f:pickle.dump(ds,f,pickle.HIGHEST_PROTOCOL)
 with (RESTRICTED/'LEGACY_CLINICAL_DESIGNS_RESTRICTED.pkl').open('wb') as f:pickle.dump(clinicalds,f,pickle.HIGHEST_PROTOCOL)
 with (RESTRICTED/'LEGACY_CASE_MASKS_RESTRICTED.pkl').open('wb') as f:pickle.dump(masks,f,pickle.HIGHEST_PROTOCOL)
 save_json(RESULTS/'LEGACY_VERIFICATION_SAP_FROZEN.json',{'before_legacy_candidate_associations':True,'Stage5_and_Emory_results_already_seen':True,
  'direct_genes':['NPTX1','KLK6'],'PDGFRB':'not_covered; PDGF-BB is ligand not receptor','peptides':peptides,
  'primary_direct_source':'L04 exact released RID/date; independent Caprion MRM longitudinal experiment; not external cohort',
  'L01_L02':'same 2014 MRM release family; no date in L01, no verified per-row technical crosswalk; descriptive same-visit corroboration only; do not double count aliases',
  'relative_test':'NULISA/Soma/CaprionMRM; all three assays use age,sex,Elecsysbatch,N,S,MRMplate,MRMrun,instrument; SUM4 CR1; strict n80/df40 gate unchanged',
  'quantitative_rule':'released normalized log transitions used directly; equal mean of two finite transitions per peptide, no second log; zero/minus4 not automatically missing without a documented rule',
  'scale_discrepancy':'dictionary log2 vs primer natural-log noted; fixed log-base scaling cannot change Wald or standardized slopes',
  'single_assay_newRID':'earliest reported MRM date per RID, exclude all Stage5 RIDs; peptide~P+age+sex+MRMtechnical+Elecsysbatch; positive direction from original clinical coefficients; separate evidence tier',
  'single_assay_inference':'SUM(plate CR1,run CR1), one-sided Student t with min cluster G-1 df; no claim this tests cross-assay disagreement; calibrate six null conditions',
  'multiplicity':'five peptide hypotheses per cohort, BY; no choice of winning fragment/peptide; no-pink sensitivity predefined',
  'source_notes':'original transitions are fragment ions of a single peptide, not independent protein indicators; source-normalized release already processed; exact date does not prove same aliquot',
  'code_sha256':{p.name:sha(p) for p in [Path(__file__),ROOT/'code/designs.py']}})
 log(f'Legacy verification frozen before legacy effects: {len(ds)} strict relative designs; {len(clinicalds)} subject-disjoint single-assay designs')
 print(pd.DataFrame(rows)[['gene','peptide','cohort','n','residual_df','eligible_strict_relative_test','reason']].to_string(index=False))

if __name__=='__main__':run()
