from config import *
from targeted_inventory import dt,rid,numeric
from designs import make
import re,pickle
from collections import defaultdict

def keys(frame):
 dates=pd.to_datetime(frame.EXAMDATE,errors='coerce').dt.strftime('%Y-%m-%d');rr=pd.to_numeric(frame.RID,errors='coerce')
 return [str(int(r))+'|'+d if np.isfinite(r) and isinstance(d,str) else '' for r,d in zip(rr,dates)]

def run():
 inventory=pd.read_csv(RESULTS/'TARGETED_PROTEIN_COVERAGE_ALL70.csv');eligible=inventory.loc[inventory.s5_design_qualified&inventory.targeted_directly_covered]
 assert set(eligible.gene)=={'CHI3L1','MDH1','NPTX2','SOD1'}
 meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
 feats=pd.read_csv(RESULTS/'TARGETED_FEATURE_VALIDITY_COVERAGE.csv');feats=feats.loc[feats.accession.isin(eligible.accession)]
 tmeta=pd.read_csv(RESTRICTED/'TARGETED_VALIDITY_AND_TIMING_RESTRICTED.csv',low_memory=False).fillna('')
 tmeta=tmeta.loc[tmeta.valid_participant&tmeta.accession.isin(eligible.accession)]
 assert tmeta.PTID_RID_match.all() and not tmeta.plate_conflict.any()
 tmeta['priority']=tmeta.source.map({'V01':0,'V02':1});tmeta=tmeta.sort_values(['key','feature_key','priority','run_date','plate','source_row'])
 selected=tmeta.drop_duplicates(['key','feature_key'],keep='first').copy()
 csvout(RESTRICTED/'TARGETED_RECORD_SELECTION_RESTRICTED.csv',selected)
 alltargets=sorted(set(selected.key))
 # N/S values are used only to determine valid cells / cache, not for targeted-P inference here.
 genes=[meta.loc[a,'nulisa_target'] for a in eligible.accession];sfeatures=[meta.loc[a,'soma_features'] for a in eligible.accession]
 assert all('|' not in f for f in sfeatures)
 n=pd.read_csv(INPUTS['I01']['path'],usecols=['RID','EXAMDATE','PTID','SampleMatrixType','SampleType','SampleQC','Target','NPQ','TargetLOD_NPQ','PlateID'],low_memory=False)
 n=n.loc[n.SampleMatrixType.eq('CSF')&n.SampleType.eq('Sample')&n.SampleQC.eq('passed')&n.Target.isin(genes)].copy();n['key']=keys(n)
 n=n.loc[n.key.isin(alltargets)].drop_duplicates();n.NPQ=pd.to_numeric(n.NPQ,errors='coerce').replace([np.inf,-np.inf],np.nan)
 n.TargetLOD_NPQ=pd.to_numeric(n.TargetLOD_NPQ,errors='coerce');n['det']=n.NPQ.ge(n.TargetLOD_NPQ)&np.isfinite(n.NPQ)&np.isfinite(n.TargetLOD_NPQ)
 # Cross-plate repeated NULISA records outside Stage5 are not assigned an arbitrary plate.
 repeated=set(n.loc[n.duplicated(['key','Target'],keep=False),'key']);n=n.loc[~n.key.isin(repeated)]
 nm=n[['key','RID','EXAMDATE','PlateID']].drop_duplicates().set_index('key');assert nm.index.is_unique
 nv=n.pivot(index='key',columns='Target',values='NPQ');nd=n.pivot(index='key',columns='Target',values='det')
 s=pd.read_csv(INPUTS['I04']['path'],usecols=['RID','EXAMDATE','PlateId']+sfeatures,low_memory=False);s['key']=keys(s);s=s.loc[s.key.isin(alltargets)].drop_duplicates()
 assert not s.key.duplicated().any();s=s.set_index('key');common=sorted(set(nm.index)&set(s.index))
 records=[];demodata=defaultdict(list);ptau=defaultdict(list)
 for _,row in enumerate(csv.DictReader(Path(INPUTS['I13']['path']).open(encoding='utf-8-sig')),2):
  r=rid(row.get('RID'))
  if not r or row.get('HAS_QC_ERROR')=='1':continue
  m=re.fullmatch(r'(\d{1,2})/(\d{4})',row.get('PTDOB','').strip());full=dt(row.get('PTDOB'));yd=dt(row.get('PTDOBYY'));yn=numeric(row.get('PTDOBYY'))
  y=full.year if full else int(m[2]) if m else int(yn) if np.isfinite(yn) else yd.year if yd else None
  mo=full.month if full else int(m[1]) if m else None;yy=int(yn) if np.isfinite(yn) else yd.year if yd else None
  demodata[r].append((y,mo,row.get('PTGENDER'),bool(y and yy and y!=yy)))
 for line,row in enumerate(csv.DictReader(Path(INPUTS['E01']['path']).open(encoding='utf-8-sig')),2):
  r=rid(row.get('RID'));day=dt(row.get('EXAMDATE'));p=numeric(row.get('PTAU'));batch=row.get('BATCH','');m=re.search(r'(\d+)$',batch)
  if r and day and 8<=p<=120:ptau[(r,day.isoformat())].append((int(m[1]) if m else 0,line,p,batch))
 scaling=json.loads((S5/'00_restricted/P_LOG_SCALING_RESTRICTED.json').read_text(encoding='utf-8'))
 old=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv').set_index('key');oldrid=set(old.RID)
 for key in common:
  r=int(key.split('|')[0]);day=dt(key.split('|')[1]);demo=demodata[r];ys={d[0] for d in demo if d[0]};mos={d[1] for d in demo if d[1]};sexs={d[2] for d in demo if d[2] in ['1','2']}
  ageok=len(ys)==len(mos)==1 and not any(d[3] for d in demo);sexok=len(sexs)==1;pp=ptau.get((r,day.isoformat()),[]);chosen=min(pp) if pp else None
  age=(day.year-next(iter(ys))+(day.month-next(iter(mos)))/12) if ageok else np.nan
  if ageok and not 0<=age<=120:age=np.nan
  pz=(np.log2(chosen[2])-scaling['mean'])/scaling['SD'] if chosen else np.nan
  records.append({'key':key,'RID':r,'EXAMDATE':day.isoformat(),'age_month_precision':age,'sex_code2':float(next(iter(sexs))=='2') if sexok else np.nan,
   'P_LOG':pz,'ptau_batch':chosen[3] if chosen else '', 'NULISA_PlateID':str(nm.loc[key,'PlateID']),'SOMASCAN_PlateId':str(s.loc[key,'PlateId']),
   'is_stage5_RID':r in oldrid,'is_stage5_exact_key':key in old.index,'TMT_plex':str(old.loc[key,'TMT_plex']) if key in old.index else ''})
 clinical=pd.DataFrame(records).sort_values(['RID','EXAMDATE']);clinical=clinical.drop_duplicates('RID',keep='first').reset_index(drop=True)
 csvout(RESTRICTED/'TARGETED_PAIR_CLINICAL_DESIGN_RESTRICTED.csv',clinical)
 nsN=np.empty((len(clinical),len(eligible)));nsS=np.empty_like(nsN);nsD=np.zeros_like(nsN,dtype=bool)
 for j,a in enumerate(eligible.accession):
  gene=meta.loc[a,'nulisa_target'];sf=meta.loc[a,'soma_features'];nsN[:,j]=nv.reindex(clinical.key)[gene].to_numpy(float)
  raw=pd.to_numeric(s.reindex(clinical.key)[sf],errors='coerce').to_numpy(float);nsS[:,j]=np.log2(np.where(np.isfinite(raw)&(raw>0),raw,np.nan));nsD[:,j]=nd.reindex(clinical.key)[gene].fillna(False).to_numpy(bool)
 np.savez_compressed(RESTRICTED/'SUPPORT_N_S_CACHE_RESTRICTED.npz',N=nsN,S=nsS,N_detect=nsD,accessions=eligible.accession.to_numpy(str),keys=clinical.key.to_numpy(str))
 rows=[];designs={};cohortsets={};chosenmetadata=selected.set_index(['key','feature_key'])
 for fk,part in feats.groupby('feature_key'):
  a=part.accession.iloc[0];gene=meta.loc[a,'nulisa_target'];j=list(eligible.accession).index(a);frame=clinical.copy();plate=[];valid=[]
  for key in frame.key:
   if (key,fk) in chosenmetadata.index:
    row=chosenmetadata.loc[(key,fk)];plate.append(str(row.source)+'_'+str(row.plate));valid.append(bool(row.positive_finite) and bool(row.plate))
   else:plate.append('');valid.append(False)
  frame['TARGETED_Plate']=plate
  base=np.isfinite(frame[['P_LOG','age_month_precision','sex_code2']]).all(axis=1).to_numpy()&np.isfinite(nsN[:,j])&np.isfinite(nsS[:,j])&np.array(valid)
  for cohort,mask in [('ALL_NS',base),('S5_SAME_DATE',base&frame.is_stage5_exact_key.to_numpy()),('NEW_RID',base&~frame.is_stage5_RID.to_numpy())]:
   d=make(frame,mask,tech=('NULISA_PlateID','SOMASCAN_PlateId','TARGETED_Plate'));reasons=[]
   if d is None:reasons.append('P_NOT_ESTIMABLE')
   else:
    if d['n']<80:reasons.append('N_LT80')
    if d['df']<40:reasons.append('DF_LT40')
    if min(len(x) for x in d['onehots'])<6:reasons.append('CLUSTERS_LT6')
   ident=cohort+'__'+a+'__'+part.peptide.iloc[0];rows.append({'design_id':ident,'cohort':cohort,'accession':a,'gene':gene,'feature_key':fk,'peptide':part.peptide.iloc[0],
    'n':int(mask.sum()),'design_eligible':not reasons,'reasons':';'.join(reasons),**(d['summary'] if d else {})})
   if not reasons:designs[ident]=d
   cohortsets[ident]={'frame':frame,'mask':mask}
 csvout(RESULTS/'ORTHOGONAL_SUPPORT_DESIGN_SCREEN.csv',rows)
 with (RESTRICTED/'ORTHOGONAL_SUPPORT_DESIGNS_RESTRICTED.pkl').open('wb') as f:pickle.dump(designs,f,pickle.HIGHEST_PROTOCOL)
 with (RESTRICTED/'ORTHOGONAL_SUPPORT_CASE_FRAMES_RESTRICTED.pkl').open('wb') as f:pickle.dump(cohortsets,f,pickle.HIGHEST_PROTOCOL)
 unique=selected.loc[selected.valid_participant].groupby('source').agg(unique_RID=('RID','nunique'),unique_keys=('key','nunique')).reset_index()
 csvout(RESULTS/'TARGETED_UNIQUE_RECORDS_AFTER_METADATA_SELECTION.csv',unique)
 validation={'date':datetime.now().astimezone().isoformat(),'Stage5_discoveries_already_known':True,'before_targeted_P_association':True,
   'three_discovery_direct_coverage':0,'three_discovery_direct_validation':'NOT_EXECUTABLE_IN_CURRENT_TARGETED_SOURCES',
   'support_accessions':sorted(eligible.accession),'support_peptides':sorted(feats.feature_key.unique()),'support_scope':'four non-discovery proteins, five peptides; no substitution for missing discoveries',
   'record_rule':'V01 before V02; then earliest run_date/plate/source_row, independently of quantitative magnitude; no averaging across releases',
   'ratio':'reported light/heavy ratio; log2 of finite positive values; REPLICATE1/2 are in-house IDs not repeat measurements',
   'technical_factors':'NULISA plate, Soma plate, source-prefixed targeted plate; Elecsys batch nuisance; not identical to original TMT conditioning',
   'inferential_eligibility':'n>=80, residual_df>=40, >=6 clusters each; all actual eligible designs require independent calibration',
   'covariance':'sum of three one-way CR1; fixed V0, no reselection','calibration_reps':4000,'null_scenarios':['C01','C02','C03','C04','C05','C06'],
   'cohorts':['ALL_NS','S5_SAME_DATE','NEW_RID'],'multiple_testing':'five fixed peptide hypotheses per cohort; BY; unavailable tests set p=1; exploratory support, not primary discovery replication',
   'source_scale':'no cross-release z normalization; source-specific plate fixed effects absorb constant log offsets; per-release descriptive results supplied',
   'QC':'sample rows from documented source; QC cohorts kept separate; no pooling unlike QC types into a single CV; published QC statistics disclosed',
   'code_sha256':{p.name:sha(p) for p in [Path(__file__),ROOT/'code/designs.py',ROOT/'code/model.py',ROOT/'code/diagnostics.py']}}
 save_json(RESULTS/'TARGETED_VERIFICATION_PLAN_FROZEN.json',validation)
 log(f'Orthogonal support plan frozen before targeted associations: {len(clinical)} unique N/S-linked clinical candidates; {len(designs)} eligible designs')
 print(pd.DataFrame(rows)[['cohort','gene','peptide','n','design_eligible','reasons']].to_string(index=False))

if __name__=='__main__':run()
