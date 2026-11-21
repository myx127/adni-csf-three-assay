from config import *
import zipfile,gc

def make_keys(frame):
 rid=pd.to_numeric(frame.RID,errors='coerce');dates=pd.to_datetime(frame.EXAMDATE,errors='coerce').dt.strftime('%Y-%m-%d')
 return pd.Series([str(int(v))+'|'+d if np.isfinite(v) and isinstance(d,str) else '' for v,d in zip(rid,dates)],index=frame.index)

def load():
 sap=json.loads((RESULTS/'MAIN_ANALYSIS_SAP_FROZEN.json').read_text(encoding='utf-8'))
 assert sap['protocol_sha256']==sha(ROOT/'STAGE5_PROTOCOL_v1.0.md')
 assert sap['real_protein_Y_read_before_freeze'] is False
 gate=RESULTS/'FIRST_PROTEIN_VALUE_ACCESS.json'
 if not gate.exists():save_json(gate,{'SAP_sha256_before_first_protein_access':sha(RESULTS/'MAIN_ANALYSIS_SAP_FROZEN.json'),
  'time':datetime.now().astimezone().isoformat(),'user_authorized_main_analysis':True})
 else:assert json.loads(gate.read_text(encoding='utf-8'))['SAP_sha256_before_first_protein_access']==sha(RESULTS/'MAIN_ANALYSIS_SAP_FROZEN.json')
 index=pd.read_csv(RESTRICTED/'CLINICAL_DESIGN_RESTRICTED.csv');keys=index.key.tolist();wanted=set(keys)
 meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv');accs=meta.accession.tolist();targetmap={}
 for m in meta.itertuples():
  for label in m.nulisa_raw_target_labels.split('|'):targetmap[label]=m.nulisa_target
 fields=['RID','EXAMDATE','SampleMatrixType','SampleType','SampleQC','Target','NPQ','TargetLOD_NPQ','PlateID']
 n=pd.read_csv(INPUTS['I01']['path'],usecols=fields,low_memory=False)
 n=n.loc[n.SampleMatrixType.eq('CSF')&n.SampleType.eq('Sample')&n.SampleQC.eq('passed')&n.Target.isin(targetmap)]
 n['key']=make_keys(n);n=n.loc[n.key.isin(wanted)].drop_duplicates();n['target']=n.Target.map(targetmap)
 assert not n.duplicated(['key','target']).any(),'Unexpected index technical duplicate; do not average silently'
 n['NPQ']=pd.to_numeric(n.NPQ,errors='coerce').replace([np.inf,-np.inf],np.nan)
 n['TargetLOD_NPQ']=pd.to_numeric(n.TargetLOD_NPQ,errors='coerce').replace([np.inf,-np.inf],np.nan)
 n['det']=np.isfinite(n.NPQ)&np.isfinite(n.TargetLOD_NPQ)&n.NPQ.ge(n.TargetLOD_NPQ)
 N=n.pivot(index='key',columns='target',values='NPQ').reindex(keys);D=n.pivot(index='key',columns='target',values='det').reindex(keys).fillna(False)
 N=N.reindex(columns=meta.nulisa_target).to_numpy(float);D=D.reindex(columns=meta.nulisa_target).to_numpy(bool)
 features=sorted({v for cell in meta.soma_features for v in cell.split('|')})
 s=pd.read_csv(INPUTS['I04']['path'],usecols=['RID','EXAMDATE','PlateId']+features,low_memory=False)
 s['key']=make_keys(s);s=s.loc[s.key.isin(wanted)].drop_duplicates();assert not s.key.duplicated().any()
 s=s.set_index('key').reindex(keys);RFU=s[features].apply(pd.to_numeric,errors='coerce').to_numpy(float)
 S=np.log2(np.where(np.isfinite(RFU)&(RFU>0),RFU,np.nan));del RFU
 tfeatures=meta.tmt_columns.tolist();assert all('|' not in f for f in tfeatures)
 with zipfile.ZipFile(INPUTS['I07']['path']) as z:
  member=next(p for p in z.namelist() if Path(p).name=='EMORY_CSF_TMT_MS.csv')
  with z.open(member) as f:t=pd.read_csv(f,usecols=['RID','EXAMDATE','ProteomicID']+tfeatures,low_memory=False)
 parsed=t.ProteomicID.str.extract(r'^adni_b(?P<plex>\d+)\.(?P<channel>\d+[NC]?)$');assert parsed.notna().all().all()
 t['plex']=parsed.plex.map(lambda x:'P'+str(int(x)).zfill(2));real=pd.to_numeric(t.RID,errors='coerce').gt(0)
 standards=t.loc[~real].set_index('plex');assert standards.index.is_unique and parsed.loc[~real,'channel'].eq('126').all()
 t['key']=make_keys(t);participants=t.loc[real&t.key.isin(wanted)].set_index('key');assert participants.index.is_unique
 participants=participants.reindex(keys)
 raw=participants[tfeatures].apply(pd.to_numeric,errors='coerce').to_numpy(float)
 gis=standards[tfeatures].apply(pd.to_numeric,errors='coerce').reindex(participants.plex).to_numpy(float)
 raw=np.where(np.isfinite(raw)&(raw>0),raw,np.nan);gis=np.where(np.isfinite(gis)&(gis>0),gis,np.nan)
 T_RAW=np.log2(raw);GIS=np.log2(gis);T_GIS=T_RAW-GIS
 np.savez_compressed(RESTRICTED/'PROTEIN_TRANSFORMED_RESTRICTED.npz',N=N,N_detect=D,S=S,T_RAW=T_RAW,T_GIS=T_GIS,GIS=GIS,
  keys=np.array(keys),accessions=np.array(accs),soma_features=np.array(features))
 masks=pd.read_parquet(S1/'00_restricted/TRIPLE_CANDIDATE_AVAILABILITY_RESTRICTED.parquet');checks=[]
 for j,m in enumerate(meta.itertuples()):
  part=masks.loc[masks.accession.eq(m.accession)].set_index('key').reindex(keys);cols=[features.index(f) for f in m.soma_features.split('|')]
  computed={'N_numeric':np.isfinite(N[:,j]),'N_detect_all_repeats':D[:,j],'S_all':np.isfinite(S[:,cols]).all(axis=1),
   'T_raw':np.isfinite(T_RAW[:,j]),'GIS_positive':np.isfinite(GIS[:,j]),'T_GIS_usable':np.isfinite(T_GIS[:,j])}
  computed['complete_GIS_S_all']=computed['N_numeric']&computed['S_all']&computed['T_GIS_usable']
  for name,actual in computed.items():checks.append({'accession':m.accession,'check':name,'passed':bool(np.array_equal(actual,part[name].to_numpy())),
   'n_reconstructed':int(actual.sum()),'n_stage1':int(part[name].sum())})
 csvout(RESULTS/'PROTEIN_RECONSTRUCTION_QA.csv',checks)
 assert all(r['passed'] for r in checks),'Reconstructed protein availability differs from sealed Stage1'
 log(f'Post-SAP protein reconstruction: {len(checks)}/{len(checks)} mask checks passed; transformed values private')
 del n,s,t,participants,raw,gis;gc.collect()

if __name__=='__main__':load()
