from config import *
import zipfile,pickle

EXTRA_ACC={'P10636','P05067','P02649'}
def tokens(v):return [s.split('-')[0] for s in re.split(r'[|;,\s]+',str(v)) if re.fullmatch(r'[A-Z][A-Z0-9]{5,9}(?:-\d+)?',s)]
def keys(d):
 rr=pd.to_numeric(d.RID,errors='coerce');dt=pd.to_datetime(d.EXAMDATE,errors='coerce').dt.strftime('%Y-%m-%d')
 return pd.Series([str(int(a))+'|'+b if np.isfinite(a) and isinstance(b,str) else '' for a,b in zip(rr,dt)],index=d.index)
def run():
 data=pickle.load((PRIVATE/'ASSEMBLED_CLINICAL_RESTRICTED.pkl').open('rb'));c=data['clinical'];wanted=c.key.tolist()
 family=pd.read_csv(S11/'results/FROZEN_FAMILY34_AND_FOCUS14.csv');excluded=set(family.accession)|EXTRA_ACC;exgenes=set(family.gene)|{'MAPT','APP','APOE','APOE2','APOE4'}
 cache=np.load(S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz');core=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
 rows=[];scores={};identities=[];duplicates=[]
 def summary(platform,arr,features,accessions,genes,allowed,det=None):
  with np.errstate(all='ignore'):
   med=np.nanmedian(arr,axis=0);quant=np.nanquantile(arr,[.25,.75],axis=0);sc=(quant[1]-quant[0])/1.349
  finite=np.isfinite(arr).mean(axis=0);detect=np.isfinite(arr).mean(axis=0) if det is None else det.mean(axis=0)
  valid=np.asarray(allowed)&(finite>=.9)&np.isfinite(sc)&(sc>1e-8)
  if platform=='NULISA':valid &= detect>=.8
  z=(arr[:,valid]-med[valid])/sc[valid];nref=int(valid.sum())
  score=np.nanmedian(z,axis=1) if nref>=10 else np.full(len(c),np.nan)
  nobs=np.isfinite(z).sum(axis=1);score[nobs<np.ceil(.5*nref)]=np.nan
  scores[platform]=score
  for j,f in enumerate(features):
   rows.append({'platform':platform,'feature':f,'accessions':str(accessions[j]),'genes':str(genes[j]),'allowed_identity':bool(allowed[j]),'finite_fraction':finite[j],'detected_fraction':detect[j],'robust_center':med[j],'robust_scale':sc[j],'included_reference':bool(valid[j]),'reason':'ELIGIBLE' if valid[j] else 'EXCLUDED_IDENTITY' if not allowed[j] else 'FINITE_LT90' if finite[j]<.9 else 'DETECTION_LT80' if platform=='NULISA' and detect[j]<.8 else 'ROBUST_SCALE_INVALID'})
  return nref
 # NULISA: use the exact index plate. Never include plasma.
 n=pd.read_csv(INPUTS['I01']['path'],low_memory=False).drop_duplicates()
 n=n.loc[n.SampleMatrixType.eq('CSF')&n.SampleType.eq('Sample')&n.SampleQC.eq('passed')].copy();n['key']=keys(n);n=n.loc[n.key.isin(wanted)]
 plates=c.set_index('key').NULISA_PlateID.astype(str);n=n.loc[n.PlateID.astype(str).eq(n.key.map(plates))].copy()
 dup=n.duplicated(['key','Target'],keep=False)
 if dup.any():
  duplicates=n.loc[dup,['key','Target','PlateID']].to_dict('records');n=n.loc[~dup]
 n.NPQ=pd.to_numeric(n.NPQ,errors='coerce').replace([np.inf,-np.inf],np.nan);lod=pd.to_numeric(n.TargetLOD_NPQ,errors='coerce');n['det']=np.isfinite(n.NPQ)&np.isfinite(lod)&n.NPQ.ge(lod)
 meta=n.groupby('Target',sort=True).agg({'UniProtID':'first','Target':'first'});nf=meta.index.tolist()
 N=n.pivot(index='key',columns='Target',values='NPQ').reindex(index=wanted,columns=nf).to_numpy(float);ND=n.pivot(index='key',columns='Target',values='det').reindex(index=wanted,columns=nf).fillna(False).to_numpy(bool)
 na=meta.UniProtID.astype(str).tolist();allowed=[not(set(tokens(a))&excluded) and f not in exgenes for f,a in zip(nf,na)]
 counts={'NULISA':summary('NULISA',N,nf,na,nf,allowed,ND)}
 for acc,gene in GENES.items():
  target=core.loc[acc,'nulisa_target'];j=list(cache['accessions']).index(acc)
  identities.append({'platform':'NULISA','gene':gene,'passed':np.allclose(N[:,nf.index(target)],cache['N'][:,j],equal_nan=True,rtol=0,atol=1e-12),'n':len(c)})
 # Soma: one deterministic first assay per single-accession human protein, before coverage filtering.
 s=pd.read_csv(INPUTS['I04']['path'],low_memory=False);sf=sorted(v for v in s if re.match(r'^X\d',v));s['key']=keys(s);s=s.loc[s.key.isin(wanted)].drop_duplicates();assert s.key.is_unique;s=s.set_index('key').reindex(wanted)
 r=s[sf].apply(pd.to_numeric,errors='coerce').to_numpy(float);S=np.log2(np.where(np.isfinite(r)&(r>0),r,np.nan))
 ann=pd.read_csv(INPUTS['I05']['path']).set_index('Analytes').reindex(sf);allowed=[];seen=set()
 for f,row in ann.iterrows():
  aa=tokens(row.UniProt);gg=set(re.split(r'[|;,\s]+',str(row.EntrezGeneSymbol)));ok=str(row.Organism).lower() in ['human','homo sapiens'] and len(aa)==1 and not(set(aa)&excluded) and not(gg&exgenes)
  if ok:
   if aa[0] in seen:ok=False
   else:seen.add(aa[0])
  allowed.append(ok)
 counts['SOMA']=summary('SOMA',S,sf,ann.UniProt.astype(str).tolist(),ann.EntrezGeneSymbol.astype(str).tolist(),allowed)
 for acc,gene in GENES.items():
  feature=core.loc[acc,'soma_features'];j=list(cache['soma_features']).index(feature)
  identities.append({'platform':'SOMA','gene':gene,'passed':np.allclose(S[:,sf.index(feature)],cache['S'][:,j],equal_nan=True,rtol=0,atol=1e-12),'n':len(c)})
 # TMT: retain the original sample/GIS definition across every ordinary protein.
 with zipfile.ZipFile(INPUTS['I07']['path']) as z:
  name=next(n for n in z.namelist() if Path(n).name=='EMORY_CSF_TMT_MS.csv');t=pd.read_csv(z.open(name),low_memory=False)
 tf=list(t.columns[7:]);parsed=t.ProteomicID.str.extract(r'^adni_b(?P<plex>\d+)\.(?P<channel>\d+[NC]?)$');assert parsed.notna().all().all()
 t['plex']=parsed.plex.map(lambda v:'P'+str(int(v)).zfill(2));real=pd.to_numeric(t.RID,errors='coerce').gt(0);std=t.loc[~real].set_index('plex');assert std.index.is_unique;t['key']=keys(t);people=t.loc[real&t.key.isin(wanted)].set_index('key');assert people.index.is_unique;people=people.reindex(wanted)
 r=people[tf].apply(pd.to_numeric,errors='coerce').to_numpy(float);g=std[tf].apply(pd.to_numeric,errors='coerce').reindex(people.plex).to_numpy(float)
 T=np.log2(np.where(np.isfinite(r)&(r>0),r,np.nan))-np.log2(np.where(np.isfinite(g)&(g>0),g,np.nan))
 ta=[v.rsplit('_',1)[-1] for v in tf];tg=[v.rsplit('_',1)[0] for v in tf];allowed=[bool(re.fullmatch(r'[A-Z][A-Z0-9]{5,9}(?:-\d+)?',a)) and not(set(tokens(a))&excluded) and gg not in exgenes for a,gg in zip(ta,tg)]
 counts['TMT']=summary('TMT',T,tf,ta,tg,allowed)
 for acc,gene in GENES.items():
  feature=core.loc[acc,'tmt_columns'];j=list(cache['accessions']).index(acc)
  identities.append({'platform':'TMT','gene':gene,'passed':np.allclose(T[:,tf.index(feature)],cache['T_GIS'][:,j],equal_nan=True,rtol=0,atol=1e-12),'n':len(c)})
 out=c[['key','RID','EXAMDATE']].copy()
 for platform,v in scores.items():out['BG_'+platform]=v
 csvout(PRIVATE/'BACKGROUND_PROXIES_RESTRICTED.csv',out);csvout(RESULTS/'BACKGROUND_FEATURE_AUDIT.csv',rows);csvout(RESULTS/'BACKGROUND_CORE_IDENTITY_CHECKS.csv',identities)
 if duplicates:csvout(PRIVATE/'BACKGROUND_NONIDENTICAL_DUPLICATE_GROUPS_RESTRICTED.csv',duplicates)
 assert all(r['passed'] for r in identities)
 for r in rows:
  if r['included_reference']:assert not(set(tokens(r['accessions']))&excluded) and not(set(re.split(r'[|;,\s]+',r['genes']))&exgenes)
 save_json(RESULTS/'BACKGROUND_CONSTRUCTION_SUMMARY.json',{'time':now(),'n_index':len(c),'reference_counts':counts,'finite_people':{k:int(np.isfinite(v).sum()) for k,v in scores.items()},'excluded_accessions':sorted(excluded),'excluded_genes':sorted(exgenes),'nonidentical_duplicate_rows_removed':len(duplicates),'new_clinical_target_coefficients_fitted':False,'proxies_are_not_measured_total_CSF_protein':True})
 log('Background proxies constructed without outcome selection: '+str(counts))
if __name__=='__main__':assert_unsealed();run()

