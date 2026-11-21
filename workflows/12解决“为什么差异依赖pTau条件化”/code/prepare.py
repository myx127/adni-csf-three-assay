from config import *
from engine import basis,nuisance,rcs,D,stats
import pickle
def construct(acc,gene,model,ids,c,X,scales,data,knots):
 f=c.iloc[ids].copy();records=[data['trajectories'][(int(i),'ADNI_MEM',5)] for i in ids];w=np.array([r['w'] for r in records]);fac=np.sqrt(w)
 Z,names,pfun,pnames=nuisance(f,model,knots);ZZ=Z*fac[:,None];xx=X[ids]*fac[:,None];y=np.array([r['slope'] for r in records])*fac
 A=[];Ks=[[] for _ in range(3)];Qall=[];ranks=[];dens=[];lookups=[];groups=[]
 scale=scales[acc]
 for col in ['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex']:
  levels=sorted(c.iloc[scale['ids']][col].astype(str).unique());lookup={v:i for i,v in enumerate(levels)};lookups.append(lookup);groups.append(np.array([lookup[str(v)] for v in f[col]]))
 Qz=basis(ZZ)
 for j in range(3):
  rr=xx[:,j]-Qz@(Qz.T@xx[:,j]);den=rr@rr;a=rr/den;Q=np.c_[Qz,rr/np.sqrt(den)];A.append(a);Qall.append(Q);ranks.append(Q.shape[1]);dens.append(den)
  for axis,col in enumerate(['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex']):
   G=len(lookups[axis]);B=np.zeros((G,len(ids)));B[groups[axis],np.arange(len(ids))]=a
   L=B-(B@Q)@Q.T;L-=L.mean(axis=0,keepdims=True);L*=np.sqrt(G/(G-1));Ks[axis].append(L)
 K=np.concatenate([np.stack(k,axis=1) for k in Ks],axis=0);A=np.array(A);sigma=scale['sigma']
 n=len(ids);df=n-max(ranks);Ah=(D*sigma[None,:])@A;cond=np.linalg.cond(Ah@Ah.T);reasons=[]
 if n<80:reasons.append('N_LT80')
 if df<40:reasons.append('INFO_DF_LT40')
 if not np.isfinite(cond) or cond>1e12:reasons.append('CONTRAST_ROW_SPACE_WEAK')
 for g in groups:
  if len(np.unique(g))<6:reasons.append('CASE_CLUSTERS_LT6')
 noise=y-Qz@(Qz.T@y);ns=np.sqrt(noise@noise/max(1,n-Qz.shape[1]))
 t={'test_id':gene+'__'+model,'gene':gene,'accession':acc,'model':model,'kind':'MODEL_DIFF','family':'SECONDARY33','ids':ids,'n_people':n,'A':A,'K':K,'sigma':sigma,'Uv':scale['Uv'],'L':D.copy(),'Ah':Ah,'y':y,'ZZ_list':[ZZ]*3,'targetx_list':[xx[:,j] for j in range(3)],'Qall':Qall,'Qcommon':Qz,'df':df,'max_hat':max(np.sum(Q*Q,axis=1).max() for Q in Qall),'cond':cond,'eligible':not reasons,'gate_reason':';'.join(reasons),'w':w,'X':X[ids],'pfun':pfun,'pnames':pnames,'Znames':names,'den':np.array(dens),'groups':groups,'full_G':[len(v) for v in lookups],'noise_scale':ns,'P':f.P_LOG.to_numpy(),'age':f.age_month_precision.to_numpy(),'visit_sd':np.ones(n),'ar_sd':np.sqrt([r['var_AR_pseudo'] for r in records]),'slope_load':fac/max(1.,np.median(fac)),'prototypes':np.column_stack([xx.mean(axis=1),xx,f.P_LOG.to_numpy()*fac,rcs(f.P_LOG.to_numpy(),knots)[:,-1]*fac])}
 return t
def comparison(acc,gene,name,ref,new):
 t={k:ref[k] for k in ['ids','n_people','y','Qcommon','w','X','groups','full_G','noise_scale','P','age','visit_sd','ar_sd','slope_load','prototypes']}
 t.update({'test_id':gene+'__'+name,'accession':acc,'gene':gene,'model':name,'kind':'CONDITIONING_CHANGE','family':'PRIMARY9' if name in PRIMARY_COMPARISONS else 'SECONDARY33','reference_model':ref['model'],'new_model':new['model'],'A':np.vstack([ref['A'],new['A']]),'K':np.concatenate([ref['K'],new['K']],axis=1),'sigma':np.tile(ref['sigma'],2),'Uv':np.tile(ref['Uv'],(1,2)),'L':np.c_[-D,D],'ZZ_list':ref['ZZ_list']+new['ZZ_list'],'targetx_list':ref['targetx_list']+new['targetx_list'],'Qall':ref['Qall']+new['Qall'],'df':min(ref['df'],new['df']),'max_hat':max(ref['max_hat'],new['max_hat'])})
 t['Ah']=(t['L']*t['sigma'][None,:])@t['A'];t['cond']=np.linalg.cond(t['Ah']@t['Ah'].T);t['eligible']=bool(ref['eligible'] and new['eligible'] and np.isfinite(t['cond']) and t['cond']<=1e12);t['gate_reason']='' if t['eligible'] else 'PARENT_OR_CONTRAST_GATE'
 return t
def run():
 data=pickle.load((PRIVATE/'ASSEMBLED_CLINICAL_RESTRICTED.pkl').open('rb'));c=data['clinical'].copy();bg=pd.read_csv(PRIVATE/'BACKGROUND_PROXIES_RESTRICTED.csv').set_index('key')
 for col in ['BG_SOMA','BG_TMT','BG_NULISA']:c[col]=c.key.map(bg[col])
 scales=pickle.load((PRIVATE/'FROZEN_SCALES_RESTRICTED.pkl').open('rb'));old=pickle.load((PRIVATE/'STAGE11_REFERENCE_TEMPLATES_RESTRICTED.pkl').open('rb'));raw=np.load(S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz');mapping=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
 knots=np.quantile(c.loc[c.base_primary.astype(bool),'P_LOG'],[.05,.35,.65,.95]);assert np.diff(knots).min()>0
 save_json(RESULTS/'RCS_KNOTS_FROZEN.json',{'time':now(),'P_LOG_knots':knots,'quantiles':[.05,.35,.65,.95],'n_reference':int(c.base_primary.astype(bool).sum()),'new_outcomes_used_to_choose_knots':False})
 models={};tests={};audit=[];repro=[];cohorts=[]
 for acc,gene in GENES.items():
  j=list(raw['accessions']).index(acc);sj=list(raw['soma_features']).index(mapping.loc[acc,'soma_features']);X=np.c_[raw['N'][:,j],raw['S'][:,sj],raw['T_GIS'][:,j]]
  original_ids=old['MEM5Y__'+acc]['ids']
  for m,oldbranch in [('NO_P','MEM_NO_PTAU'),('LOG_P','MEM5Y')]:
   t=construct(acc,gene,m,original_ids,c,X,scales,data,knots);v=stats(t['y'],t);z=old[oldbranch+'__'+acc];oldb=z['A']@z['y'];olds=oldb*z['sigma']
   repro.append({'gene':gene,'model':m,'n_people':len(original_ids),'max_abs_gamma_error':float(abs(v['gamma'][0]-olds).max()),'passed':bool(np.allclose(v['gamma'][0],olds,rtol=1e-8,atol=1e-10))})
  available=[col for col in ['BG_SOMA','BG_TMT','BG_NULISA'] if c[col].notna().sum()>0]
  ids=np.array([i for i in original_ids if np.isfinite(c.loc[i,available].to_numpy(float)).all()],int)
  cohorts.append({'gene':gene,'original_stage11_n':len(original_ids),'common_stage12_n':len(ids),'n_lost_for_background':len(original_ids)-len(ids)})
  for m in MODELS:
   required=['BG_SOMA','BG_TMT','BG_NULISA'] if m.startswith('BG3') else ['BG_SOMA','BG_TMT'] if m.startswith('BG2') else []
   if any(col not in available for col in required):continue
   t=construct(acc,gene,m,ids,c,X,scales,data,knots);models[t['test_id']]=t;tests[t['test_id']]=t
  for name,(a,b) in COMPARISONS.items():
   if gene+'__'+a not in models or gene+'__'+b not in models:continue
   t=comparison(acc,gene,name,models[gene+'__'+a],models[gene+'__'+b]);tests[t['test_id']]=t
  for name in MODELS+list(COMPARISONS):
   key=gene+'__'+name;t=tests.get(key)
   audit.append({'test_id':key,'gene':gene,'model':name,'kind':'MODEL_DIFF' if name in MODELS else 'CONDITIONING_CHANGE','family':'PRIMARY9' if name in PRIMARY_COMPARISONS else 'SECONDARY33','n_people':len(ids),'information_df':t['df'] if t else np.nan,'rank_max':len(ids)-t['df'] if t else np.nan,'max_hat':t['max_hat'] if t else np.nan,'contrast_Gram_condition':t['cond'] if t else np.nan,'eligible':t['eligible'] if t else False,'gate_reason':t['gate_reason'] if t else 'BACKGROUND_NOT_AVAILABLE'})
 assert len(audit)==42 and sum(r['family']=='PRIMARY9' for r in audit)==9 and all(r['passed'] for r in repro)
 with (PRIVATE/'MODELS_RESTRICTED.pkl').open('wb') as f:pickle.dump(models,f,pickle.HIGHEST_PROTOCOL)
 with (PRIVATE/'TEST_TEMPLATES_RESTRICTED.pkl').open('wb') as f:pickle.dump(tests,f,pickle.HIGHEST_PROTOCOL)
 csvout(PRIVATE/'COMMON_CLINICAL_RESTRICTED.csv',c);csvout(RESULTS/'ACTUAL_DESIGN_AUDIT.csv',audit);csvout(RESULTS/'COMMON_CASE_COHORT_AUDIT.csv',cohorts);csvout(RESULTS/'STAGE11_BASELINE_REPRODUCTION.csv',repro)
 save_json(RESULTS/'TEMPLATE_FREEZE_RECORD.json',{'time':now(),'protocol_sha256':sha(ROOT/'STAGE12_PROTOCOL_v1.0.md'),'templates_sha256':sha(PRIVATE/'TEST_TEMPLATES_RESTRICTED.pkl'),'n_tests':len(audit),'geometry_passed':sum(r['eligible'] for r in audit),'new_extension_target_results_fitted':False,'previously_known_stage11_models_refitted_only':True})
 log(f'Prepared {len(audit)} tests, {sum(r["eligible"] for r in audit)} geometrically eligible; old known coefficients reproduced. No new extension target fitted.')
 print(pd.DataFrame(cohorts).to_string(index=False));print(pd.DataFrame(audit)[['gene','model','information_df','eligible','gate_reason']].to_string(index=False))
if __name__=='__main__':assert_unsealed();run()

