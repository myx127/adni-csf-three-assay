from config import *
from engine import basis,nuisance
import pickle

def run():
 data=pickle.load((PRIVATE/'ASSEMBLED_CLINICAL_RESTRICTED.pkl').open('rb'));clinical=data['clinical'];traj=data['trajectories'];pet=data['PET']
 raw=np.load(S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz');meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
 templates={};audit=[];scaleaudit=[];registry=[];scales={};private_y={};schedule_audit=[]
 for acc,gene in CANDIDATES.items():
  j=list(raw['accessions']).index(acc);sj=list(raw['soma_features']).index(meta.loc[acc,'soma_features'])
  X=np.c_[raw['N'][:,j],raw['S'][:,sj],raw['T_GIS'][:,j]]
  mask=clinical.base_primary.astype(bool).to_numpy()&np.isfinite(X).all(axis=1)&np.isfinite(clinical[['education','apoe4_carrier']]).all(axis=1).to_numpy()
  scale_ids=np.flatnonzero(mask);base=clinical.iloc[scale_ids];XX=X[scale_ids]
  Zs,_=nuisance(base,include_ptau=False);Qs=basis(Zs);E=XX-Qs@(Qs.T@XX);df=len(base)-Qs.shape[1];v=(E*E).sum(axis=0)/df;IFv=(E*E-(E*E).mean(axis=0))/df
  source_labels=[];group_lookup={};gv=[]
  for c in ['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex']:
   label,levels=pd.factorize(base[c].astype(str),sort=True);source_labels.append(label);group_lookup[c]={value:i for i,value in enumerate(levels)}
   summed=np.zeros((len(levels),3));np.add.at(summed,label,IFv);summed-=summed.mean(axis=0);gv.append(summed*np.sqrt(len(levels)/(len(levels)-1)))
  Uv=np.concatenate(gv);sigma=np.sqrt(v)
  scales[acc]={'ids':scale_ids,'X':XX,'v':v,'sigma':sigma,'IFv':IFv,'Z':Zs,'Q':Qs,'Uv':Uv,'group_labels':source_labels}
  for k,assay in enumerate(['N','S','T']):scaleaudit.append({'accession':acc,'gene':gene,'assay':assay,'n_scale_cohort':len(base),'scale_nuisance_rank':Qs.shape[1],'scale_df':df,'residual_SD_fixed_unit':sigma[k],'includes_ptau_in_scale_model':False,'new_outcomes_used_for_scaling':False})
  branches=FORMAL_BRANCHES+SENSITIVITIES
  for branch in branches:
   if branch=='MEM_APOE_INTERACTION' and gene!='PDGFRB':continue
   formal=branch in FORMAL_BRANCHES
   marker='ADNI_EF' if branch=='EF5Y' else 'ADNI_MEM';horizon=3 if branch=='MEM3Y' else 5
   is_long=branch not in ['MEM_BASE','AMY_CL'];is_pet=branch in ['AMY_CL','MEM_AMY_INTERACTION','MEM_PET_MAIN']
   cases=[];records=[]
   for i in scale_ids:
    tr=traj.get((i,marker,horizon))
    eligible=(tr is not None and tr.get('eligible',False)) if is_long else (tr is not None if branch=='MEM_BASE' else i in pet.index)
    if not eligible:continue
    if is_pet and i not in pet.index:continue
    if branch=='MEM_HIGHINFO' and not tr['highinfo']:continue
    if branch=='MEM_BATCH9' and str(clinical.loc[i,'ptau_batch'])!='UPENNBIOMK9':continue
    if branch=='MEM_BASE_ADJUSTED' and not np.isfinite(tr['baseline']):continue
    cases.append(i);records.append(tr)
   ids=np.array(cases,int);f=clinical.iloc[ids].copy();n=len(ids)
   row={'design_id':branch+'__'+acc,'branch':branch,'accession':acc,'gene':gene,'formal':formal,'n_people':n,'n_scale_cohort':len(base),'nuisance_exposure_scaling':'Full original complete-case cohort; common demographic+technical, no pTau/outcome'}
   reasons=[]
   if n<80:reasons.append('N_LT80')
   if n<5:
    row.update({'eligible':False,'gate_reason':';'.join(reasons)});audit.append(row);continue
   modifier=None
   if branch=='MEM_AMY_INTERACTION':modifier=pet.loc[ids,'CL'].to_numpy()/50
   elif branch=='MEM_APOE_INTERACTION':modifier=f.apoe4_carrier.to_numpy(float)
   if is_pet:f['PET_tracer']=pet.loc[ids,'tracer'].to_numpy()
   extra_base=[tr['baseline'] for tr in records] if branch=='MEM_BASE_ADJUSTED' else None
   Z,znames=nuisance(f,include_ptau=branch!='MEM_NO_PTAU',ptau='P_RAW' if branch=='MEM_RAW_P' else 'P_LOG',tracer=is_pet,base_score=extra_base,modifier=modifier)
   if branch=='MEM_PET_MAIN':Z=np.c_[Z,pet.loc[ids,'CL'].to_numpy()/50]
   x=X[ids]
   if is_long:
    w=np.array([tr['w'] for tr in records]);fac=np.ones(n) if branch=='MEM_EQUAL_PERSON' else np.sqrt(w)
    y=np.array([tr['slope'] for tr in records])*fac
    ZZ=Z*fac[:,None];targetx=x*fac[:,None]
    people=ids.copy();case_local=np.arange(n);slope_load=fac/max(1.,np.median(np.sqrt(w)))
    visit_sd=fac/np.sqrt(w);ar_sd=fac*np.sqrt(np.array([tr['var_AR_pseudo'] for tr in records])/w)
    ar_cross=np.zeros(n);ar_second=np.zeros(n)
    if branch=='MEM_QUADRATIC':
     h=np.array([tr['h'] for tr in records]);q=np.array([tr['q'] for tr in records]);orth=np.maximum(q-h*h/w,0);root=np.sqrt(orth);good=root>1e-6
     y2=np.zeros((n,2));y2[:,0]=y
     sy=np.array([tr['sqy'] for tr in records]);sty=np.array([tr['sty'] for tr in records]);y2[good,1]=(sy[good]-h[good]/w[good]*sty[good])/root[good]
     z2=np.zeros((2*n,Z.shape[1]+1));z2[::2,:-1]=ZZ;z2[::2,-1]=h/np.sqrt(w);z2[1::2,-1]=root;ZZ=z2
     x2=np.zeros((2*n,3));x2[::2]=targetx;targetx=x2;y=y2.ravel()
     people=np.repeat(ids,2);case_local=np.repeat(np.arange(n),2)
     visit_sd=np.column_stack([np.ones(n),good.astype(float)]).ravel();slope_load=np.column_stack([slope_load,np.zeros(n)]).ravel()
     for ii,tr in enumerate(records):
      tt=tr['t'];tc=tt-tt.mean();qq=tt*tt-np.mean(tt*tt);v1=tc/np.sqrt(w[ii]);v2=(qq-h[ii]/w[ii]*tc)/root[ii] if good[ii] else np.zeros(len(tt));RR=np.exp(-.6*abs(tt[:,None]-tt[None,:]));ar_cross[ii]=v1@RR@v2;ar_second[ii]=v2@RR@v2
     ar_sd=np.column_stack([np.array([tr['var_AR_pseudo'] for tr in records])**.5,np.maximum(ar_second,0)**.5]).ravel()
    for tr,i0 in zip(records,ids):
     schedule_audit.append({'design_id':row['design_id'],'RID':clinical.loc[i0,'RID'],'index':i0,'n_visits':tr['n_visits'],'span_years':tr['span'],'s_tt':tr['w'],'highinfo':tr['highinfo'],'baseline_MEM':tr['baseline']})
   else:
    y=np.array([tr['baseline'] for tr in records]) if branch=='MEM_BASE' else pet.loc[ids,'CL'].to_numpy(float)
    ZZ=Z;targetx=x;people=ids.copy();case_local=np.arange(n);visit_sd=np.ones(n);ar_sd=np.ones(n);ar_cross=np.zeros(n);ar_second=np.zeros(n);slope_load=np.ones(n)
   Qcommon=basis(ZZ);ranks=[];A=[];Kparts=[[] for _ in range(3)];fullQ=[];targetden=[]
   for k in range(3):
    nuisance_j=ZZ.copy()
    target=targetx[:,k].copy()
    if modifier is not None:
     main=target.copy();nuisance_j=np.c_[nuisance_j,main]
     target=main*modifier[case_local]
    Qz=basis(nuisance_j);rt=target-Qz@(Qz.T@target);den=rt@rt
    if den<1e-10:reasons.append('TARGET_NOT_ESTIMABLE');den=max(den,1e-30)
    a=rt/den;Qj=np.c_[Qz,rt/np.sqrt(den)];A.append(a);ranks.append(Qj.shape[1]);fullQ.append(Qj);targetden.append(den)
    for axis,c in enumerate(['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex']):
     G=len(group_lookup[c]);lab=np.array([group_lookup[c][str(clinical.loc[i,c])] for i in people])
     B=np.zeros((G,len(people)));B[lab,np.arange(len(people))]=a
     L=B-(B@Qj)@Qj.T;L-=L.mean(axis=0,keepdims=True);L*=np.sqrt(G/(G-1));Kparts[axis].append(L)
   A=np.array(A);K=np.concatenate([np.stack(one,axis=1) for one in Kparts],axis=0)
   rankmax=max(ranks);df_info=n-rankmax
   if df_info<40:reasons.append('INFORMATION_DF_LT40')
   for axis,c in enumerate(['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex']):
    G=f[c].astype(str).nunique();row[c+'_G_cases']=G
    if G<6:reasons.append(c+'_G_LT6')
   condition=np.linalg.cond((A*sigma[:,None])@(A*sigma[:,None]).T)
   if not np.isfinite(condition) or condition>1e12:reasons.append('JOINT_TARGET_ROW_SPACE_WEAK')
   noise=y-Qcommon@(Qcommon.T@y);noise_scale=np.sqrt((noise@noise)/max(1,len(y)-Qcommon.shape[1]))
   noise_scale=max(float(noise_scale),1e-8)
   target_unit=noise_scale/max(1.,np.median(np.sqrt(w))) if is_long else noise_scale
   row.update({'n_pseudo_rows':len(y),'rank_max':rankmax,'information_df':df_info,'n_common_nuisance_rank':Qcommon.shape[1],'max_hat':max(float(np.sum(q*q,axis=1).max()) for q in fullQ),'target_row_Gram_condition':condition,'eligible':not reasons,'gate_reason':';'.join(sorted(set(reasons))),'outcome_noise_scale_nuisance_only':noise_scale,'target_calibration_effect_unit':target_unit,'longitudinal':is_long,'quadratic':branch=='MEM_QUADRATIC'})
   if is_long:row.update({'median_visits':float(np.median([tr['n_visits'] for tr in records])),'min_span_years':min(tr['span'] for tr in records),'median_span_years':float(np.median([tr['span'] for tr in records])),'max_span_years':max(tr['span'] for tr in records),'effective_time_weight_people':float(w.sum()**2/(w@w)),'max_time_weight_fraction':float(w.max()/w.sum())})
   audit.append(row)
   t={'design_id':row['design_id'],'accession':acc,'gene':gene,'branch':branch,'formal':formal,'eligible':not reasons,'ids':ids,'people':people,'n_people':n,'A':A,'K':K,'Uv':Uv,'sigma':sigma,'groups':source_labels,'fullQ':fullQ,'Qcommon':Qcommon,'y':y,'Xcases':x,'ZZ':ZZ,'targetx':targetx,'modifier':modifier,'noise_scale':noise_scale,'target_unit':target_unit,'is_long':is_long,'is_quadratic':branch=='MEM_QUADRATIC','case_local':case_local,'visit_sd':visit_sd,'ar_sd':ar_sd,'ar_cross':ar_cross,'ar_second':ar_second,'slope_load':slope_load,'age':f.age_month_precision.to_numpy(float),'P':f.P_LOG.to_numpy(float),'case_groups':[pd.factorize(f[c].astype(str),sort=True)[0] for c in ['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex']],'summary':row}
   templates[row['design_id']]=t;private_y[row['design_id']]=y
 # Fixed32 formal hypotheses, even gated designs retain P1.
 for acc,gene in CANDIDATES.items():
  for branch in FORMAL_BRANCHES:
   if branch=='MEM_APOE_INTERACTION' and gene!='PDGFRB':continue
   interaction='INTERACTION' in branch
   for test in (['ANY','DIFF'] if interaction else ['SAME','DIFF']):
    registry.append({'hypothesis_id':branch+'__'+acc+'__'+test,'design_id':branch+'__'+acc,'accession':acc,'gene':gene,'branch':branch,'test':test,'family':'PRIMARY6' if branch==PRIMARY else 'SECONDARY26'})
 assert len(registry)==32 and sum(r['family']=='PRIMARY6' for r in registry)==6
 with (PRIVATE/'ANALYSIS_TEMPLATES_RESTRICTED.pkl').open('wb') as f:pickle.dump(templates,f,pickle.HIGHEST_PROTOCOL)
 with (PRIVATE/'EXPOSURE_SCALING_RESTRICTED.pkl').open('wb') as f:pickle.dump(scales,f,pickle.HIGHEST_PROTOCOL)
 csvout(RESULTS/'EXPOSURE_SCALING_AUDIT.csv',scaleaudit);csvout(RESULTS/'ACTUAL_DESIGN_AUDIT.csv',audit)
 csvout(RESULTS/'FIXED_FORMAL_HYPOTHESIS_REGISTRY.csv',registry)
 csvout(PRIVATE/'MODEL_FOLLOWUP_SCHEDULE_AUDIT_RESTRICTED.csv',schedule_audit)
 save_json(RESULTS/'TEMPLATE_FREEZE_RECORD.json',{'time':now(),'protocol_sha256':sha(ROOT/'STAGE10_PROTOCOL_v1.0.md'),'templates_sha256':sha(PRIVATE/'ANALYSIS_TEMPLATES_RESTRICTED.pkl'),'n_designs':len(templates),'n_formal_designs':sum(t['formal'] for t in templates.values()),'eligible':sum(t['eligible'] for t in templates.values()),'new_target_coefficients_or_P_values_fitted':False,'code_sha256':{p.name:sha(p) for p in [Path(__file__),ROOT/'code/engine.py',ROOT/'code/config.py']}})
 log(f'{len(templates)} three-assay designs prepared;{sum(t["eligible"] for t in templates.values())} pass geometry;32 hypotheses fixed before new target coefficients')
 print(pd.DataFrame(audit)[['gene','branch','n_people','rank_max','information_df','eligible','gate_reason']].to_string(index=False))
if __name__=='__main__':
 assert_unsealed();run()
