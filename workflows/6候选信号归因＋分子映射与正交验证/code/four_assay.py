from config import *
from designs import make
from model import bh,wilson
from scipy.stats import norm,chi2
from concurrent.futures import ProcessPoolExecutor,as_completed
import pickle

PAIRS=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]

def moments(Y,d):
 if Y.ndim==2:Y=Y[None,:,:]
 B,n,m=Y.shape;flat=Y.transpose(1,0,2).reshape(n,-1);Q=d['Q'];E=(flat-Q@(Q.T@flat)).reshape(n,B,m).transpose(1,0,2)
 prod=np.stack([E[:,:,i]*E[:,:,j] for i,j in PAIRS],axis=2);c=prod.sum(axis=1)/d['df'];IF=(prod-prod.mean(axis=1,keepdims=True))/d['df']
 V=np.zeros((B,6,6));flat=IF.transpose(1,0,2).reshape(n,-1)
 for one in d['onehots']:
  g=len(one);s=(one@flat).reshape(g,B,6).transpose(1,0,2);s-=s.mean(axis=1,keepdims=True);V+=np.einsum('bgi,bgj->bij',s,s)*g/(g-1)
 h=c[:,1]*c[:,4]-c[:,2]*c[:,3];J=np.zeros((B,6));J[:,1]=c[:,4];J[:,4]=c[:,1];J[:,2]=-c[:,3];J[:,3]=-c[:,2]
 variance=np.einsum('bi,bij,bj->b',J,V,J);stat=h*h/variance;p=chi2.sf(stat,1)
 p[~np.isfinite(variance)|(variance<=0)]=np.nan
 return p,h,c,V,E

def generate(d,scenario,B,rng):
 n=d['n'];sig=np.array([.6,.8,1.2,.7]);u=rng.normal(size=(B,n));z=rng.normal(size=(B,n,4))
 r1=.4 if scenario in ['F02','F03','F04','F05','F06'] else 0.;r2=.6 if r1 else 0.
 if scenario=='F03':r1=-.4;r2=.4
 z[:,:,1]=r1*z[:,:,0]+np.sqrt(1-r1*r1)*z[:,:,1];z[:,:,3]=r2*z[:,:,2]+np.sqrt(1-r2*r2)*z[:,:,3]
 if scenario=='F06':
  z[:,:,:2]/=np.sqrt(rng.chisquare(5,size=(B,n,1))/3);z[:,:,2:]/=np.sqrt(rng.chisquare(5,size=(B,n,1))/3);u=rng.standard_t(5,size=(B,n))*np.sqrt(3/5)
 e=z*sig
 if scenario=='F05':
  scale=np.exp(np.clip(.3*d['P'],-3,3));e*=scale[None,:,None]/np.sqrt(np.mean(scale**2))
 if scenario in ['F04','F06']:
  for axis in range(4):
   lab=d['groups'][:,axis];one=d['onehots'][axis];age=d['q_age'];q=age-(one@age/one.sum(axis=1))[lab]
   coef=rng.normal(scale=.35,size=(B,len(one),5));u+=coef[:,lab,0]*q[None,:]
   # Independent across affinity/MS blocks; within-block errors may be correlated.
   coef[:,:,2]=.4*coef[:,:,1]+np.sqrt(.84)*coef[:,:,2];coef[:,:,4]=.6*coef[:,:,3]+.8*coef[:,:,4]
   e+=coef[:,lab,1:]*q[None,:,None]
 return (.55*d['P'][None,:,None]+u[:,:,None]+e)/np.sqrt(1+sig**2)

def cell(task):
 key,d,scenario=task;seed=int(np.random.SeedSequence([2026102906,int(hashlib.sha256(key.encode()).hexdigest()[:8],16),int(scenario[1:])]).generate_state(1,dtype=np.uint64)[0])
 rng=np.random.default_rng(seed);p=[]
 for i in range(0,4000,64):p.extend(moments(generate(d,scenario,min(64,4000-i),rng),d)[0])
 p=np.array(p);valid=np.isfinite(p);k=int((valid&(p<=.05)).sum());lo,hi=wilson(k,len(p));file='FOUR__'+key+'__'+scenario+'.parquet'
 pd.DataFrame({'replicate':np.arange(len(p)),'p_tetrad':p}).to_parquet(SIM/file,index=False)
 return {'design_id':key,'scenario':scenario,'reps':len(p),'n_failed':int((~valid).sum()),'n_rejected05':k,'rate':k/len(p),'Wilson95_low':lo,'Wilson95_high':hi,'seed':str(seed),'file':file}

def run():
 frames=pickle.load((RESTRICTED/'ORTHOGONAL_SUPPORT_CASE_FRAMES_RESTRICTED.pkl').open('rb'));screen=pd.read_csv(RESULTS/'ORTHOGONAL_SUPPORT_DESIGN_SCREEN.csv')
 old=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv').set_index('key');raw=np.load(S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz');nsc=np.load(RESTRICTED/'SUPPORT_N_S_CACHE_RESTRICTED.npz');targ=np.load(RESTRICTED/'TARGETED_LOG_RATIOS_RESTRICTED.npz')
 ds={};ys={};info={};rows=[];corr_rows=[]
 for r in screen.loc[screen.cohort.eq('S5_SAME_DATE')].itertuples():
  frame=frames[r.design_id]['frame'].copy();frame['TMT_plex']=[str(old.loc[k,'TMT_plex']) if k in old.index else '' for k in frame.key]
  j=list(nsc['accessions']).index(r.accession);jr=list(raw['accessions']).index(r.accession);k=list(targ['feature_keys']).index(r.feature_key);idx={key:i for i,key in enumerate(raw['keys'])}
  T=np.array([raw['T_GIS'][idx[key],jr] if key in idx else np.nan for key in frame.key]);Y=np.c_[nsc['N'][:,j],nsc['S'][:,j],T,targ['Y'][:,k]]
  mask=frames[r.design_id]['mask']&np.isfinite(Y).all(axis=1);d=make(frame,mask,tech=('NULISA_PlateID','SOMASCAN_PlateId','TMT_plex','TARGETED_Plate'))
  ident=r.accession+'__'+r.peptide;ds[ident]=d;ys[ident]=Y[d['indices']];info[ident]=r
 save_json(RESULTS/'FOUR_ASSAY_PLAN_FROZEN.json',{'post_targeted_triplet_results_extension':True,'all_five_annotation_selected_peptides':True,
   'model':'rank1 affinity-versus-MS cross-covariance block; correlated errors allowed within affinity pair and within MS pair',
   'before_four_assay_statistics':True,'equivalence_ratio_factors':[1.25,1.5],'conditional_cross_block_error_assumption_not_verified':True,
   'code_sha256':sha(Path(__file__)),'independent_null_reps_per_design_scenario':4000})
 with (RESTRICTED/'FOUR_ASSAY_DESIGNS_RESTRICTED.pkl').open('wb') as f:pickle.dump(ds,f,pickle.HIGHEST_PROTOCOL)
 cal=[];tasks=[(key,d,s) for key,d in ds.items() for s in ['F01','F02','F03','F04','F05','F06']]
 with ProcessPoolExecutor(max_workers=6) as pool:
  for f in as_completed([pool.submit(cell,t) for t in tasks]):cal.append(f.result())
 cal=pd.DataFrame(cal).sort_values(['design_id','scenario']);csvout(RESULTS/'FOUR_ASSAY_NULL_CALIBRATION.csv',cal)
 for ident,d in ds.items():
  p,h,c,V,E=moments(ys[ident],d);c=c[0];V=V[0];E=E[0];r=info[ident];calq=cal.loc[cal.design_id.eq(ident)];passed=bool((calq.Wilson95_high<=.06).all() and (calq.n_failed/calq.reps<=.005).all())
  C=E.T@E/d['df'];R=C/np.sqrt(np.outer(np.diag(C),np.diag(C)))
  base={'design_id':ident,'accession':r.accession,'gene':r.gene,'peptide':r.peptide,'n':d['n'],'rank_X':d['rank'],'residual_df':d['df'],
   'calibration_passed':passed,'tetrad':float(h[0]),'p_tetrad':float(p[0]),'p_for_family':float(p[0]) if passed else 1.}
  for i,j in PAIRS:corr_rows.append({'accession':r.accession,'gene':r.gene,'peptide':r.peptide,'n':d['n'],'assay_i':['NULISA','Soma','TMT','targetMS'][i],
    'assay_j':['NULISA','Soma','TMT','targetMS'][j],'residual_correlation_after_same_four_assay_design':R[i,j],'point_diagnostic_not_independent_replication':True})
  cross=[1,2,3,4];strong=bool(all(c[j]-1.959964*np.sqrt(V[j,j])>0 for j in cross));base['cross_block_cov_CI95_all_positive']=strong
  if all(c[j]>0 for j in cross):
   L=np.log(c[1])+np.log(c[4])-np.log(c[2])-np.log(c[3]);g=np.zeros(6);g[1]=1/c[1];g[4]=1/c[4];g[2]=-1/c[2];g[3]=-1/c[3];se=np.sqrt(g@V@g)
   base.update({'log_direction_ratio':L,'direction_ratio':np.exp(L),'direction_ratio_CI95_low':np.exp(L-1.959964*se),'direction_ratio_CI95_high':np.exp(L+1.959964*se)})
   for tol in [1.25,1.5]:base['p_equiv_factor_'+str(tol)]=max(norm.sf((L+np.log(tol))/se),norm.sf((np.log(tol)-L)/se)) if strong and passed else 1.
  else:
   for tol in [1.25,1.5]:base['p_equiv_factor_'+str(tol)]=1.
  rows.append(base)
 out=pd.DataFrame(rows);out['q_BY5_tetrad_exploratory']=bh(out.p_for_family.to_numpy(),True)
 for tol in [1.25,1.5]:out['q_BY5_equiv_factor_'+str(tol)]=bh(out['p_equiv_factor_'+str(tol)].to_numpy(),True)
 csvout(RESULTS/'FOUR_ASSAY_CROSSBLOCK_RESULTS.csv',out);csvout(RESULTS/'FOUR_ASSAY_RESIDUAL_CORRELATIONS.csv',corr_rows)
 save_json(RESULTS/'FOUR_ASSAY_RUN_SUMMARY.json',{'independent_replicates':len(tasks)*4000,'designs':len(ds),'new_primary_discovery_validation':False,'postdiscovery_exploratory_extension':True})
 log('Four-assay cross-block rank/working-equivalence exploration complete; no Stage5 results changed')
 print(out[['gene','peptide','n','calibration_passed','p_tetrad','q_BY5_tetrad_exploratory','direction_ratio','cross_block_cov_CI95_all_positive']].to_string(index=False))

if __name__=='__main__':run()
