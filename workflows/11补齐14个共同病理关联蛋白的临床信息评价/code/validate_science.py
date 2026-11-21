from config import *
from scipy.stats import norm,chi2
from scipy.linalg import lstsq
from datetime import datetime
from simulation_dgp import innovations,patterns
import pickle

def run():
 checks=[]
 def add(name,ok,detail=''):checks.append({'check':name,'passed':bool(ok),'detail':str(detail)})
 def close(name,a,b,rtol=1e-7,atol=1e-10):
  add(name,np.allclose(a,b,rtol=rtol,atol=atol,equal_nan=True),float(np.nanmax(abs(np.asarray(a,dtype=float)-np.asarray(b,dtype=float)))))
 T=pickle.load((PRIVATE/'ANALYSIS_TEMPLATES_RESTRICTED.pkl').open('rb'));scales=pickle.load((PRIVATE/'EXPOSURE_SCALING_RESTRICTED.pkl').open('rb'))
 clinical=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv')
 real=pd.read_csv(RESULTS/'ALL144_DESIGN_CLINICAL_RESULTS.csv').set_index('design_id');effects=pd.read_csv(RESULTS/'PLATFORM_CLINICAL_EFFECTS_ALL144.csv')
 for r in pd.read_csv(RESULTS/'ENGINE_INDEPENDENT_CHECKS.csv').itertuples():add('SYNTHETIC/'+r.check,r.passed,r.max_abs_error)
 D=np.array([[-1.,1.,0.],[-1.,0.,1.]])
 for key,t in T.items():
  bb=[];parts=[]
  for j in range(3):
   X=np.c_[t['ZZ'],t['targetx'][:,j]];fit=lstsq(X,t['y'],cond=max(X.shape)*np.finfo(float).eps,lapack_driver='gelsy')[0]
   res=t['y']-X@fit;bb.append(fit[-1]);scores=[]
   for c in ['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex']:
    levels=sorted(clinical.loc[scales[t['accession']]['ids'],c].astype(str).unique());lookup={v:i for i,v in enumerate(levels)}
    labels=np.array([lookup[str(clinical.loc[i,c])] for i in t['people']]);sg=np.bincount(labels,weights=t['A'][j]*res,minlength=len(levels));sg-=sg.mean();sg*=np.sqrt(len(levels)/(len(levels)-1));scores.append(sg)
   parts.append(np.concatenate(scores));ref=effects.loc[effects.design_id.eq(key)&effects.assay.eq(['N','S','T'][j])].iloc[0]
   close('REAL_INDEPENDENT_QR/'+key+str(j),fit[-1],ref.b_original_scale)
  b=np.array(bb);U=np.array(parts).T;V=U.T@U;gamma=b*t['sigma']
  Us=U*t['sigma'][None,:]+b[None,:]*t['Uv']/(2*t['sigma'][None,:]);Vs=Us.T@Us;z=b/np.sqrt(np.diag(V))
  ps=min(1.,2*min(norm.sf(z).max(),norm.cdf(z).max()));delta=D@gamma;Vd=D@Vs@D.T;pd_=chi2.sf(delta@np.linalg.solve(Vd,delta),2)
  close('REAL_IUT/'+key,ps,real.loc[key,'p_same_computed'])
  close('REAL_PAIRED_WALD/'+key,pd_,real.loc[key,'p_diff_computed'])
  for j,a in enumerate(['N','S','T']):
   rr=effects.loc[effects.design_id.eq(key)&effects.assay.eq(a)].iloc[0]
   close('REAL_FULL_JOINT_SE/'+key+a,np.sqrt(Vs[j,j]),rr.standardized_SE_full_joint)
 formal=pd.read_csv(RESULTS/'FORMAL204_CLINICAL_RESULTS.csv');dec=pd.read_csv(RESULTS/'CALIBRATION_QUALIFICATION_DECISIONS.csv').set_index('design_id')
 add('FIXED34',set(formal.accession)==set(FAMILY34));add('FIXED14',set(formal.loc[formal.focus14,'gene'])==set(FOCUS14))
 add('PRIMARY68',int(formal.family.eq('PRIMARY68').sum())==68);add('SECONDARY136',int(formal.family.eq('SECONDARY136').sum())==136);add('FORMAL204',len(formal)==204)
 for r in formal.itertuples():add('GATED_P1/'+r.hypothesis_id,bool(r.calibration_qualified and r.math_valid) or r.p_for_family==1.)
 def independent_by(v):
  v=np.asarray(v);m=len(v);ix=np.argsort(v);q=np.minimum.accumulate((v[ix]*m*sum(1/i for i in range(1,m+1))/np.arange(1,m+1))[::-1])[::-1];out=np.ones(m);out[ix]=np.minimum(1,q);return out
 for family,d in formal.groupby('family'):close('INDEPENDENT_BY/'+family,independent_by(d.p_for_family),d.q_BY_within_family)
 close('INDEPENDENT_BY204',independent_by(formal.p_for_family),formal.q_BY_combined204)
 sens=pd.read_csv(RESULTS/'PRESPECIFIED_SENSITIVITY_RESULTS.csv')
 add('SENSITIVITY84',len(sens)==84)
 for br,d in sens.groupby('branch'):close('SENSITIVITY_BY28/'+br,independent_by(d.p_for_diagnostic_family),d.q_BY28_diagnostic)
 focus=pd.read_csv(RESULTS/'FOCUS14_PRIMARY28_FROM_FULL_FAMILY.csv')
 add('FOCUS14_28',len(focus)==28)
 for r in focus.itertuples():
  rr=formal.loc[formal.hypothesis_id.eq(r.hypothesis_id)].iloc[0];close('FOCUS_FULL_FAMILY_Q/'+r.hypothesis_id,r.q_BY_within_family,rr.q_BY_within_family)
 for r in pd.read_csv(RESULTS/'STAGE10_OVERLAPPING_DESIGN_REPRODUCTION.csv').itertuples():add('STAGE10_REPRODUCTION/'+r.design_id,r.passed,max(r.max_abs_P_error,r.max_abs_effect_error))
 aux=pd.read_csv(RESULTS/'FOCUS14_PTAU_CONDITIONING_DECOMPOSITION.csv')
 for r in aux.itertuples():
  add('CONDITIONING_IDENTITY/'+r.gene+r.assay,r.identity_abs_error<1e-8)
  mainref=effects.loc[effects.gene.eq(r.gene)&effects.assay.eq(r.assay)&effects.branch.eq('MEM5Y')].iloc[0]
  nopref=effects.loc[effects.gene.eq(r.gene)&effects.assay.eq(r.assay)&effects.branch.eq('MEM_NO_PTAU')].iloc[0]
  close('CONDITIONING_MAIN_MATCH/'+r.gene+r.assay,r.standardized_with_pTau,mainref.standardized_effect)
  close('CONDITIONING_NOP_MATCH/'+r.gene+r.assay,r.standardized_without_pTau,nopref.standardized_effect)
 inv=pd.read_csv(RESULTS/'SIMULATION_FILE_INVENTORY.csv',dtype={'seed':str});recounts=[];replay=[];total=0
 zcrit=norm.ppf(.975)
 def upper(k,n):
  pp=k/n;den=1+zcrit*zcrit/n
  return (pp+zcrit*zcrit/(2*n)+zcrit*np.sqrt(pp*(1-pp)/n+zcrit*zcrit/(4*n*n)))/den
 for count,r in enumerate(inv.itertuples(),1):
  f=SIM/r.file;d=pd.read_parquet(f);total+=len(d);meta=json.loads((SIM/(r.design_id+'__'+r.scenario+'.json')).read_text(encoding='utf-8'))
  add('SIM_SHA/'+r.design_id+r.scenario,sha(f)==r.file_sha256)
  add('SIM_REPLICATES/'+r.design_id+r.scenario,len(d)==r.reps and np.array_equal(d.replicate,np.arange(r.reps)))
  allcounts=True;allci=True
  for row in meta['summary']:
   k=int(np.count_nonzero(d[row['test']+'_'+row['pattern']].to_numpy()<=row['alpha']));u=upper(k,len(d));ok=k==row['rejected'] and np.isclose(u,row['Wilson95_high'],rtol=1e-12,atol=1e-12)
   recounts.append({'design_id':r.design_id,'scenario':r.scenario,'pattern':row['pattern'],'test':row['test'],'alpha_name':row['alpha_name'],'rejected_recounted':k,'Wilson95_upper_recounted':u,'passed':ok});allcounts&=k==row['rejected'];allci&=np.isclose(u,row['Wilson95_high'],rtol=1e-12,atol=1e-12)
  add('SIM_COUNTS/'+r.design_id+r.scenario,allcounts);add('SIM_WILSON/'+r.design_id+r.scenario,allci)
  t=T[r.design_id];rng=np.random.default_rng(int(r.seed));Y=innovations(rng,t,r.scenario,min(BATCH,r.reps))[:,:3]
  pats=patterns(t);As=t['sigma'][:,None]*t['A'];names=['ZERO','COMMON_0.35_1','ONE_ZERO_0_1',next(pp['name'] for pp in pats if pp['name'].startswith('POWER_BALANCED') and pp['name'].endswith('_1'))]
  for pp in pats:
   if pp['name'] not in names:continue
   mu=As.T@np.linalg.solve(As@As.T,pp['effect']);YY=Y+mu[:,None]
   for j in range(3):
    bb=t['A']@YY[:,j];ss=np.einsum('gkn,n->gk',t['K'],YY[:,j],optimize=True);vv=ss.T@ss;zz=bb/np.sqrt(np.diag(vv));ps=min(1.,2*min(norm.sf(zz).max(),norm.cdf(zz).max()))
    gg=bb*t['sigma'];uu=ss*t['sigma'][None,:]+bb[None,:]*t['Uv']/(2*t['sigma'][None,:]);dd=D@gg;vvv=D@(uu.T@uu)@D.T;pd_=chi2.sf(dd@np.linalg.solve(vvv,dd),2)
    for typ,value in [('SAME',ps),('DIFF',pd_)]:
     expected=d.loc[j,typ+'_'+pp['name']];ok=np.isclose(value,expected,rtol=1e-6,atol=1e-10)
     replay.append({'design_id':r.design_id,'scenario':r.scenario,'pattern':pp['name'],'replicate':j,'test':typ,'P_replayed':value,'P_saved':expected,'passed':ok})
  if count%40==0:log(f'Independent saved-P recount/seed replay {count}/{len(inv)}')
 csvout(RESULTS/'INDEPENDENT_SIMULATION_RECOUNTS.csv',recounts);csvout(RESULTS/'INDEPENDENT_SEED_REPLAY_CHECKS.csv',replay)
 add('ALL_RECOUNTS',all(r['passed'] for r in recounts),len(recounts));add('ALL_REPLAY',all(r['passed'] for r in replay),len(replay));add('NOISE_COUNT',total==53100000,total)
 full=pd.read_csv(RESULTS/'CALIBRATION_AND_POWER_FULL.csv');comp=pd.read_csv(RESULTS/'COMPONENT_NULL_CALIBRATION.csv')
 for key,t in T.items():
  part=full.loc[full.design_id.eq(key)&full.true_null];nom=part.loc[part.alpha_name.eq('NOMINAL')];tail=part.loc[~part.alpha_name.eq('NOMINAL')];cc=comp.loc[comp.design_id.eq(key)];ii=inv.loc[inv.design_id.eq(key)]
  qualified=bool(t['eligible'] and nom.Wilson95_high.le(.06).all() and tail.Wilson95_high.le(tail.alpha).all() and cc.Wilson95_high.le(.035).all() and (ii.n_failed/ii.reps).max()<=.005)
  add('QUALIFICATION_RECOUNT/'+key,qualified==bool(dec.loc[key,'qualified']))
  add('COMPONENT_WILSON/'+key,all(np.isclose(upper(int(r.rejected025),int(r.reps)),r.Wilson95_high,rtol=1e-12,atol=1e-12) for r in cc.itertuples()))
 for acc,g in FAMILY34.items():
  ts=[T[br+'__'+acc] for br in FORMAL_BRANCHES]
  if g in FOCUS14:ts += [T[br+'__'+acc] for br in SENSITIVITIES]
  for t in ts:close('FIXED_SCALE/'+t['design_id'],t['sigma'],ts[0]['sigma'],rtol=0,atol=0)
 for r in pd.read_csv(RESULTS/'INPUT_SOURCE_MANIFEST.csv').itertuples():add('INPUT_UNCHANGED/'+r.path,sha(r.path)==r.sha256)
 prior=pd.read_csv(RESULTS/'PRIOR_TEN_STAGE_SEAL_CHECKS_INITIAL.csv');cache={};old=[]
 for r in prior.itertuples():
  if r.path not in cache:cache[r.path]=sha(r.path) if Path(r.path).exists() else ''
  old.append({'path':r.path,'expected_sha256':r.expected_sha256,'actual_sha256':cache[r.path],'passed':cache[r.path]==r.expected_sha256})
 csvout(RESULTS/'PRIOR_TEN_STAGE_SEAL_CHECKS_FINAL.csv',old);add('ALL_PRIOR_SEALED_UNCHANGED',all(r['passed'] for r in old),len(old))
 access=json.loads((RESULTS/'ACCESS_AND_PROTOCOL_RECORD.json').read_text(encoding='utf-8'));freeze=json.loads((RESULTS/'CALIBRATION_FREEZE_BEFORE_RUN.json').read_text(encoding='utf-8'));realaccess=json.loads((RESULTS/'FIRST_REAL_TARGET_ASSOCIATION_ACCESS.json').read_text(encoding='utf-8'));cal=json.loads((RESULTS/'CALIBRATION_RUN_SUMMARY.json').read_text(encoding='utf-8'))
 add('PROTOCOL_FIXED',access['protocol_sha256']==sha(ROOT/'STAGE11_PROTOCOL_v1.0.md')==realaccess['protocol_sha256'])
 add('CALIBRATION_BEFORE_NEW_REAL_TARGET',datetime.fromisoformat(cal['time'])<=datetime.fromisoformat(realaccess['time']))
 for n,expected in freeze['code_sha256'].items():add('FROZEN_ENGINE/'+n,sha(ROOT/'code'/n)==expected)
 csvout(RESULTS/'SCIENTIFIC_QA_CHECKS.csv',checks)
 out={'time':now(),'checks':len(checks),'passed':sum(r['passed'] for r in checks),'failed':sum(not r['passed'] for r in checks),'independent_real_designs':len(T),'saved_P_recount_rows':len(recounts),'seed_replay_P_values':len(replay),'independent_noise_draws':total,'prior_sealed_entries':len(old),'prior_unique_paths':len(cache),'marginal_component_P_scope':'All aggregate counts with independent Wilson checking; deterministic replay included; full marginal arrays not saved','independent_external_validation':False}
 save_json(RESULTS/'SCIENTIFIC_QA_SUMMARY.json',out)
 assert out['failed']==0,[r for r in checks if not r['passed']][:10]
 log(f'Scientific QA complete: {len(checks)}/{len(checks)} passed.')
if __name__=='__main__':assert_unsealed();run()

