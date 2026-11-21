from config import *
from engine import stats
from calibrate import patterns,innovations,qualification
from scipy.linalg import lstsq
from scipy.stats import norm,chi2
from datetime import datetime
import pickle
def run():
 checks=[]
 def add(n,ok,d=''):checks.append({'check':n,'passed':bool(ok),'detail':str(d)})
 def close(n,a,b):
  add(n,np.allclose(a,b,rtol=1e-7,atol=1e-9,equal_nan=True),float(np.nanmax(abs(np.asarray(a,float)-np.asarray(b,float)))))
 for r in pd.read_csv(RESULTS/'ENGINE_INDEPENDENT_CHECKS.csv').itertuples():add('SYNTHETIC/'+r.check,r.passed,r.max_abs_error)
 T=pickle.load((PRIVATE/'TEST_TEMPLATES_RESTRICTED.pkl').open('rb'));c=pd.read_csv(PRIVATE/'COMMON_CLINICAL_RESTRICTED.csv');scale=pickle.load((PRIVATE/'FROZEN_SCALES_RESTRICTED.pkl').open('rb'));real=pd.read_csv(RESULTS/'ALL42_DIRECTIONAL_EXPLANATION_RESULTS.csv').set_index('test_id');repro=[]
 oldr=pd.read_csv(S11/'results/ALL144_DESIGN_CLINICAL_RESULTS.csv').set_index('design_id');olde=pd.read_csv(S11/'results/PLATFORM_CLINICAL_EFFECTS_ALL144.csv')
 for key,t in T.items():
  b=[];parts=[]
  for j in range(t['A'].shape[0]):
   X=np.c_[t['ZZ_list'][j],t['targetx_list'][j]];bb=lstsq(X,t['y'],cond=max(X.shape)*np.finfo(float).eps,lapack_driver='gelsy')[0];e=t['y']-X@bb;b.append(bb[-1]);ss=[]
   for axis,col in enumerate(['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex']):
    labels=t['groups'][axis];G=t['full_G'][axis];s=np.bincount(labels,weights=t['A'][j]*e,minlength=G);s-=s.mean();s*=np.sqrt(G/(G-1));ss.append(s)
   parts.append(np.concatenate(ss))
  b=np.array(b);U=np.array(parts).T;g=b*t['sigma'];Us=U*t['sigma'][None,:]+b[None,:]*t['Uv']/(2*t['sigma'][None,:]);h=t['L']@g;Uh=Us@t['L'].T;V=Uh.T@Uh;p=chi2.sf(h@np.linalg.solve(V,h),2);r=real.loc[key]
  close('REAL_P/'+key,p,r.p_computed)
  for j,a in enumerate(['S_MINUS_N','T_MINUS_N']):
   close('REAL_CONTRAST/'+key+a,h[j],r[a]);close('REAL_SE/'+key+a,np.sqrt(V[j,j]),r[a+'_SE'])
  if t['model'] in ['NO_P','LOG_P']:
   oldkey=('MEM5Y' if t['model']=='LOG_P' else 'MEM_NO_PTAU')+'__'+t['accession'];oo=oldr.loc[oldkey];ee=olde.loc[olde.design_id.eq(oldkey)].set_index('assay').loc[['N','S','T']]
   errors=[float(abs(p-oo.p_diff_computed)),float(abs(g-ee.standardized_effect.to_numpy()).max()),float(abs(np.sqrt(np.diag(Us.T@Us))-ee.standardized_SE_full_joint.to_numpy()).max())]
   repro.append({'test_id':key,'same_n':int(oo.n_people)==t['n_people'],'max_abs_error':max(errors),'passed':max(errors)<1e-8 and int(oo.n_people)==t['n_people']})
 csvout(RESULTS/'STAGE11_FULL_STATISTIC_REPRODUCTION.csv',repro)
 for r in repro:add('STAGE11/'+r['test_id'],r['passed'],r['max_abs_error'])
 f=real.reset_index();add('REGISTRY42',len(f)==42);add('PRIMARY9',f.family.eq('PRIMARY9').sum()==9);add('SECONDARY33',f.family.eq('SECONDARY33').sum()==33)
 def independent_by(a):
  a=np.array(a);m=len(a);ix=np.argsort(a);h=sum(1/j for j in range(1,m+1));v=np.minimum.accumulate((a[ix]*m*h/np.arange(1,m+1))[::-1])[::-1];z=np.ones(m);z[ix]=np.minimum(1,v);return z
 for fam,d in f.groupby('family'):close('BY/'+fam,independent_by(d.p_for_family),d.q_BY_family)
 close('BY42',independent_by(f.p_for_family),f.q_BY_combined42)
 for r in f.itertuples():add('GATED_P1/'+r.test_id,bool(r.calibration_qualified and r.math_valid) or r.p_for_family==1.)
 inv=pd.read_csv(RESULTS/'SIMULATION_FILE_INVENTORY.csv',dtype={'seed':str});recounts=[];replay=[];nnoise=0;zc=1.959963984540054
 def upper(k,n):
  p=k/n;d=1+zc*zc/n
  return (p+zc*zc/(2*n)+zc*np.sqrt(p*(1-p)/n+zc*zc/(4*n*n)))/d
 for i,r in enumerate(inv.itertuples(),1):
  file=SIM/r.file;d=pd.read_parquet(file);meta=json.loads(file.with_suffix('.json').read_text(encoding='utf-8'));nnoise+=len(d);add('SIM_SHA/'+r.test_id+r.scenario,sha(file)==r.file_sha256);add('SIM_N/'+r.test_id+r.scenario,len(d)==r.reps and np.array_equal(d.replicate,np.arange(r.reps)))
  allgood=True
  for row in meta['summary']:
   k=int((d[row['pattern']].to_numpy()<=row['alpha']).sum());u=upper(k,len(d));ok=k==row['rejected'] and np.isclose(u,row['Wilson95_high'],rtol=1e-12,atol=1e-12);allgood &=ok;recounts.append({'test_id':r.test_id,'scenario':r.scenario,'pattern':row['pattern'],'alpha_name':row['alpha_name'],'rejected_recounted':k,'Wilson_upper_recounted':u,'passed':ok})
  add('SIM_ALL_COUNTS/'+r.test_id+r.scenario,allgood)
  t=T[r.test_id];rng=np.random.default_rng(int(r.seed));Y=innovations(rng,t,r.scenario,min(BATCH,r.reps))[:,:3];pp=patterns(t)
  for item in [pp[0],pp[1],pp[6],pp[7]]:
   YY=Y+item['mu'][:,None]
   for j in range(3):
    b=t['A']@YY[:,j];U=np.einsum('gpn,n->gp',t['K'],YY[:,j],optimize=True);Us=U*t['sigma'][None,:]+b[None,:]*t['Uv']/(2*t['sigma'][None,:]);h=t['L']@(b*t['sigma']);v=(Us@t['L'].T).T@(Us@t['L'].T);pv=chi2.sf(h@np.linalg.solve(v,h),2);expected=d.loc[j,item['name']];ok=np.isclose(pv,expected,rtol=1e-6,atol=1e-10);replay.append({'test_id':r.test_id,'scenario':r.scenario,'pattern':item['name'],'replicate':j,'P_replayed':pv,'P_saved':expected,'passed':ok})
  if i%30==0:log(f'Independent saved-P recount and replay {i}/{len(inv)}')
 csvout(RESULTS/'INDEPENDENT_SIMULATION_RECOUNTS.csv',recounts);csvout(RESULTS/'INDEPENDENT_SEED_REPLAY_CHECKS.csv',replay);add('ALL_RECOUNTS',all(r['passed'] for r in recounts),len(recounts));add('ALL_REPLAY',all(r['passed'] for r in replay),len(replay));add('NOISE_COUNT',nnoise==12600000,nnoise)
 full=pd.read_csv(RESULTS/'CALIBRATION_AND_POWER_FULL.csv',keep_default_na=False);dec=pd.read_csv(RESULTS/'CALIBRATION_QUALIFICATION_DECISIONS.csv').set_index('test_id');audit=pd.read_csv(RESULTS/'ACTUAL_DESIGN_AUDIT.csv');ind=qualification(audit,full,inv).set_index('test_id')
 add('CSV_NULL_MEAN_LABEL_PRESERVED',full.mean_kind.eq('NULL').sum()==42*5*7*4)
 for key in ind.index:add('QUALIFICATION/'+key,ind.loc[key,'qualified']==dec.loc[key,'qualified'])
 for name in ['OVB_ASSAY_DECOMPOSITION.csv','OVB_PAIRED_CHANGE_RECONSTRUCTION.csv']:
  d=pd.read_csv(RESULTS/name);add('EXACT_ALGEBRA/'+name,d.identity_abs_error.max()<1e-8,d.identity_abs_error.max())
 add('FWL_MATCH',pd.read_csv(RESULTS/'FWL_WEIGHTED_MOMENT_EXPLANATION.csv').joint_formula_matches_full_model.all())
 add('BACKGROUND_CORE_IDENTITY',pd.read_csv(RESULTS/'BACKGROUND_CORE_IDENTITY_CHECKS.csv').passed.all())
 bg=pd.read_csv(RESULTS/'BACKGROUND_FEATURE_AUDIT.csv');kept=bg.loc[bg.included_reference];add('BACKGROUND90',kept.finite_fraction.ge(.9).all());add('NULISA_LOD80',kept.loc[kept.platform.eq('NULISA'),'detected_fraction'].ge(.8).all());add('BACKGROUND_ALL_IDENTITIES_ALLOWED',kept.allowed_identity.all())
 excluded=set(pd.read_csv(S11/'results/FROZEN_FAMILY34_AND_FOCUS14.csv').accession)|{'P10636','P05067','P02649'}
 for p,dd in kept.groupby('platform'):add('NO_TARGETS_IN_BACKGROUND/'+p,not dd.accessions.astype(str).str.split(r'[|;,]').map(lambda aa:bool({a.split('-')[0] for a in aa}&excluded)).any(),len(dd))
 toy=pd.read_csv(RESULTS/'SINGLE_TRUE_SIGNAL_SYNTHETIC_REPRODUCTION.csv');add('TOY_POPULATION_REPRODUCTION',toy.within_6_MC_SE.all())
 # Exact visit-level check using a new nonlinear/background model and synthetic visits.
 models=pickle.load((PRIVATE/'MODELS_RESTRICTED.pkl').open('rb'));data=pickle.load((PRIVATE/'ASSEMBLED_CLINICAL_RESTRICTED.pkl').open('rb'));rng=np.random.default_rng(120101)
 for acc,gene in GENES.items():
  t=models[gene+'__BG2_RCS_P'];Y=[];pseudo=[];Xdesign=[[],[],[]]
  for k,i in enumerate(t['ids']):
   tr=data['trajectories'][(int(i),'ADNI_MEM',5)];tt=tr['t']-np.mean(tr['t']);yy=rng.normal(size=len(tt));yy-=yy.mean();Y.extend(yy);pseudo.append(tt@yy/np.sqrt(tr['w']))
   for j in range(3):Xdesign[j].append(tt[:,None]*np.r_[t['ZZ_list'][j][k]/np.sqrt(tr['w']),t['X'][k,j]][None,:])
  for j in range(3):
   X=np.concatenate(Xdesign[j]);bb=lstsq(X,np.array(Y),cond=max(X.shape)*np.finfo(float).eps,lapack_driver='gelsy')[0][-1];close('EXACT_VISIT_FE/'+gene+str(j),bb,t['A'][j]@pseudo)
 for r in pd.read_csv(RESULTS/'INPUT_SOURCE_MANIFEST.csv').itertuples():add('INPUT_UNCHANGED/'+r.path,sha(r.path)==r.sha256)
 prior=pd.read_csv(RESULTS/'PRIOR_ELEVEN_STAGE_SEAL_CHECKS_INITIAL.csv');cache={};old=[]
 for r in prior.itertuples():
  if r.path not in cache:cache[r.path]=sha(r.path) if Path(r.path).exists() else ''
  old.append({'path':r.path,'expected_sha256':r.expected_sha256,'actual_sha256':cache[r.path],'passed':cache[r.path]==r.expected_sha256})
 csvout(RESULTS/'PRIOR_ELEVEN_STAGE_SEAL_CHECKS_FINAL.csv',old);add('ALL_PREVIOUS_UNCHANGED',all(r['passed'] for r in old),len(old))
 access=json.loads((RESULTS/'PROTOCOL_ACCESS_RECORD.json').read_text(encoding='utf-8'));new=json.loads((RESULTS/'FIRST_NEW_EXTENSION_TARGET_ACCESS.json').read_text(encoding='utf-8'));cal=json.loads((RESULTS/'CALIBRATION_RUN_SUMMARY.json').read_text(encoding='utf-8'));freeze=json.loads((RESULTS/'CALIBRATION_FREEZE_BEFORE_RUN.json').read_text(encoding='utf-8'))
 add('PROTOCOL_UNCHANGED',access['protocol_sha256']==new['protocol_sha256']==sha(ROOT/'STAGE12_PROTOCOL_v1.0.md'));add('BOUNDARY_RECORDED',new['inference_boundary_sha256']==sha(ROOT/'EXECUTION_AND_INFERENCE_BOUNDARY_v1.1.md'));add('CALIBRATION_BEFORE_NEW_TARGET',datetime.fromisoformat(cal['time'])<=datetime.fromisoformat(new['time']))
 for n,h in freeze['code_sha256'].items():add('FROZEN_ENGINE/'+n,sha(ROOT/'code'/n)==h)
 csvout(RESULTS/'SCIENTIFIC_QA_CHECKS.csv',checks);out={'time':now(),'checks':len(checks),'passed':sum(r['passed'] for r in checks),'failed':sum(not r['passed'] for r in checks),'real_tests_independently_refit':len(T),'saved_P_recount_rows':len(recounts),'seed_replay_P_values':len(replay),'independent_noise_draws':nnoise,'prior_sealed_entries':len(old),'prior_unique_paths':len(cache),'post_selection_confirmatory_FDR_claimed':False,'causal_mechanism_identified':False};save_json(RESULTS/'SCIENTIFIC_QA_SUMMARY.json',out)
 assert out['failed']==0,[r for r in checks if not r['passed']][:15]
 log(f'Scientific QA {len(checks)}/{len(checks)} passed.')
if __name__=='__main__':assert_unsealed();run()

