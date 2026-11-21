from config import *
from scipy.linalg import lstsq
from scipy.stats import norm,chi2
from run_calibration import innovations,patterns
from datetime import datetime
import pickle

def independent_wilson(k,n):
 z=1.959963984540054;p=k/n;den=1+z*z/n
 return (p+z*z/(2*n)+z*np.sqrt(p*(1-p)/n+z*z/(4*n*n)))/den

def run():
 checks=[]
 def add(name,ok,detail=''):checks.append({'check':name,'passed':bool(ok),'detail':str(detail)})
 templates=pickle.load((PRIVATE/'ANALYSIS_TEMPLATES_RESTRICTED.pkl').open('rb'));scales=pickle.load((PRIVATE/'EXPOSURE_SCALING_RESTRICTED.pkl').open('rb'))
 clinical=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv');real=pd.read_csv(RESULTS/'ALL43_DESIGN_CLINICAL_RESULTS.csv').set_index('design_id');effects=pd.read_csv(RESULTS/'PLATFORM_CLINICAL_EFFECTS_ALL43.csv')
 original=pd.read_csv(RESULTS/'ENGINE_INDEPENDENT_CHECKS.csv')
 for r in original.itertuples():add('ENGINE/'+r.check,r.passed,r.max_abs_error)
 D=np.array([[-1.,1.,0.],[-1.,0.,1.]])
 for key,t in templates.items():
  y=t['y'];bb=[];scores=[]
  for j in range(3):
   T=np.c_[t['ZZ'],t['targetx'][:,j]] if t['modifier'] is None else np.c_[t['ZZ'],t['targetx'][:,j],t['targetx'][:,j]*t['modifier'][t['case_local']]]
   fit=lstsq(T,y,cond=max(T.shape)*np.finfo(float).eps,lapack_driver='gelsy')[0];res=y-T@fit;b=fit[-1];bb.append(b);ss=[]
   for c in ['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex']:
    source=clinical.loc[scales[t['accession']]['ids'],c].astype(str);levels=sorted(source.unique());lookup={v:i for i,v in enumerate(levels)};labels=np.array([lookup[str(clinical.loc[i,c])] for i in t['people']])
    sg=np.bincount(labels,weights=t['A'][j]*res,minlength=len(levels));sg-=sg.mean();sg*=np.sqrt(len(levels)/(len(levels)-1));ss.append(sg)
   scores.append(np.concatenate(ss));ref=effects.loc[effects.design_id.eq(key)&effects.assay.eq(['N','S','T'][j])].iloc[0]
   add('REAL_BETA/'+key+str(j),np.isclose(b,ref.b_original_scale,rtol=1e-7,atol=1e-10))
  b=np.array(bb);S=np.array(scores).T;V=S.T@S;gamma=b*t['sigma'];scaledscore=S*t['sigma'][None,:]+b[None,:]*t['Uv']/(2*t['sigma'][None,:]);Vs=scaledscore.T@scaledscore;z=b/np.sqrt(np.diag(V))
  p_same=min(1.,2*min(norm.sf(z).max(),norm.cdf(z).max()));d=D@gamma;Vd=D@Vs@D.T;p_diff=chi2.sf(d@np.linalg.solve(Vd,d),2);p_any=chi2.sf(b@np.linalg.solve(V,b),3)
  for name,p in [('same',p_same),('diff',p_diff),('any',p_any)]:add('REAL_P/'+key+name,np.isclose(p,real.loc[key,'p_'+name+'_computed'],rtol=1e-6,atol=1e-10))
  for j in range(3):
   rr=effects.loc[effects.design_id.eq(key)&effects.assay.eq(['N','S','T'][j])].iloc[0]
   add('REAL_SCALED_SE/'+key+str(j),np.isclose(np.sqrt(Vs[j,j]),rr.standardized_SE_full_joint,rtol=1e-7,atol=1e-10))
 formal=pd.read_csv(RESULTS/'FORMAL32_CLINICAL_RESULTS.csv');decision=pd.read_csv(RESULTS/'CALIBRATION_QUALIFICATION_DECISIONS.csv').set_index('design_id')
 add('FORMAL32_FIXED',len(formal)==32);add('PRIMARY6_FIXED',int(formal.family.eq('PRIMARY6').sum())==6);add('SECONDARY26_FIXED',int(formal.family.eq('SECONDARY26').sum())==26)
 for r in formal.itertuples():
  qual=bool(decision.loc[r.design_id,'qualified']);add('QUALIFICATION_P1/'+r.hypothesis_id,qual or r.p_for_family==1.)
 for family,part in formal.groupby('family'):
  order=part.p_for_family.sort_values();m=len(part);harm=sum(1/i for i in range(1,m+1));vals=np.minimum.accumulate((order.to_numpy()*m*harm/np.arange(1,m+1))[::-1])[::-1]
  for ix,q in zip(order.index,np.minimum(1,vals)):add('BY_FAMILY/'+formal.loc[ix,'hypothesis_id'],np.isclose(q,formal.loc[ix,'q_BY_within_family']))
 order=formal.p_for_family.sort_values();vals=np.minimum.accumulate((order.to_numpy()*32*sum(1/i for i in range(1,33))/np.arange(1,33))[::-1])[::-1]
 for ix,q in zip(order.index,np.minimum(1,vals)):add('BY32/'+formal.loc[ix,'hypothesis_id'],np.isclose(q,formal.loc[ix,'q_BY_combined32']))
 # Full independent recount of every saved test P value and its Wilson upper.
 inventory=pd.read_csv(RESULTS/'SIMULATION_FILE_INVENTORY.csv',dtype={'seed':str});recounts=[];replay=[]
 totalnoise=0
 for count,r in enumerate(inventory.itertuples(),1):
  file=SIM/r.file;df=pd.read_parquet(file);meta=json.loads((SIM/(r.design_id+'__'+r.scenario+'.json')).read_text(encoding='utf-8'));totalnoise+=len(df)
  add('SIM_HASH/'+r.design_id+r.scenario,sha(file)==r.file_sha256)
  add('SIM_REPS/'+r.design_id+r.scenario,len(df)==r.reps and np.array_equal(df.replicate,np.arange(r.reps)))
  allcounts=True;allintervals=True
  for row in meta['summary']:
   values=df[row['test']+'_'+row['pattern']].to_numpy();k=int(np.count_nonzero(values<=row['alpha']));upper=independent_wilson(k,len(values));passcount=k==row['rejected'];passinterval=np.isclose(upper,row['Wilson95_high'],rtol=1e-12,atol=1e-12);allcounts&=passcount;allintervals&=passinterval
   recounts.append({'design_id':r.design_id,'scenario':r.scenario,'pattern':row['pattern'],'test':row['test'],'alpha_name':row['alpha_name'],'rejected_recounted':k,'rejected_recorded':row['rejected'],'Wilson95_upper_recounted':upper,'passed':passcount and passinterval})
  add('SIM_ALL_COUNTS/'+r.design_id+r.scenario,allcounts);add('SIM_ALL_WILSON/'+r.design_id+r.scenario,allintervals)
  # Replay first RNG block exactly; independently calculate P from scores for three draws and four patterns.
  t=templates[r.design_id];rng=np.random.default_rng(int(r.seed));Y=innovations(rng,t,r.scenario,min(BATCH,r.reps))[:,:3];ps=patterns(t)
  names=['ZERO','COMMON_0.35_1','ONE_ZERO_0_1',next(p['name'] for p in ps if p['name'].startswith('POWER_BALANCED') and p['name'].endswith('_1'))]
  As=t['sigma'][:,None]*t['A'];Gram=As@As.T
  for p in ps:
   if p['name'] not in names:continue
   mu=As.T@np.linalg.solve(Gram,p['effect']);YY=Y+mu[:,None]
   for j in range(3):
    b=t['A']@YY[:,j];score=np.einsum('gkn,n->gk',t['K'],YY[:,j],optimize=True);V=score.T@score;z=b/np.sqrt(np.diag(V));p_same=min(1.,2*min(norm.sf(z).max(),norm.cdf(z).max()))
    gamma=b*t['sigma'];ss=score*t['sigma'][None,:]+b[None,:]*t['Uv']/(2*t['sigma'][None,:]);Vs=ss.T@ss;d=D@gamma;Vd=D@Vs@D.T;p_diff=chi2.sf(d@np.linalg.solve(Vd,d),2);p_any=chi2.sf(b@np.linalg.solve(V,b),3)
    for typ,value in [('SAME',p_same),('DIFF',p_diff),('ANY',p_any)]:
     expected=df.loc[j,typ+'_'+p['name']];ok=np.isclose(value,expected,rtol=1e-6,atol=1e-10)
     replay.append({'design_id':r.design_id,'scenario':r.scenario,'pattern':p['name'],'replicate':j,'test':typ,'P_replayed':value,'P_saved':expected,'passed':ok})
  if count%40==0:log(f'Independent simulation recount/replay {count}/{len(inventory)}')
 csvout(RESULTS/'INDEPENDENT_SIMULATION_RECOUNTS.csv',recounts);csvout(RESULTS/'INDEPENDENT_SEED_REPLAY_CHECKS.csv',replay)
 add('ALL_SEED_REPLAY_P_VALUES',all(r['passed'] for r in replay),len(replay));add('TOTAL_NOISE_DRAWS',totalnoise==2410000,totalnoise)
 # Qualification is recomputed from full saved test recounts; component zero counts use saved aggregate plus seed checks.
 full=pd.read_csv(RESULTS/'CALIBRATION_AND_POWER_FULL.csv');comp=pd.read_csv(RESULTS/'COMPONENT_NULL_CALIBRATION.csv')
 for key,t in templates.items():
  rel=['ANY','DIFF'] if 'INTERACTION' in t['branch'] else ['SAME','DIFF'];part=full.loc[full.design_id.eq(key)&full.true_null&full.test.isin(rel)];nom=part.loc[part.alpha_name.eq('NOMINAL')];tail=part.loc[~part.alpha_name.eq('NOMINAL')];cc=comp.loc[comp.design_id.eq(key)]
  inv=inventory.loc[inventory.design_id.eq(key)];qualified=bool(t['eligible'] and nom.Wilson95_high.le(.06).all() and tail.Wilson95_high.le(tail.alpha).all() and (cc.Wilson95_high.le(.035).all() if 'SAME' in rel else True) and (inv.n_failed/inv.reps).max()<=.005)
  add('QUALIFICATION_RECOUNT/'+key,qualified==decision.loc[key,'qualified'])
 # Exact data selection, scale invariance and auxiliary identities.
 cov=pd.read_csv(RESULTS/'CLINICAL_COVERAGE_RECOUNT.csv').set_index(['marker','horizon_years'])
 add('STAGE4_MEM_BASE568',cov.loc[('ADNI_MEM',5),'baseline']==568);add('STAGE4_MEM_HIGH420',cov.loc[('ADNI_MEM',5),'highinfo']==420);add('STAGE4_EF_HIGH417',cov.loc[('ADNI_EF',5),'highinfo']==417)
 for r in pd.read_csv(RESULTS/'PET_SOURCE_REPRODUCTION_CHECKS.csv').itertuples():add('PET_SOURCE/'+str(r.index_number),r.same_stage5_CL and r.same_tracer)
 for acc,t0 in [(a,templates['MEM5Y__'+a]) for a in CANDIDATES]:
  for branch in SENSITIVITIES:
   add('FIXED_EXPOSURE_UNIT/'+acc+branch,np.array_equal(t0['sigma'],templates[branch+'__'+acc]['sigma']))
 aux=pd.read_csv(RESULTS/'PTAU_CONDITIONING_ALGEBRAIC_DECOMPOSITION.csv')
 for r in aux.itertuples():add('PTAU_ALGEBRA/'+r.gene+r.assay,r.identity_abs_error<1e-8)
 source=pd.read_csv(RESULTS/'INPUT_SOURCE_MANIFEST.csv')
 for r in source.itertuples():add('INPUT_UNCHANGED/'+r.path,sha(r.path)==r.sha256)
 access=json.loads((RESULTS/'ACCESS_AND_PROTOCOL_RECORD.json').read_text(encoding='utf-8'));first=json.loads((RESULTS/'FIRST_NEW_CLINICAL_VALUE_ACCESS.json').read_text(encoding='utf-8'));realf=json.loads((RESULTS/'FIRST_REAL_TARGET_ASSOCIATION_ACCESS.json').read_text(encoding='utf-8'));cal=json.loads((RESULTS/'CALIBRATION_RUN_SUMMARY.json').read_text(encoding='utf-8'))
 add('PROTOCOL_UNCHANGED',access['protocol_sha256']==sha(ROOT/'STAGE10_PROTOCOL_v1.0.md')==first['protocol_sha256'])
 add('CALIBRATION_PRECEDES_REAL_TARGET',datetime.fromisoformat(cal['time'])<datetime.fromisoformat(realf['time']))
 freeze=json.loads((RESULTS/'CALIBRATION_FREEZE_BEFORE_RUN.json').read_text(encoding='utf-8'))
 for name,expected in freeze['code_sha256'].items():add('FROZEN_CALIBRATION_ENGINE/'+name,sha(ROOT/'code'/name)==expected)
 prior=pd.read_csv(S9/'results/PRIOR_ALL_STAGE_SEAL_CHECKS_FINAL.csv')[['path','expected_sha256']].to_dict('records')
 for r in pd.read_csv(S9/'SHA256_MANIFEST.csv').itertuples():prior.append({'path':str(S9/r.relative_path),'expected_sha256':r.sha256})
 cache={};old=[]
 for r in prior:
  p=Path(r['path']);key=str(p.resolve())
  if key not in cache:cache[key]=sha(p) if p.exists() else ''
  old.append({**r,'actual_sha256':cache[key],'passed':cache[key]==r['expected_sha256']})
 csvout(RESULTS/'PRIOR_ALL_NINE_STAGES_SEAL_CHECKS_FINAL.csv',old);add('PREVIOUS_ALL_SEALED_UNCHANGED',all(r['passed'] for r in old),f'{len(old)} entries/{len(cache)} unique paths')
 csvout(RESULTS/'SCIENTIFIC_QA_CHECKS.csv',checks)
 result={'analysis_date':'2026-10-29','time':now(),'checks':len(checks),'passed':sum(r['passed'] for r in checks),'failed':sum(not r['passed'] for r in checks),'simulation_P_recount_rows':len(recounts),'seed_replay_P_values':len(replay),'independent_noise_draws':totalnoise,'component_null_full_recount_scope':'Aggregate counts/Wilson plus deterministic initial-block replay; full marginal P arrays not saved','prior_sealed_entries':len(old),'prior_unique_paths':len(cache),'real_designs_independently_refit':len(templates),'preexisting_scientific_results_modified':False}
 save_json(RESULTS/'SCIENTIFIC_QA_SUMMARY.json',result)
 assert all(r['passed'] for r in checks),[r for r in checks if not r['passed']][:10]
 log(f'Scientific QA:{len(checks)}/{len(checks)} passed;{len(recounts)} simulation-count rows and{len(replay)} replay P values checked')
if __name__=='__main__':
 assert_unsealed();run()
