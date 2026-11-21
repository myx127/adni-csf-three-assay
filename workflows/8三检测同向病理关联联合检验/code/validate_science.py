from config import *
from bootstrap import verify_prior
from freeze_engine import verify
from scipy.linalg import lstsq
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from concurrent.futures import ProcessPoolExecutor,as_completed
import pickle,math

def interval(k,n):
 z=1.959963984540054;p=k/n;den=1+z*z/n;mid=(p+z*z/(2*n))/den;half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 return max(0.,mid-half),min(1.,mid+half)

def recount(task):
 file,records,components,expected_sha=task;r=pd.read_parquet(SIM/file)
 ok=len(r)==REPS and np.array_equal(r.replicate,np.arange(REPS)) and sha(SIM/file)==expected_sha
 noise=r[[x+'_noise_over_templateSD' for x in ['N','S','T']]].to_numpy();se=r[[x+'_SE_over_templateSD' for x in ['N','S','T']]].to_numpy();valid=r.valid.to_numpy(bool)
 for pattern in PATTERNS:
  z=(noise+np.array(pattern['effect']))/se;low=z.min(axis=1);high=z.max(axis=1);strength=np.maximum(low,-high)
  p=np.minimum(1.,2*norm.sf(strength));p[~valid]=np.nan;direction=np.where(low>0,1,-1)
  for row in records[pattern['pattern']]:
   hit=valid&(p<=row['alpha']);count=int(hit.sum());lo,hi=interval(count,REPS)
   ok&=count==row['success_count'] and np.isclose(lo,row['Wilson95_low'],atol=1e-13) and np.isclose(hi,row['Wilson95_high'],atol=1e-13)
   ok&=int((hit&(direction==1)).sum())==row['positive_rejection_count'] and int((hit&(direction==-1)).sum())==row['negative_rejection_count']
   if pattern['kind']=='POWER':ok&=int((hit&(direction!=pattern['sign'])).sum())==row['wrong_direction_count']
   ok&=row['failure_count']==int((~valid).sum())
 for row in components:
  j=['N','S','T'].index(row['assay']);p=norm.sf(row['direction']*noise[:,j]/se[:,j]);count=int((valid&(p<=.025)).sum());lo,hi=interval(count,REPS)
  ok&=count==row['reject_count'] and np.isclose(hi,row['Wilson95_high'])
 return {'file':file,'reps':len(r),'pattern_threshold_evaluations':len(PATTERNS)*len(ALPHAS),'passed':bool(ok)}

def run():
 verify();checks=[]
 def check(name,ok,detail=''):checks.append({'check':name,'passed':bool(ok),'detail':str(detail)})
 plan=json.loads((RESULTS/'STAGE8_PLAN_FROZEN.json').read_text(encoding='utf-8'))
 check('original_protocol_hash',sha(ROOT/'STAGE8_PROTOCOL_v1.0.md')==plan['protocol_sha256'])
 first=json.loads((RESULTS/'FIRST_NEW_SHARED_STATISTIC_ACCESS.json').read_text(encoding='utf-8'))
 check('new_statistics_after_frozen_complete_calibration',first['after_complete_frozen_calibration'] and first['calibration_decisions_sha256']==sha(RESULTS/'CALIBRATION_DECISIONS.csv'))
 engine=pd.read_csv(RESULTS/'ENGINE_INDEPENDENT_CHECKS.csv');check('1651_engine_checks',len(engine)==1651 and engine.passed.all())
 for row in pd.read_csv(RESULTS/'REFERENCE_CODE_ORIGIN.csv').itertuples():check('reference_code_'+row.file,sha(ROOT/'code'/row.file)==row.sha256==sha(row.source))
 for row in pd.read_csv(RESULTS/'SOURCE_MANIFEST.csv').itertuples():check('source_'+row.source_id,sha(row.path)==row.sha256)
 templates=pickle.load((PRIVATE/'SIMULATION_TEMPLATES_RESTRICTED.pkl').open('rb'));memo=pickle.load((PRIVATE/'REAL_MOMENTS_RESTRICTED.pkl').open('rb'))
 cache=np.load(S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz');meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
 result=pd.read_csv(RESULTS/'ALL35_THREE_MODEL_SHARED_RESULTS.csv').set_index('design_id')
 for key,t in templates.items():
  d=t['d'];acc=t['accession'];j=list(cache['accessions']).index(acc);k=list(cache['soma_features']).index(meta.loc[acc,'soma_features'])
  Y=np.c_[cache['N'][:,j],cache['S'][:,k],cache['T_GIS'][:,j]][d['indices']];X=np.c_[d['Z'],d['P']]
  coef,_,rank,_=lstsq(X,Y,cond=max(X.shape)*np.finfo(float).eps);E=Y-X@coef;rp=d['P']-d['Z']@lstsq(d['Z'],d['P'],cond=max(d['Z'].shape)*np.finfo(float).eps)[0]
  a=rp/(rp@rp);df=len(Y)-rank;prod=np.c_[E[:,0]**2,E[:,1]**2,E[:,2]**2,E[:,0]*E[:,1],E[:,0]*E[:,2],E[:,1]*E[:,2]]
  theta=np.r_[coef[-1],prod.sum(axis=0)/df];IF=np.c_[a[:,None]*E,(prod-prod.mean(axis=0))/df];V=np.zeros((9,9))
  for axis in range(3):
   labels=d['groups'][:,axis];s=np.vstack([IF[labels==g].sum(axis=0) for g in np.unique(labels)]);s-=s.mean(axis=0);V+=s.T@s*len(s)/(len(s)-1)
  err=max(np.max(abs(theta-memo[key]['theta'])),np.max(abs(V-memo[key]['V'])))
  check('independent_SVD_theta_SUM3_'+key,err<1e-8,err)
  z=theta[:3]/np.sqrt(V.diagonal()[:3]);strength=max(z.min(),-z.max());p=min(1.,2*norm.sf(strength))
  check('independent_joint_formula_'+key,np.isclose(p,result.loc[key,'p_shared_computed'],rtol=1e-11,atol=1e-14))
  grad=np.zeros((3,9))
  for j in range(3):grad[j,j]=1/np.sqrt(theta[3+j]);grad[j,3+j]=-.5*theta[j]/theta[3+j]**1.5
  for j in range(9):
   th=theta.astype(complex);th[j]+=1e-20j;numeric=np.imag(th[:3]/np.sqrt(th[3:6]))/1e-20
   check('scaled_effect_complex_gradient_'+key+'_'+str(j),np.allclose(numeric,grad[:,j],atol=1e-10,rtol=1e-9))
 check('105_old_coefficients_reproduced',pd.read_csv(RESULTS/'PRIOR_PARAMETER_REPRODUCTION.csv').passed.all())
 for model in MODELS:
  part=result.loc[result.model.eq(model)&result.in_family34];check('family34_'+model,len(part)==34)
  ids=part.index;p=part.p_shared_for_BY.to_numpy()
  for method,col in [('fdr_by','q_BY_shared34'),('fdr_bh','q_BH_shared34'),('holm','p_Holm_shared34')]:
   check('statsmodels_'+model+'_'+method,np.allclose(multipletests(p,method=method)[1],part[col],atol=1e-13))
  combined=multipletests(np.r_[p,part.prior_H_p],method='fdr_by')[1]
  check('combined68_'+model,np.allclose(combined[:34],part.q_BY_combined68_shared) and np.allclose(combined[34:],part.q_BY_combined68_H))
  dirs=multipletests(np.r_[part.p_plus_for_directional_BY,part.p_minus_for_directional_BY],method='fdr_by')[1]
  check('directional68_'+model,np.allclose(dirs[:34],part.q_BY_direction_plus68) and np.allclose(dirs[34:],part.q_BY_direction_minus68))
  qualified=part.calibration_qualified&part.math_valid;check('failed_gates_never_shrink_family_'+model,part.loc[~qualified,'p_shared_for_BY'].eq(1).all())
 check('no_factor_exclusion',len(result)==105 and result.loc[result.model.eq('M0')&result.in_family34].prior_M0_factor_structure_compatible.sum()==15)
 check('MAPT_separate',result.loc[~result.in_family34,'q_BY_shared34'].isna().all() and len(result.loc[~result.in_family34])==3)
 check('supported_shared_all_observed_signs_agree',result.loc[result.shared_BY34_supported,'all_three_observed_signs_agree'].all())
 sensitivity=pd.read_csv(RESULTS/'APOE_SENSITIVITY_COMPARISONS.csv')
 check('sensitivity_decompositions',len(sensitivity)==68 and sensitivity.q_identity_error.max()<1e-12 and sensitivity.P_identity_error.max()<1e-12)
 null=pd.read_csv(RESULTS/'NULL_CALIBRATION_SUMMARY.csv',dtype={'seed':str});power=pd.read_csv(RESULTS/'POWER_SCENARIO_SUMMARY.csv',dtype={'seed':str});allsummary=pd.concat([null,power],ignore_index=True)
 components=pd.read_csv(RESULTS/'COMPONENT_NULL_CALIBRATION.csv',dtype={'seed':str});inventory=pd.read_csv(RESULTS/'SIMULATION_FILE_INVENTORY.csv',dtype={'seed':str})
 check('525_unique_noise_files_seeds',len(inventory)==525 and inventory.file.is_unique and inventory.seed.is_unique)
 tasks=[]
 for file,group in allsummary.groupby('file'):
  records={name:g.to_dict('records') for name,g in group.groupby('pattern')}
  inv=inventory.loc[inventory.file.eq(file)].iloc[0];tasks.append((file,records,components.loc[components.file.eq(file)].to_dict('records'),inv.sha256))
 recounted=[]
 with ProcessPoolExecutor(max_workers=WORKERS) as pool:
  for i,future in enumerate(as_completed([pool.submit(recount,t) for t in tasks]),1):
   r=future.result();recounted.append(r);check('independent_recount_'+r['file'],r['passed'])
   if i%35==0 or i==len(tasks):log(f'Independent full-grid recount {i}/{len(tasks)}')
 csvout(RESULTS/'INDEPENDENT_SIMULATION_RECOUNTS.csv',recounted)
 check('21000000_unique_repetitions',sum(r['reps'] for r in recounted)==21000000)
 decision=pd.read_csv(RESULTS/'CALIBRATION_DECISIONS.csv')
 for row in decision.itertuples():
  n=null.loc[null.design_id.eq(row.design_id)];cp=components.loc[components.design_id.eq(row.design_id)];nom=n.loc[n.threshold.eq('NOMINAL')];tail=n.loc[~n.threshold.eq('NOMINAL')]
  qualified=nom.Wilson95_high.le(.06).all() and cp.Wilson95_high.le(.03).all() and tail.Wilson95_high.le(tail.alpha).all() and (n.failure_count/n.reps).le(.005).all()
  check('independent_gate_'+row.design_id,bool(qualified)==row.qualified)
 for p in RESULTS.glob('*.csv'):
  if p.name in ['INDEPENDENT_QA_CHECKS.csv','DELIVERY_QA_CHECKS.csv','WORKBOOK_SHEET_INDEX.csv']:continue
  check('aggregate_columns_'+p.stem,not set(pd.read_csv(p,nrows=0).columns)&{'RID','PTID','key','EXAMDATE','source_row','GENOTYPE','PTDOB'})
 prior=verify_prior();csvout(RESULTS/'PRIOR_SEALED_FILE_CHECKS_FINAL.csv',prior);check('all4323_prior_entries_unchanged',len(prior)==4323 and all(r['passed'] for r in prior))
 out=pd.DataFrame(checks);csvout(RESULTS/'INDEPENDENT_QA_CHECKS.csv',out)
 summary={'checks':len(out),'passed':int(out.passed.sum()),'failed':int((~out.passed).sum()),'independent_noise_repetitions_recounted':sum(r['reps'] for r in recounted),'noise_cells':len(recounted),'real_designs_independently_refit':105,'previous_sealed_entries_verified':len(prior)}
 save_json(RESULTS/'INDEPENDENT_QA_SUMMARY.json',summary);print(summary)
 if not out.passed.all():print(out.loc[~out.passed].to_string(index=False));raise SystemExit(1)
if __name__=='__main__':run()
