"""Independent algebra, reference multiplicity, and retained-replicate recounts."""
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import numpy as np,pandas as pd,json,hashlib,pickle
from scipy.linalg import lstsq
from statsmodels.stats.multitest import multipletests
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results';D=ROOT/'00_restricted';checks=[]
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def check(name,ok,detail=''):checks.append({'check':name,'passed':bool(ok),'detail':str(detail)})
sap=json.loads((R/'MAIN_ANALYSIS_SAP_FROZEN.json').read_text(encoding='utf-8'));v=sap['version'];family=sap['primary_family_accessions']
first=json.loads((R/'FIRST_PROTEIN_VALUE_ACCESS.json').read_text(encoding='utf-8'))
check('SAP_linked_before_protein_access',first['SAP_sha256_before_first_protein_access']==sha(R/'MAIN_ANALYSIS_SAP_FROZEN.json'))
check('protocol_unchanged',sap['protocol_sha256']==sha(ROOT/'STAGE5_PROTOCOL_v1.0.md'))
for name,digest in json.loads((R/'ENGINE_FREEZE_BEFORE_SIMULATION.json').read_text(encoding='utf-8'))['code'].items():check('frozen_engine_'+name,sha(ROOT/'code'/name)==digest)
for name in ['model.py','diagnostics.py','load_proteins.py','analyze_main.py']:check('frozen_analysis_'+name,sha(ROOT/'code'/name)==sap['code_sha256'][name])
sources=pd.read_csv(R/'SOURCE_HASH_CHECKS.csv')
for row in sources.itertuples():check('source_unchanged_'+row.id,sha(row.path)==row.expected_sha256)
check('protein_reconstruction_all490',pd.read_csv(R/'PROTEIN_RECONSTRUCTION_QA.csv').passed.all())
ds=pickle.load((D/'DESIGNS_MAIN_RESTRICTED.pkl').open('rb'));cache=np.load(D/'PROTEIN_TRANSFORMED_RESTRICTED.npz')
meta=pd.read_csv(ROOT.parent/'1候选逐蛋白可评价性审计/results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
main=pd.read_csv(R/'PRIMARY_ALL43_RESULTS.csv').set_index('accession');numeric=[]
for acc,d in ds.items():
 j=list(cache['accessions']).index(acc);k=list(cache['soma_features']).index(meta.loc[acc,'soma_features'])
 Y=np.c_[cache['N'][:,j],cache['S'][:,k],cache['T_GIS'][:,j]][d['indices']]
 X=np.c_[d['Z'],d['P']];coef,_,rank,_=lstsq(X,Y);E=Y-X@coef;b=coef[-1];df=len(Y)-rank
 cov=E.T@E/df;c12,c13,c23=cov[0,1],cov[0,2],cov[1,2]
 h=np.array([b[1]*c13-b[0]*c23,b[2]*c12-b[0]*c23])
 pcoef,_,_,_=lstsq(d['Z'],d['P']);rP=d['P']-d['Z']@pcoef;a=rP/(rP@rP)
 prod=np.c_[E[:,0]*E[:,1],E[:,0]*E[:,2],E[:,1]*E[:,2]];ci=(prod-prod.mean(axis=0))/df;bi=a[:,None]*E
 IF=np.c_[-c23*bi[:,0]+c13*bi[:,1]+b[1]*ci[:,1]-b[0]*ci[:,2],-c23*bi[:,0]+c12*bi[:,2]+b[2]*ci[:,0]-b[0]*ci[:,2]]
 V=np.zeros((2,2))
 for axis in range(3):
  labels=d['groups'][:,axis];s=np.vstack([IF[labels==g].sum(axis=0) for g in sorted(set(labels))]);s-=s.mean(axis=0);V+=s.T@s*len(s)/(len(s)-1)
 factor=[1.,(len(Y)-1)/df,2*(len(Y)-1)/df][v];V*=factor;Q=float(h@np.linalg.solve(V,h));p=float(np.exp(-Q/2));r=main.loc[acc]
 ok=np.allclose(Q,r.Q,rtol=2e-7,atol=1e-9) and np.allclose(p,r.p_value,rtol=2e-7,atol=1e-12)
 check('independent_full_OLS_H_Wald_'+acc,ok)
 numeric.append({'accession':acc,'Q_independent':Q,'Q_main':r.Q,'p_independent':p,'p_main':r.p_value,'passed':ok})
pd.DataFrame(numeric).to_csv(R/'INDEPENDENT_MAIN_NUMERIC_VALIDATION.csv',index=False,encoding='utf-8-sig')
check('43_candidates_retained',len(main)==43)
check('34_primary_no_MAPT',len(family)==34 and 'P10636' not in family and main.loc[family].role.eq('PRIMARY_DISCOVERY_CANDIDATE').all())
for name,frame in [('PRIMARY',main.reset_index()),('VARIANTS',pd.read_csv(R/'ALL_FIXED_VARIANT_RESULTS.csv'))]:
 groups=[(name,frame)] if name=='PRIMARY' else list(frame.groupby('variant'))
 for label,part in groups:
  part=part.loc[part.in_frozen_family];p=part.p_value.fillna(1).to_numpy()
  for column,method in [('q_BY','fdr_by'),('q_BH','fdr_bh'),('p_Holm','holm')]:check('reference_'+str(label)+'_'+column,np.allclose(part[column],multipletests(p,method=method)[1],atol=1e-10))
profiles=pd.read_csv(R/'COVARIANCE_OFFSET_GRID125.csv')
for acc,part in profiles.groupby('accession'):
 zero=part.loc[part.rho_NS.eq(0)&part.rho_NT.eq(0)&part.rho_ST.eq(0)].iloc[0]
 check('rho_grid125_zero_anchor_'+acc,len(part)==125 and np.allclose(zero.p,main.loc[acc,'p_value'],rtol=1e-8))
check('GIS_FE_same_subset_equivalence',pd.read_csv(R/'GIS_FIXED_EFFECT_EQUIVALENCE_QA.csv').passed.all())
# Recount every retained replicate file once and verify all summary rates/failure counts.
repchecks=[]
for table in ['DEV_SIMULATION_SUMMARY.csv','CONFIRM_SIMULATION_SUMMARY.csv','VARIANT_CALIBRATION_SUMMARY.csv','POWER_SIMULATION_SUMMARY.csv']:
 tableframe=pd.read_csv(R/table)
 for filename,part in tableframe.groupby('file'):
  q=pd.read_parquet(R/'simulation_replicates'/filename).Q0.to_numpy();ident=part.design_id.iloc[0]
  if ident in ds:des=ds[ident]
  else:
   if 'extra' not in locals():extra=pickle.load((D/'DESIGNS_VARIANTS_RESTRICTED.pkl').open('rb'))
   des=extra[ident]
  factors=[1.,(des['n']-1)/des['df'],2*(des['n']-1)/des['df']];allok=True
  for row in part.itertuples():
   p=np.exp(-q/(2*factors[row.version]));ok=np.isfinite(p);count=int(np.sum(ok&(p<=row.alpha)))
   allok&=len(q)==row.reps and int((~ok).sum())==row.n_failed and count==row.n_rejected and abs(count/len(q)-row.rate)<1e-12
  repchecks.append({'summary_table':table,'file':filename,'reps':len(q),'passed':bool(allok)})
stress=pd.read_csv(R/'IDENTIFICATION_STRESS_RESULTS.csv')
for row in stress.itertuples():
 filename=f'STRESS__{row.accession}__{row.scenario}.parquet';raw=pd.read_parquet(R/'simulation_replicates'/filename);p=raw.p_value.to_numpy();valid=np.isfinite(p)
 ok=len(raw)==row.reps and int((~valid).sum())==row.n_failed and int((valid&(p<=.05)).sum())==row.n_rejected05
 if v==0:ok=ok and np.allclose(p,np.exp(-raw.Q0.to_numpy()/2),equal_nan=True)
 repchecks.append({'summary_table':'IDENTIFICATION_STRESS_RESULTS.csv','file':filename,'reps':len(raw),'passed':bool(ok)})
check('all_simulation_recount_files',all(r['passed'] for r in repchecks),len(repchecks))
pd.DataFrame(repchecks).to_csv(R/'INDEPENDENT_REPLICATE_RECOUNTS.csv',index=False,encoding='utf-8-sig')
for table in ['INDEPENDENT_CONFIRMATION_DECISIONS.csv','VARIANT_CALIBRATION_DECISIONS.csv']:
 frame=pd.read_csv(R/table);flag='passed_independent_confirmation' if table.startswith('INDEPENDENT') else 'calibration_passed'
 passed=frame.loc[frame[flag]];check('calibration_threshold_'+table,(passed.max_null_wilson95_upper<=.06).all() and (passed.max_failure_rate<=.005).all())
power=pd.read_csv(R/'PER_PROTEIN_CONDITIONAL_POWER.csv')
check('power_all_35_not_observed_effect',power.accession.nunique()==35 and power.observed_effect_not_used.all())
check('power_6_effects_2targets_3thresholds',len(power)==35*6*2*3 and power.n_replicates.eq(2000).all())
check('power_probability_bounds',power.power.between(0,1).all() and (power.wilson95_low<=power.power).all() and (power.wilson95_high>=power.power).all())
for path in R.glob('*.csv'):
 cols=pd.read_csv(path,nrows=0).columns
 check('public_no_person_ID_'+path.stem,not any(c in cols for c in ['RID','PTID','key','EXAMDATE','GENOTYPE','PTDOB','PTDOBYY']))
out=pd.DataFrame(checks);out.to_csv(R/'INDEPENDENT_QA_CHECKS.csv',index=False,encoding='utf-8-sig')
summary={'checks':len(out),'passed':int(out.passed.sum()),'failed':int((~out.passed).sum()),'simulation_files_recounted':len(repchecks),
 'recounted_independent_replicates':sum(r['reps'] for r in repchecks),'actual_main_genes_independently_refit':len(ds),'no_primary_family_changes':True}
(R/'INDEPENDENT_QA_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(summary)
if not out.passed.all():print(out.loc[~out.passed].to_string(index=False));raise SystemExit(1)
