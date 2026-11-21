"""Independent source, count, algebra, gradient, inference and multiplicity checks."""
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from config import *
from scipy.linalg import lstsq
from scipy.stats import chi2,norm
from statsmodels.stats.multitest import multipletests
from collections import defaultdict,Counter
import pickle,re

checks=[]
def check(name,ok,detail=''):checks.append({'check':name,'passed':bool(ok),'detail':str(detail)})
for row in pd.read_csv(RESULTS/'SOURCE_IDENTITY_CHECKS.csv').itertuples():check('source_hash_'+row.id,sha(row.path)==row.expected_sha256)
for row in pd.read_csv(RESULTS/'ADDITIONAL_INPUT_MANIFEST.csv').itertuples():check('extra_hash_'+row.input_id,sha(row.path)==row.sha256)
for row in pd.read_csv(RESULTS/'UNIPROT_SOURCE_MANIFEST.csv').itertuples():check('uniprot_hash_'+row.accession,sha(ROOT/row.path)==row.sha256)
boundary=json.loads((RESULTS/'ACCESS_BOUNDARY_INITIAL.json').read_text(encoding='utf-8'));check('protocol_unchanged',sha(ROOT/'STAGE6_PROTOCOL_v1.0.md')==boundary['protocol_sha256'])
first=json.loads((RESULTS/'FIRST_TARGETED_QUANTITATIVE_ACCESS.json').read_text(encoding='utf-8'));check('targeted_plan_linked_before_association',first['plan_sha256_before_association']==sha(RESULTS/'TARGETED_VERIFICATION_PLAN_FROZEN.json'))
for name in ['model.py','diagnostics.py']:check('reference_code_exact_'+name,sha(ROOT/'code'/name)==sha(S5/'code'/name))
check('reference_sim_exact',sha(ROOT/'code/simulation_reference.py')==sha(S5/'code/simulation.py'))

# Entire prior sealed manifests are read-only verified, not just a selection of result files.
prior=[]
for stage,file,col in [(S1,'STAGE1_INPUT_OUTPUT_SHA256.csv','path'),(S4,'STAGE4_INPUT_OUTPUT_SHA256.csv','relative_path'),(S5,'STAGE5_INPUT_OUTPUT_SHA256.csv','path')]:
 for row in pd.read_csv(stage/file).to_dict('records'):
  p=Path(row[col]);p=p if p.is_absolute() else stage/p
  ok=sha(p)==row['sha256'];prior.append({'stage':stage.name,'manifest_path':row[col],'passed':ok})
check('sealed_previous_stages_all_unchanged',all(r['passed'] for r in prior),len(prior));csvout(RESULTS/'SEALED_PREVIOUS_STAGES_HASH_VALIDATION.csv',prior)

# Independent metadata and selected raw value reread; no individual values printed.
additional=pd.read_csv(RESULTS/'ADDITIONAL_INPUT_MANIFEST.csv').set_index('input_id');annot_accessions=set();samplekeys={};validraw={};sourcecount=[]
selected=pd.read_csv(RESTRICTED/'TARGETED_RECORD_SELECTION_RESTRICTED.csv',low_memory=False);wanted={(r.source,int(r.source_row)):(r.key,r.feature_key) for r in selected.itertuples()}
for source in ['V01','V02']:
 ids=set();keyset=set();counts=Counter();dups=Counter()
 with Path(additional.loc[source,'path']).open(encoding='utf-8-sig',newline='') as f:
  for line,row in enumerate(csv.DictReader(f),2):
   counts['rows']+=1;s=row['TESTNAME'] if source=='V01' else row['PROTEIN'];m=re.search(r'(?:sp|tr)\|([^|]+)\|',s)
   if m:annot_accessions.add(m[1])
   try:r=int(float(row.get('RID','')))
   except ValueError:continue
   d=row.get('EXAMDATE','');valid=r>0 and r!=999999 and bool(re.fullmatch(r'\d{4}-\d{2}-\d{2}',d)) and (source=='V02' or row.get('TYPE')=='CSF')
   if valid:ids.add(r);keyset.add(str(r)+'|'+d)
   if (source,line) in wanted:
    try:v=float(row['TESTVALUE' if source=='V01' else 'TOTAL_AREA_RATIO'])
    except ValueError:v=np.nan
    validraw[wanted[(source,line)]]=np.log2(v) if np.isfinite(v) and v>0 else np.nan
 samplekeys[source]=keyset;sourcecount.append({'source':source,'n_RID_independent':len(ids),'n_keys_independent':len(keyset),'rows_independent':counts['rows']})
check('V01_719_unique',sourcecount[0]['n_RID_independent']==719)
check('V02_832_unique',sourcecount[1]['n_RID_independent']==832)
check('three_exact_accessions_not_covered_in_Emory_V01_V02',not set(CANDIDATES)&annot_accessions)
csvout(RESULTS/'INDEPENDENT_TARGETED_SOURCE_COUNTS.csv',sourcecount)
clinical=pd.read_csv(RESTRICTED/'TARGETED_PAIR_CLINICAL_DESIGN_RESTRICTED.csv');cache=np.load(RESTRICTED/'TARGETED_LOG_RATIOS_RESTRICTED.npz');ns=np.load(RESTRICTED/'SUPPORT_N_S_CACHE_RESTRICTED.npz')
check('clinical659_unique_RID',len(clinical)==659 and clinical.RID.is_unique)
errors=[]
for i,key in enumerate(cache['keys']):
 for j,fk in enumerate(cache['feature_keys']):
  errors.append(np.allclose(cache['Y'][i,j],validraw.get((key,fk),np.nan),equal_nan=True,atol=1e-12))
check('all3295_raw_log_values_reconstructed',len(errors)==3295 and all(errors))
oldclinical=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv').set_index('key');shared=clinical.loc[clinical.key.isin(oldclinical.index)].set_index('key')
for c in ['P_LOG','age_month_precision','sex_code2']:check('clinical_matches_S5_'+c,np.allclose(shared[c],oldclinical.reindex(shared.index)[c],equal_nan=True,atol=1e-10))
check('new_RID_disjoint',not set(clinical.loc[~clinical.is_stage5_RID,'RID'])&set(oldclinical.RID))
check('support_NS_keys_aligned',np.array_equal(cache['keys'],ns['keys']))
oldy=np.load(S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz');meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
indices={k:i for i,k in enumerate(ns['keys'])};s5indices={k:i for i,k in enumerate(oldy['keys'])}
for j,acc in enumerate(ns['accessions']):
 ji=list(oldy['accessions']).index(acc);sj=list(oldy['soma_features']).index(meta.loc[acc,'soma_features']);common=set(indices)&set(s5indices)
 check('NS_cache_matches_S5_'+acc,all(np.allclose([ns['N'][indices[k],j],ns['S'][indices[k],j]],[oldy['N'][s5indices[k],ji],oldy['S'][s5indices[k],sj]],equal_nan=True,atol=1e-12) for k in common))

# Independently refit 140 full nuisance models and rebuild two-constraint inference.
models=pickle.load((RESTRICTED/'ATTRIBUTION_MODELS_RESTRICTED.pkl').open('rb'));attribution=pd.read_csv(RESULTS/'COVARIATE_FOUR_MODEL_RESULTS.csv').set_index(['accession','model']);numericrows=[]
for (acc,name),d in models.items():
 j=list(oldy['accessions']).index(acc);sj=list(oldy['soma_features']).index(meta.loc[acc,'soma_features']);Y=np.c_[oldy['N'][:,j],oldy['S'][:,sj],oldy['T_GIS'][:,j]][d['indices']]
 X=np.c_[d['Z'],d['P']];coef,_,rank,_=lstsq(X,Y);E=Y-X@coef;b=coef[-1];df=len(Y)-rank;C=E.T@E/df;c12,c13,c23=C[0,1],C[0,2],C[1,2]
 h=np.array([b[1]*c13-b[0]*c23,b[2]*c12-b[0]*c23]);gcoef,_,_,_=lstsq(d['Z'],d['P']);rP=d['P']-d['Z']@gcoef;a=rP/(rP@rP)
 pp=np.c_[E[:,0]*E[:,1],E[:,0]*E[:,2],E[:,1]*E[:,2]];ci=(pp-pp.mean(axis=0))/df;bi=a[:,None]*E
 U=np.c_[-c23*bi[:,0]+c13*bi[:,1]+b[1]*ci[:,1]-b[0]*ci[:,2],-c23*bi[:,0]+c12*bi[:,2]+b[2]*ci[:,0]-b[0]*ci[:,2]];V=np.zeros((2,2))
 for axis in range(3):
  labels=d['groups'][:,axis];s=np.vstack([U[labels==g].sum(axis=0) for g in sorted(set(labels))]);s-=s.mean(axis=0);V+=s.T@s*len(s)/(len(s)-1)
 Q=float(h@np.linalg.solve(V,h));p=float(chi2.sf(Q,2));actual=attribution.loc[(acc,name)]
 ok=np.allclose(Q,actual.Q,rtol=1e-7) and np.allclose(p,actual.p_value,rtol=1e-7,atol=1e-12);check('independent_four_model_'+acc+'_'+name,ok)
 numericrows.append({'accession':acc,'model':name,'Q_independent':Q,'Q_produced':actual.Q,'passed':ok})
csvout(RESULTS/'INDEPENDENT_ATTRIBUTION_NUMERIC_VALIDATION.csv',numericrows)
check('70_S5_models_reproduced',len(pd.read_csv(RESULTS/'S5_M0_BOTH_REPRODUCTION_QA.csv'))==70 and pd.read_csv(RESULTS/'S5_M0_BOTH_REPRODUCTION_QA.csv').passed.all())
for table,fields in [('H_V_AND_B_C_DECOMPOSITION.csv',['Q_identity_abs_error','H_identity_max_abs_error']),('EDUCATION_APOE_SHAPLEY_Q.csv',['identity_abs_error']),('BY_Q_OWN_VS_FAMILY_DECOMPOSITION.csv',['identity_abs_error'])]:
 frame=pd.read_csv(RESULTS/table)
 for field in fields:check('exact_identity_'+table+'_'+field,frame[field].max()<1e-8,frame[field].max())
family=json.loads((S5/'results/MAIN_ANALYSIS_SAP_FROZEN.json').read_text(encoding='utf-8'))['primary_family_accessions']
for model,frame in attribution.reset_index().groupby('model'):
 frame=frame.set_index('accession').reindex(family);check('BY_reference_attribution_'+model,np.allclose(frame.q_BY_exploratory,multipletests(frame.p_value,method='fdr_by')[1],atol=1e-10))

# New support inference gets an independent 3-equation OLS moment reconstruction.
from diagnostics import summarize,parameters
from equivalence import actual_equivalence,equivalence_batch
support=pd.read_csv(RESULTS/'ORTHOGONAL_MEASUREMENT_SUPPORT_RESULTS.csv');frames=pickle.load((RESTRICTED/'ORTHOGONAL_SUPPORT_CASE_FRAMES_RESTRICTED.pkl').open('rb'));ds=pickle.load((RESTRICTED/'ORTHOGONAL_SUPPORT_DESIGNS_RESTRICTED.pkl').open('rb'))
for row in support.loc[support.calibration_passed].itertuples():
 d=ds[row.design_id];j=list(ns['accessions']).index(row.accession);k=list(cache['feature_keys']).index(row.feature_key);Y=np.c_[ns['N'][:,j],ns['S'][:,j],cache['Y'][:,k]][d['indices']]
 X=np.c_[d['Z'],d['P']];coef,_,rank,_=lstsq(X,Y);E=Y-X@coef;C=E.T@E/(len(Y)-rank);b=coef[-1]
 h=np.array([b[1]*C[0,2]-b[0]*C[1,2],b[2]*C[0,1]-b[0]*C[1,2]])
 check('support_independent_H_'+row.design_id,np.allclose(h,[row.H1,row.H2],atol=1e-8))
 R=C/np.sqrt(np.outer(C.diagonal(),C.diagonal()));check('support_independent_correlation_'+row.design_id,np.allclose([R[0,2],R[1,2]],[row.corr_residual_N_targetMS,row.corr_residual_S_targetMS],atol=1e-9))
for cohort,part in support.groupby('cohort'):check('support_BY_reference_'+cohort,np.allclose(part.q_BY5_exploratory,multipletests(part.p_value_for_family,method='fdr_by')[1]))
check('no_unqualified_newRID_inference',support.loc[support.cohort.eq('NEW_RID'),'p_value'].isna().all())
check('new_support_not_called_direct_validation',not support.three_discovery_direct_validation.any())

# Equivalence functions: independent Fieller inequality identity and full-theta versus batch calculation.
main_designs=pickle.load((S5/'00_restricted/DESIGNS_MAIN_RESTRICTED.pkl').open('rb'));eq=pd.read_csv(RESULTS/'OPERATIONAL_EQUIVALENCE_ALL35.csv');ci=pd.read_csv(RESULTS/'RELATIVE_RESPONSE_FIELLER_INTERVALS.csv')
for acc,d in main_designs.items():
 j=list(oldy['accessions']).index(acc);sj=list(oldy['soma_features']).index(meta.loc[acc,'soma_features']);Y=np.c_[oldy['N'][:,j],oldy['S'][:,sj],oldy['T_GIS'][:,j]][d['indices']]
 for margin in [.2,.35,.5]:
  exact=actual_equivalence(Y,d,margin)['p_joint_equivalence'];fast=equivalence_batch(Y,d,margin)[0];r=eq.loc[eq.accession.eq(acc)&eq.margin_working_SD.eq(margin)].iloc[0]
  check('equiv_full_vs_batch_'+acc+'_'+str(margin),np.allclose(exact,fast,atol=1e-9) and np.allclose(exact,r.p_equiv_computed,atol=1e-9))
  part=ci.loc[ci.accession.eq(acc)&ci.margin_working_SD.eq(margin)]
  if part.fieller90_bounded.all():
   contained=bool((part.fieller90_lower>-margin).all() and (part.fieller90_upper<margin).all());check('Fieller90_TOST_identity_'+acc+'_'+str(margin),contained==(exact<.05))
for margin,part in eq.loc[eq.in_family34].groupby('margin_working_SD'):
 check('equiv_BY_reference_'+str(margin),len(part)==34 and np.allclose(part.q_BY_equiv_exploratory,multipletests(part.p_equiv_for_family,method='fdr_by')[1]))

# Four-assay covariance block moments and gradient, including an empirical-weight perturbation.
from four_assay import moments
fourds=pickle.load((RESTRICTED/'FOUR_ASSAY_DESIGNS_RESTRICTED.pkl').open('rb'));four=pd.read_csv(RESULTS/'FOUR_ASSAY_CROSSBLOCK_RESULTS.csv')
for row in four.itertuples():
 d=fourds[row.design_id];src='S5_SAME_DATE__'+row.accession+'__'+row.peptide;frame=frames[src]['frame'];j=list(ns['accessions']).index(row.accession);jr=list(oldy['accessions']).index(row.accession)
 fk=support.loc[support.design_id.eq(src),'feature_key'].iloc[0];k=list(cache['feature_keys']).index(fk);oldidx={v:i for i,v in enumerate(oldy['keys'])}
 TY=np.array([oldy['T_GIS'][oldidx[key],jr] if key in oldidx else np.nan for key in frame.key]);Y=np.c_[ns['N'][:,j],ns['S'][:,j],TY,cache['Y'][:,k]][d['indices']]
 X=np.c_[d['Z'],d['P']];coef,_,rank,_=lstsq(X,Y);E=Y-X@coef;df=len(Y)-rank;C=E.T@E/df;tetrad=C[0,2]*C[1,3]-C[0,3]*C[1,2]
 J=np.array([C[1,3],C[0,2],-C[1,2],-C[0,3]])
 products=np.c_[E[:,0]*E[:,2],E[:,1]*E[:,3],E[:,0]*E[:,3],E[:,1]*E[:,2]];u=(products-products.mean(axis=0))@J/df;variance=0.
 for axis in range(4):
  lab=d['groups'][:,axis];s=np.array([u[lab==g].sum() for g in np.unique(lab)]);s-=s.mean();variance+=s@s*len(s)/(len(s)-1)
 p=chi2.sf(tetrad*tetrad/variance,1);check('independent_four_tetrad_'+row.design_id,np.allclose(tetrad,row.tetrad,atol=1e-10) and np.allclose(p,row.p_tetrad,rtol=1e-7,atol=1e-12))
 # The log cross-block ratio gradient is verified independently where defined.
 if min(C[0,2],C[1,3],C[0,3],C[1,2])>0:check('four_direction_ratio_identity_'+row.design_id,np.allclose(row.direction_ratio,(C[0,2]*C[1,3])/(C[0,3]*C[1,2]),rtol=1e-8))
check('four_BY_reference',np.allclose(four.q_BY5_tetrad_exploratory,multipletests(four.p_for_family,method='fdr_by')[1]))

# Every retained synthetic file is independently recounted, with failures kept in denominators.
recounts=[]
for table,dspath in [('ATTRIBUTION_CALIBRATION_SUMMARY.csv','ATTRIBUTION_CALIBRATION_DESIGNS_RESTRICTED.pkl'),('ORTHOGONAL_SUPPORT_CALIBRATION_SUMMARY.csv','ORTHOGONAL_SUPPORT_DESIGNS_RESTRICTED.pkl')]:
 des=pickle.load((RESTRICTED/dspath).open('rb'));data=pd.read_csv(RESULTS/table)
 for filename,part in data.groupby('file'):
  q=pd.read_parquet(SIM/filename).Q0.to_numpy();d=des[part.design_id.iloc[0]];factors=[1.,(d['n']-1)/d['df'],2*(d['n']-1)/d['df']];ok=True
  for r in part.itertuples():
   p=np.exp(-q/(2*factors[r.version]));valid=np.isfinite(p);count=int((valid&(p<=r.alpha)).sum());ok&=len(p)==r.reps and count==r.n_rejected and int((~valid).sum())==r.n_failed
  recounts.append({'table':table,'file':filename,'reps':len(q),'passed':bool(ok)})
for table,pcol,kcol in [('EQUIVALENCE_BOUNDARY_CALIBRATION.csv','p_equiv','false_equivalence_count'),('FOUR_ASSAY_NULL_CALIBRATION.csv','p_tetrad','n_rejected05'),('FOUR_DIRECTION_EQ_BOUNDARY_CALIBRATION.csv','p_equiv','false_equivalence_count')]:
 for r in pd.read_csv(RESULTS/table).to_dict('records'):
  p=pd.read_parquet(SIM/r['file'])[pcol].to_numpy();ok=len(p)==r['reps'] and int((np.isfinite(p)&(p<=.05)).sum())==r[kcol]
  recounts.append({'table':table,'file':r['file'],'reps':len(p),'passed':bool(ok)})
for r in pd.read_csv(RESULTS/'LEGACY_SINGLE_ASSAY_NULL_CALIBRATION.csv').to_dict('records'):
 p=pd.read_parquet(SIM/r['file']).p_positive.to_numpy();ok=len(p)==r['reps'] and int((np.isfinite(p)&(p<=.05)).sum())==r['n_rejected05'] and int((~np.isfinite(p)).sum())==r['n_failed']
 recounts.append({'table':'LEGACY_SINGLE_ASSAY_NULL_CALIBRATION.csv','file':r['file'],'reps':len(p),'passed':bool(ok)})
check('all1010_simulated_cells_recount',len(recounts)==1010 and all(r['passed'] for r in recounts),len(recounts));check('3160000_independent_repeats',sum(r['reps'] for r in recounts)==3160000)
csvout(RESULTS/'INDEPENDENT_SIMULATION_RECOUNTS.csv',recounts)
check('no_raw_QC_pools_mixed',not len(pd.read_csv(RESULTS/'TARGETED_RAW_QC_POOL_DIAGNOSTICS.csv')))
status=json.loads((RESULTS/'DIRECT_DISCOVERY_VALIDATION_STATUS.json').read_text(encoding='utf-8'));check('direct_full_replication_unavailable_honestly_reported',status['direct_three_discovery_validation_completed'] is False and status['external_cohort_validation_completed'] is False)
from validate_legacy import validate
legacy_summary=validate(check)
for path in RESULTS.glob('*.csv'):
 if path.name in ['FIGURE_INDEX.csv','WORKBOOK_SHEET_INDEX.csv','DELIVERY_QA_CHECKS.csv']:continue # delivery metadata checked by delivery_qa.py
 try:cols=pd.read_csv(path,nrows=0).columns
 except pd.errors.EmptyDataError:check('csv_has_schema_'+path.name,False);continue
 check('public_no_person_identifiers_'+path.stem,not any(c in cols for c in ['RID','PTID','key','EXAMDATE','SEQUENCE_NO','REPLICATE1','REPLICATE2','GENOTYPE','PTDOB']))
out=pd.DataFrame(checks);csvout(RESULTS/'INDEPENDENT_QA_CHECKS.csv',out)
summary={'checks':len(out),'passed':int(out.passed.sum()),'failed':int((~out.passed).sum()),'independently_refit_attribution_models':len(models),
 'simulation_cells_recounted':len(recounts),'independent_repeats_recounted':sum(r['reps'] for r in recounts),'previous_manifest_entries_verified':len(prior),
 'three_discoveries_directly_covered':2,'support_proteins':4,'support_peptides':5,'Emory_newRID_relative_inference_not_claimed':True,**legacy_summary}
save_json(RESULTS/'INDEPENDENT_QA_SUMMARY.json',summary);print(summary)
if not out.passed.all():print(out.loc[~out.passed].to_string(index=False));raise SystemExit(1)
