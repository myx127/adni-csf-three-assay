from config import *
from scipy.linalg import lstsq
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from openpyxl import load_workbook
from information_audit import load,y_for,h_and_den,p_iut
from model import wilson
import pickle,math

def run():
 checks=[]
 def check(name,ok,detail=''):checks.append({'check':name,'passed':bool(ok),'detail':str(detail)})
 plan=json.loads((RESULTS/'STAGE7_PLAN_FROZEN.json').read_text(encoding='utf-8'));check('original_protocol_unchanged',sha(ROOT/'STAGE7_PROTOCOL_v1.0.md')==plan['protocol_sha256'])
 freeze=json.loads((RESULTS/'POWER_ENGINE_FROZEN.json').read_text(encoding='utf-8'))
 for key,p in [('engine_sha256',ROOT/'code/power_engine.py'),('runner_sha256',ROOT/'code/run_power.py'),('template_sha256',PRIVATE/'POWER_TEMPLATES_RESTRICTED.pkl'),('execution_supplement_sha256',ROOT/'POWER_EXECUTION_DETAIL_v1.1.md')]:check('frozen_'+key,sha(p)==freeze[key])
 for r in pd.read_csv(RESULTS/'REFERENCE_CODE_ORIGIN.csv').itertuples():check('reference_code_'+r.file,sha(ROOT/'code'/r.file)==r.sha256==sha(r.source))
 for r in pd.read_csv(RESULTS/'INPUT_SOURCE_CHECKS.csv').itertuples():check('raw_source_'+r.source_id,sha(r.path)==r.sha256)
 for r in pd.read_csv(RESULTS/'ADDITIONAL_METADATA_INPUTS.csv').itertuples():check('metadata_source_'+r.source_id,sha(r.path)==r.sha256)
 ex=json.loads((RESULTS/'EXTERNAL_ANNOTATION_SOURCE_MANIFEST.json').read_text(encoding='utf-8'));check('external_original_workbook_hash',sha(ROOT/ex['path'])==ex['sha256'])
 prior=[]
 for stage,manifest,col in [(S1,'STAGE1_INPUT_OUTPUT_SHA256.csv','path'),(S4,'STAGE4_INPUT_OUTPUT_SHA256.csv','relative_path'),(S5,'STAGE5_INPUT_OUTPUT_SHA256.csv','path'),(S6,'STAGE6_INPUT_OUTPUT_SHA256.csv','path')]:
  for r in pd.read_csv(stage/manifest).to_dict('records'):
   p=Path(r[col]);p=p if p.is_absolute() else stage/p;prior.append({'stage':stage.name,'manifest_path':r[col],'passed':sha(p)==r['sha256']})
 check('all3469_prior_sealed_entries_unchanged',len(prior)==3469 and all(r['passed'] for r in prior));csvout(RESULTS/'PRIOR_SEALED_FILE_CHECKS_FINAL.csv',prior)
 engine=pd.read_csv(RESULTS/'POWER_ENGINE_INDEPENDENT_CHECKS.csv');check('all147_engine_checks_pass',len(engine)==147 and engine.passed.all())
 ds,raw,meta=load();memo=pickle.load((PRIVATE/'REAL_MOMENTS_RESTRICTED.pkl').open('rb'));info=pd.read_csv(RESULTS/'CONSENSUS_INFORMATION_BOTTLENECKS.csv');ci=pd.read_csv(RESULTS/'FIELLER_PRECISION_AND_DENOMINATOR_INFORMATION.csv');precision=pd.read_csv(RESULTS/'CONDITIONAL_PRECISION_EXPANSION_SCENARIOS.csv')
 for acc,d in ds.items():
  Y=y_for(raw,meta,acc,d);X=np.c_[d['Z'],d['P']];beta,_,rank,_=lstsq(X,Y,cond=max(X.shape)*np.finfo(float).eps);E=Y-X@beta;df=len(Y)-rank
  coefp,_,_,_=lstsq(d['Z'],d['P'],cond=max(d['Z'].shape)*np.finfo(float).eps);rp=d['P']-d['Z']@coefp;a=rp/(rp@rp)
  prod=np.c_[E[:,0]**2,E[:,1]**2,E[:,2]**2,E[:,0]*E[:,1],E[:,0]*E[:,2],E[:,1]*E[:,2]];theta=np.r_[beta[-1],prod.sum(axis=0)/df];IF=np.c_[a[:,None]*E,(prod-prod.mean(axis=0))/df];V=np.zeros((9,9))
  for axis in range(3):
   labels=d['groups'][:,axis];s=np.vstack([IF[labels==g].sum(axis=0) for g in np.unique(labels)]);s-=s.mean(axis=0);V+=s.T@s*len(s)/(len(s)-1)
  check('independent_actual_theta_SUM3_V_'+acc,np.allclose(theta,memo[acc]['theta'],atol=1e-8) and np.allclose(V,memo[acc]['V'],atol=1e-8))
  A,JA,B,JB=h_and_den(theta)
  for j in range(9):
   step=1e-5*max(abs(theta[j]),.001);lo=theta.copy();hi=theta.copy();lo[j]-=step;hi[j]+=step;aa,_,bb,_=h_and_den(lo);ac,_,bc,_=h_and_den(hi)
   check('independent_numeric_AB_gradient_'+acc+'_'+str(j),np.allclose((ac-aa)/(2*step),JA[:,j],atol=1e-7,rtol=1e-5) and np.allclose((bc-bb)/(2*step),JB[:,j],atol=1e-7,rtol=1e-5))
  for r in precision.loc[precision.accession.eq(acc)&precision.fixed_cluster_variance_fraction_scenario_not_estimated.eq(0)].itertuples():
   k=r.required_full_independent_information_multiplier
   if np.isfinite(k):check('precision_multiplier_meets_IUT_'+acc+'_'+str(r.margin)+'_'+r.threshold,p_iut(A,B,JA,JB,V/(k*(1+1e-8)),r.margin)<=r.alpha+1e-7)
   else:check('precision_infinite_has_explicit_reason_'+acc+'_'+str(r.margin)+'_'+r.threshold,r.limitation in ['STRUCTURE_INCOMPATIBLE','POINT_OUTSIDE_MARGIN'])
  for r in ci.loc[ci.accession.eq(acc)&ci.bounded].itertuples():
   j=0 if r.contrast=='S_relative_to_N' else 1;va=JA[j]@V@JA[j];vb=JB[j]@V@JB[j];cab=JA[j]@V@JB[j];z=norm.isf(r.one_sided_alpha)
   error=max(abs((A[j]-x*B[j])**2-z*z*(va-2*x*cab+x*x*vb)) for x in [r.lower,r.upper]);check('independent_Fieller_endpoint_'+acc+'_'+r.contrast+'_'+r.threshold,error<1e-7,error)
 check('105_equivalence_reproduced',len(pd.read_csv(RESULTS/'STAGE6_EQUIVALENCE_REPRODUCTION.csv'))==105 and pd.read_csv(RESULTS/'STAGE6_EQUIVALENCE_REPRODUCTION.csv').passed.all())
 for m,part in info.loc[info.in_family34].groupby('margin'):
  check('fixed_family34_'+str(m),len(part)==34 and part.equivalence_eligible_Stage6.sum()==11 and part.frozen_exact_zero_deviation_BY.sum()==3 and part.BY_equivalence_Stage6.sum()==0)
 wide=info.loc[info.in_family34&info.margin.eq(.5)];check('six_nominal_not_BY_at_widest_margin',wide.nominal_equiv_eligible.sum()==6)
 rng=np.random.default_rng(2026103009)
 for i in range(20):
  p=rng.uniform(0,1,34);p[0]=ALPHA_BY_SUFFICIENT*.999;check('BY_sufficient_threshold_identity_'+str(i),multipletests(p,method='fdr_by')[1][0]<=.05)
 # Independent repeat counts, including nonfinite outcomes in total denominators.
 recount=[];seeds=[]
 for kind in ['POWER','BOUNDARY']:
  frame=pd.read_csv(RESULTS/(kind+'_SIMULATION_SUMMARY.csv'),dtype={'seed':str})
  for file,part in frame.groupby('file'):
   r=pd.read_parquet(SIM/file);ok=len(r)==2000 and np.array_equal(r.replicate,np.arange(2000))
   for x in part.itertuples():
    pp=r['p_equiv_d'+str(int(round(x.margin*100)))].to_numpy();success=np.isfinite(pp)&(pp<=x.alpha);ng=int((success&r.eligible.to_numpy(bool)).sum());nr=int(success.sum())
    low,high=wilson(ng,len(r));ok&=ng==x.gate_success_count and nr==x.raw_success_count and np.isclose(low,x.Wilson95_low) and np.isclose(high,x.Wilson95_high) and int((~r.valid).sum())==x.failure_count
   recount.append({'kind':kind,'file':file,'reps':len(r),'passed':bool(ok)});seeds.append(part.seed.iloc[0])
 check('all720_cells_independently_recounted',len(recount)==720 and all(r['passed'] for r in recount));check('1440000_independent_new_repetitions',sum(r['reps'] for r in recount)==1440000);check('720_unique_cell_seeds',len(set(seeds))==720)
 csvout(RESULTS/'INDEPENDENT_SIMULATION_RECOUNTS.csv',recount)
 template=pickle.load((PRIVATE/'POWER_TEMPLATES_RESTRICTED.pkl').open('rb'));check('15_templates_all_in_original_family',len(template)==15 and set(template)<=set(FAMILY))
 for acc,t in template.items():
  C=np.outer(t['lambda'],t['lambda'])+np.diag(t['psi'])
  for gamma in GAMMAS:
   b=t['lambda']*gamma;h=np.array([b[1]*C[0,2]-b[0]*C[1,2],b[2]*C[0,1]-b[0]*C[1,2]]);check('true_equivalence_template_'+acc+'_'+str(gamma),np.allclose(h,0,atol=1e-12))
  for m in MARGINS:
   for j in [1,2]:
    for sign in [-1,1]:
     b=t['lambda']*.55;b[j]+=sign*m*np.sqrt(C[j,j]);h=np.array([b[1]*C[0,2]-b[0]*C[1,2],b[2]*C[0,1]-b[0]*C[1,2]]);den=np.array([C[0,2]*np.sqrt(C[1,1]),C[0,1]*np.sqrt(C[2,2])]);truth=h/den
     check('true_boundary_exact_'+acc+'_'+str(m)+'_'+str(j)+'_'+str(sign),np.isclose(truth[j-1],sign*m) and np.isclose(truth[2-j],0))
 boundary=pd.read_csv(RESULTS/'REAL_TEMPLATE_BOUNDARY_DECISIONS.csv');check('90_endpoint_margin_boundaries_pass',len(boundary)==90 and boundary.raw_reference_pass.all())
 # Whole-vector source identity and all recovered technical labels independently re-read.
 source=pd.read_csv(S6/'results/LEGACY_INPUT_MANIFEST.csv').set_index('input_id');a=pd.read_csv(source.loc['L01','path']);resolved=pd.read_csv(PRIVATE/'MRM2014_RESOLVED_TECHNICAL_CROSSWALK_RESTRICTED.csv')
 book=load_workbook(source.loc['L02','path'],read_only=True,data_only=True);e=list(book['CSV-Export'].iter_rows(values_only=True));metadata=list(book['Sample Info'].iter_rows(values_only=True));book.close()
 eh=list(e[0]);cols=[(j,c.replace('-','_',1)) for j,c in enumerate(eh) if isinstance(c,str) and c.replace('-','_',1) in a]
 mh=list(metadata[0]);refcol=mh.index('External_Reference');injections={str(row[refcol]).upper():row for row in metadata[1:] if row[refcol] is not None}
 allok=[]
 for r in resolved.itertuples():
  ar=a.iloc[int(r.source_row)-2];er=e[int(r.export_row)+1];ok=int(ar.RID)==int(r.RID) and str(ar.VISCODE)==r.VISCODE
  for j,col in cols:
   av=float(ar[col]);bv=float(er[j]);resolution=10.**(math.floor(math.log10(abs(av)))-5) if av else 1e-5;ok&=abs(av-bv)<=resolution+1e-10
  mr=injections[str(r.External_Reference).upper()];ok&=str(mr[mh.index('LIMS_ID')])==str(r.LIMS_ID) and float(mr[mh.index('Processing_Day')])==r.Processing_Day and float(mr[mh.index('Injection_Order')])==r.Injection_Order
  allok.append(ok)
 check('303_source_rows_full320_vector_and_metadata_independent',len(allok)==303 and len(cols)==320 and all(allok))
 check('no_duplicates_promoted_to_new_people',resolved.RID.nunique()==287)
 dates=pd.read_csv(PRIVATE/'MRM2014_CANDIDATE_DATE_CROSSWALK_RESTRICTED.csv');check('287_unique_visit_date_candidates_not_verified',dates.n_unique_master_visit_drawdates.eq(1).sum()==287 and dates.candidate_equals_Stage5_date.sum()==149 and not dates.is_verified_MRM_specimen_collection_date.any())
 screen=pd.read_csv(RESULTS/'MRM2014_INFERRED_VISIT_DESIGN_SCREEN.csv');check('conditional2014_geometry_without_effect_fit',len(screen)==5 and sorted(screen.n.unique())==[141,144] and sorted(screen.residual_df.unique())==[75,78] and not screen.strict_same_day_replication_eligible.any() and not screen.fresh_candidate_association_fitted.any())
 # External platform membership is annotation, not per-person complete measurement or an effect.
 book=load_workbook(ROOT/ex['path'],read_only=True,data_only=True);coverage=pd.read_csv(RESULTS/'EXTERNAL_EMORY2022_FILTERED_GENE_PLATFORM_COVERAGE.csv')
 for r in coverage.itertuples():
  ws=book[r.sheet];values=[row for row in ws.iter_rows(values_only=True) if row[0]==r.gene];check('external_gene_membership_'+r.gene+'_'+r.sheet,len(values)==1 and list(map(bool,values[0][1:4]))==[r.SomaScan_filtered_gene_present,r.Olink_filtered_gene_present,r.MS_filtered_gene_present])
 book.close();check('external_coverage_three_genes_two_MS_methods',len(coverage)==6 and not coverage.NULISA_in_this_external_study.any() and coverage.per_gene_complete_case_n_not_verified.all())
 for p in RESULTS.glob('*.csv'):
  if p.name in ['FIGURE_INDEX.csv','WORKBOOK_SHEET_INDEX.csv','DELIVERY_QA_CHECKS.csv','INDEPENDENT_QA_CHECKS.csv']:continue
  columns=pd.read_csv(p,nrows=0).columns;check('public_no_individual_columns_'+p.stem,not set(columns)&{'RID','PTID','key','EXAMDATE','LIMS_ID','External_Reference','PTDOB','GENOTYPE','source_row'})
 output=pd.DataFrame(checks);csvout(RESULTS/'INDEPENDENT_QA_CHECKS.csv',output)
 summary={'checks':len(output),'passed':int(output.passed.sum()),'failed':int((~output.passed).sum()),'previous_sealed_entries_verified':len(prior),
  'new_simulation_cells_recounted':len(recount),'new_independent_repetitions_recounted':sum(r['reps'] for r in recount),'real_designs_refit':len(ds),'original_equivalence_results_reproduced':105,
  'recovered2014_technical_rows':len(resolved),'verified2014_MRM_collection_dates':0,'fresh_candidate_pathology_associations_fitted':0}
 save_json(RESULTS/'INDEPENDENT_QA_SUMMARY.json',summary);print(summary)
 if not output.passed.all():print(output.loc[~output.passed].to_string(index=False));raise SystemExit(1)

if __name__=='__main__':run()
