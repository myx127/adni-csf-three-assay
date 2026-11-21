from config import *
from model import make_design,bh,factors
from analyze_main import actual_Y
from diagnostics import parameters,constraints_jacobian,summarize
import pickle

sap=json.loads((RESULTS/'MAIN_ANALYSIS_SAP_FROZEN.json').read_text(encoding='utf-8'));version=sap['version']
main=pd.read_csv(RESULTS/'PRIMARY_ALL43_RESULTS.csv');rows=[]
for r in main.loc[main.analysis_status.ne('EXCLUDED_PRE_EFFECT')].itertuples():
 corr=np.array([r.corr_residual_NS,r.corr_residual_NT,r.corr_residual_ST]);tri=float(np.prod(corr));load=np.full(3,np.nan)
 if tri>0 and np.min(abs(corr))>1e-12:load=np.array([corr[0]*corr[1]/corr[2],corr[0]*corr[2]/corr[1],corr[1]*corr[2]/corr[0]])
 uniq=1-load
 rows.append({'accession':r.accession,'gene':r.gene,'cross_cov_product_positive':tri>0,'all_cross_cov_positive':bool((corr>0).all()),
  'one_factor_unique_variance_nonnegative_at_point':bool(np.isfinite(uniq).all() and np.min(uniq)>=0),
  **{label+'_unique_variance_fraction_point':v for label,v in zip(['N','S','T'],uniq)},
  'point_diagnostic_not_a_new_filter':True,'no_primary_family_or_pvalue_changed':True})
csvout(RESULTS/'COMMON_COVARIANCE_COMPATIBILITY_DIAGNOSTICS.csv',rows)
clinical=pd.read_csv(RESTRICTED/'CLINICAL_DESIGN_RESTRICTED.csv');cache=np.load(RESTRICTED/'PROTEIN_TRANSFORMED_RESTRICTED.npz');meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
ds=pickle.load((RESTRICTED/'DESIGNS_MAIN_RESTRICTED.pkl').open('rb'));loco=[];gischecks=[]
for acc,d in ds.items():
 Yall=actual_Y(cache,meta,acc);Y=Yall[d['indices']];theta,IF,V,E=parameters(Y,d);h,J=constraints_jacobian(theta)
 # GIS log reference is a protein-specific plex intercept; common-subset residuals are equal after plex FE.
 Yr=actual_Y(cache,meta,acc,'RAW_TMT')[d['indices']];er=Yr-d['Q']@(d['Q'].T@Yr)
 gischecks.append({'accession':acc,'check':'RAW_vs_GIS_same_rows_after_full_plex_FE','max_residual_absolute_difference':float(abs(er-E).max()),'passed':bool(np.allclose(er,E,atol=1e-9))})
 for axis,col in enumerate(['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex']):
  one=d['onehots'][axis];s=one@(IF@J.T);ener=(s*s/np.maximum(np.diag(J@V@J.T),1e-100)).sum(axis=1);g=int(np.argmax(ener))
  remove=d['groups'][:,axis]==g;keep=d['indices'][~remove];mask=np.zeros(len(clinical),bool);mask[keep]=True;newd=make_design(clinical,mask)
  batch=clinical.loc[d['indices'][remove],col].iloc[0]
  valid=newd is not None and newd['df']>=100
  stats=summarize(Yall[newd['indices']],newd,version)[0] if valid else {'p_value':np.nan,'Q':np.nan}
  loco.append({'accession':acc,'gene':meta.loc[acc,'nulisa_target'],'axis':col,'omitted_batch_label':str(batch),'n_omitted':int(remove.sum()),
   'remaining_n':int((~remove).sum()),'p_value':stats['p_value'],'Q':stats['Q'],'nuisance_refit':True,'selection':'largest observed H influence energy; diagnostic only',
   'not_independently_recalibrated_after_omission':True})
csvout(RESULTS/'GIS_FIXED_EFFECT_EQUIVALENCE_QA.csv',gischecks);assert all(r['passed'] for r in gischecks)
loco=pd.DataFrame(loco);loco['q_BY_diagnostic']=np.nan
for axis,part in loco.groupby('axis'):
 keep=part.accession.isin(sap['primary_family_accessions']);pp=part.loc[keep,'p_value'].fillna(1).to_numpy()
 loco.loc[part.index[keep],'q_BY_diagnostic']=bh(pp,True)
csvout(RESULTS/'TOP_INFLUENCE_CLUSTER_OMISSION_DIAGNOSTICS.csv',loco)
save_json(RESULTS/'DIAGNOSTIC_EXTENSION_NOTE.json',{'added_after_first_main_result':True,
 'scope':'Common covariance admissibility and top-influence-cluster diagnostics; no changed test, family, transform or discovery definition',
 'rationale':'Statistical rejection must be distinguished from molecular common-factor interpretation; plex FE makes raw/GIS equivalent on same rows',
 'all_qualified_proteins_included_not_only_significant':True,'no_primary_p_value_changed':True})
log('All-qualified covariance compatibility, GIS-FE equivalence and influence omission diagnostics completed; primary family unchanged')
