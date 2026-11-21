from config import *
from equivalence import actual_equivalence,MARGINS
from attribution import y_for
from model import bh
import pickle

cal=pd.read_csv(RESULTS/'EQUIVALENCE_BOUNDARY_CALIBRATION.csv');ds=pickle.load((S5/'00_restricted/DESIGNS_MAIN_RESTRICTED.pkl').open('rb'))
meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession');cache=np.load(S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz')
main=pd.read_csv(S5/'results/PRIMARY_ALL43_RESULTS.csv').set_index('accession');compat=pd.read_csv(S5/'results/COMMON_COVARIANCE_COMPATIBILITY_DIAGNOSTICS.csv').set_index('accession')
sap=json.loads((S5/'results/MAIN_ANALYSIS_SAP_FROZEN.json').read_text(encoding='utf-8'));family=sap['primary_family_accessions'];rows=[];contrasts=[]
for acc,d in ds.items():
 Y=y_for(cache,meta,acc)[d['indices']];struct=bool(compat.loc[acc,'all_cross_cov_positive'] and compat.loc[acc,'one_factor_unique_variance_nonnegative_at_point'])
 strong=bool(main.loc[acc,'positive_cov_CI_all_three'])
 for margin in MARGINS:
  check=cal.loc[cal.accession.eq(acc)&cal.margin.eq(margin)];cal_pass=bool(len(check)==4 and check.Wilson95_high.le(.06).all())
  result=actual_equivalence(Y,d,margin);bounded=all(r['fieller90_bounded'] for r in result['two_contrasts'])
  eligible=struct and strong and bounded and cal_pass
  rows.append({'accession':acc,'gene':meta.loc[acc,'nulisa_target'],'margin_working_SD':margin,'n':d['n'],'in_family34':acc in family,
   'simple_cov_structure_point_compatible':struct,'all_three_cov_CI95_positive':strong,'Fieller90_both_bounded':bounded,'boundary_calibration_passed':cal_pass,
   'eligible_for_operational_equivalence':eligible,'p_equiv_computed':result['p_joint_equivalence'],'p_equiv_for_family':result['p_joint_equivalence'] if eligible else 1.,
   'S5_BY_discovery':bool(main.loc[acc,'BY_discovery']),'not_true_abundance_or_individual_consensus_score':True})
  for r in result['two_contrasts']:contrasts.append({'accession':acc,'gene':meta.loc[acc,'nulisa_target'],'margin_working_SD':margin,**r})
out=pd.DataFrame(rows);out['q_BY_equiv_exploratory']=np.nan
for margin,part in out.groupby('margin_working_SD'):
 ids=part.index[part.in_family34];out.loc[ids,'q_BY_equiv_exploratory']=bh(out.loc[ids,'p_equiv_for_family'].to_numpy(),True)
out['operational_equivalence_supported_BY']=out.eligible_for_operational_equivalence&out.q_BY_equiv_exploratory.le(.05)
csvout(RESULTS/'OPERATIONAL_EQUIVALENCE_ALL35.csv',out);csvout(RESULTS/'RELATIVE_RESPONSE_FIELLER_INTERVALS.csv',contrasts)
summary=out.loc[out.in_family34].groupby('margin_working_SD').agg(n_family=('accession','size'),n_eligible=('eligible_for_operational_equivalence','sum'),n_BY_equiv=('operational_equivalence_supported_BY','sum')).reset_index()
csvout(RESULTS/'EQUIVALENCE_WORKING_MARGIN_SUMMARY.csv',summary);print(summary.to_string(index=False));log('Operational relative-response equivalence analyzed; untested/non-significant proteins not automatically labelled consensus')
