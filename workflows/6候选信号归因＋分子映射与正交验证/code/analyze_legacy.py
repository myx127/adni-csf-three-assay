from config import *
from designs import make
from legacy_assay import single_summary
from analyze_support import correlations
from model import bh
import pickle

def run():
 plan=json.loads((RESULTS/'LEGACY_VERIFICATION_SAP_FROZEN.json').read_text(encoding='utf-8'));first=json.loads((RESULTS/'FIRST_LEGACY_QUANTITATIVE_ACCESS.json').read_text(encoding='utf-8'))
 assert first['plan_sha256_before_association']==sha(RESULTS/'LEGACY_VERIFICATION_SAP_FROZEN.json')
 clinical=pd.read_csv(RESTRICTED/'LEGACY_L04_CLINICAL_DESIGN_RESTRICTED.csv');cache=np.load(RESTRICTED/'LEGACY_PEPTIDE_CACHE_RESTRICTED.npz')
 assert np.array_equal(clinical.key,cache['keys']);masks=pickle.load((RESTRICTED/'LEGACY_CASE_MASKS_RESTRICTED.pkl').open('rb'))
 ds=pickle.load((RESTRICTED/'LEGACY_CLINICAL_DESIGNS_RESTRICTED.pkl').open('rb'));cal=pd.read_csv(RESULTS/'LEGACY_SINGLE_ASSAY_CALIBRATION_DECISIONS.csv').set_index('design_id')
 screen=pd.read_csv(RESULTS/'LEGACY_DIRECT_DESIGN_SCREEN.csv');old=np.load(S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz');mapping=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
 relative=[];assoc=[]
 for row in screen.itertuples():
  mask=masks[row.design_id];k=list(cache['peptides']).index(row.peptide);y=cache['Y'][:,k]
  base={'design_id':row.design_id,'cohort':row.cohort,'accession':row.accession,'gene':row.gene,'peptide':row.peptide,'n':row.n,
    'same_ADNI_parent_cohort_not_external':True,'exact_reported_date_not_proof_same_aliquot':True,'relative_disagreement_validated':False}
  if row.cohort.startswith('NEW_RID'):
   if row.design_id in ds:
    d=ds[row.design_id];passed=bool(cal.loc[row.design_id,'passed']);base.update(single_summary(y[d['indices']],d));base.update(d['summary'])
    base.update({'calibration_passed':passed,'p_for_family':base['p_positive_direction'] if passed else 1.,
      'analysis_status':'SUBJECT_DISJOINT_PEPTIDE_PTAU_ASSOCIATION_ONLY' if passed else 'CALIBRATION_FAILED_DESCRIPTIVE_ONLY'})
   else:base.update({'calibration_passed':False,'p_for_family':1.,'analysis_status':'DESIGN_INSUFFICIENT_DESCRIPTIVE_ONLY'})
   assoc.append(base)
  else:
   d=make(clinical,mask,tech=('NULISA_PlateID','SOMASCAN_PlateId','MRM_PLATE','MRM_RUN'),additional_factors=('ptau_batch','MRM_INSTRUMENT'))
   j=list(old['accessions']).index(row.accession);sj=list(old['soma_features']).index(mapping.loc[row.accession,'soma_features']);idx=clinical.loc[mask,'S5_index'].to_numpy(int)
   Y=np.c_[old['N'][idx,j],old['S'][idx,sj],y[mask]]
   base.update({'formal_relative_p_value':np.nan,'analysis_status':'DF_BELOW_FROZEN_40_GATE_DESCRIPTIVE_ONLY'})
   if d is not None and d['df']>0:
    b=d['a']@Y;E=Y-d['Q']@(d['Q'].T@Y);sd=np.sqrt((E*E).sum(axis=0)/d['df'])
    base.update({**d['summary'],**correlations(Y,d),'N_P_descriptive':b[0],'S_P_descriptive':b[1],'targetMS_P_descriptive':b[2],
     'N_P_per_residual_SD_descriptive':b[0]/sd[0],'S_P_per_residual_SD_descriptive':b[1]/sd[1],'targetMS_P_per_residual_SD_descriptive':b[2]/sd[2]})
   relative.append(base)
 out=pd.DataFrame(assoc);out['q_BY5_single_assay']=np.nan
 for cohort,part in out.groupby('cohort'):
  assert len(part)==5;out.loc[part.index,'q_BY5_single_assay']=bh(part.p_for_family.to_numpy(),True)
 out['positive_association_BY']=out.calibration_passed&out.q_BY5_single_assay.le(.05)
 csvout(RESULTS/'LEGACY_NEW_RID_PEPTIDE_PTAU_ASSOCIATION.csv',out);csvout(RESULTS/'LEGACY_SAME_DATE_RELATIVE_DESCRIPTIVE.csv',relative)
 log('Direct legacy peptide results complete; subject-disjoint pathology association kept separate from unperformed relative-pattern validation')
 print(out[['gene','peptide','cohort','n','b_P_per_residual_SD','p_positive_direction','q_BY5_single_assay','positive_association_BY']].to_string(index=False))

if __name__=='__main__':run()
