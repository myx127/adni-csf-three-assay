from config import *
from scipy.linalg import lstsq
import pickle

def solve(T,y):return lstsq(T,y,cond=max(T.shape)*np.finfo(float).eps,lapack_driver='gelsy')[0]
def run():
 ts=pickle.load((PRIVATE/'ANALYSIS_TEMPLATES_RESTRICTED.pkl').open('rb'));c=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv');effects=pd.read_csv(RESULTS/'PLATFORM_CLINICAL_EFFECTS_ALL43.csv')
 save_json(RESULTS/'DESCRIPTIVE_EXTENSION_ACCESS_RECORD.json',{'time':now(),'original_real_results_already_known':True,'extension_sha256':sha(ROOT/'DESCRIPTIVE_CONDITIONING_EXTENSION_v1.2.md'),'new_hypothesis_tests':0,'not_causal_mediation':True})
 decomposition=[];pet=[]
 for acc,gene in CANDIDATES.items():
  main=ts['MEM5Y__'+acc];nop=ts['MEM_NO_PTAU__'+acc];subcl=ts['MEM_PET_MAIN__'+acc]
  assert np.array_equal(main['ids'],nop['ids'])
  assembled=pickle.load((PRIVATE/'ASSEMBLED_CLINICAL_RESTRICTED.pkl').open('rb'))
  weights=np.array([assembled['trajectories'][(i,'ADNI_MEM',5)]['w'] for i in main['ids']]);Plift=c.loc[main['ids'],'P_LOG'].to_numpy()*np.sqrt(weights)
  lookup={i:j for j,i in enumerate(main['ids'])};sel=np.array([lookup[i] for i in subcl['ids']])
  for j,assay in enumerate(['N','S','T']):
   T0=np.c_[nop['ZZ'],nop['targetx'][:,j]];T1=np.c_[T0,Plift]
   fit0=solve(T0,main['y']);fit1=solve(T1,main['y']);pfit=solve(T0,Plift)
   b0=fit0[-1];b1=fit1[-2];delta=fit1[-1];pi=pfit[-1];sigma=main['sigma'][j]
   decomposition.append({'gene':gene,'assay':assay,'n_people':main['n_people'],'standardized_b_without_ptau':b0*sigma,'standardized_b_with_ptau':b1*sigma,'change_without_minus_with':(b0-b1)*sigma,'ptau_outcome_slope_coefficient':delta,'ptau_per_proteinSD_conditional':pi*sigma,'product_delta_pi':delta*pi*sigma,'identity_abs_error':abs(b0-b1-delta*pi)*sigma,'not_causal_mediation_or_signal_fraction':True})
   coeff_sub=solve(np.c_[main['ZZ'][sel],main['targetx'][sel,j]],main['y'][sel])[-1]*sigma
   bfull=effects.loc[effects.design_id.eq(main['design_id'])&effects.assay.eq(assay),'standardized_effect'].iloc[0];bcl=effects.loc[effects.design_id.eq(subcl['design_id'])&effects.assay.eq(assay),'standardized_effect'].iloc[0]
   pet.append({'gene':gene,'assay':assay,'n_full':main['n_people'],'n_PET_subset':subcl['n_people'],'full_main_standardized_b':bfull,'PET_subset_same_covariates_b':coeff_sub,'PET_subset_plus_CL_b':bcl,'point_change_due_to_subset':coeff_sub-bfull,'point_change_adding_CL_same_subset':bcl-coeff_sub,'descriptive_only_no_new_P':True})
 csvout(RESULTS/'PTAU_CONDITIONING_ALGEBRAIC_DECOMPOSITION.csv',decomposition);csvout(RESULTS/'PET_SUBSET_VS_CL_POINT_DECOMPOSITION.csv',pet)
 assert max(r['identity_abs_error'] for r in decomposition)<1e-8
 log('Result-informed descriptive conditioning decomposition completed;no new tests or main-result changes')
 print(pd.DataFrame(decomposition)[['gene','assay','standardized_b_without_ptau','standardized_b_with_ptau','product_delta_pi','identity_abs_error']].to_string(index=False))
 print(pd.DataFrame(pet).to_string(index=False))
if __name__=='__main__':
 assert_unsealed();run()
