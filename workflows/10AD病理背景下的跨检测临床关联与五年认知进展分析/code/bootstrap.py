from config import *
def run():
 wanted=[Path(INPUTS[k]['path']) for k in ['I09','I11']]
 wanted += [S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv',S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz',S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv',S4/'00_restricted/AMYLOID_PET_FLAGS_RESTRICTED.csv',S4/'00_restricted/COGNITION_DATE_FLAGS_RESTRICTED.csv',S4/'00_restricted/COGNITION_TIMING_FLAGS_RESTRICTED.csv',S9/'SHA256_MANIFEST.csv']
 source=[]
 for p in wanted:
  digest=sha(p);expect=next((INPUTS[k]['sha256'] for k in ['I09','I11'] if Path(INPUTS[k]['path'])==p),None)
  source.append({'path':str(p),'bytes':p.stat().st_size,'sha256':digest,'original_expected_sha256':expect,'passed':expect is None or expect==digest})
 assert all(r['passed'] for r in source);csvout(RESULTS/'INPUT_SOURCE_MANIFEST.csv',source)
 old=[]
 for r in pd.read_csv(S9/'SHA256_MANIFEST.csv').itertuples():
  p=S9/r.relative_path;actual=sha(p);old.append({'path':str(p),'expected_sha256':r.sha256,'actual_sha256':actual,'passed':actual==r.sha256})
 assert all(r['passed'] for r in old);csvout(RESULTS/'STAGE9_SEAL_INITIAL_CHECKS.csv',old)
 save_json(RESULTS/'ACCESS_AND_PROTOCOL_RECORD.json',{'analysis_date':'2026-10-29','recorded_time':now(),'protocol_sha256':sha(ROOT/'STAGE10_PROTOCOL_v1.0.md'),'prior_stages_known':True,'new_cognitive_associations_not_yet_fitted':True,'authorization':'Explicit user request to complete Stage10 clinical association and five-year cognition analyses supersedes initial manifest POST_FREEZE_CLINICAL access labels; Stage10 protocol now fixed before reading new cognition values','original_files_remain_read_only':True})
 log(f'Stage10 protocol/access fixed;{len(source)} inputs registered;130 Stage9 sealed files checked')
if __name__=='__main__':
 assert_unsealed();run()
