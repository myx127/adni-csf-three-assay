from config import *
from model import make_design
import pickle

clinical=pd.read_csv(RESTRICTED/'CLINICAL_DESIGN_RESTRICTED.csv')
masks=pd.read_parquet(S1/'00_restricted/TRIPLE_CANDIDATE_AVAILABILITY_RESTRICTED.parquet')
screen=pd.read_csv(RESULTS/'PRE_EFFECT_PROTEIN_SCREEN.csv');wanted=screen.loc[screen.design_status.eq('DESIGN_QUALIFIED')]
rows=[];designs={};links=[];fingerprints={}
main=pickle.load((RESTRICTED/'DESIGNS_MAIN_RESTRICTED.pkl').open('rb'))
def fingerprint(d):
 h=hashlib.sha256()
 for k in ['indices','P','Z','groups']:h.update(np.ascontiguousarray(d[k]).tobytes())
 return h.hexdigest()
for acc,d in main.items():fingerprints[fingerprint(d)]='MAIN__'+acc
for variant in ['RAW_P','RAW_TMT','LOD','BATCH9','FULLCOV','AMY']:
 for m in wanted.itertuples():
  part=masks.loc[masks.accession.eq(m.accession)].set_index('key').reindex(clinical.key)
  col='complete_raw_S_all' if variant=='RAW_TMT' else 'complete_GIS_detect_all_N_S_all' if variant=='LOD' else 'complete_GIS_S_all'
  good=part[col].to_numpy()&clinical.base_primary.to_numpy()
  if variant=='AMY':good=part[col].to_numpy()&np.isfinite(clinical[['P_AMY','age_month_precision','sex_code2']]).all(axis=1).to_numpy()
  if variant=='BATCH9':good&=clinical.ptau_batch.eq('UPENNBIOMK9').to_numpy()
  if variant=='FULLCOV':good&=np.isfinite(clinical[['education','apoe4_carrier']]).all(axis=1).to_numpy()
  des=make_design(clinical,good,variant);reasons=[]
  if des is None:reasons.append('P_NOT_ESTIMABLE')
  else:
   if des['n']<200:reasons.append('N_LT200')
   if des['df']<100:reasons.append('RESIDUAL_DF_LT100')
   if min(des['summary'][p+'_G'] for p in ['N','S','T'])<10:reasons.append('CLUSTERS_LT10')
  ident=variant+'__'+m.accession;link=''
  if not reasons:
   fp=fingerprint(des)
   if fp in fingerprints:link=fingerprints[fp]
   else:fingerprints[fp]=ident;designs[ident]=des;link=ident
  rows.append({'variant':variant,'accession':m.accession,'gene':m.gene,'role':m.role,'n':int(good.sum()),'status':'DESIGN_QUALIFIED' if not reasons else 'EXCLUDED_VARIANT_DESIGN',
   'reasons':';'.join(reasons),'calibration_design_id':link,**(des['summary'] if des else {})})
  links.append({'variant':variant,'accession':m.accession,'unique_design_id':link})
csvout(RESULTS/'VARIANT_DESIGN_SCREEN.csv',rows)
with (RESTRICTED/'DESIGNS_VARIANTS_RESTRICTED.pkl').open('wb') as f:pickle.dump(designs,f,pickle.HIGHEST_PROTOCOL)
save_json(RESULTS/'VARIANT_DESIGN_FREEZE.json',{'code_sha256':sha(Path(__file__)),'unique_additional_designs':len(designs),'rows':len(rows),'before_real_protein_Y_connection':True,
 'calibration_reps':4000,'scenarios':['C01','C02','C03','C04','C05','C06'],'same_primary_selected_variance_no_reselection':True})
log(f'Variant designs frozen: {len(rows)} rows; {len(designs)} unique additional designs after exact deduplication')
