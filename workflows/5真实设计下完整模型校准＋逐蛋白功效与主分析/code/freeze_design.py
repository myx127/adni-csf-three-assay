from config import *
from model import make_design

def build():
 clinical=pd.read_csv(RESTRICTED/'CLINICAL_DESIGN_RESTRICTED.csv')
 masks=pd.read_parquet(S1/'00_restricted/TRIPLE_CANDIDATE_AVAILABILITY_RESTRICTED.parquet')
 meta=pd.read_csv(S1/'results/PROTEIN_EVALUABILITY_70.csv').set_index('accession')
 designs={};rows=[];lookup={}
 for acc,m in meta.iterrows():
  mask=masks.loc[masks.accession.eq(acc)].set_index('key').reindex(clinical.key)
  good=mask.complete_GIS_S_all.to_numpy()&clinical.base_primary.to_numpy()
  des=make_design(clinical,good)
  n=int(good.sum());det=float(mask.N_detect_all_repeats.mean());detsub=float(mask.loc[good,'N_detect_all_repeats'].mean())
  reasons=[]
  if not m.annotation_strict_candidate:reasons.append('MULTI_FEATURE_OR_COMPOSITE')
  if n<200:reasons.append('N_LT200')
  if des is None or des['df']<100:reasons.append('RESIDUAL_DF_LT100')
  if des is not None and min(des['summary'][p+'_G'] for p in ['N','S','T'])<10:reasons.append('FEWER_THAN10_CLUSTERS')
  if det<.8 or detsub<.8:reasons.append('NULISA_DETECTION_LT80PCT')
  status='DESIGN_QUALIFIED' if not reasons else 'EXCLUDED_PRE_EFFECT'
  role='MAPT_CONSTRUCT_CONTROL' if m.nulisa_target=='MAPT' else 'PRIMARY_DISCOVERY_CANDIDATE'
  record={'accession':acc,'gene':m.nulisa_target,'single_annotation':bool(m.annotation_strict_candidate),'role':role,'design_status':status,
   'exclusion_reasons':';'.join(reasons),'N_detection_fraction_571':det,'N_detection_fraction_analysis':detsub,**(des['summary'] if des else {'n':n})}
  rows.append(record)
  if not reasons:designs[acc]=des;lookup[acc]=record
 csvout(RESULTS/'PRE_EFFECT_PROTEIN_SCREEN.csv',rows)
 # Designs contain clinical values and identities and remain restricted.
 import pickle
 with (RESTRICTED/'DESIGNS_MAIN_RESTRICTED.pkl').open('wb') as f:pickle.dump(designs,f,pickle.HIGHEST_PROTOCOL)
 save_json(RESULTS/'SIMULATION_CONFIG_FROZEN.json',{'protocol_sha256':sha(ROOT/'STAGE5_PROTOCOL_v1.0.md'),
  'code_sha256':{p.name:sha(p) for p in (ROOT/'code').glob('*.py')},'seeds':SEEDS,'n_designs':len(designs),
  'qualified_accessions':sorted(designs),'screen_sha256':sha(RESULTS/'PRE_EFFECT_PROTEIN_SCREEN.csv'),
  'clinical_design_sha256':sha(RESTRICTED/'CLINICAL_DESIGN_RESTRICTED.csv'),'individual_protein_values_opened':False,
  'dev_reps':1000,'confirm_reps':4000,'power_reps':2000,'scenarios':['C01','C02','C03','C04','C05','C06'],
  'selection':'first variance version V0,V1,V2 passing every design/scenario development cell; confirmation independent'})
 log(f'Pre-effect screen frozen: {len(designs)} qualified including MAPT control; no protein Y read')
 print(pd.DataFrame(rows).query('single_annotation')[['gene','design_status','exclusion_reasons','n']].to_string(index=False))

if __name__=='__main__':build()
