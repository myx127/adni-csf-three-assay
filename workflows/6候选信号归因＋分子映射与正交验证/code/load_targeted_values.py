from config import *
from targeted_inventory import numeric,parse_feature
from collections import defaultdict
import re

def run():
 plan=RESULTS/'TARGETED_VERIFICATION_PLAN_FROZEN.json';assert plan.exists()
 first=RESULTS/'FIRST_TARGETED_QUANTITATIVE_ACCESS.json'
 if not first.exists():save_json(first,{'time':datetime.now().astimezone().isoformat(),'plan_sha256_before_association':sha(plan),
    'support_features_chosen_by_annotation_not_effect':True,'three_discoveries_not_covered':True})
 else:assert json.loads(first.read_text(encoding='utf-8'))['plan_sha256_before_association']==sha(plan)
 selected=pd.read_csv(RESTRICTED/'TARGETED_RECORD_SELECTION_RESTRICTED.csv',low_memory=False);lookup={(r.source,int(r.source_row)):(r.key,r.feature_key) for r in selected.itertuples()}
 clinical=pd.read_csv(RESTRICTED/'TARGETED_PAIR_CLINICAL_DESIGN_RESTRICTED.csv');featurelist=sorted(selected.feature_key.unique());values=np.full((len(clinical),len(featurelist)),np.nan)
 ridx={k:i for i,k in enumerate(clinical.key)};fidx={k:i for i,k in enumerate(featurelist)};qc=defaultdict(list);duplicates=defaultdict(list);nflags=0
 extras=pd.read_csv(RESULTS/'ADDITIONAL_INPUT_MANIFEST.csv').set_index('input_id');chosen_acc=set(selected.accession)
 for source in ['V01','V02']:
  with Path(extras.loc[source,'path']).open(encoding='utf-8-sig',newline='') as f:
   for line,row in enumerate(csv.DictReader(f),2):
    feature=parse_feature(row,source)
    if feature['accession'] not in chosen_acc:continue
    v=numeric(row.get('TESTVALUE' if source=='V01' else 'TOTAL_AREA_RATIO'))
    if (source,line) in lookup:
     key,fk=lookup[(source,line)]
     if key in ridx and v>0:values[ridx[key],fidx[fk]]=np.log2(v);nflags+=1
    rr=row.get('RID','');person=numeric(rr);isqc=person==999999 or (source=='V01' and row.get('TYPE')!='CSF')
    if isqc and np.isfinite(v) and v>0:
     text=' '.join(str(row.get(c,'')) for c in ['TYPE','PLATEPOSITION','REPLICATE1','REPLICATE2'])
     pool='ATpos' if re.search(r'ATpos|(?:^|_)AD(?:_|$)',text,re.I) else 'ATneg' if re.search(r'ATneg|CTL',text,re.I) else 'QCpool' if re.search('QCpool',text,re.I) else 'UNSPECIFIED_DO_NOT_POOL'
     if pool!='UNSPECIFIED_DO_NOT_POOL':qc[(source,feature['feature_key'],pool)].append(v)
    if np.isfinite(person) and person!=999999 and np.isfinite(v) and v>0:
     duplicates[(source,row.get('RID'),row.get('EXAMDATE'),feature['feature_key'])].append(np.log2(v))
 np.savez_compressed(RESTRICTED/'TARGETED_LOG_RATIOS_RESTRICTED.npz',Y=values,keys=clinical.key.to_numpy(str),feature_keys=np.array(featurelist))
 qrows=[]
 for (source,fk,pool),vv in qc.items():
  vals=np.array(vv);qrows.append({'source':source,'feature_key':fk,'QC_pool':pool,'n':len(vals),'mean_ratio_QC_only':float(vals.mean()),
   'CV_QC_only':float(vals.std(ddof=1)/vals.mean()) if len(vals)>1 else np.nan,'biological_participant_values_not_exported':True})
 csvout(RESULTS/'TARGETED_RAW_QC_POOL_DIAGNOSTICS.csv',pd.DataFrame(qrows,columns=['source','feature_key','QC_pool','n','mean_ratio_QC_only','CV_QC_only','biological_participant_values_not_exported']))
 csvout(RESULTS/'QC_POOL_LABEL_AVAILABILITY.csv',[{'source':s,'pool_type_recovered_for_selected_features':any(r['source']==s for r in qrows),
   'individual_QC_pool_CV_recomputable':any(r['source']==s for r in qrows),
   'rule':'Do not mix ATpos/ATneg/QCpool when exported in-house IDs lack pool type; retain published method QC statistics'} for s in ['V01','V02']])
 repeat=[]
 for source in ['V01','V02']:
  for fk in featurelist:
   groups=[v for (s,r,d,f),v in duplicates.items() if s==source and f==fk and len(v)>1]
   diffs=[max(v)-min(v) for v in groups]
   repeat.append({'source':source,'feature_key':fk,'same_RID_date_multiple_record_groups':len(groups),
    'median_abs_log2_range':float(np.median(diffs)) if diffs else np.nan,'max_abs_log2_range':max(diffs) if diffs else np.nan,
    'not_assumed_independent_technical_replicates':True,'not_used_to_calibrate_cross_assay_error_correlation':True})
 csvout(RESULTS/'TARGETED_DUPLICATE_VALUE_DIAGNOSTICS.csv',repeat)
 log(f'Targeted log ratios cached privately after plan freeze; {nflags} selected positive cells; no source re-aggregation by outcome')

if __name__=='__main__':run()
