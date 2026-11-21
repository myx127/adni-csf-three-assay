from config import *
import zipfile

def run():
 sap=RESULTS/'LEGACY_VERIFICATION_SAP_FROZEN.json';first=RESULTS/'FIRST_LEGACY_QUANTITATIVE_ACCESS.json'
 if first.exists():
  assert json.loads(first.read_text(encoding='utf-8'))['plan_sha256_before_association']==sha(sap)
 else:
  save_json(first,{'plan_sha256_before_association':sha(sap),'time':datetime.now().astimezone().isoformat(),
    'purpose':'Fixed equal mean of two released normalized log transitions per peptide; no second log',
    'legacy_association_results_not_previously_seen':True})
 plan=json.loads(sap.read_text(encoding='utf-8'));manifest=pd.read_csv(RESULTS/'LEGACY_INPUT_MANIFEST.csv').set_index('input_id')
 assert sha(manifest.loc['L04','path'])==manifest.loc['L04','sha256']
 fields=[x for p in plan['peptides'] for x in p['transition_fields']]
 with zipfile.ZipFile(manifest.loc['L04','path']) as z:
  with z.open('csfproteomics/adni_csfproteomics.csv') as f:raw=pd.read_csv(f,usecols=fields)
 frame=pd.read_csv(RESTRICTED/'LEGACY_L04_CLINICAL_DESIGN_RESTRICTED.csv');selected=raw.iloc[frame.source_row.to_numpy(int)-2]
 cols=[];summ=[]
 for p in plan['peptides']:
  two=selected[p['transition_fields']].apply(pd.to_numeric,errors='coerce').to_numpy(float)
  y=two.mean(axis=1);y[~np.isfinite(two).all(axis=1)]=np.nan;cols.append(y)
  pair=two[np.isfinite(two).all(axis=1)]
  summ.append({'accession':p['accession'],'gene':p['gene'],'peptide':p['peptide'],'n_selected_keys':len(frame),
    'n_both_finite':len(pair),'transition_difference_SD_descriptive':float(np.std(pair[:,0]-pair[:,1],ddof=1)),
    'transition_difference_median_descriptive':float(np.median(pair[:,0]-pair[:,1])),
    'are_two_fragment_ions_of_same_peptide_not_independent_protein_indicators':True,
    'raw_LH_ratio_or_QC_pool_not_reconstructed_from_processed_log_release':True})
 np.savez_compressed(RESTRICTED/'LEGACY_PEPTIDE_CACHE_RESTRICTED.npz',keys=frame.key.to_numpy(str),
   peptides=np.array([p['peptide'] for p in plan['peptides']]),accessions=np.array([p['accession'] for p in plan['peptides']]),Y=np.column_stack(cols))
 csvout(RESULTS/'LEGACY_TRANSITION_PAIR_DIAGNOSTICS.csv',summ)
 log('Legacy normalized-log peptide cache loaded under frozen rules; association results not yet fitted')

if __name__=='__main__':run()
