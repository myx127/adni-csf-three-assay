from config import *
from designs import make
from collections import Counter
from profile_release_precision import significant
import re

def run():
 sources=pd.read_csv(S6/'results/LEGACY_INPUT_MANIFEST.csv').set_index('input_id');rawrows=list(csv.DictReader(Path(sources.loc['L01','path']).open(encoding='utf-8-sig')));src=pd.DataFrame(rawrows);exp=pd.read_csv(PRIVATE/'MRM2014_CSV_EXPORT_RESTRICTED.csv');link=pd.read_csv(PRIVATE/'MRM2014_EXPORT_SAMPLEINFO_LINK_RESTRICTED.csv');info=pd.read_csv(PRIVATE/'MRM2014_SAMPLEINFO_RESTRICTED.csv')
 cols=[c for c in exp if c.replace('-','_',1) in src];decimals=Counter()
 for c in cols:
  for v in src[c.replace('-','_',1)]:
   m=re.fullmatch(r'[+-]?\d+(?:\.(\d+))?',v)
   if m:decimals[len(m[1] or '')]+=1
 quantized=src[[c.replace('-','_',1) for c in cols]].apply(pd.to_numeric,errors='coerce').to_numpy(float)
 export_values=exp[cols].apply(pd.to_numeric,errors='coerce').to_numpy(float)
 rounding_grid_error=float(np.nanmax(abs(quantized-significant(quantized,6))));export_grid_error=float(np.nanmax(abs(export_values-significant(export_values,6))))
 assert rounding_grid_error<1e-10
 save_json(RESULTS/'MRM2014_IDENTITY_ROUNDING_RULE.json',{'before_any_new_candidate_pathology_fit':True,'source_lexical_decimal_counts':dict(decimals),
  'rule':'Require all320 peptide values within one unit of the last6-significant-digit position, preserving missingness, sameRID/visit and exactly one matching export record. No affine fit used for final matching.',
  'L01_max_distance_to_6significant_digit_grid':rounding_grid_error,'L02_export_max_distance_to_6significant_digit_grid':export_grid_error,
  'reason':'L01 uses6 significant digits. Exact6-digit rounding still differs in0-6/320 features per row; observed raw discrepancies reach0.541 last-place units, compatible with a small additional numeric/storage difference. Final provenance tolerance is one published last-place unit, not tuned to pathology effects. Not claimed to be byte/numerically exact equality.',
  'no_new_candidate_association_tests':True})
 numeric=src.copy();numeric.RID=pd.to_numeric(numeric.RID);A=numeric[[c.replace('-','_',1) for c in cols]].apply(pd.to_numeric,errors='coerce').to_numpy(float);B=exp[cols].apply(pd.to_numeric,errors='coerce').to_numpy(float)
 rows=[]
 for i,r in numeric.iterrows():
  candidates=link.loc[link.RID.eq(r.RID)&link.VISCODE.eq(str(r.VISCODE))];matches=[]
  for c in candidates.itertuples():
   y=B[c.export_row];finite=np.isfinite(A[i])&np.isfinite(y);same_missing=np.array_equal(np.isfinite(A[i]),np.isfinite(y))
   units=10.**(np.floor(np.log10(np.maximum(abs(A[i,finite]),1e-300)))-5)
   ok=same_missing and bool((abs(A[i,finite]-y[finite])<=units+1e-10).all())
   if ok:matches.append(c)
  assert len(matches)==1,(i,len(matches));m=matches[0];part=info.loc[info.External_Reference.astype(str).str.upper().eq(m.lab_reference.upper())]
  assert len(part)==1
  s=part.iloc[0];rows.append({'source_row':i+2,'RID':r.RID,'VISCODE':r.VISCODE,'export_row':m.export_row,'LIMS_ID':s.LIMS_ID,'External_Reference':s.External_Reference,
    'Processing_Day':s.Processing_Day,'Depletion_Column':s.Depletion_Column,'Injection_Order':s.Injection_Order,'Sample_Duplicated':s['Sample Duplicated'],
    'all320_within_one_source_precision_unit':True,'has_source_MRM_collection_date':False})
 matched=pd.DataFrame(rows);csvout(PRIVATE/'MRM2014_RESOLVED_TECHNICAL_CROSSWALK_RESTRICTED.csv',matched)
 csvout(RESULTS/'MRM2014_FINAL_IDENTITY_STATUS.csv',[{'item':'L01 records uniquely matched to L02 all320-feature vector within source precision','count':len(matched),'status':'CLOSED_WITH_DISCLOSED_PRECISION_TOLERANCE'},
  {'item':'L01 records linked to one LIMS_ID and processing/injection record','count':len(matched),'status':'CLOSED'},
  {'item':'distinct Processing_Day','count':matched.Processing_Day.nunique(),'status':'METADATA_RECOVERED'},
  {'item':'distinct Depletion_Column','count':matched.Depletion_Column.nunique(),'status':'METADATA_RECOVERED'},
  {'item':'verified source-specific MRM specimen dates','count':0,'status':'NOT_CLOSED'}])
 status=pd.read_csv(RESULTS/'MRM2014_LINKAGE_STATUS.csv')
 status.loc[status.layer.eq('Release family'),['status','evidence']]=['CLOSED_WITH_PUBLICATION_PRECISION_TOLERANCE','303 source rows, all320 features within one6-significant-digit last-place unit, unique within sameRID/visit; not exact raw-value equality']
 status.loc[status.layer.eq('Metadata identity'),['status','evidence','remaining']]=['CLOSED_FOR303_SOURCE_ROWS','One export lab reference to one LIMS/RID/visit technical record after full-vector identity check','No new independent experiment created by two releases']
 status.loc[status.layer.eq('Technical assignment'),['status','evidence','remaining']]=['CLOSED_FOR303_SOURCE_ROWS','Processing day, depletion column and injection order recovered per source row','Date bridge and actual-design calibration still required; depletion-column dependence not validated']
 csvout(RESULTS/'MRM2014_LINKAGE_STATUS.csv',status)
 # Conditional design geometry only, without reading legacy candidate values for effects.
 old=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv');matched=matched.sort_values(['Processing_Day','Injection_Order','source_row']).drop_duplicates(['RID','VISCODE'])
 base=old.merge(matched,on='RID',suffixes=('','_MRM'),how='inner');base=base.loc[base.NULISA_VISCODE2.eq(base.VISCODE)]
 base['MRM_DAY']=base.Processing_Day.astype(str);base['MRM_DEPLETION']=base.Depletion_Column.astype(str)
 rawcache=np.load(S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz');index={k:i for i,k in enumerate(rawcache['keys'])};meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession');screen=[]
 for gene,acc,peps in [('NPTX1','Q15818',['FQLTFPLR','LENLEQYSR']),('KLK6','Q92876',['ESSQEQSSVVR','LSELIQPLPLER','YTNWIQK'])]:
  for pep in peps:
   f=gene+'_'+pep;value_by_source=pd.to_numeric(src[f],errors='coerce').to_numpy();j=list(rawcache['accessions']).index(acc);sj=list(rawcache['soma_features']).index(meta.loc[acc,'soma_features']);ix=np.array([index[k] for k in base.key])
   mask=np.isfinite(value_by_source[base.source_row.to_numpy(int)-2])&np.isfinite(rawcache['N'][ix,j])&np.isfinite(rawcache['S'][ix,sj])&np.isfinite(base[['P_LOG','age_month_precision','sex_code2']]).all(axis=1).to_numpy()
   d=make(base,mask,tech=('NULISA_PlateID','SOMASCAN_PlateId','MRM_DAY'),additional_factors=('ptau_batch','MRM_DEPLETION'))
   screen.append({'accession':acc,'gene':gene,'peptide':pep,'date_status':'INFERRED_FROM_OTHER_ASSAY_SAME_VISIT_NOT_VERIFIED_MRM_DATE','n_same_visit_before_validity':len(base),
     **(d['summary'] if d else {'n':int(mask.sum())}),'geometry_passes_S6_n80_df40_G6':bool(d is not None and d['n']>=80 and d['df']>=40 and min(len(x) for x in d['onehots'])>=6),
     'strict_same_day_replication_eligible':False,'fresh_candidate_association_fitted':False,'remaining':'Source specimen-date bridge and actual2014-design inferential calibration; depletion-column dependence unresolved'})
 csvout(RESULTS/'MRM2014_INFERRED_VISIT_DESIGN_SCREEN.csv',screen)
 csvout(PRIVATE/'MRM2014_INFERRED_VISIT_DESIGN_RESTRICTED.csv',base)
 log('All3032014 rows uniquely linked using all320 values within one6-digit last-place unit; source numeric discrepancy disclosed; technical mapping closed, exact specimen dates unresolved; no effects fitted')
 print(pd.DataFrame(screen)[['gene','peptide','n','rank_X','residual_df','geometry_passes_S6_n80_df40_G6']].to_string(index=False))

if __name__=='__main__':run()
