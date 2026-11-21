from config import *
from collections import defaultdict
from openpyxl import load_workbook
import re,tarfile

def normid(v):return str(v).strip().upper() if pd.notna(v) else ''
def run():
 legacy=pd.read_csv(S6/'results/LEGACY_INPUT_MANIFEST.csv').set_index('input_id');src=pd.read_csv(legacy.loc['L01','path']);exp=pd.read_csv(PRIVATE/'MRM2014_CSV_EXPORT_RESTRICTED.csv');info=pd.read_csv(PRIVATE/'MRM2014_SAMPLEINFO_RESTRICTED.csv')
 info['RID_NUM']=pd.to_numeric(info.RID,errors='coerce');info['_ref']=info.External_Reference.map(normid);exp['_ref']=exp.iloc[:,0].map(normid)
 valid=info.loc[info.RID_NUM.gt(0)&info.Study.eq('ADNI1')].copy();groups={k:g for k,g in valid.groupby('_ref')};records=[]
 for i,row in exp.iterrows():
  if not row['_ref']:continue
  part=groups.get(row['_ref'],pd.DataFrame());rids=set(pd.to_numeric(part.get('RID_NUM',pd.Series(dtype=float)),errors='coerce').dropna().astype(int));vis=set(part.get('Viscode',pd.Series(dtype=str)).dropna().astype(str));tech=part[['Processing_Day','Depletion_Column','Injection_Order']].drop_duplicates() if len(part) else []
  records.append({'export_row':i,'lab_reference':row['_ref'],'RID':next(iter(rids)) if len(rids)==1 else np.nan,'VISCODE':next(iter(vis)) if len(vis)==1 else '',
   'n_metadata_matches':len(part),'RID_unambiguous':len(rids)==1,'visit_unambiguous':len(vis)==1,'technical_assignment_unambiguous':len(tech)==1,
   'Processing_Day':part.Processing_Day.iloc[0] if len(tech)==1 else np.nan,'Depletion_Column':part.Depletion_Column.iloc[0] if len(tech)==1 else np.nan})
 links=pd.DataFrame(records);csvout(PRIVATE/'MRM2014_EXPORT_SAMPLEINFO_LINK_RESTRICTED.csv',links)
 # Whole-measurement vectors verify release identity, not selection by candidate associations.
 columns=[c for c in exp.columns if c!='_ref' and c.replace('-','_',1) in src.columns];assert len(columns)>=300
 matrix=exp[columns].apply(pd.to_numeric,errors='coerce').to_numpy(float);original=src[[c.replace('-','_',1) for c in columns]].apply(pd.to_numeric,errors='coerce').to_numpy(float)
 identity=[]
 for i,row in src.iterrows():
  cand=links.loc[links.RID.eq(row.RID)&links.VISCODE.eq(str(row.VISCODE)),'export_row'].astype(int).tolist();matches=[];best=np.inf
  for j in cand:
   a=original[i];b=matrix[j];validpair=np.isfinite(a)&np.isfinite(b);same_missing=np.array_equal(np.isfinite(a),np.isfinite(b))
   distance=np.max(np.abs(a[validpair]-b[validpair])) if validpair.any() else np.inf;best=min(best,distance)
   if same_missing and validpair.sum()>=200 and distance<=1e-5:matches.append(j)
  identity.append({'source_row':i+2,'RID':row.RID,'VISCODE':row.VISCODE,'n_export_same_RID_visit_candidates':len(cand),'n_whole_vector_matches':len(matches),
    'export_row_match':matches[0] if len(matches)==1 else np.nan,'best_whole_vector_max_abs_difference':best})
 idf=pd.DataFrame(identity);csvout(PRIVATE/'MRM2014_RELEASE_IDENTITY_RESTRICTED.csv',idf)
 # Dates on a different measurement release are visit-level candidates only.
 master=pd.read_csv(INPUTS['I10']['path'],usecols=['RID','VISCODE','DRAWDTE'],dtype=str);master['parsed_date']=pd.to_datetime(master.DRAWDTE,errors='coerce',format='mixed').dt.strftime('%Y-%m-%d')
 md=master.groupby(['RID','VISCODE']).parsed_date.apply(lambda x:sorted(set(x.dropna()))).to_dict();old=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv');oldmap=old.set_index('RID');dates=[]
 for i,r in src.drop_duplicates(['RID','VISCODE']).iterrows():
  candidates=md.get((str(int(r.RID)),str(r.VISCODE)),[]);stage5=oldmap.loc[r.RID] if r.RID in oldmap.index else None
  dates.append({'RID':r.RID,'VISCODE':r.VISCODE,'n_unique_master_visit_drawdates':len(candidates),'candidate_drawdate':candidates[0] if len(candidates)==1 else '',
    'in_Stage5':stage5 is not None,'Stage5_visit_matches':bool(stage5 is not None and str(stage5.NULISA_VISCODE2)==str(r.VISCODE)),
    'candidate_equals_Stage5_date':bool(stage5 is not None and len(candidates)==1 and str(stage5.EXAMDATE)==candidates[0]),
    'is_verified_MRM_specimen_collection_date':False})
 datef=pd.DataFrame(dates);csvout(PRIVATE/'MRM2014_CANDIDATE_DATE_CROSSWALK_RESTRICTED.csv',datef)
 book=load_workbook(legacy.loc['L02','path'],read_only=True,data_only=True);headers=list(next(book['log(Peptide Intensity)'].iter_rows(values_only=True)));book.close()
 refset=set(info._ref);limsset=set(info.LIMS_ID.dropna().map(normid));external_matches=sum(normid(h) in refset for h in headers);lims_matches=sum(normid(h) in limsset for h in headers)
 public=[
  ('L01_release_rows',len(src)),('L01_unique_RID',src.RID.nunique()),('L01_baseline_rows',int(src.VISCODE.eq('bl').sum())),
  ('L02_metadata_ADNI1_rows',len(valid)),('L02_export_nonempty_rows',len(links)),('L02_export_RID_unambiguous_rows',int(links.RID_unambiguous.sum())),
  ('L02_export_visit_unambiguous_rows',int(links.visit_unambiguous.sum())),('L02_export_technical_unambiguous_rows',int(links.technical_assignment_unambiguous.sum())),
  ('whole_vector_features_compared',len(columns)),('L01_rows_with_at_least_one_identical_vector',int(idf.n_whole_vector_matches.gt(0).sum())),
  ('L01_rows_with_unique_identical_vector',int(idf.n_whole_vector_matches.eq(1).sum())),('L01_rows_with_multiple_identical_vectors',int(idf.n_whole_vector_matches.gt(1).sum())),
  ('peptide_matrix_headers_matching_External_Reference',external_matches),('peptide_matrix_headers_matching_LIMS_ID',lims_matches),
  ('RID_visit_with_unique_MASTER_drawdate',int(datef.n_unique_master_visit_drawdates.eq(1).sum())),('RID_visit_with_multiple_MASTER_drawdates',int(datef.n_unique_master_visit_drawdates.gt(1).sum())),
  ('RID_in_Stage5',int(datef.in_Stage5.sum())),('same_visit_in_Stage5',int(datef.Stage5_visit_matches.sum())),
  ('unique_visit_candidate_drawdate_equals_Stage5',int(datef.candidate_equals_Stage5_date.sum())),('verified_MRM_specimen_collection_dates',0)]
 csvout(RESULTS/'MRM2014_CROSSWALK_AUDIT_COUNTS.csv',[{'audit_item':k,'count':int(v)} for k,v in public])
 csvout(RESULTS/'MRM2014_LINKAGE_STATUS.csv',[
  {'layer':'Release family','status':'PARTLY_OR_FULLY_MATCHED_BY_VECTORS_SEE_COUNTS','evidence':'All available shared peptide features compared within identity/visit; no candidate-P effects fitted','remaining':'L01/L02 do not count as separate experiments'},
  {'layer':'Metadata identity','status':'RESOLVED_WHERE_LAB_REFERENCE_HAS_UNIQUE_RID_VISIT','evidence':'CSV-Export first column matched Sample Info External_Reference','remaining':'Some repeated lab references may map to multiple technical injections'},
  {'layer':'Technical assignment','status':'ROW_SPECIFIC_AMBIGUITY_RETAINED','evidence':f'{int(links.technical_assignment_unambiguous.sum())}/{len(links)} export rows have one processing/injection combination','remaining':'Need per-export-row or LIMS/injection mapping for ambiguous duplicates'},
  {'layer':'Exact specimen date','status':'NOT_CLOSED','evidence':'No collection date in 2014 source or specimen-ID bridge to a dated CSF record','remaining':'MASTER DRAWDTE via RID/visit is an inferred visit-level candidate, not verified MRM aliquot date'},
  {'layer':'Formal relative-pattern replication','status':'NOT_EXECUTED_BY_THIS_AUDIT','evidence':'This stage audits metadata, not new candidate associations','remaining':'Resolve exact dates/technical assignment then re-audit frozen design gates before any new test'}])
 log('2014 release/technical/date crosswalk audited; inferred same-visit dates not promoted to exact MRM specimen dates')
 print(pd.DataFrame(public,columns=['item','count']).to_string(index=False))

if __name__=='__main__':run()
