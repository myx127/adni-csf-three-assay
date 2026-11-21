from config import *

def run():
    candidate=pd.read_csv(S7/'00_restricted/MRM2014_CANDIDATE_DATE_CROSSWALK_RESTRICTED.csv')
    dictionary=RAW/'Data dictionary/DATADIC_15May2026.csv';dic=pd.read_csv(dictionary,dtype=str,low_memory=False)
    selected=dic.loc[dic.TBLNAME.str.upper().isin(['LOCLAB','CSFMETH'])&dic.FLDNAME.str.upper().isin(['SAMPLECOLL','DONE','VISDATE'])]
    csvout(RESULTS/'COLLECTION_CRF_FIELD_DEFINITIONS.csv',selected)
    details=[];summary=[]
    for name,flag in [('LOCLAB_15May2026.csv','SAMPLECOLL'),('CSFMETH_15May2026.csv','DONE')]:
        p=RAW/'CSF'/name;source=pd.read_csv(p,dtype=str)
        source.RID=pd.to_numeric(source.RID,errors='coerce');source['VISCODE2']=source.VISCODE2.astype(str).str.lower()
        source['VISDATE_PARSED']=pd.to_datetime(source.VISDATE,errors='coerce').dt.strftime('%Y-%m-%d')
        for row in candidate.itertuples():
            sub=source.loc[source.RID.eq(row.RID)&source.VISCODE2.eq(str(row.VISCODE).lower())]
            positive=sub.loc[pd.to_numeric(sub[flag],errors='coerce').eq(1)]
            alldates=set(sub.VISDATE_PARSED.dropna());yesdates=set(positive.VISDATE_PARSED.dropna());target=str(row.candidate_drawdate)[:10]
            details.append({'RID':row.RID,'VISCODE':row.VISCODE,'table':name,'n_matching_visit_rows':len(sub),'n_yes_collection_rows':len(positive),'n_blank_or_unavailable_collection_flag':int(pd.to_numeric(sub[flag],errors='coerce').isna().sum()),'n_explicit_no_collection_rows':int(pd.to_numeric(sub[flag],errors='coerce').eq(0).sum()),'n_unique_VISDATE_all':len(alldates),'n_unique_VISDATE_yes':len(yesdates),'MASTER_candidate_date':target,'all_VISDATEs':'|'.join(sorted(alldates)),'yes_VISDATEs':'|'.join(sorted(yesdates)),'MASTER_matches_any_visit_date':target in alldates,'MASTER_matches_yes_VISDATE':target in yesdates,'source_field_may_use_registry_fallback':True,'not_verified_MRM_specimen_date':True})
    d=pd.DataFrame(details);csvout(PRIVATE/'MRM2014_COLLECTION_CRF_TRIANGULATION_RESTRICTED.csv',d)
    for table,part in d.groupby('table'):
        summary.append({'table':table,'n_MRM_RID_visit':len(part),'n_with_any_matching_CRF_row':int(part.n_matching_visit_rows.gt(0).sum()),'n_with_yes_collection_row':int(part.n_yes_collection_rows.gt(0).sum()),'n_with_blank_or_unavailable_collection_flag':int(part.n_blank_or_unavailable_collection_flag.gt(0).sum()),'n_with_explicit_no_collection_row':int(part.n_explicit_no_collection_rows.gt(0).sum()),'n_MASTER_matches_any_VISDATE':int(part.MASTER_matches_any_visit_date.sum()),'n_MASTER_matches_yes_VISDATE':int(part.MASTER_matches_yes_VISDATE.sum()),'n_with_matching_rows_but_no_MASTER_date_match':int((part.n_matching_visit_rows.gt(0)&~part.MASTER_matches_any_visit_date).sum()),'date_definition':'Assessment date or registry fallback; not guaranteed specimen collection time','missing_flag_does_not_mean_no_collection':True,'specimen_ID_link_present':False,'verified_MRM_date_bridge':False})
    csvout(RESULTS/'MRM2014_COLLECTION_CRF_DATE_TRIANGULATION.csv',summary)
    save_json(RESULTS/'DATE_BRIDGE_ADJUDICATION.json',{'time':now(),'exact_ID_hits_outside_own2014_workbook':0,'verified_MRM_specimen_dates':0,'decision':'NOT_CLOSED','reason':'Only own2014 Sample Info matches target LIMS/ExternalReference; 2019 namespaces do not overlap. Additional clinical VISDATE may fall back to registry date, so RID/visit triangulation cannot establish MRM aliquot collection date.','2014_pathology_associations_fitted':0,'additional_inputs':[{'path':str(p),'sha256':sha(p)} for p in [dictionary,RAW/'CSF/LOCLAB_15May2026.csv',RAW/'CSF/CSFMETH_15May2026.csv']]})
    log('Collection-CRF triangulation complete; dictionary confirms VISDATE may fall back to registry. No specimen-level MRM date bridge established.')
    print(pd.DataFrame(summary).to_string(index=False))

if __name__=='__main__':
    assert_unsealed()
    run()
