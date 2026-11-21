from config import *
from independent_analysis import design
import re,zipfile,io
from openpyxl import load_workbook

ID_PATTERN=re.compile(r'LIMS|SAMP|SPEC|ALIQUOT|BARCODE|EXTERNAL|REFERENCE|ACCNO|LAB.?ID|PROTEOMIC|GUSPEC',re.I)
DATE_PATTERN=re.compile(r'DATE|DTE|COLLECT',re.I)
def norm_id(x):return str(x).strip().upper()

def run():
    cross=pd.read_csv(S7/'00_restricted/MRM2014_RESOLVED_TECHNICAL_CROSSWALK_RESTRICTED.csv',dtype=str)
    ids={field:set(cross[field].map(norm_id)) for field in ['LIMS_ID','External_Reference']}
    tables=[];hits=[];inputs=[];schemas=[]
    def scan(df,label,path,member=''):
        idcols=[c for c in df if ID_PATTERN.search(str(c))]
        dates=[c for c in df if DATE_PATTERN.search(str(c))]
        schemas.append({'table':label,'path':str(path),'member':member,'n_rows':len(df),'id_fields':'|'.join(map(str,idcols)),'date_fields':'|'.join(map(str,dates)),'all_fields_json':json.dumps(list(map(str,df.columns)),ensure_ascii=False)})
        for c in idcols:
            values=df[c].map(norm_id)
            for field,tokens in ids.items():
                mask=values.isin(tokens)
                tables.append({'table':label,'field':c,'target_field':field,'n_matched_rows':int(mask.sum()),'n_unique_matched_ids':int(values[mask].nunique()),'date_fields':'|'.join(dates),'specimen_bridge_verified':False})
                for index in df.index[mask]:
                    row=df.loc[index];hit={'table':label,'row_index':int(index),'field':c,'target_field':field,'matched_id':norm_id(row[c]),'path':str(path),'member':member}
                    for k in list(dict.fromkeys([*dates,*[x for x in df if re.fullmatch('RID|PTID|VISCODE|VISCODE2',str(x),re.I)]])):hit[str(k)]=row[k]
                    hits.append(hit)
    for p in sorted(RAW.rglob('*.csv')):
        if '_audit' in str(p):continue
        try:
            cols=pd.read_csv(p,nrows=0).columns.tolist()
            selected=[c for c in cols if ID_PATTERN.search(c) or DATE_PATTERN.search(c) or c.upper() in ['RID','PTID','VISCODE','VISCODE2']]
            if not any(ID_PATTERN.search(c) for c in cols):continue
            df=pd.read_csv(p,usecols=selected,dtype=str,keep_default_na=False,low_memory=False)
            scan(df,str(p.relative_to(RAW)),p)
            inputs.append({'path':str(p),'member':'','bytes':p.stat().st_size,'sha256':sha(p)})
        except Exception as e:
            raise RuntimeError(f'Failed metadata scan: {p}: {e}') from e
    manifests=pd.read_csv(S6/'results/LEGACY_INPUT_MANIFEST.csv').set_index('input_id')
    xlsx=Path(manifests.loc['L02','path']);book=load_workbook(xlsx,read_only=True,data_only=True)
    for ws in book:
        header=next(ws.iter_rows(values_only=True));names=[str(c) for c in header]
        # Only metadata sheets. Quantitative matrices are not candidate-date evidence.
        if any(ID_PATTERN.search(c) for c in names) and ws.max_column<50:
            df=pd.DataFrame(ws.iter_rows(min_row=2,values_only=True),columns=names).fillna('')
            scan(df,'2014_WORKBOOK/'+ws.title,xlsx,ws.title)
    book.close();inputs.append({'path':str(xlsx),'member':'','bytes':xlsx.stat().st_size,'sha256':sha(xlsx)})
    l04=Path(manifests.loc['L04','path'])
    with zipfile.ZipFile(l04) as z:
        # Both supporting workbooks are new, source-level metadata opportunities.
        for member in z.namelist():
            if member.endswith('.xlsx'):
                book=load_workbook(io.BytesIO(z.read(member)),read_only=True,data_only=True)
                for ws in book:
                    header=next(ws.iter_rows(values_only=True));names=[str(c) for c in header]
                    if any(ID_PATTERN.search(c) for c in names):
                        indices=[i for i,c in enumerate(names) if ID_PATTERN.search(c) or DATE_PATTERN.search(c) or c.upper() in ['RID','PTID','VISCODE','VISCODE2']]
                        data=[[row[i] if i<len(row) else '' for i in indices] for row in ws.iter_rows(min_row=2,values_only=True)]
                        df=pd.DataFrame(data,columns=[names[i] for i in indices]).fillna('')
                        scan(df,'2019_ARCHIVE/'+Path(member).name+'/'+ws.title,l04,member+'::'+ws.title)
                book.close()
        a=pd.read_csv(z.open('csfproteomics/adni_csfproteomics.csv'),usecols=['RID','SAMPLE.NAME','EXAMDATE'])
        # A shared-looking CSF serial is not a specimen identity. Assess namespace collisions explicitly.
        a['base_name']=a['SAMPLE.NAME'].astype(str).str.replace(r'-\d+$','',regex=True)
        overlap=set(a.base_name)&ids['LIMS_ID']
        test=a.merge(cross[['RID','LIMS_ID']].assign(RID=lambda d:pd.to_numeric(d.RID)),left_on='base_name',right_on='LIMS_ID',suffixes=('_2019','_2014'))
        namespace={'2014_LIMS_overlap_2019_sample_base':len(overlap),'joined_rows':len(test),'RID_agree_rows':int(test.RID_2019.eq(test.RID_2014).sum()) if len(test) else 0,'status':'NO_IDENTITY_INFERENCE_FROM_PREFIX_SERIAL'}
    csvout(RESULTS/'METADATA_ID_SCAN_RESULTS.csv',tables)
    csvout(RESULTS/'METADATA_TABLE_SCHEMAS.csv',schemas)
    csvout(RESULTS/'METADATA_SOURCE_MANIFEST.csv',inputs)
    csvout(PRIVATE/'METADATA_EXACT_ID_HITS_RESTRICTED.csv',hits)
    save_json(RESULTS/'SAMPLE_ID_NAMESPACE_AUDIT.json',namespace)
    # New, separate intersection/rank audit for three and four assays.
    base=pd.read_csv(S7/'00_restricted/MRM2014_INFERRED_VISIT_DESIGN_RESTRICTED.csv')
    raw=np.load(S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz')
    ixmap={k:i for i,k in enumerate(raw['keys'])};ix=np.array([ixmap[k] for k in base.key])
    meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
    source=pd.read_csv(manifests.loc['L01','path'])
    geometry=[];selection=[]
    for acc,gene in CANDIDATES.items():
        if gene not in PEPTIDES:continue
        j=list(raw['accessions']).index(acc);sj=list(raw['soma_features']).index(meta.loc[acc,'soma_features'])
        for pep in PEPTIDES[gene]:
            yy=pd.to_numeric(source[gene+'_'+pep],errors='coerce').to_numpy()[base.source_row.to_numpy(int)-2]
            for assays,tech in [('N_S_MRM',('NULISA_PlateID','SOMASCAN_PlateId','MRM_DAY')),('N_S_TMT_MRM',('NULISA_PlateID','SOMASCAN_PlateId','TMT_plex','MRM_DAY'))]:
                valid=np.isfinite(yy)&np.isfinite(raw['N'][ix,j])&np.isfinite(raw['S'][ix,sj])&np.isfinite(base[['P_LOG','age_month_precision','sex_code2']]).all(axis=1).to_numpy()
                if assays=='N_S_TMT_MRM':valid&=np.isfinite(raw['T_GIS'][ix,j])
                f=base.loc[valid];d=design(f,tech=tech,factors=('ptau_batch','MRM_DEPLETION'))
                row={'accession':acc,'gene':gene,'peptide':pep,'assay_design':assays,'n_same_visit_before_validity':len(base),'n':len(f),'date_status':'VISIT_LEVEL_INFERENCE_NOT_MRM_SPECIMEN_VERIFIED','strict_validation_eligible':False,'candidate_pathology_effect_fitted':False,'depletion_column_G':f.MRM_DEPLETION.nunique(),'fixed_depletion_FE_does_not_remove_arbitrary_column_dependence':True}
                if d:
                    hat=np.sum(d['Q']**2,axis=1);partial=d['a']**2/(d['a']@d['a'])
                    row.update({'rank_X':d['rank'],'residual_df':d['df'],'P_residual_SS':d['ss'],'max_hat_leverage':hat.max(),'n_hat_ge_0999':int((hat>=.999).sum()),'max_P_row_information_share':partial.max()})
                    for c,lab in zip(tech,d['groups']):
                        shares=np.bincount(lab,weights=partial)
                        row[c+'_G']=int(lab.max())+1;row[c+'_P_effective_clusters']=1/float(shares@shares);row[c+'_max_P_cluster_share']=shares.max()
                    row['geometry_n80_df40_G6']=bool(d['n']>=80 and d['df']>=40 and all(int(a.max())+1>=6 for a in d['groups']))
                    row['inference_calibrated_for_2014_design']=False
                geometry.append(row)
                selection.extend([{'gene':gene,'peptide':pep,'assay_design':assays,'key':k,'source_row':sr} for k,sr in zip(f.key,f.source_row)])
    csvout(RESULTS/'MRM2014_THREE_FOUR_ASSAY_DESIGN_AUDIT.csv',geometry)
    csvout(PRIVATE/'MRM2014_CONDITIONAL_DESIGN_SELECTION_RESTRICTED.csv',selection)
    counts=[{'item':'metadata_tables_scanned','count':len(schemas)},{'item':'ID_field_by_target_comparisons','count':len(tables)},{'item':'exact_identifier_hit_rows_including_source_redundancy','count':len(hits)},{'item':'verified_MRM2014_specimen_dates','count':0},{'item':'conditional_three_four_assay_designs','count':len(geometry)},{'item':'new_2014_pathology_associations_fitted','count':0}]
    csvout(RESULTS/'MRM2014_STAGE9_DATE_BRIDGE_STATUS.csv',counts)
    log(f'Metadata audit: {len(schemas)} tables; {len(hits)} exact ID-hit rows retained privately; 10 separate three/four-assay geometries audited. Direct date bridges require adjudication.')
    print(pd.DataFrame(geometry)[['gene','peptide','assay_design','n','rank_X','residual_df','max_hat_leverage']].to_string(index=False))
    print('Matched tables:',pd.DataFrame(tables).loc[lambda d:d.n_matched_rows.gt(0),['table','field','target_field','n_matched_rows','date_fields']].to_string(index=False))

if __name__=='__main__':
    assert_unsealed()
    run()
