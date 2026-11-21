from config import *
from scipy.stats import spearmanr
import zipfile,io,re,pickle

def pair_metrics(x,y,seed):
    ok=np.isfinite(x)&np.isfinite(y);x=np.asarray(x)[ok];y=np.asarray(y)[ok];n=len(x)
    if n<3:return {'n_complete_pairs':n}
    def corr(a,b):return float(np.corrcoef(a,b)[0,1]) if np.std(a)>0 and np.std(b)>0 else np.nan
    r=corr(x,y);vx=np.var(x,ddof=1);vy=np.var(y,ddof=1);cv=np.cov(x,y,ddof=1)[0,1]
    ccc=2*cv/(vx+vy+(x.mean()-y.mean())**2);diff=y-x
    rng=np.random.default_rng(seed);boot=[]
    for i in range(4000):
        pick=rng.integers(0,n,n);boot.append(corr(x[pick],y[pick]))
    jack=[corr(np.delete(x,i),np.delete(y,i)) for i in range(n)]
    lo,hi=np.nanquantile(boot,[.025,.975])
    return {'n_complete_pairs':n,'pearson_r':r,'pearson_pair_bootstrap95_low':lo,'pearson_pair_bootstrap95_high':hi,'bootstrap_reps':4000,'bootstrap_seed':str(seed),
            'spearman_r':spearmanr(x,y).statistic,'Lin_CCC_descriptive':ccc,'mean_second_minus_first':diff.mean(),'SD_pair_difference':diff.std(ddof=1),
            'RMS_pair_difference':np.sqrt(np.mean(diff**2)),'RMS_pair_difference_per_pooled_SD':np.sqrt(np.mean(diff**2))/np.std(np.r_[x,y],ddof=1),
            'repeat_error_SD_equal_independent_error_working':diff.std(ddof=1)/np.sqrt(2),'pearson_leave_one_pair_min':min(jack),'pearson_leave_one_pair_max':max(jack)}

def run():
    manifests=pd.read_csv(S6/'results/LEGACY_INPUT_MANIFEST.csv').set_index('input_id')
    cross=pd.read_csv(S7/'00_restricted/MRM2014_RESOLVED_TECHNICAL_CROSSWALK_RESTRICTED.csv')
    raw14=pd.read_csv(manifests.loc['L01','path']);repeat=[];pair_rows=[];identity=[]
    groups14=[g.sort_values(['Processing_Day','Injection_Order','source_row']) for _,g in cross.groupby(['RID','VISCODE']) if len(g)==2]
    assert len(groups14)==16 and all(g.Sample_Duplicated.astype(str).eq('Yes').all() for g in groups14)
    exp=pd.read_csv(S7/'00_restricted/MRM2014_CSV_EXPORT_RESTRICTED.csv')
    columns=[c for c in exp if c.replace('-','_',1) in raw14];assert len(columns)==320
    for row in cross.itertuples():
        a=raw14.loc[row.source_row-2,[c.replace('-','_',1) for c in columns]].apply(pd.to_numeric,errors='coerce').to_numpy(float)
        b=exp.loc[row.export_row,columns].apply(pd.to_numeric,errors='coerce').to_numpy(float)
        finite=np.isfinite(a);same=np.array_equal(finite,np.isfinite(b));unit=10**(np.floor(np.log10(np.maximum(abs(a[finite]),1e-300)))-5)
        passed=same and bool((abs(a[finite]-b[finite])<=unit+1e-10).all())
        identity.append({'source_record_number':int(row.source_row-1),'n_features':320,'same_missingness':same,'within_one_6significant_digit_unit':passed,'max_last_digit_units':float(np.max(abs(a[finite]-b[finite])/unit))})
    assert all(x['within_one_6significant_digit_unit'] for x in identity)
    csvout(RESULTS/'INDEPENDENT_2014_RELEASE_IDENTITY_CHECKS.csv',identity)
    with zipfile.ZipFile(manifests.loc['L04','path']) as z:
        data19=pd.read_csv(z.open('csfproteomics/adni_csfproteomics.csv'),low_memory=False)
        member=next(n for n in z.namelist() if 'Supplementary' in n and n.endswith('.xlsx'))
        payload=z.read(member);book=pd.ExcelFile(io.BytesIO(payload),engine='openpyxl')
        ann=pd.read_excel(book,'Sample Annotations');trans=pd.read_excel(book,'Relative Quant Transitions')
        selected=trans.loc[trans['Peptide.Modified.Sequence'].isin(sum(PEPTIDES.values(),[]))].copy();assert len(selected)==10
        codecol='TransitionCode(Protein-Peptide-Fragment-Exclude)'
        code_to_field={str(r[codecol]):str(r['Transition']).strip().replace('-','.').upper().replace('+','.') for r in selected.to_dict('records')}
        assert all(field in data19 for field in code_to_field.values())
        assert ann['Sample Name'].is_unique
        merged=data19[['RID','EXAMDATE','SAMPLE.NAME']].merge(ann,left_on='SAMPLE.NAME',right_on='Sample Name',how='left',validate='one_to_one')
        assert merged['SampleCode.Replicate'].notna().all()
        csvout(PRIVATE/'MRM2019_SAMPLE_ANNOTATION_LINK_RESTRICTED.csv',merged)
        dupe=[g for _,g in merged.groupby(['RID','EXAMDATE']) if len(g)==2];pair_meta=[]
        for g in dupe:
            codes=g['SampleCode.Replicate'].astype(str);bases=codes.str.rsplit('.',n=1).str[0];rep=codes.str.rsplit('.',n=1).str[-1]
            verified=bases.nunique()==1 and set(rep)=={'01','02'} and g['Sample Group'].str.upper().eq('RARC').all()
            pair_meta.append({'pair_number':len(pair_meta)+1,'same_sample_code_base':bases.nunique()==1,'replicate_labels_01_02':set(rep)=={'01','02'},'both_RARC':g['Sample Group'].str.upper().eq('RARC').all(),'verified_replicate_pair':verified,'different_instrument':g.Instrument.nunique()>1,'different_run':g['MS_Run.SubRun'].nunique()>1,'different_MS_plate':g['MS PlateIndex'].nunique()>1})
        assert len(dupe)==20 and all(x['verified_replicate_pair'] for x in pair_meta)
        csvout(RESULTS/'MRM2019_REPLICATE_PAIR_METADATA_AUDIT.csv',pair_meta)
        tables={}
        for name in ['Unfiltered Light PA','Unfiltered Heavy PA','Sample-Transition-Filter-Flags','Filtered-log2(ratio(light+1))','Normalized-log2(ratio(light+1))']:
            tables[name]=pd.read_excel(book,name,usecols=lambda c:c in list(code_to_field)+['Sample Name','Sort'])
            assert tables[name]['Sample Name'].is_unique
            tables[name]=tables[name].set_index('Sample Name')
        legend=pd.read_excel(book,'Legend',header=None)
        (EVIDENCE/'MRM2019_SUPPLEMENT_LEGEND.txt').write_text('\n'.join(' | '.join(str(x) for x in row if pd.notna(x)) for row in legend.to_numpy()),encoding='utf-8')
        book.close()
    csvout(RESULTS/'MRM2019_TRANSITION_CODE_MAPPING.csv',selected)
    for name,df in tables.items():csvout(PRIVATE/('MRM2019_'+re.sub('[^A-Za-z0-9]+','_',name)+'_TARGETS_RESTRICTED.csv'),df.reset_index())
    transreg=pd.read_csv(S6/'results/LEGACY_CANDIDATE_TRANSITION_REGISTRY.csv')
    for gene,peps in PEPTIDES.items():
        for pep in peps:
            a14=[];b14=[];a19=[];b19=[]
            for num,g in enumerate(groups14,1):
                i,j=g.source_row.to_numpy(int)-2;x=float(pd.to_numeric(raw14.loc[i,gene+'_'+pep],errors='coerce'));y=float(pd.to_numeric(raw14.loc[j,gene+'_'+pep],errors='coerce'))
                a14.append(x);b14.append(y);pair_rows.append({'release':'MRM2014','gene':gene,'peptide':pep,'pair_number':num,'value_first':x,'value_second':y,'pair_identity_verified':True})
            fields=transreg.loc[transreg.peptide.eq(pep),'transition_field'].tolist()
            for num,g in enumerate(dupe,1):
                g=g.sort_values('SampleCode.Replicate');names=g['SAMPLE.NAME'].tolist()
                pair=data19.set_index('SAMPLE.NAME').loc[names,fields].apply(pd.to_numeric,errors='coerce').mean(axis=1,skipna=False).to_numpy()
                a19.append(pair[0]);b19.append(pair[1]);pair_rows.append({'release':'MRM2019','gene':gene,'peptide':pep,'pair_number':num,'value_first':pair[0],'value_second':pair[1],'pair_identity_verified':True})
            for release,x,y in [('MRM2014',a14,b14),('MRM2019',a19,b19)]:
                seed=202610280900+int(hashlib.sha256((release+pep).encode()).hexdigest()[:8],16)
                repeat.append({'release':release,'gene':gene,'peptide':pep,'n_identified_pairs':len(x),**pair_metrics(np.array(x),np.array(y),seed),'scale':'RELEASED_NORMALIZED_LOG_PEPTIDE','evidence':'TECHNICAL_REPEATABILITY_ONLY','no_pTau_or_new_pathology_test':True})
    csvout(RESULTS/'CORE_PEPTIDE_TECHNICAL_REPEATABILITY.csv',repeat)
    csvout(PRIVATE/'CORE_PEPTIDE_REPLICATE_VALUES_RESTRICTED.csv',pair_rows)
    lf=pd.read_csv(S6/'00_restricted/LEGACY_L04_CLINICAL_DESIGN_RESTRICTED.csv')
    designs=pickle.load((S6/'00_restricted/LEGACY_CLINICAL_DESIGNS_RESTRICTED.pkl').open('rb'))
    first=next(v for k,v in designs.items() if k.startswith('NEW_RID_CLINICAL__'))
    new_names=set(data19.iloc[lf.loc[first['indices'],'source_row'].to_numpy(int)-2]['SAMPLE.NAME'])
    pqc_names=ann.loc[ann['Sample Group'].str.upper().eq('PQC'),'Sample Name'].tolist()
    assert len(new_names)==95 and len(pqc_names)==24
    qc=[];agreement=[];pqc=[]
    for code,field in code_to_field.items():
        gene,pep=field.split('.')[:2]
        for group,names in [('ALL_PUBLISHED',data19['SAMPLE.NAME'].tolist()),('NEW_RID95',sorted(new_names)),('PQC24',pqc_names)]:
            L=pd.to_numeric(tables['Unfiltered Light PA'].loc[names,code],errors='coerce').to_numpy(float)
            H=pd.to_numeric(tables['Unfiltered Heavy PA'].loc[names,code],errors='coerce').to_numpy(float)
            flag=tables['Sample-Transition-Filter-Flags'].loc[names,code]
            filt=pd.to_numeric(tables['Filtered-log2(ratio(light+1))'].loc[names,code],errors='coerce').to_numpy(float)
            normalized=pd.to_numeric(tables['Normalized-log2(ratio(light+1))'].loc[names,code],errors='coerce').to_numpy(float)
            valid=(H>0)&(L>=0)&np.isfinite(L)&np.isfinite(H);rawlog=np.full(len(L),np.nan);rawlog[valid]=np.log2((L[valid]+1)/H[valid])
            difference=filt-rawlog;finite=np.isfinite(difference)
            qc.append({'gene':gene.replace('_HUMAN',''),'peptide':pep,'transition':field,'code':code,'group':group,'n_rows':len(names),'n_positive_heavy_and_nonnegative_light':int(valid.sum()),'flag_counts_json':json.dumps(flag.fillna('BLANK').astype(str).value_counts().to_dict()),'n_filtered_log_finite':int(np.isfinite(filt).sum()),'n_normalized_log_finite':int(np.isfinite(normalized).sum()),'rawlog_formula':'log2((Light+1)/Heavy), per source methods','max_filtered_log_reconstruction_error_finite':float(np.max(abs(difference[finite]))) if finite.any() else np.nan})
            if group=='PQC24':
                v=L[valid]/H[valid];good=v>0
                pqc.append({'gene':gene.replace('_HUMAN',''),'peptide':pep,'transition':field,'n_PQC':24,'n_positive_finite_raw_LH':int(good.sum()),'raw_LH_CV_percent_all_positive_no_flag_exclusion':float(100*np.std(v[good],ddof=1)/np.mean(v[good])) if good.sum()>1 else np.nan,'flag_counts_json':json.dumps(flag.fillna('BLANK').astype(str).value_counts().to_dict()),'purpose':'SOURCE_QC_DIAGNOSTIC_NOT_ASSAY_CERTIFICATION','no_phenotype_selection':True})
            if group=='ALL_PUBLISHED':
                released=pd.to_numeric(data19.set_index('SAMPLE.NAME').loc[names,field],errors='coerce').to_numpy(float)
                both=np.isfinite(released)&np.isfinite(normalized);same=np.array_equal(np.isfinite(released),np.isfinite(normalized))
                err=float(np.max(abs(released[both]-normalized[both]))) if both.any() else np.nan
                agreement.append({'gene':gene.replace('_HUMAN',''),'peptide':pep,'transition':field,'n_pairs_compared':int(both.sum()),'same_missingness':same,'max_abs_release_workbook_difference':err,'within_1e_minus5':bool(same and err<=1e-5)})
    csvout(RESULTS/'MRM2019_TARGET_TRANSITION_RAW_QC.csv',qc)
    csvout(RESULTS/'MRM2019_PQC24_RAW_LH_VARIATION.csv',pqc)
    csvout(RESULTS/'MRM2019_WORKBOOK_RELEASE_VALUE_AGREEMENT.csv',agreement)
    save_json(RESULTS/'MRM2019_SUPPLEMENT_SOURCE.json',{'archive':str(manifests.loc['L04','path']),'archive_sha256':sha(manifests.loc['L04','path']),'member':member,'member_sha256':hashlib.sha256(payload).hexdigest(),'member_bytes':len(payload),'processed_release_names_linked':len(merged),'annotated_sample_rows':len(ann),'PQC_recovered':24,'verified_RARC_pairs':20,'raw_QC_is_new':True,'raw_intensities_not_original_chromatograms_or_S_N':True})
    log('2019 raw supplementary source linked:20 verified RARC pairs,24 PQC;2014 16 flagged pairs;10 peptide-repeatability summaries;target raw transition QC recovered')
    print(pd.DataFrame(repeat)[['release','gene','peptide','n_complete_pairs','pearson_r','Lin_CCC_descriptive','RMS_pair_difference_per_pooled_SD']].to_string(index=False))
    print(pd.DataFrame(agreement).to_string(index=False))

if __name__=='__main__':
    assert_unsealed()
    run()
