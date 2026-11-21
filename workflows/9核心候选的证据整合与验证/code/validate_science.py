from config import *
from scipy.stats import pearsonr,norm
from itertools import combinations
import pickle

def run():
    checks=[]
    def add(name,ok,detail=''):
        checks.append({'check':name,'passed':bool(ok),'detail':str(detail)})
    original=pd.read_csv(RESULTS/'INDEPENDENT_NUMERIC_CHECKS.csv')
    for r in original.itertuples():add('NUMERIC/'+r.check,r.passed,r.max_absolute_difference)
    molecular=pd.read_csv(RESULTS/'INDEPENDENT_MOLECULAR_CHECKS.csv')
    for r in molecular.itertuples():add('MOLECULAR/'+r.check,r.passed,f'{r.actual}/{r.expected}')
    new=pd.read_csv(RESULTS/'INDEPENDENT_PLATFORM_EFFECTS.csv')
    old=pd.read_csv(S8/'results/PLATFORM_EFFECTS_ALL_MODELS.csv')
    for r in new.itertuples():
        ref=old.loc[old.accession.eq(r.accession)&old.model.eq(r.model)&old.assay.eq(r.assay)].iloc[0]
        for f,g in [('scaled_b','scaled_b'),('scaled_SE_full_joint','scaled_SE_full_joint_uncertainty'),('scaled_CI95_low','scaled_CI95_low'),('scaled_CI95_high','scaled_CI95_high')]:
            actual=getattr(r,f)
            if f in ['scaled_CI95_low','scaled_CI95_high']:
                sign=-1 if f.endswith('low') else 1
                actual=r.scaled_b+sign*1.959964*r.scaled_SE_full_joint
                add('EXACT_NORMAL_CI/'+r.accession+'_'+r.model+'_'+r.assay+'_'+f,np.isclose(getattr(r,f),r.scaled_b+sign*norm.ppf(.975)*r.scaled_SE_full_joint,atol=1e-12))
            add('SCALED/'+r.accession+'_'+r.model+'_'+r.assay+'_'+f,np.isclose(actual,ref[g],rtol=2e-7,atol=1e-10),'Legacy CI uses rounded1.959964; new CI uses exact normal .975 quantile' if 'CI95' in f else '')
    family=pd.read_csv(RESULTS/'FAMILY34_INTEGRATED_EVIDENCE_ALL_MODELS.csv')
    for model,part in family.groupby('model'):
        add('FAMILY34/'+model,len(part)==34)
        add('SHARED14/'+model,int(part.shared_supported.sum())==14)
        oldgene=set(family.loc[family.model.eq('M0')&family.shared_supported,'gene'])
        add('SHARED_IDENTICAL_GENES/'+model,set(part.loc[part.shared_supported,'gene'])==oldgene)
    m0=family.loc[family.model.eq('M0')]
    add('M0_H_exact_three',set(m0.loc[m0.H_deviation_supported,'gene'])==set(CANDIDATES.values()))
    add('M0_both_NPTX1_KLK6',set(m0.loc[m0.shared_supported&m0.H_deviation_supported,'gene'])=={'NPTX1','KLK6'})
    add('NO_EQ_DISCOVERY',not family.M0_BY_equivalence.any())
    p=m0.loc[m0.gene.eq('PDGFRB')].iloc[0]
    add('PDGFRB_no_rounding_upgrade',p.q_BY_shared34>.05 and not p.shared_supported,p.q_BY_shared34)
    previous_assoc=pd.read_csv(S8/'results/ALL35_THREE_MODEL_SHARED_RESULTS.csv')
    for model in MODELS:
        part=family.loc[family.model.eq(model)]
        for row in part.itertuples():
            add('INTEGRATED_H_Q/'+row.accession+'_'+model,np.isclose(row.q_BY_H34,previous_assoc.loc[previous_assoc.accession.eq(row.accession)&previous_assoc.model.eq(model),'prior_H_q_BY_by_model'].iloc[0]))
    pairs=pd.read_csv(PRIVATE/'CORE_PEPTIDE_REPLICATE_VALUES_RESTRICTED.csv')
    reps=pd.read_csv(RESULTS/'CORE_PEPTIDE_TECHNICAL_REPEATABILITY.csv',dtype={'bootstrap_seed':str})
    for r in reps.itertuples():
        part=pairs.loc[pairs.release.eq(r.release)&pairs.peptide.eq(r.peptide)]
        x=part.value_first.to_numpy();y=part.value_second.to_numpy();ok=np.isfinite(x)&np.isfinite(y);x=x[ok];y=y[ok];n=len(x)
        add('REPEAT/'+r.release+r.peptide+'_n',n==r.n_complete_pairs)
        add('REPEAT/'+r.release+r.peptide+'_pearson',np.isclose(pearsonr(x,y).statistic,r.pearson_r,atol=1e-12))
        # Vectorized independent recount of all4000 bootstrap draws.
        rng=np.random.default_rng(int(r.bootstrap_seed));pick=rng.integers(0,n,size=(4000,n));xx=x[pick];yy=y[pick]
        xx-=xx.mean(axis=1,keepdims=True);yy-=yy.mean(axis=1,keepdims=True)
        values=np.sum(xx*yy,axis=1)/np.sqrt(np.sum(xx**2,axis=1)*np.sum(yy**2,axis=1));lo,hi=np.quantile(values,[.025,.975])
        add('REPEAT/'+r.release+r.peptide+'_bootstrap_interval',np.allclose([lo,hi],[r.pearson_pair_bootstrap95_low,r.pearson_pair_bootstrap95_high],atol=1e-12))
        add('REPEAT/'+r.release+r.peptide+'_pair_identity',part.pair_identity_verified.all())
    meta=pd.read_csv(RESULTS/'MRM2019_REPLICATE_PAIR_METADATA_AUDIT.csv')
    add('RARC20_pairs_verified',len(meta)==20 and meta.verified_replicate_pair.all())
    # Read independently exported source target tables and verify raw QC summaries.
    light=pd.read_csv(PRIVATE/'MRM2019_Unfiltered_Light_PA_TARGETS_RESTRICTED.csv').set_index('Sample Name')
    heavy=pd.read_csv(PRIVATE/'MRM2019_Unfiltered_Heavy_PA_TARGETS_RESTRICTED.csv').set_index('Sample Name')
    flags=pd.read_csv(PRIVATE/'MRM2019_Sample_Transition_Filter_Flags_TARGETS_RESTRICTED.csv').set_index('Sample Name')
    filtered=pd.read_csv(PRIVATE/'MRM2019_Filtered_log2_ratio_light_1__TARGETS_RESTRICTED.csv').set_index('Sample Name')
    ann=pd.read_csv(PRIVATE/'MRM2019_SAMPLE_ANNOTATION_LINK_RESTRICTED.csv')
    # PQC rows are absent from clinical750 annotations, so recover names from source table labels R999.
    pqc_names=[k for k in light.index if str(k).startswith('PQC')]
    if len(pqc_names)!=24:
        # The frozen supplement has Caprion CSF-style names; metadata source cache retains the full annotation.
        import zipfile,io
        mf=pd.read_csv(S6/'results/LEGACY_INPUT_MANIFEST.csv').set_index('input_id')
        with zipfile.ZipFile(mf.loc['L04','path']) as z:
            mem=json.loads((RESULTS/'MRM2019_SUPPLEMENT_SOURCE.json').read_text(encoding='utf-8'))['member']
            fullann=pd.read_excel(io.BytesIO(z.read(mem)),sheet_name='Sample Annotations')
        pqc_names=fullann.loc[fullann['Sample Group'].str.upper().eq('PQC'),'Sample Name'].tolist()
    add('PQC24_metadata_count',len(pqc_names)==24)
    pqc=pd.read_csv(RESULTS/'MRM2019_PQC24_RAW_LH_VARIATION.csv')
    rawqc=pd.read_csv(RESULTS/'MRM2019_TARGET_TRANSITION_RAW_QC.csv')
    mapping=pd.read_csv(RESULTS/'MRM2019_TRANSITION_CODE_MAPPING.csv')
    codecol='TransitionCode(Protein-Peptide-Fragment-Exclude)'
    for r in mapping.to_dict('records'):
        code=r[codecol];L=pd.to_numeric(light.loc[pqc_names,code],errors='coerce').to_numpy();H=pd.to_numeric(heavy.loc[pqc_names,code],errors='coerce').to_numpy();ratio=L/H
        cv=100*ratio.std(ddof=1)/ratio.mean();field=r['Transition'].strip().replace('-','.').upper().replace('+','.')
        target=pqc.loc[pqc.transition.eq(field)].iloc[0]
        add('PQC_CV/'+code,np.isclose(cv,target.raw_LH_CV_percent_all_positive_no_flag_exclusion,rtol=1e-12))
        add('PQC_QC0/'+code,pd.to_numeric(flags.loc[pqc_names,code],errors='coerce').eq(0).all())
        Ll=pd.to_numeric(light[code],errors='coerce').to_numpy();Hh=pd.to_numeric(heavy[code],errors='coerce').to_numpy();F=pd.to_numeric(filtered[code],errors='coerce').to_numpy()
        expected=np.log2((Ll+1)/Hh);finite=np.isfinite(expected)&np.isfinite(F)
        add('RAW_LOG_FORMULA/'+code,np.max(abs(F[finite]-expected[finite]))<1e-10)
        rr=rawqc.loc[rawqc.code.eq(code)&rawqc.group.eq('NEW_RID95')].iloc[0]
        add('NEW95_ALL_QC0/'+code,json.loads(rr.flag_counts_json)=={'0':95})
    agreement=pd.read_csv(RESULTS/'MRM2019_WORKBOOK_RELEASE_VALUE_AGREEMENT.csv')
    for r in agreement.itertuples():add('SOURCE_WORKBOOK_RELEASE/'+r.transition,r.same_missingness and r.max_abs_release_workbook_difference<1e-10)
    identity=pd.read_csv(RESULTS/'INDEPENDENT_2014_RELEASE_IDENTITY_CHECKS.csv')
    for r in identity.itertuples():add('SOURCE2014_VECTOR/'+str(r.source_record_number),r.n_features==320 and r.same_missingness and r.within_one_6significant_digit_unit)
    assoc=pd.read_csv(RESULTS/'INDEPENDENT_2019_NEW_RID_PEPTIDE_ASSOCIATIONS.csv')
    add('2019_FIXED_BY5_NO_DISCOVERY',len(assoc)==5 and not assoc.q_BY5.le(.05).any())
    influence=pd.read_csv(RESULTS/'NEW_RID_LEAVE_ONE_BATCH_DESCRIPTIVE.csv')
    add('BATCH_DELETIONS85_ALL_ESTIMABLE_POSITIVE',len(influence)==85 and influence.status.eq('DESCRIPTIVE_ONLY').all() and influence.coefficient_positive.all())
    geom=pd.read_csv(RESULTS/'MRM2014_THREE_FOUR_ASSAY_DESIGN_AUDIT.csv')
    for r in geom.itertuples():
        add('DESIGN/'+r.gene+r.peptide+r.assay_design+'_rank_df',r.n==r.rank_X+r.residual_df)
        add('DESIGN/'+r.gene+r.peptide+r.assay_design+'_no_unverified_fit',not r.strict_validation_eligible and not r.candidate_pathology_effect_fitted)
        if r.assay_design=='N_S_TMT_MRM':add('DESIGN/'+r.gene+r.peptide+'_four_df_lt40',r.residual_df<40)
    scan=pd.read_csv(RESULTS/'METADATA_ID_SCAN_RESULTS.csv');outside=scan.loc[scan.n_matched_rows.gt(0)&~scan.table.str.startswith('2014_WORKBOOK/')]
    add('NO_EXTERNAL_SPECIMEN_ID_MATCH',len(outside)==0)
    adjud=json.loads((RESULTS/'DATE_BRIDGE_ADJUDICATION.json').read_text(encoding='utf-8'))
    add('DATE_REMAINS_UNVERIFIED',adjud['verified_MRM_specimen_dates']==0 and adjud['2014_pathology_associations_fitted']==0)
    overlap=pd.read_csv(RESULTS/'SUBJECT_OVERLAP_AND_INDEPENDENCE_AUDIT.csv')
    z=overlap.loc[overlap.cohort_A.eq('S5_INDEX571')&overlap.cohort_B.eq('MRM2019_NEW95')].iloc[0]
    add('SUBJECT_DISJOINT_95_VS571',z.RID_overlap==0 and z.RID_B==95)
    # Check every frozen prior-stage entry, caching repeated paths to avoid rereading large sources.
    pri=pd.read_csv(S8/'results/PRIOR_SEALED_FILE_CHECKS_FINAL.csv')
    newest=pd.read_csv(S8/'SHA256_MANIFEST.csv')
    items=pri[['stage','path','expected_sha256']].to_dict('records')
    items.extend({'stage':S8.name,'path':str(S8/r.relative_path),'expected_sha256':r.sha256} for r in newest.itertuples())
    cache={};prior=[]
    for row in items:
        p=Path(row['path']);k=str(p.resolve())
        if k not in cache:cache[k]=sha(p) if p.exists() else ''
        actual=cache[k];prior.append({**row,'actual_sha256':actual,'passed':actual==row['expected_sha256']})
    csvout(RESULTS/'PRIOR_ALL_STAGE_SEAL_CHECKS_FINAL.csv',prior)
    add('ALL_PRIOR_SEALED_ENTRIES_UNCHANGED',all(x['passed'] for x in prior),f'{len(prior)} entries/{len(cache)} unique paths')
    # Aggregate every local source actually used by this stage.
    paths=set()
    for filename,pathcol in [('INPUT_MANIFEST.csv','path'),('METADATA_SOURCE_MANIFEST.csv','path'),('ORIGINAL_ASSAY_METHOD_SOURCES.csv','path')]:
        paths.update(pd.read_csv(RESULTS/filename)[pathcol])
    paths.update(pd.read_csv(S6/'results/LEGACY_INPUT_MANIFEST.csv').path)
    paths.update(str(S5/'results'/n) for n in ['ALL_FIXED_VARIANT_RESULTS.csv','COVARIANCE_OFFSET_SENSITIVITY_SUMMARY.csv'])
    paths.update(str(S7/'00_restricted'/n) for n in ['MRM2014_CSV_EXPORT_RESTRICTED.csv'])
    paths.update(str(RAW/'CSF'/n) for n in ['LOCLAB_15May2026.csv','CSFMETH_15May2026.csv'])
    paths.add(str(RAW/'Data dictionary/DATADIC_15May2026.csv'))
    sources=[]
    for path in sorted(paths):
        p=Path(path);actual=cache.get(str(p.resolve())) or sha(p)
        sources.append({'path':str(p),'bytes':p.stat().st_size,'sha256':actual,'restricted_person_level':('00_restricted' in str(p)) or (p.suffix.lower() in ['.csv','.xlsx','.npz','.pkl'] and str(RAW) in str(p))})
    csvout(RESULTS/'FINAL_LOCAL_SOURCE_MANIFEST.csv',sources)
    csvout(RESULTS/'SCIENTIFIC_QA_CHECKS.csv',checks)
    save_json(RESULTS/'SCIENTIFIC_QA_SUMMARY.json',{'time':now(),'checks':len(checks),'passed':sum(x['passed'] for x in checks),'failed':sum(not x['passed'] for x in checks),'prior_sealed_entries':len(prior),'prior_unique_paths':len(cache),'independent_three_assay_designs':105,'new2019_peptides_independently_refit':5,'new_pathology_simulations':0,'paired_bootstrap_draws_for_descriptive_intervals':40000,'external_individual_replication_completed':False})
    assert all(x['passed'] for x in checks),[x for x in checks if not x['passed']][:10]
    log(f'Scientific QA:{len(checks)}/{len(checks)} passed;{len(prior)} prior sealed entries unchanged across{len(cache)} unique paths')

if __name__=='__main__':
    assert_unsealed()
    run()
