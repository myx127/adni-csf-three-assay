from config import *
import zipfile,pickle

def run():
    independent=pd.read_csv(RESULTS/'INDEPENDENT_ALL35_THREE_MODEL_REPRODUCTION.csv')
    old=pd.read_csv(S8/'results/ALL35_THREE_MODEL_SHARED_RESULTS.csv')
    eq=pd.read_csv(S7/'results/CONSENSUS_INFORMATION_BOTTLENECKS.csv')
    matrix=independent.loc[independent.in_family34].copy()
    eq05=eq.loc[eq.margin.eq(.5)&eq.in_family34].set_index('accession')
    prior_eq=eq05[['equivalence_eligible_Stage6','q_BY_Stage6','BY_equivalence_Stage6','point_structure_compatible','gate_reasons']].reset_index()
    matrix=matrix.merge(prior_eq,on='accession',how='left',validate='many_to_one').rename(columns={'q_BY_Stage6':'M0_q_equiv_working_margin05','equivalence_eligible_Stage6':'M0_equivalence_eligible','BY_equivalence_Stage6':'M0_BY_equivalence','point_structure_compatible':'M0_factor_point_compatible','gate_reasons':'M0_equivalence_gate_reasons'})
    matrix['shared_supported']=matrix.q_BY_shared34.le(.05)
    matrix['H_deviation_supported']=matrix.q_BY_H34.le(.05)
    matrix['evidence_state']=np.select([matrix.shared_supported&matrix.H_deviation_supported,matrix.shared_supported,matrix.H_deviation_supported],['SHARED_AND_H_DEVIATION','SHARED_ONLY_H_UNRESOLVED','H_ONLY_SHARED_UNRESOLVED'],default='BOTH_UNRESOLVED')
    matrix['equivalence_evidence_is_M0_only']=True
    csvout(RESULTS/'FAMILY34_INTEGRATED_EVIDENCE_ALL_MODELS.csv',matrix)
    csvout(RESULTS/'CORE_CANDIDATE_MODEL_EVIDENCE.csv',matrix.loc[matrix.accession.isin(CANDIDATES)])
    counts=[]
    for model,part in matrix.groupby('model',sort=False):
        for state,g in part.groupby('evidence_state'):
            counts.append({'model':model,'evidence_state':state,'n':len(g),'genes':'|'.join(sorted(g.gene))})
    csvout(RESULTS/'INTEGRATED_EVIDENCE_STATE_COUNTS.csv',counts)
    variants=pd.read_csv(S5/'results/ALL_FIXED_VARIANT_RESULTS.csv')
    variants=variants.loc[variants.accession.isin(CANDIDATES)]
    csvout(RESULTS/'CORE_PRIOR_FIXED_SENSITIVITY_RESULTS.csv',variants)
    offsets=pd.read_csv(S5/'results/COVARIANCE_OFFSET_SENSITIVITY_SUMMARY.csv')
    offsets=offsets.loc[offsets.accession.isin(CANDIDATES)]
    offsets['gene']=offsets.accession.map(CANDIDATES)
    csvout(RESULTS/'CORE_PRIOR_ERROR_COVARIANCE_SENSITIVITY.csv',offsets)
    association=pd.read_csv(RESULTS/'INDEPENDENT_2019_NEW_RID_PEPTIDE_ASSOCIATIONS.csv')
    rep=pd.read_csv(RESULTS/'CORE_PEPTIDE_TECHNICAL_REPEATABILITY.csv')
    qc=pd.read_csv(RESULTS/'MRM2019_TARGET_TRANSITION_RAW_QC.csv')
    pqc=pd.read_csv(RESULTS/'MRM2019_PQC24_RAW_LH_VARIATION.csv')
    pep=pd.read_csv(RESULTS/'CORE_PEPTIDE_SEQUENCE_REGION_SPECIFICITY.csv')
    summary=[];tiers=[]
    molecular=pd.read_csv(S6/'results/THREE_DISCOVERY_MOLECULAR_EVIDENCE_MATRIX.csv').set_index('accession')
    for acc,gene in CANDIDATES.items():
        part=matrix.loc[matrix.accession.eq(acc)].set_index('model');row={'accession':acc,'gene':gene,'candidate_fixed_from_original_H_discovery':True}
        for model in MODELS:
            for key in ['q_BY_shared34','q_BY_H34','q_BY_combined68_shared','q_BY_combined68_H','p_shared','p_H','n','residual_df','evidence_state']:
                row[model+'_'+key]=part.loc[model,key]
        row.update({'M0_eq_working_margin05_q':part.loc['M0','M0_q_equiv_working_margin05'],'M0_eq_eligible':part.loc['M0','M0_equivalence_eligible'],'M0_factor_point_compatible':part.loc['M0','M0_factor_point_compatible']})
        a=association.loc[association.gene.eq(gene)];p=pep.loc[pep.gene.eq(gene)];q=qc.loc[qc.gene.eq(gene)&qc.group.eq('NEW_RID95')];qp=pqc.loc[pqc.gene.eq(gene)]
        row.update({'direct_MRM_peptides':len(p),'peptide_names':'|'.join(p.peptide),'human_reference_gene_specific_all_target_peptides':bool(len(p)>0 and p.EXACT_gene_specific_in_frozen_reference.all() and p.IL_EQUIVALENT_gene_specific_in_frozen_reference.all()),'new_RID_association_n':int(a.n.min()) if len(a) else 0,'new_RID_association_best_BY5_q':a.q_BY5.min(),'new_RID_association_BY5_supported_peptides':int(a.q_BY5.le(.05).sum()),'new_RID_all_target_QC0':bool(len(q)>0 and all(json.loads(v)=={'0':95} for v in q.flag_counts_json)),'PQC24_raw_LH_CV_min_percent':qp.raw_LH_CV_percent_all_positive_no_flag_exclusion.min(),'PQC24_raw_LH_CV_max_percent':qp.raw_LH_CV_percent_all_positive_no_flag_exclusion.max()})
        for release in ['MRM2014','MRM2019']:
            rp=rep.loc[rep.gene.eq(gene)&rep.release.eq(release)]
            row[release+'_replicate_r_min']=rp.pearson_r.min();row[release+'_replicate_r_max']=rp.pearson_r.max();row[release+'_complete_pair_n_min']=rp.n_complete_pairs.min() if len(rp) else 0
        row.update({'source_specific_2014_dates_verified':0,'strict_relative_pattern_validation_completed':False,'independent_external_cohort_validation_completed':False,'specific_molecular_form_or_mechanism_verified':False,'binding_epitope_available':False,'original_TMT_peptide_matrix_obtained':False})
        for sourcekey in ['soma_features','nulisa_label','TMT_aggregate_column']:row[sourcekey]=molecular.loc[acc,sourcekey]
        raw=variants.loc[variants.accession.eq(acc)&variants.variant.eq('RAW_P')]
        row['prior_RAW_P_H_BY_q']=raw.q_BY.iloc[0] if len(raw) else np.nan
        summary.append(row)
        tiers.extend([
            {'gene':gene,'evidence_layer':'Original H primary discovery','status':'SUPPORTED_IN_M0','support':'Original BY34; APOE sensitivity retained'},
            {'gene':gene,'evidence_layer':'Three-assay same-direction pTau association','status':'SUPPORTED_ALL_THREE_MODELS' if part.q_BY_shared34.le(.05).all() else 'NOT_BY34_SUPPORTED','support':'Same ADNI participants; not independent replication'},
            {'gene':gene,'evidence_layer':'Peptide identity and source QC','status':'DIRECT_GENE_SPECIFIC_PEPTIDES_AND_QC' if len(p) else 'DIRECT_RECEPTOR_PEPTIDES_NOT_COVERED','support':'No intact-species/activity claim'},
            {'gene':gene,'evidence_layer':'Subject-disjoint same-ADNI peptide association','status':'POSITIVE_DIRECTION_BY5_NOT_PASSED' if len(a) else 'NOT_COVERED','support':'Five peptide family retained; no pooled/selected-peptide rescue'},
            {'gene':gene,'evidence_layer':'Strict relative response pattern replication','status':'NOT_ESTABLISHED','support':'2019 low df;2014 date/geometry/inference gates'},
            {'gene':gene,'evidence_layer':'Independent external cohort replication','status':'NOT_COMPLETED','support':'Emory public metadata available; file download requires authorized access'},
            {'gene':gene,'evidence_layer':'Specific mechanism/molecular form','status':'NOT_IDENTIFIED','support':'No binding-region/original TMT peptide evidence'}
        ])
    csvout(RESULTS/'THREE_CORE_CANDIDATE_EVIDENCE_MATRIX.csv',summary)
    csvout(RESULTS/'CANDIDATE_VALIDATION_STATUS_BY_LAYER.csv',tiers)
    # Audit subject overlap across releases and technical-repeatability samples.
    s5=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv');c14=pd.read_csv(S7/'00_restricted/MRM2014_RESOLVED_TECHNICAL_CROSSWALK_RESTRICTED.csv')
    mf=pd.read_csv(S6/'results/LEGACY_INPUT_MANIFEST.csv').set_index('input_id')
    with zipfile.ZipFile(mf.loc['L04','path']) as z:l04=pd.read_csv(z.open('csfproteomics/adni_csfproteomics.csv'),usecols=['RID','EXAMDATE'])
    valid=l04.loc[pd.to_numeric(l04.RID,errors='coerce').notna()&pd.to_numeric(l04.RID,errors='coerce').ne(999999)]
    lf=pd.read_csv(S6/'00_restricted/LEGACY_L04_CLINICAL_DESIGN_RESTRICTED.csv')
    ds=pickle.load((S6/'00_restricted/LEGACY_CLINICAL_DESIGNS_RESTRICTED.pkl').open('rb'))
    di=next(v['indices'] for k,v in ds.items() if k.startswith('NEW_RID_CLINICAL__'))
    cohorts={'S5_INDEX571':set(s5.RID.astype(int)),'MRM2014_ALL':set(c14.RID.astype(int)),'MRM2019_VALID_ALL':set(valid.RID.astype(int)),'MRM2019_NEW95':set(lf.loc[di,'RID'].astype(int))}
    overlaps=[]
    for i,(name,a) in enumerate(cohorts.items()):
        for name2,b in list(cohorts.items())[i:]:
            overlaps.append({'cohort_A':name,'cohort_B':name2,'RID_A':len(a),'RID_B':len(b),'RID_overlap':len(a&b),'counts_not_independent_experiment_numbers':True})
    assert len(cohorts['S5_INDEX571']&cohorts['MRM2019_NEW95'])==0
    csvout(RESULTS/'SUBJECT_OVERLAP_AND_INDEPENDENCE_AUDIT.csv',overlaps)
    registry=[
        ['E01','Original N/S/TMT and Stage8 shared test','Same subjects,assays,pathology','MAIN_OBSERVATIONAL_EVIDENCE','Different questions; not independent validation'],
        ['E02','Independent Stage9 implementation','Same data and same estimands','COMPUTATIONAL_VERIFICATION','Code consistency cannot remove measurement confounding'],
        ['E03','MRM2014 CSV and original workbook','Same experiment and release aliases','SOURCE_PROVENANCE_AND_TECHNICAL_REPEATABILITY','Not two replications; date remains open'],
        ['E04','MRM2019 same-date intersection','Shares participants,N/S and pTau with main','SUPPLEMENTAL_MEASUREMENT_OPPORTUNITY','Residual df20-24; original relative test not executed'],
        ['E05','MRM2019 new95','Zero RID overlap with all571; same ADNI parent','SUBJECT_DISJOINT_SINGLE_ASSAY_ASSOCIATION','BY5 not passed; different hypothesis from H'],
        ['E06','2014 sixteen pairs and2019 twenty pairs','Technical replicates from same releases','TECHNICAL_REPEATABILITY','No independent biology; stable errors not excluded'],
        ['E07','2019 PQC24','Repeated pooled material; not24 distinct donors','SOURCE_PROCESS_QC','CV applies to target raw L/H and these pools'],
        ['E08','UniProt reference and source annotations','Independent annotation; no participant data','MOLECULAR_ANNOTATION','Not validation of actual epitope/species'],
        ['E09','External Emory2022','Separate parent cohort; Soma/Olink/TMT','PUBLIC_METADATA_LEAD_ONLY','Individual files unauthenticated403; n35 cap below frozen80']
    ]
    csvout(RESULTS/'EVIDENCE_SOURCE_INDEPENDENCE_REGISTRY.csv',pd.DataFrame(registry,columns=['evidence_id','source','dependence','evidence_tier','claim_boundary']))
    request=c14[['source_row','RID','VISCODE','LIMS_ID','External_Reference']].copy()
    for field in ['specimen_id','parent_specimen_id','aliquot_id','actual_collection_date','actual_collection_time','collection_date_source','biobank_or_lab_verified','same_parent_as_NULISA','same_parent_as_Soma','same_parent_as_TMT']:request[field]=''
    csvout(PRIVATE/'MRM2014_SPECIMEN_DATE_LINK_REQUEST_TEMPLATE_RESTRICTED.csv',request)
    gaps=[
        ['V01','2014 exact specimen date','NOT_CLOSED','LIMS_ID and External_Reference to specimen/aliquot IDs and actual collection date/time','Required before paired pathology inference'],
        ['V02','2014 four-assay design','CURRENT_DESIGN_BELOW_DF40','NPTX1 df19;KLK6 df17. Need more independent paired samples/batches or separately justified and calibrated design','Date correction alone cannot establish four-assay eligibility'],
        ['V03','2014 three-assay inference','DATE_AND_CALIBRATION_OPEN','N/S/MRM n144/141,df78/75 conditionally; day/column dependence and high leverage require actual-design calibration','Do not equate nominal geometry with inferential qualification'],
        ['V04','PDGFRB direct receptor coverage','LOCAL_GAP','P09619 receptor-derived peptides plus sequence and sample map; PDGF-BB is ligand','Required for direct receptor measurement support'],
        ['V05','Original TMT peptide and grouping','ACCESS_IDENTIFIED_NOT_OBTAINED','Method points to syn59804727; anonymous access403. Need peptide/PSM matrix,UniProt mapping,uniqueness,modifications,plex/channel','Required to locate domains/forms or diagnose protein rollup'],
        ['V06','NULISA/Soma binding regions','NOT_AVAILABLE','Capture/detection antibody and aptamer binding regions,cross-reactivity,complex/matrix responses','Required to attribute molecular recognition differences'],
        ['V07','External individual-level paired data','PUBLIC_ENTITY_FOUND_FILES_ACCESS403','Emory syn30549894/syn30550193 plus subject/pTau/batch metadata; frozen n80 remains unmet at n35 cap','External descriptive support possible after authorized access; no formal frozen replication now'],
        ['V08','Full error covariance/causality','NOT_IDENTIFIED','Matched aliquots,randomized technical repeats across all platforms and form-specific perturbations','MRM repeats alone do not identify latent abundance or causal pathway']
    ]
    csvout(RESULTS/'FINAL_VALIDATION_GAP_AND_ACTION_REGISTRY.csv',pd.DataFrame(gaps,columns=['gap_id','issue','status','minimum_evidence_or_action','decision_impact']))
    log('34-protein integrated evidence,three-candidate validation matrix,independence audit,and explicit unresolved-validation registry completed')
    print(pd.DataFrame(summary)[['gene','M0_q_BY_shared34','M0_q_BY_H34','APOE_q_BY_H34','new_RID_association_best_BY5_q','strict_relative_pattern_validation_completed']].to_string(index=False))
    print(pd.DataFrame(overlaps).to_string(index=False))

if __name__=='__main__':
    assert_unsealed()
    run()
