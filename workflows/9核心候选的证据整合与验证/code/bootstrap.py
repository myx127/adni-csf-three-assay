from config import *
import shutil

def run():
    rows=[]
    wanted={
        S1:['results/EXTENDED_CORE_MAPPING_DRAFT.csv'],
        S5:['00_restricted/CLINICAL_DESIGN_RESTRICTED.csv','00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz','results/MAIN_ANALYSIS_SAP_FROZEN.json'],
        S6:['results/COVARIATE_FOUR_MODEL_RESULTS.csv','results/LEGACY_INPUT_MANIFEST.csv','results/LEGACY_NEW_RID_PEPTIDE_PTAU_ASSOCIATION.csv','00_restricted/LEGACY_L04_CLINICAL_DESIGN_RESTRICTED.csv','00_restricted/LEGACY_CLINICAL_DESIGNS_RESTRICTED.pkl','00_restricted/LEGACY_CASE_MASKS_RESTRICTED.pkl','results/LEGACY_CANDIDATE_TRANSITION_REGISTRY.csv','results/LEGACY_DIRECT_DESIGN_SCREEN.csv','results/RELATIVE_RESPONSE_FIELLER_INTERVALS.csv','results/THREE_DISCOVERY_MOLECULAR_EVIDENCE_MATRIX.csv','results/HUMAN_REFERENCE_PEPTIDE_SPECIFICITY.csv','results/LEGACY_VERIFICATION_SAP_FROZEN.json'],
        S7:['00_restricted/MRM2014_RESOLVED_TECHNICAL_CROSSWALK_RESTRICTED.csv','00_restricted/MRM2014_INFERRED_VISIT_DESIGN_RESTRICTED.csv','00_restricted/MRM2014_CANDIDATE_DATE_CROSSWALK_RESTRICTED.csv','results/CONSENSUS_INFORMATION_BOTTLENECKS.csv','results/VALIDATION_GAP_REGISTRY.csv','results/EXTERNAL_VALIDATION_OPPORTUNITIES.csv'],
        S8:['results/ALL35_THREE_MODEL_SHARED_RESULTS.csv','results/PLATFORM_EFFECTS_ALL_MODELS.csv','results/CALIBRATION_DECISIONS.csv','STAGE8_PROTOCOL_v1.0.md']
    }
    for stage,rels in wanted.items():
        for rel in rels:
            p=stage/rel
            rows.append({'source_id':f'I{len(rows)+1:03d}','path':str(p),'bytes':p.stat().st_size,'sha256':sha(p),'scope':'PERSON_LEVEL_RESTRICTED' if '00_restricted' in rel else 'RESULT_OR_PROTOCOL'})
    for p in list((S6/'evidence').glob('*_UNIPROT.json'))+list((S6/'evidence').glob('*.fasta.gz')):
        rows.append({'source_id':f'I{len(rows)+1:03d}','path':str(p),'bytes':p.stat().st_size,'sha256':sha(p),'scope':'PUBLIC_REFERENCE'})
    csvout(RESULTS/'INPUT_MANIFEST.csv',rows)
    prior=[]
    for stage in [S8]:
        man=pd.read_csv(stage/'SHA256_MANIFEST.csv')
        log('Prior manifest fields: '+str(list(man)))
        for row in man.to_dict('records'):
            rel=row.get('relative_path',row.get('path',row.get('file')));p=stage/str(rel)
            expected=row.get('sha256',row.get('SHA256'));actual=sha(p) if p.exists() else ''
            prior.append({'stage':stage.name,'relative_path':rel,'bytes':p.stat().st_size if p.exists() else -1,'expected_sha256':expected,'actual_sha256':actual,'passed':actual==expected})
    csvout(RESULTS/'PRIOR_STAGE8_SEAL_CHECKS_INITIAL.csv',prior)
    assert all(x['passed'] for x in prior)
    save_json(RESULTS/'EXECUTION_PLAN_RECORD.json',{'time':now(),'protocol_sha256':sha(ROOT/'STAGE9_PROTOCOL_v1.0.md'),'all_prior_stage_results_already_known':True,'new_stage_only':True,'initial_prior_seal_entries':len(prior)})
    log(f'Stage9 initialized; {len(rows)} input records; {len(prior)} previous Stage8 sealed entries verified')

if __name__=='__main__':
    assert_unsealed()
    run()
