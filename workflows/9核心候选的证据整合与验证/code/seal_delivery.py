from config import *

def run():
    science=json.loads((RESULTS/'SCIENTIFIC_QA_SUMMARY.json').read_text(encoding='utf-8'))
    delivery=json.loads((RESULTS/'DELIVERY_QA_SUMMARY.json').read_text(encoding='utf-8'))
    assert science['failed']==0 and delivery['failed']==0
    log('All scientific and delivery checks passed; stage execution complete, strict/external validation gaps explicitly retained. Sealing immutable delivery.')
    status={'stage':9,'status':'COMPLETE','date':'2026-10-28','timezone':'Asia/Shanghai','known_prior_results_informed_integration':True,'scientific_QA':science,'delivery_QA':delivery,
            'same_direction_association_discoveries':14,'core_candidates':['NPTX1','PDGFRB','KLK6'],'new_RID_peptide_BY5_discoveries':0,
            'verified2014_technical_replicate_pairs':16,'verified2019_technical_replicate_pairs':20,'recovered2019_process_QC_samples':24,
            'new95_target_transition_measurements_QC0':950,'2014_specimen_dates_verified':0,'new2014_pathology_associations_fitted':0,
            'conditional_four_assay_df':{'NPTX1':19,'KLK6':17},'strict_relative_pattern_validation_completed':False,'external_individual_replication_completed':False,'specific_mechanism_identified':False,
            'original_BY_equivalence_discoveries_unchanged':0,'prior_stage_files_modified':False,'created':now()}
    save_json(ROOT/'STAGE9_COMPLETION_STATUS.json',status)
    rows=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or p.name=='SHA256_MANIFEST.csv' or '__pycache__' in p.parts or p.suffix in ['.pyc','.tmp'] or p.name.startswith('~$'):continue
        rel=str(p.relative_to(ROOT)).replace('\\','/')
        rows.append({'relative_path':rel,'bytes':p.stat().st_size,'sha256':sha(p),'restricted':p.is_relative_to(PRIVATE),'category':p.parts[len(ROOT.parts)] if len(p.parts)>len(ROOT.parts)+1 else 'DELIVERY_DOCUMENT'})
    csvout(ROOT/'SHA256_MANIFEST.csv',rows)
    print(json.dumps({'status':'SEALED','files':len(rows),'bytes':sum(r['bytes'] for r in rows),'strict_validation':'NOT_ESTABLISHED','external_replication':'NOT_COMPLETED'},ensure_ascii=False))

if __name__=='__main__':
    assert_unsealed()
    run()
