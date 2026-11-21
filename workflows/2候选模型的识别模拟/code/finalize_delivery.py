from __future__ import annotations

import csv as csvlib
import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from sim_config import *


def mean_influence_check():
    # Independent weighted empirical contamination derivative for the 2-df test.
    rng=np.random.default_rng(270991);n=120;g=rng.normal(size=n)
    y=rng.normal(size=(n,3))+(g+rng.normal(size=n))[:,None]
    def weighted_h(w):
        w=w/w.sum();gc=g-w@g;mean=w@y
        b=(w*gc)@y/np.sum(w*gc*gc);e=y-mean-gc[:,None]*b
        cov=(e*w[:,None]).T@e;c12,c13,c23=cov[0,1],cov[0,2],cov[1,2]
        return np.array([b[1]*c13-b[0]*c23,b[2]*c12-b[0]*c23])
    gc=g-g.mean();b=gc@y/(gc@gc);e=y-y.mean(axis=0)-gc[:,None]*b
    products=np.c_[e[:,0]*e[:,1],e[:,0]*e[:,2],e[:,1]*e[:,2]]
    c12,c13,c23=products.mean(axis=0)
    j=np.array([[-c23,c13,0.,0.,b[1],-b[0]],[-c23,0.,c12,b[2],0.,-b[0]]])
    contributions=np.c_[gc[:,None]*e/(gc@gc),(products-products.mean(axis=0))/n]@j.T
    error=0.
    for i in [0,2,15,78]:
        unit=np.zeros(n);unit[i]=1.;base=np.ones(n)/n;step=1e-6
        numeric=(weighted_h((1-step)*base+step*unit)-weighted_h((1+step)*base-step*unit))/(2*step)
        error=max(error,float(np.max(abs(numeric-n*contributions[i]))))
    return error


def main():
    checks=[]
    def add(name,passed,detail):checks.append({'check':name,'passed':bool(passed),'detail':detail})
    run=json.loads((RESULTS/'SIMULATION_RUN_SUMMARY.json').read_text(encoding='utf-8'))
    validation=json.loads((RESULTS/'SIMULATION_VALIDATION_QA.json').read_text(encoding='utf-8'))
    access=json.loads((RESULTS/'ACCESS_AND_PROTOCOL_RECORD.json').read_text(encoding='utf-8'))
    add('all_planned_simulations_complete',run['independent_replicates']==67000 and run['method_result_rows']==603000 and run['cells']==49,run)
    add('independent_and_replay_validation',validation['status']=='PASS' and validation['checks_passed']==99,{'passed':validation['checks_passed'],'total':validation['checks_total']})
    add('preanalysis_protocol_unchanged',sha(ROOT/'SIMULATION_PROTOCOL_v1.0.md')==access['protocol_sha256_before_simulation'],access['protocol_sha256_before_simulation'])
    add('clinical_protein_restricted_boundary',all(access[k]==0 for k in ['clinical_files_opened','participant_level_protein_files_opened','restricted_files_opened']),
        'Synthetic-only generation; code/input-boundary evidence, not OS forensic monitoring')
    bound=read_input=pd.read_csv(RESULTS/'INPUT_AGGREGATE_HASH_CHECKS.csv')
    add('stage1_aggregates_unchanged',all(sha(r.path)==r.expected_sha256 for r in bound.itertuples()),len(bound))
    error=mean_influence_check()
    add('independent_contamination_derivative_for_mean_test',error<1e-6,error)
    summary=pd.read_csv(RESULTS/'MONTE_CARLO_SUMMARY.csv')
    denominator_checks=[]
    for cell in design_grid(scenarios()):
        raw=pd.read_parquet(RESULTS/'replicates'/(cell['cell_id']+'.parquet'))
        for method,q in raw.groupby('method',observed=True):
            out=summary.loc[summary.cell_id.eq(cell['cell_id'])&summary.method.eq(method)].iloc[0]
            rejected=int((q.status.eq('OK')&q.p_value.lt(.05)).sum());failed=int(q.status.ne('OK').sum())
            denominator_checks.append(out.n_replicates==len(q) and out.n_rejected==rejected and out.n_failed==failed and abs(out.rejection_rate_all_replicates-rejected/len(q))<1e-12)
    add('all441_summary_counts_from_raw_replicates',len(denominator_checks)==441 and all(denominator_checks),len(denominator_checks))
    add('failure_denominators_retained',summary.n_replicates.eq(summary.n_successful+summary.n_failed).all(),int(summary.n_failed.sum()))
    recovery=pd.read_csv(RESULTS/'PARAMETER_RECOVERY_REPLAY.csv')
    add('8000_replays_not_new_simulations',len(recovery)==8000 and recovery.groupby('scenario_id').size().eq(2000).all(),'Four existing main-experiment seed sets replayed')
    family=pd.read_csv(RESULTS/'OBSERVATIONALLY_EQUIVALENT_PARAMETER_FAMILY.csv')
    add('all_equivalent_family_models_positive_definite_and_equal',family.minimum_error_covariance_eigenvalue.gt(0).all() and
        family.maximum_covariance_reconstruction_error.lt(1e-12).all() and family.maximum_slope_reconstruction_error.lt(1e-12).all(),len(family))
    schema=pd.read_csv(RESULTS/'WORKBOOK_SCHEMA.csv');wb=load_workbook(ROOT/'候选模型识别模拟总表.xlsx',read_only=True,data_only=False)
    rows=[]
    for row in schema.itertuples(index=False):
        ws=wb[row.sheet];rows.append({'sheet':row.sheet,'rows':ws.max_row-1,'columns':ws.max_column,
            'passed':ws.max_row==row.rows+1 and ws.max_column==row.columns})
    add('all18_workbook_dimensions',len(wb.sheetnames)==18 and all(r['passed'] for r in rows),rows)
    ws=wb['实验与种子'];headers=list(next(ws.values));seedcol=headers.index('seed')
    text_seeds=[row[seedcol] for row in list(ws.values)[1:]]
    add('64bit_seeds_stored_as_text',all(isinstance(x,str) for x in text_seeds),'Prevents Excel 15-digit numeric truncation')
    add('workbook_seed_values_match_plan',text_seeds==[c['seed'] for c in design_grid(scenarios())],len(text_seeds))
    wb.close()
    required=['README_阅读顺序.md','MODEL_IDENTIFICATION_SIMULATION_REPORT_ZH.md','MODEL_RECOMMENDATION_AND_FREEZE_GATES.md',
        'OBSERVABLE_CONSTRAINT_DERIVATION.md','PARAMETER_DIAGNOSTIC_SCOPE.md','候选模型识别模拟总表.xlsx','code/run_all.ps1']
    for stem in ['01_MAIN_REJECTION_AND_FAILURE','02_OFFSET_EFFECT_CURVE','03_OBSERVATIONAL_EQUIVALENCE','04_PARAMETER_RECOVERY']:
        required += [f'figures/{stem}.png',f'figures/{stem}.pdf']
    add('required_files_present',all((ROOT/f).is_file() and (ROOT/f).stat().st_size for f in required),required)
    missing=[]
    for path in ROOT.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',path.read_text(encoding='utf-8')):
            if link.startswith('http') or link in ['SIMULATION_INPUT_OUTPUT_SHA256.csv','results/DELIVERY_COMPLETION_QA.json']:continue
            if not (ROOT/link).exists():missing.append({'file':path.name,'link':link})
    add('report_links_resolve',not missing,missing)
    passed=all(r['passed'] for r in checks)
    result={'status':'SIMULATION_COMPLETE_WITH_IDENTIFICATION_LIMITS' if passed else 'DELIVERY_QA_FAILED',
        'completed_at_local':datetime.now().isoformat(),'completion_checks_passed':sum(r['passed'] for r in checks),'completion_checks_total':len(checks),
        'independent_synthetic_replicates':67000,'method_result_rows':603000,'clinical_files_opened':0,
        'participant_level_protein_files_opened':0,'formal_mechanistic_SAP_authorized':False,'final_protein_set_frozen':False,
        'checks':checks,'interpretation':'Delivery and implementation validation passed. Scientific identification restrictions and candidate-method failures remain explicit.'}
    save_json(RESULTS/'DELIVERY_COMPLETION_QA.json',result)
    if not passed:raise RuntimeError('Delivery QA failed; inspect DELIVERY_COMPLETION_QA.json')
    (LOGS/'FINALIZATION.log').write_text(datetime.now().isoformat()+f' {len(checks)}/{len(checks)} completion checks passed; manifest written next.\n',encoding='utf-8')
    manifest=[]
    for row in bound.itertuples():
        path=Path(row.path);manifest.append({'type':'STAGE1_AGGREGATE_INPUT','path':str(path),'bytes':path.stat().st_size,'sha256':sha(path)})
    stage1_manifest=STAGE1/'STAGE1_INPUT_OUTPUT_SHA256.csv'
    manifest.append({'type':'STAGE1_HASH_BINDING_MANIFEST','path':str(stage1_manifest),'bytes':stage1_manifest.stat().st_size,'sha256':sha(stage1_manifest)})
    target=ROOT/'SIMULATION_INPUT_OUTPUT_SHA256.csv'
    for path in sorted(ROOT.rglob('*')):
        if not path.is_file() or path==target or '__pycache__' in path.parts:continue
        manifest.append({'type':'SIMULATION_OUTPUT','path':path.relative_to(ROOT).as_posix(),'bytes':path.stat().st_size,'sha256':sha(path)})
    csv(target,manifest)
    with target.open(encoding='utf-8-sig',newline='') as f:
        written=list(csvlib.DictReader(f))
    for row in written:
        path=ROOT/row['path'] if row['type']=='SIMULATION_OUTPUT' else Path(row['path'])
        if path.stat().st_size!=int(row['bytes']) or sha(path)!=row['sha256']:raise RuntimeError('Output hash verification failed')
    print(json.dumps({'status':result['status'],'checks':f'{len(checks)}/{len(checks)}','manifest_entries_verified':len(written),'method_failures_retained':int(summary.n_failed.sum())},ensure_ascii=False))


if __name__=='__main__':main()
