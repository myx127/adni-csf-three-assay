from __future__ import annotations

import os
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='1'

import importlib.metadata
import json
import platform
import sys
import time
import traceback

import numpy as np
import pandas as pd

from sim_config import *
from synthetic_data import generate
from candidate_models import run_methods


def bind_inputs():
    manifest=pd.read_csv(STAGE1/'STAGE1_INPUT_OUTPUT_SHA256.csv')
    rows=[]
    for name in ['AUDIT_SUMMARY.json','TRIPLE_BATCH_SAMPLE_SIZES.csv']:
        relative='results/'+name
        expected=manifest.loc[manifest.path.eq(relative)].iloc[0]
        path=STAGE1/relative
        observed=sha(path)
        passed=observed==expected.sha256 and path.stat().st_size==expected.bytes
        rows.append({'path':str(path),'sha256':observed,'expected_sha256':expected.sha256,'passed':passed})
    csv(RESULTS/'INPUT_AGGREGATE_HASH_CHECKS.csv',rows)
    if not all(r['passed'] for r in rows): raise RuntimeError('Stage1 aggregate input identity changed')
    sizes=pd.read_csv(STAGE1/'results/TRIPLE_BATCH_SAMPLE_SIZES.csv')
    audited={p:d.n_triple_participants.to_list() for p,d in sizes.groupby('platform')}
    assert {p:len(v) for p,v in audited.items()}=={'NULISA':19,'SOMASCAN':36,'TMT_MS':65}
    assert all(sum(v)==571 for v in audited.values())
    save_json(RESULTS/'AUDITED_BATCH_SIZES.json',audited)
    access={'clinical_files_opened':0,'participant_level_protein_files_opened':0,'restricted_files_opened':0,
            'aggregate_inputs':rows,'stage1_manifest_read_for_hash_binding_only':True,
            'protocol_sha256_before_simulation':sha(ROOT/'SIMULATION_PROTOCOL_v1.0.md'),
            'seed':SEED,'all_individual_data':'SYNTHETIC','platform_precision_ranking_is_synthetic_only':True}
    save_json(RESULTS/'ACCESS_AND_PROTOCOL_RECORD.json',access)
    return audited


def equivalence_proofs():
    # Two different mechanisms imply exactly the same conditional Gaussian (Y | P).
    lamA=np.ones(3); phiA=1.; psiA=np.diag([.4,.6,1.]); alphaA=.55; deltaA=np.zeros(3)
    lamB=np.array([1.,1.,1.25]); phiB=.8
    psiB=np.array([[.6,.2,0.],[.2,.8,0.],[0.,0.,.75]])
    alphaB=.55; deltaB=np.array([0.,0.,-.1375])
    covA=phiA*np.outer(lamA,lamA)+psiA
    covB=phiB*np.outer(lamB,lamB)+psiB
    betaA=alphaA*lamA+deltaA; betaB=alphaB*lamB+deltaB
    assert np.allclose(covA,covB,atol=1e-14) and np.allclose(betaA,betaB,atol=1e-14)
    # Known lambda: 3 observable pathology slopes cannot identify alpha + 3 direct effects.
    jac=np.column_stack([lamA,np.eye(3)])
    nullvector=np.r_[1.,-lamA]
    proof={'observational_equivalence':{
        'model_A':{'lambda':lamA,'latent_residual_variance':phiA,'error_covariance':psiA,'alpha':alphaA,'direct_effect':deltaA},
        'model_B':{'lambda':lamB,'latent_residual_variance':phiB,'error_covariance':psiB,'alpha':alphaB,'direct_effect':deltaB},
        'same_observed_pathology_slope':betaA,'same_conditional_covariance':covA,
        'maximum_covariance_difference':float(np.max(np.abs(covA-covB))),
        'maximum_slope_difference':float(np.max(np.abs(betaA-betaB))),
        'model_B_error_correlation_12':float(.2/np.sqrt(.6*.8)),
        'minimum_error_covariance_eigenvalue_A':float(np.linalg.eigvalsh(psiA).min()),
        'minimum_error_covariance_eigenvalue_B':float(np.linalg.eigvalsh(psiB).min()),
        'implication':'Same complete conditional Gaussian distribution for every P. No observed-data algorithm can distinguish these two models without additional assumptions or information.'},
        'common_vs_direct_alias':{'jacobian':jac,'rank':int(np.linalg.matrix_rank(jac)),
            'unknown_parameters':4,'nullvector':nullvector,'maximum_J_nullvector':float(np.max(np.abs(jac@nullvector))),
            'reparameterization':'alpha*=alpha+c; delta*=delta-c*lambda leaves beta=lambda*alpha+delta unchanged'},
        'clinical_files_opened':0}
    save_json(RESULTS/'EXACT_IDENTIFICATION_AND_EQUIVALENCE.json',proof)


def main():
    log('START synthetic identification simulation; no clinical or participant-level protein inputs')
    audited=bind_inputs()
    equivalence_proofs()
    registry=scenarios(); grid=design_grid(registry)
    registry_table=pd.DataFrame(registry)
    registry_table['settings']=registry_table.settings.map(lambda x:json.dumps(x,ensure_ascii=False,sort_keys=True))
    csv(RESULTS/'SCENARIO_REGISTRY.csv',registry_table)
    csv(RESULTS/'EXPERIMENT_GRID.csv',grid)
    save_json(RESULTS/'SOFTWARE_ENVIRONMENT.json',{'python':sys.version,'platform':platform.platform(),
        'packages':{key:importlib.metadata.version(key) for key in ['numpy','scipy','pandas','pyarrow','matplotlib','openpyxl']},
        'threads_per_process':1})
    rawdir=RESULTS/'replicates'; rawdir.mkdir(exist_ok=True)
    lookup={r['scenario_id']:r for r in registry}
    start=time.time(); total=0
    for cell in grid:
        row=lookup[cell['scenario_id']]
        results=[]
        for rep in range(cell['repetitions']):
            rng=np.random.default_rng(np.random.SeedSequence([int(cell['seed']),rep]))
            data=generate(cell['n'],rng,row['settings'],audited,cell['continuous_P'],cell['offset_effect_override'])
            estimates=run_methods(data)
            for method,result in estimates.items():
                results.append({'replicate':rep,'method':method,**result})
        frame=pd.DataFrame(results)
        for column in ['method','status']: frame[column]=frame[column].astype('category')
        frame.to_parquet(rawdir/(cell['cell_id']+'.parquet'),index=False)
        total+=cell['repetitions']
        log(f"Completed {cell['cell_id']}: R={cell['repetitions']}, total={total}/67000, elapsed={time.time()-start:.1f}s")
    save_json(RESULTS/'SIMULATION_RUN_SUMMARY.json',{'status':'NUMERICAL_SIMULATIONS_COMPLETE','cells':len(grid),
        'independent_replicates':total,'method_result_rows':total*len(METHODS),'methods':METHODS,
        'elapsed_seconds':time.time()-start,'clinical_files_opened':0,'final_SAP_frozen':False,
        'protocol_sha256_after_simulation':sha(ROOT/'SIMULATION_PROTOCOL_v1.0.md')})
    log('COMPLETE all planned simulations')


if __name__=='__main__':
    try: main()
    except Exception:
        (LOGS/'SIMULATION_FAILURE_TRACEBACK.txt').write_text(traceback.format_exc(),encoding='utf-8')
        raise
