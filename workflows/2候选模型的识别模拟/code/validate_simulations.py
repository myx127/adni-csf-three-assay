from __future__ import annotations

import os
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']: os.environ[key]='1'
import argparse
import json
import numpy as np
import pandas as pd
from scipy.stats import chi2
from sim_config import ROOT,RESULTS,METHODS,save_json,scenarios,design_grid
from candidate_models import invariance_constraints,invariance_moments,ols_joint,multiway_covariance,run_methods
from synthetic_data import generate


def preflight():
    checks=[]
    def add(name,passed,detail):
        checks.append({'check':name,'passed':bool(passed),'detail':detail})
    lam=np.array([1.,.8,1.2]);means=np.array([-.5,.1,.8]);variances=np.array([.6,1.,1.4])
    moment=np.zeros((3,6))
    moment[:,:3]=np.array([.2,-.1,.3])+means[:,None]*lam
    for g in range(3):
        cov=variances[g]*np.outer(lam,lam)+np.diag([.3+g*.1,.7,1.2])
        moment[g,3:]=cov[0,1],cov[0,2],cov[1,2]
    h,jac=invariance_constraints(moment,np.array([.3,.3,.4]))
    add('population_null_constraints_exact',np.max(np.abs(h))<1e-12,float(np.max(np.abs(h))))
    add('population_constraint_jacobian_rank8',np.linalg.matrix_rank(jac)==8,int(np.linalg.matrix_rank(jac)))
    maxerror=0.
    rng=np.random.default_rng(27110)
    for _ in range(20):
        m=moment+rng.normal(scale=.1,size=moment.shape)
        w=np.array([.3,.3,.4]);h,j=invariance_constraints(m,w)
        numeric=np.zeros_like(j)
        for col in range(18):
            plus=m.copy().reshape(-1);minus=plus.copy();plus[col]+=1e-6;minus[col]-=1e-6
            numeric[:,col]=(invariance_constraints(plus.reshape(3,6),w)[0]-invariance_constraints(minus.reshape(3,6),w)[0])/2e-6
        maxerror=max(maxerror,float(np.max(np.abs(numeric-j))))
    add('20_constraint_jacobians_vs_central_difference',maxerror<1e-7,maxerror)
    X=np.column_stack([np.ones(150),rng.normal(size=(150,3))]);Y=rng.normal(size=(150,2))
    result=ols_joint(Y,X,[1])
    b=np.linalg.lstsq(X,Y,rcond=None)[0];e=Y-X@b
    bread=np.linalg.solve(X.T@X,np.eye(X.shape[1]));lever=np.diag(X@bread@X.T)
    variance=np.zeros((2,2))
    for i in range(len(Y)):
        v=float((bread@X[i])[1])/(1-lever[i])*e[i]
        variance+=np.outer(v,v)
    stat=float(b[1]@np.linalg.solve(variance,b[1]))
    add('OLS_HC3_joint_Wald_independent_loop',abs(stat-result['statistic'])<1e-9,{'expected':stat,'observed':result['statistic']})
    add('Wald_p_vs_scipy_chi2',abs(chi2.sf(stat,2)-result['p_value'])<1e-12,result['p_value'])
    U=rng.normal(size=(40,4));clusters=np.column_stack([rng.integers(0,k,40) for k in [6,8,10]])
    independent=np.zeros((4,4))
    import itertools
    for size in range(1,4):
        for columns in itertools.combinations(range(3),size):
            sums={}
            for value,label in zip(U,clusters[:,columns]):
                key=tuple(np.atleast_1d(label));sums[key]=sums.get(key,np.zeros(4))+value
            independent+=(-1)**(size+1)*sum(np.outer(v,v) for v in sums.values())*len(sums)/(len(sums)-1)
    err=float(np.max(np.abs(independent-multiway_covariance(U,clusters))))
    add('multiway_sandwich_independent_dictionary_loop',err<1e-9,err)
    # Independent list/matrix calculation of the observational equivalence.
    covariance_A=np.ones((3,3))+np.diag([.4,.6,1.])
    covariance_B=.8*np.outer([1.,1.,1.25],[1.,1.,1.25])+np.array([[.6,.2,0],[.2,.8,0],[0,0,.75]])
    slope_B=.55*np.array([1.,1.,1.25])+[0,0,-.1375]
    add('equivalent_gaussian_models_same_observed_distribution',np.allclose(covariance_A,covariance_B) and np.allclose(slope_B,.55),
        {'covariance_max_error':float(np.max(abs(covariance_A-covariance_B))),'slope_max_error':float(np.max(abs(slope_B-.55)))})
    return checks


def main(preflight_only=False):
    checks=preflight()
    if not preflight_only:
        grid=design_grid(scenarios());lookup={r['scenario_id']:r for r in scenarios()}
        audited=json.loads((RESULTS/'AUDITED_BATCH_SIZES.json').read_text(encoding='utf-8'))
        total_rows=0;total_reps=0;failures={};cell_checks=[]
        for cell in grid:
            frame=pd.read_parquet(RESULTS/'replicates'/(cell['cell_id']+'.parquet'))
            valid=frame.status.eq('OK')
            passed=(len(frame)==cell['repetitions']*len(METHODS) and
                    frame.groupby('method',observed=True).size().eq(cell['repetitions']).all() and
                    not frame.duplicated(['replicate','method']).any() and
                    frame.loc[valid,'p_value'].between(0,1).all() and frame.loc[~valid,'p_value'].isna().all())
            cell_checks.append({'cell_id':cell['cell_id'],'passed':bool(passed),'rows':len(frame)})
            total_rows+=len(frame);total_reps+=cell['repetitions']
            for key,value in frame.loc[~valid,'status'].value_counts().items():
                if value: failures[str(key)]=failures.get(str(key),0)+int(value)
            if cell['phase']=='MAIN' and cell['scenario_id'] in ['N03','N07','N09','N17','A01']:
                for rep in [0,3]:
                    data=generate(cell['n'],np.random.default_rng(np.random.SeedSequence([int(cell['seed']),rep])),lookup[cell['scenario_id']]['settings'],audited,
                                  cell['continuous_P'],cell['offset_effect_override'])
                    rerun=run_methods(data)
                    stored=frame.loc[frame.replicate.eq(rep)].set_index('method')
                    for method in METHODS:
                        a,b=stored.loc[method],rerun[method]
                        same=a.status==b['status'] and (np.isnan(a.p_value) and np.isnan(b['p_value']) or abs(a.p_value-b['p_value'])<1e-12)
                        checks.append({'check':f'seed_replay_{cell["scenario_id"]}_{rep}_{method}','passed':bool(same),'detail':'Full generator and method pipeline replay'})
        checks.append({'check':'all49_replicate_cells_complete','passed':all(r['passed'] for r in cell_checks),'detail':cell_checks})
        checks.append({'check':'67000_repeats_603000_method_rows','passed':total_reps==67000 and total_rows==603000,'detail':{'replicates':total_reps,'rows':total_rows}})
    else:
        failures={}
    result={'status':'PASS' if all(r['passed'] for r in checks) else 'FAIL','checks_passed':sum(r['passed'] for r in checks),'checks_total':len(checks),
            'clinical_files_opened':0,'expected_explicit_method_failures':failures,'checks':checks,
            'scope':'Independent algebra/Jacobian/OLS/cluster checks; seed replay shares production generator and is a reproducibility check, not independent scientific validation'}
    save_json(RESULTS/('PREFLIGHT_QA.json' if preflight_only else 'SIMULATION_VALIDATION_QA.json'),result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['checks']},ensure_ascii=False))
    if result['status']!='PASS': raise RuntimeError('Simulation validation failed')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--preflight',action='store_true');args=parser.parse_args()
    main(args.preflight)
