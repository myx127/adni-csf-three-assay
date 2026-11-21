from __future__ import annotations

import os
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[key]='1'
import json
import numpy as np
import pandas as pd
from sim_config import *
from synthetic_data import generate


def main():
    audited=json.loads((RESULTS/'AUDITED_BATCH_SIZES.json').read_text(encoding='utf-8'))
    lookup={r['scenario_id']:r for r in scenarios()};rows=[]
    for cell in design_grid(scenarios()):
        if cell['phase']!='MAIN' or cell['scenario_id'] not in ['N03','N05','N07','N19']:continue
        for rep in range(cell['repetitions']):
            d=generate(cell['n'],np.random.default_rng(np.random.SeedSequence([int(cell['seed']),rep])),lookup[cell['scenario_id']]['settings'],audited)
            P,Y=d['P'],d['Y'];g=P-P.mean();b=g@Y/(g@g);e=Y-Y.mean(axis=0)-g[:,None]*b
            covariance=e.T@e/len(e);c12,c13,c23=covariance[0,1],covariance[0,2],covariance[1,2]
            products=np.column_stack([e[:,0]*e[:,1],e[:,0]*e[:,2],e[:,1]*e[:,2]])
            inf=(products-products.mean(axis=0))/len(e)
            l2=c23/c13;l3=c23/c12;phi=c12*c13/c23
            grad2=np.array([0.,-c23/c13**2,1/c13]);grad3=np.array([-c23/c12**2,0.,1/c12])
            se2=np.sqrt(np.sum((inf@grad2)**2));se3=np.sqrt(np.sum((inf@grad3)**2))
            uniqueness=np.diag(covariance)-phi*np.array([1.,l2,l3])**2
            rows.append({'scenario_id':cell['scenario_id'],'replicate':rep,'relative_loading2_estimate':l2,'relative_loading3_estimate':l3,
                         'loading2_se_delta':se2,'loading3_se_delta':se3,'loading2_interval_covers_true1':abs(l2-1)<=1.96*se2,
                         'loading3_interval_covers_true1':abs(l3-1)<=1.96*se3,'latent_residual_variance_estimate':phi,
                         'any_negative_uniqueness':bool(np.any(uniqueness<0)),
                         'nonpositive_cross_covariance':bool(min(c12,c13,c23)<=0),
                         'minimum_uniqueness_estimate':float(uniqueness.min())})
    data=pd.DataFrame(rows);csv(RESULTS/'PARAMETER_RECOVERY_REPLAY.csv',data)
    summary=[]
    for sid,q in data.groupby('scenario_id'):
        for field in ['relative_loading2_estimate','relative_loading3_estimate']:
            values=q[field];prefix='loading2' if 'loading2' in field else 'loading3'
            summary.append({'scenario_id':sid,'parameter':field,'synthetic_true_value':1.,'n_replays':len(q),
                            'mean_estimate':float(values.mean()),'median_estimate':float(values.median()),'mean_bias':float(values.mean()-1),
                            'RMSE':float(np.sqrt(np.mean((values-1)**2))),'p01':float(values.quantile(.01)),'p99':float(values.quantile(.99)),
                            'delta_interval_coverage':float(q[prefix+'_interval_covers_true1'].mean()),
                            'negative_uniqueness_rate':float(q.any_negative_uniqueness.mean()),
                            'nonpositive_cross_covariance_rate':float(q.nonpositive_cross_covariance.mean()),
                            'caution':'Conditional-independence estimator; N07 violates that assumption. Weak-signal ratios can be unstable; no trimming.'})
    csv(RESULTS/'PARAMETER_RECOVERY_SUMMARY.csv',summary)
    family=[]
    target_cov=np.ones((3,3))+np.diag([.4,.6,1.])
    for phi in np.linspace(.75,1.10,36):
        lam=np.array([1.,1.,1/phi]);psi=target_cov-phi*np.outer(lam,lam)
        delta=np.repeat(.55,3)-.55*lam
        family.append({'latent_residual_variance':phi,'relative_loading3':lam[2],
                       'error_correlation12':psi[0,1]/np.sqrt(psi[0,0]*psi[1,1]),'direct_pathology_effect3':delta[2],
                       'minimum_error_covariance_eigenvalue':float(np.linalg.eigvalsh(psi).min()),
                       'maximum_covariance_reconstruction_error':float(np.max(abs(phi*np.outer(lam,lam)+psi-target_cov))),
                       'maximum_slope_reconstruction_error':float(np.max(abs(.55*lam+delta-.55)))})
    csv(RESULTS/'OBSERVATIONALLY_EQUIVALENT_PARAMETER_FAMILY.csv',family)
    print(pd.DataFrame(summary)[['scenario_id','parameter','median_estimate','RMSE','negative_uniqueness_rate','nonpositive_cross_covariance_rate']].to_string(index=False))


if __name__=='__main__':main()
