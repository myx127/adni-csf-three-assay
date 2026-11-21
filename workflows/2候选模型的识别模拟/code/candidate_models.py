from __future__ import annotations

import itertools

import numpy as np
from scipy.special import gammaincc

from sim_config import METHODS


def failure(status):
    return {'p_value':np.nan,'statistic':np.nan,'df':0,'status':status}


def wald(h, covariance):
    h=np.atleast_1d(np.asarray(h,dtype=float))
    covariance=(covariance+covariance.T)/2
    scale=np.sqrt(np.maximum(np.diag(covariance),0))
    if not np.isfinite(covariance).all() or np.any(scale<=1e-15): return failure('SINGULAR_SANDWICH')
    correlation=covariance/np.outer(scale,scale)
    eigen=np.linalg.eigvalsh(correlation)
    if eigen[0]<-1e-8: return failure('INDEFINITE_SANDWICH')
    if eigen[0]<=1e-9: return failure('SINGULAR_SANDWICH')
    z=h/scale
    statistic=float(z@np.linalg.solve(correlation,z))
    return {'p_value':float(gammaincc(len(h)/2,max(0.,statistic)/2)),
            'statistic':max(0.,statistic),'df':len(h),'status':'OK'}


def grouped_covariance(contributions, labels):
    _, inverse=np.unique(labels,axis=0,return_inverse=True)
    k=int(inverse.max())+1
    if k<2: return np.full((contributions.shape[1],)*2,np.nan)
    sums=np.zeros((k,contributions.shape[1]))
    np.add.at(sums,inverse,contributions)
    return sums.T@sums*k/(k-1)


def multiway_covariance(contributions, batches):
    covariance=np.zeros((contributions.shape[1],)*2)
    for size in [1,2,3]:
        for columns in itertools.combinations(range(3),size):
            covariance+=(-1)**(size+1)*grouped_covariance(contributions,batches[:,columns])
    return covariance


def ols_joint(Y,X, tested, hc3=True):
    try:
        inverse=np.linalg.inv(X.T@X)
    except np.linalg.LinAlgError:
        return failure('SINGULAR_DESIGN')
    coefficients=inverse@X.T@Y
    residual=Y-X@coefficients
    weights=X@inverse
    leverage=np.sum(weights*X,axis=1)
    if np.any(leverage>.999999): return failure('EXTREME_LEVERAGE')
    influence=np.concatenate([weights[:,j,None]*residual for j in tested],axis=1)
    if hc3: influence/=1-leverage[:,None]
    h=np.concatenate([coefficients[j] for j in tested])
    return wald(h,influence.T@influence)


def pca_basis(covariance):
    values,vectors=np.linalg.eigh(covariance)
    loading=vectors[:,-1]
    if loading.sum()<0: loading=-loading
    # Deterministic Gram-Schmidt avoids fold-wise arbitrary rotation of D1/D2.
    for raw in [np.array([1.,-1.,0.]),np.array([1.,0.,-1.])]:
        q1=raw-loading*(raw@loading)
        if np.linalg.norm(q1)>1e-8: break
    q1/=np.linalg.norm(q1)
    q2=np.cross(loading,q1)
    return np.column_stack([loading,q1,q2])


def pca_test(Y,P,crossfit=False):
    n=len(Y)
    if not crossfit:
        sd=Y.std(axis=0,ddof=1)
        if np.any(sd<=1e-12): return failure('ZERO_SD')
        Z=(Y-Y.mean(axis=0))/sd
        scores=Z@pca_basis(Z.T@Z/n)
        X=np.column_stack([np.ones(n),P,scores[:,0]])
    else:
        # Input ordering is independently shuffled by generator. Fixed alternating folds.
        fold=np.arange(n)%2
        scores=np.zeros_like(Y)
        for level in [0,1]:
            train=Y[fold!=level]
            sd=train.std(axis=0,ddof=1)
            if np.any(sd<=1e-12): return failure('ZERO_TRAINING_SD')
            trainZ=(train-train.mean(axis=0))/sd
            basis=pca_basis(trainZ.T@trainZ/len(trainZ))
            scores[fold==level]=((Y[fold==level]-train.mean(axis=0))/sd)@basis
        X=np.column_stack([np.ones(n),P,scores[:,0],fold])
    return ols_joint(scores[:,1:],X,[1])


def common_direction_test(Y,P):
    n=len(Y)
    gc=P-P.mean()
    sxx=float(gc@gc)
    if sxx<=1e-12: return failure('NO_PATHOLOGY_VARIATION')
    slope=gc@Y/sxx
    residual=Y-Y.mean(axis=0)-gc[:,None]*slope
    products=np.column_stack([residual[:,0]*residual[:,1],residual[:,0]*residual[:,2],residual[:,1]*residual[:,2]])
    c12,c13,c23=products.mean(axis=0)
    b1,b2,b3=slope
    h=np.array([b2*c13-b1*c23,b3*c12-b1*c23])
    jac=np.array([[-c23,c13,0.,0.,b2,-b1],[-c23,0.,c12,b3,0.,-b1]])
    contributions=np.column_stack([gc[:,None]*residual/sxx,(products-products.mean(axis=0))/n])
    influence=contributions@jac.T
    return wald(h,influence.T@influence*n/(n-2))


def invariance_constraints(moment,weights):
    # moment shape (3 groups, 6 fields): mean1/2/3, covariance12/13/23.
    h=np.zeros(8)
    jac=np.zeros((8,18))
    pooled=weights@moment[:,3:]
    c12,c13,c23=pooled
    for g in [1,2]:
        d=moment[g,:3]-moment[0,:3]
        r=2*(g-1)
        h[r]=d[1]*c13-d[0]*c23
        h[r+1]=d[2]*c12-d[0]*c23
        for group,sign in [(g,1.),(0,-1.)]:
            jac[r,6*group:6*group+3]=sign*np.array([-c23,c13,0.])
            jac[r+1,6*group:6*group+3]=sign*np.array([-c23,0.,c12])
        for group,w in enumerate(weights):
            jac[r,6*group+3:6*group+6]+=w*np.array([0.,d[1],-d[0]])
            jac[r+1,6*group+3:6*group+6]+=w*np.array([d[2],0.,-d[0]])
        a=moment[g,3:]; b=moment[0,3:]
        r=4+2*(g-1)
        h[r]=a[0]*b[1]-a[1]*b[0]
        h[r+1]=a[0]*b[2]-a[2]*b[0]
        jac[r,6*g+3:6*g+6]=[b[1],-b[0],0.]
        jac[r,3:6]=[-a[1],a[0],0.]
        jac[r+1,6*g+3:6*g+6]=[b[2],0.,-b[0]]
        jac[r+1,3:6]=[-a[2],0.,a[0]]
    return h,jac


def invariance_moments(Y,groups):
    n=len(Y)
    moments=np.zeros((3,6))
    contributions=np.zeros((n,18))
    counts=np.bincount(groups,minlength=3)
    for g in range(3):
        mask=groups==g
        mean=Y[mask].mean(axis=0)
        residual=Y[mask]-mean
        products=np.column_stack([residual[:,0]*residual[:,1],residual[:,0]*residual[:,2],residual[:,1]*residual[:,2]])
        moments[g]=np.r_[mean,products.mean(axis=0)]
        contributions[mask,6*g:6*g+6]=np.column_stack([residual,products-products.mean(axis=0)])/counts[g]
    return moments,contributions,counts


def invariance_tests(Y,groups,batches):
    moments,contributions,counts=invariance_moments(Y,groups)
    h,jac=invariance_constraints(moments,counts/len(Y))
    influence=contributions@jac.T
    covariance=influence.T@influence
    out={'INVARIANCE_JOINT_IID':wald(h,covariance),
         'INVARIANCE_MEAN_IID':wald(h[:4],covariance[:4,:4]),
         'INVARIANCE_COV_IID':wald(h[4:],covariance[4:,4:]),
         'INVARIANCE_JOINT_TMT_CLUSTER':wald(h,grouped_covariance(influence,batches[:,[2]])),
         'INVARIANCE_JOINT_3WAY_CLUSTER':wald(h,multiway_covariance(influence,batches))}
    pooled=(counts/len(Y))@moments[:,3:]
    return out,bool(np.any(pooled<=0))


def run_methods(data):
    good=np.isfinite(data['Y']).all(axis=1)
    Y,P,groups,batches=(data[key][good] for key in ['Y','P','groups','batches'])
    n=len(Y)
    counts=np.bincount(groups,minlength=3)
    out={name:failure('NOT_RUN') for name in METHODS}
    oracle_X=np.column_stack([np.ones(len(data['P'])),data['P'],data['P']*data['T']])
    out['ORACLE_COMPLETE_TRUTH']=ols_joint(data['oracle_residual'],oracle_X,[1,2])
    sign_flag=False
    if n<30 or counts.min()<10:
        for name in METHODS:
            if name!='ORACLE_COMPLETE_TRUTH': out[name]=failure('INSUFFICIENT_COMPLETE_GROUPS')
    else:
        out['PCA_IN_SAMPLE']=pca_test(Y,P)
        out['PCA_CROSSFIT_2FOLD']=pca_test(Y,P,True)
        out['COMMON_LOADING_MEAN_WALD']=common_direction_test(Y,P)
        group_results,sign_flag=invariance_tests(Y,groups,batches)
        out.update(group_results)
    for name,result in out.items():
        result['n_complete']=len(data['P']) if name=='ORACLE_COMPLETE_TRUTH' else n
        result['smallest_group_n']=int(counts.min())
        result['nonpositive_pooled_cross_covariance']=sign_flag if 'INVARIANCE' in name else False
    return out
