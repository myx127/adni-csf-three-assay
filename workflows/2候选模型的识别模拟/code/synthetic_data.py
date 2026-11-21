from __future__ import annotations

import numpy as np
from scipy.special import expit, ndtr, ndtri


def batch_labels(n, rng, audited_sizes):
    result=[]
    for platform in ['NULISA','SOMASCAN','TMT_MS']:
        sizes=np.asarray(audited_sizes[platform],dtype=int)
        if n==int(sizes.sum()):
            labels=np.repeat(np.arange(len(sizes)),sizes)
        else:
            scaled=sizes*n/sizes.sum()
            counts=np.floor(scaled).astype(int)
            counts[np.argsort(-(scaled-counts))[:n-counts.sum()]]+=1
            labels=np.repeat(np.arange(len(sizes)),counts)
        rng.shuffle(labels)
        result.append(labels)
    return np.column_stack(result)


def generate(n, rng, settings, audited_sizes, continuous=False, effect_override=None):
    cfg={'beta':.55,'sd':[.25,.8,2.0],'uvar':[1.,1.,1.],**settings}
    batches=batch_labels(n,rng,audited_sizes)
    if continuous:
        P=rng.normal(size=n)
        groups=np.digitize(P,ndtri([1/3,2/3]))
    else:
        groups=np.arange(n)%3
        rng.shuffle(groups)
        P=groups.astype(float)-1.
    batch_shift=rng.normal(scale=.75,size=int(batches[:,2].max())+1)
    if cfg.get('batch') in ['raw_confounded','gis_confounded']:
        score=batch_shift[batches[:,2]]+rng.normal(scale=.25,size=n)
        order=np.argsort(score)
        groups=np.empty(n,dtype=int)
        groups[order]=np.minimum(2,np.arange(n)*3//n)
        P=groups.astype(float)-1.
    T=cfg['beta']*P+np.sqrt(np.asarray(cfg['uvar'])[groups])*rng.normal(size=n)
    z=rng.normal(size=(n,3))
    rho=cfg.get('rho',0.)
    rho=np.asarray(rho)[groups] if isinstance(rho,list) else rho
    z[:,1]=rho*z[:,0]+np.sqrt(1-np.asarray(rho)**2)*z[:,1]
    eps=z*np.asarray(cfg['sd'])
    eps[:,2]*=np.exp(cfg.get('hetero',0.)*P)
    response=np.repeat(T[:,None],3,axis=1)
    if cfg.get('response')=='quadratic': response[:,2]=T+.35*T*T
    if cfg.get('response')=='saturation': response[:,2]=1.8*np.tanh(T/1.8)
    delta=np.asarray(cfg.get('delta',[0.,0.,0.]),dtype=float)
    if effect_override is not None: delta=np.array([0.,effect_override,0.])
    load=np.asarray(cfg.get('loading_delta',[0.,0.,0.]),dtype=float)
    truth_effect=P[:,None]*delta+(P*T)[:,None]*load
    # Oracle residual has full unselected data and exactly known fixed response/technical effects.
    oracle_residual=eps+truth_effect
    Y=response+oracle_residual
    batch=cfg.get('batch','none')
    if batch=='raw_confounded': Y[:,2]+=batch_shift[batches[:,2]]
    if batch in ['gis_balanced','gis_confounded']:
        gis_noise=rng.normal(scale=.6,size=len(batch_shift))
        Y[:,2]-=gis_noise[batches[:,2]]
    if batch=='crossed':
        for j in range(3):
            shift=rng.normal(scale=.45,size=int(batches[:,j].max())+1)
            Y[:,j]+=shift[batches[:,j]]
    missing=cfg.get('missing','none')
    if missing=='mcar': Y[rng.uniform(size=Y.shape)<.20]=np.nan
    if missing in ['lod30','lod70','floor30']:
        q=.70 if missing=='lod70' else .30
        # A fixed theoretical marginal threshold, not a sample-selected quantile.
        vp=1. if continuous else 2/3
        threshold=np.sqrt(cfg['beta']**2*vp+float(np.mean(cfg['uvar']))+cfg['sd'][0]**2)*ndtri(q)
        if missing=='floor30': Y[:,0]=np.maximum(Y[:,0],threshold)
        else: Y[Y[:,0]<threshold,0]=np.nan
    if missing=='selective':
        observed=rng.uniform(size=n)<expit(1.6-.7*Y[:,2]+.6*P)
        Y[~observed,2]=np.nan
    return {'Y':Y,'P':P,'groups':groups,'batches':batches,'T':T,'oracle_residual':oracle_residual}
