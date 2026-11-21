from config import *
from scipy.linalg import svd
from scipy.stats import norm,t
import zipfile,pickle

def design(frame,extra=(),tech=('NULISA_PlateID','SOMASCAN_PlateId','TMT_plex'),factors=('ptau_batch',),pcol='P_LOG'):
    n=len(frame);blocks=[np.ones((n,1))]
    for c in ['age_month_precision','sex_code2',*extra]:
        v=frame[c].to_numpy(float);s=v.std();blocks.append(((v-v.mean())/(s if s else 1))[:,None])
    for c in [*tech,*factors]:
        blocks.append(pd.get_dummies(frame[c].astype(str),drop_first=True,dtype=float).to_numpy())
    Z=np.column_stack(blocks);u,s,_=svd(Z,full_matrices=False,check_finite=True)
    rank=int((s>s[0]*max(Z.shape)*np.finfo(float).eps).sum());Qz=u[:,:rank]
    P=frame[pcol].to_numpy(float);rp=P-Qz@(Qz.T@P);ss=float(rp@rp)
    if ss<1e-10 or n<=rank+1:return None
    Q=np.column_stack([Qz,rp/np.sqrt(ss)]);groups=[pd.factorize(frame[c].astype(str),sort=True)[0] for c in tech]
    return {'n':n,'df':n-rank-1,'rank':rank+1,'Q':Q,'a':rp/ss,'groups':groups,'P':P,'ss':ss,'frame_index':frame.index.to_numpy()}

def moments(Y,d):
    E=Y-d['Q']@(d['Q'].T@Y);b=d['a']@Y;df=d['df']
    products=np.column_stack([E[:,0]**2,E[:,1]**2,E[:,2]**2,E[:,0]*E[:,1],E[:,0]*E[:,2],E[:,1]*E[:,2]])
    theta=np.r_[b,products.sum(axis=0)/df]
    IF=np.column_stack([d['a'][:,None]*E,(products-products.mean(axis=0))/df])
    V=np.zeros((9,9))
    for labels in d['groups']:
        g=int(labels.max())+1;sums=np.zeros((g,9));np.add.at(sums,labels,IF);sums-=sums.mean(axis=0)
        V+=sums.T@sums*g/(g-1)
    cNS,cNT,cST=theta[6:];h=np.array([b[1]*cNT-b[0]*cST,b[2]*cNS-b[0]*cST])
    J=np.array([[-cST,cNT,0,0,0,0,0,b[1],-b[0]],[-cST,0,cNS,0,0,0,b[2],0,-b[0]]])
    Vh=J@V@J.T;Q=float(h@np.linalg.solve(Vh,h));se=np.sqrt(np.diag(V)[:3]);z=b/se
    pplus=float(norm.sf(z).max());pminus=float(norm.cdf(z).max())
    return theta,V,E,{'Q':Q,'p_H':np.exp(-Q/2),'p_plus':pplus,'p_minus':pminus,'p_shared':min(1.,2*min(pplus,pminus))}

def simple_assay(y,d):
    E=y-d['Q']@(d['Q'].T@y);b=float(d['a']@y);IF=d['a']*E;v=0.
    for labels in d['groups']:
        g=int(labels.max())+1;sums=np.bincount(labels,weights=IF,minlength=g);sums-=sums.mean();v+=float(sums@sums)*g/(g-1)
    se=np.sqrt(v);tdf=min(int(a.max())+1 for a in d['groups'])-1;sd=np.sqrt(float(E@E)/d['df'])
    return {'b_P':b,'se_SUM2_CR1':se,'p_positive_direction':float(t.sf(b/se,tdf)),'b_P_per_residual_SD':b/sd,'t_df_min_Gminus1':tdf,'residual_SD':sd}

def run():
    clinical=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv')
    raw=np.load(S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz')
    assert np.array_equal(raw['keys'],clinical.key)
    mapping=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
    prior=pd.read_csv(S8/'results/ALL35_THREE_MODEL_SHARED_RESULTS.csv')
    hs=pd.read_csv(S6/'results/COVARIATE_FOUR_MODEL_RESULTS.csv')
    checks=[];rows=[];effects=[]
    def check(name,actual,expected,rtol=2e-7,atol=2e-11):
        ok=bool(np.allclose(actual,expected,rtol=rtol,atol=atol,equal_nan=True))
        checks.append({'check':name,'passed':ok,'max_absolute_difference':float(np.nanmax(abs(np.asarray(actual)-np.asarray(expected))))})
        return ok
    for acc in prior.accession.unique():
        j=list(raw['accessions']).index(acc);sj=list(raw['soma_features']).index(mapping.loc[acc,'soma_features'])
        Yall=np.column_stack([raw['N'][:,j],raw['S'][:,sj],raw['T_GIS'][:,j]])
        mask=clinical.base_primary.astype(bool).to_numpy()&np.isfinite(Yall).all(axis=1)
        Y=Yall[mask];f=clinical.loc[mask]
        for model in MODELS:
            extra=() if model=='M0' else ('apoe4_carrier',) if model=='APOE' else ('education','apoe4_carrier')
            d=design(f,extra=extra);theta,V,E,stats=moments(Y,d)
            ref=prior.loc[prior.accession.eq(acc)&prior.model.eq(model)].iloc[0]
            href=hs.loc[hs.accession.eq(acc)&hs.model.eq(model)].iloc[0]
            check(acc+'_'+model+'_n',d['n'],ref.n)
            check(acc+'_'+model+'_rank',d['rank'],ref.rank_X)
            for key,old in [('p_shared',ref.p_shared_computed),('p_plus',ref.p_plus_computed),('p_minus',ref.p_minus_computed),('Q',href.Q),('p_H',href.p_value)]:
                check(acc+'_'+model+'_'+key,stats[key],old)
            for k,label in enumerate(['N','S','T']):
                check(acc+'_'+model+'_'+label+'_b',theta[k],ref[label+'_b'])
                check(acc+'_'+model+'_'+label+'_se',np.sqrt(V[k,k]),ref[label+'_SE'])
                v=theta[3+k];gradient=np.zeros(9);gradient[k]=1/np.sqrt(v);gradient[3+k]=-.5*theta[k]/v**1.5
                ses=np.sqrt(gradient@V@gradient);effects.append({'accession':acc,'gene':ref.gene,'model':model,'assay':label,'n':d['n'],'b':theta[k],'SE':np.sqrt(V[k,k]),'scaled_b':theta[k]/np.sqrt(v),'scaled_SE_full_joint':ses,'scaled_CI95_low':theta[k]/np.sqrt(v)-norm.ppf(.975)*ses,'scaled_CI95_high':theta[k]/np.sqrt(v)+norm.ppf(.975)*ses})
            rows.append({'accession':acc,'gene':ref.gene,'model':model,'n':d['n'],'rank_X':d['rank'],'residual_df':d['df'],'in_family34':acc in FAMILY,**stats})
    out=pd.DataFrame(rows)
    for col in ['q_BY_shared34','q_BY_H34','q_BY_combined68_shared','q_BY_combined68_H']:out[col]=np.nan
    for model in MODELS:
        ids=out.index[out.model.eq(model)&out.in_family34];part=out.loc[ids];assert len(part)==34
        out.loc[ids,'q_BY_shared34']=by(part.p_shared)
        out.loc[ids,'q_BY_H34']=by(part.p_H)
        both=by(np.r_[part.p_shared.to_numpy(),part.p_H.to_numpy()])
        out.loc[ids,'q_BY_combined68_shared']=both[:34];out.loc[ids,'q_BY_combined68_H']=both[34:]
        for i in ids:
            rr=out.loc[i];pr=prior.loc[prior.accession.eq(rr.accession)&prior.model.eq(model)].iloc[0]
            for col,oldcol in [('q_BY_shared34','q_BY_shared34'),('q_BY_H34','prior_H_q_BY_by_model'),('q_BY_combined68_shared','q_BY_combined68_shared'),('q_BY_combined68_H','q_BY_combined68_H')]:
                check(rr.accession+'_'+model+'_'+col,rr[col],pr[oldcol])
    csvout(RESULTS/'INDEPENDENT_ALL35_THREE_MODEL_REPRODUCTION.csv',out)
    csvout(RESULTS/'INDEPENDENT_PLATFORM_EFFECTS.csv',effects)
    # Rebuild peptide values from raw2019, without reusing processed peptide cache.
    manifests=pd.read_csv(S6/'results/LEGACY_INPUT_MANIFEST.csv').set_index('input_id')
    transitions=pd.read_csv(S6/'results/LEGACY_CANDIDATE_TRANSITION_REGISTRY.csv')
    lframe=pd.read_csv(S6/'00_restricted/LEGACY_L04_CLINICAL_DESIGN_RESTRICTED.csv')
    with zipfile.ZipFile(manifests.loc['L04','path']) as z:
        source=pd.read_csv(z.open('csfproteomics/adni_csfproteomics.csv'),low_memory=False)
    prev=pd.read_csv(S6/'results/LEGACY_NEW_RID_PEPTIDE_PTAU_ASSOCIATION.csv');assoc=[];influence=[]
    ds=pickle.load((S6/'00_restricted/LEGACY_CLINICAL_DESIGNS_RESTRICTED.pkl').open('rb'))
    for (acc,pep),part in transitions.groupby(['accession','peptide'],sort=False):
        cols=part.transition_field.tolist();selected=source.iloc[lframe.source_row.to_numpy(int)-2]
        yy=selected[cols].apply(pd.to_numeric,errors='coerce').mean(axis=1,skipna=False).to_numpy()
        ident='NEW_RID_CLINICAL__'+acc+'__'+pep;ix=ds[ident]['indices'];sub=lframe.loc[ix].copy();y=yy[ix]
        d=design(sub,tech=('MRM_PLATE','MRM_RUN'),factors=('ptau_batch','MRM_INSTRUMENT'))
        stat=simple_assay(y,d);pr=prev.loc[prev.design_id.eq(ident)].iloc[0]
        check('2019_'+pep+'_n',len(ix),pr.n)
        for key in ['b_P','se_SUM2_CR1','p_positive_direction','b_P_per_residual_SD']:
            check('2019_'+pep+'_'+key,stat[key],pr[key])
        assoc.append({'accession':acc,'gene':pr.gene,'peptide':pep,'n':len(ix),**stat,'evidence':'SAME_ADNI_SUBJECT_DISJOINT_ASSOCIATION_ONLY'})
        for axis in ['MRM_PLATE','MRM_RUN']:
            for number,label in enumerate(sorted(sub[axis].astype(str).unique()),1):
                keep=sub[axis].astype(str).ne(label).to_numpy();dd=design(sub.loc[keep],tech=('MRM_PLATE','MRM_RUN'),factors=('ptau_batch','MRM_INSTRUMENT'))
                if dd is None:
                    influence.append({'gene':pr.gene,'peptide':pep,'deleted_axis':axis,'deleted_level_index':number,'n':int(keep.sum()),'status':'P_NOT_ESTIMABLE'})
                else:
                    a=simple_assay(y[keep],dd)
                    influence.append({'gene':pr.gene,'peptide':pep,'deleted_axis':axis,'deleted_level_index':number,'n':dd['n'],'rank_X':dd['rank'],'residual_df':dd['df'],'status':'DESCRIPTIVE_ONLY','b_P':a['b_P'],'scaled_b':a['b_P_per_residual_SD'],'coefficient_positive':a['b_P']>0,'no_selected_P_value':True})
    assoc=pd.DataFrame(assoc);assoc['q_BY5']=by(assoc.p_positive_direction)
    for rr in assoc.itertuples():
        old=prev.loc[prev.cohort.eq('NEW_RID_CLINICAL')&prev.peptide.eq(rr.peptide)].iloc[0]
        check('2019_'+rr.peptide+'_q_BY5',rr.q_BY5,old.q_BY5_single_assay)
    csvout(RESULTS/'INDEPENDENT_2019_NEW_RID_PEPTIDE_ASSOCIATIONS.csv',assoc)
    csvout(RESULTS/'NEW_RID_LEAVE_ONE_BATCH_DESCRIPTIVE.csv',influence)
    inf=pd.DataFrame(influence);summ=[]
    for (gene,pep,axis),part in inf.groupby(['gene','peptide','deleted_axis']):
        valid=part.loc[part.status.eq('DESCRIPTIVE_ONLY')]
        summ.append({'gene':gene,'peptide':pep,'deleted_axis':axis,'n_deletions':len(part),'n_estimable':len(valid),'n_positive':int(valid.coefficient_positive.sum()),'scaled_b_min':valid.scaled_b.min(),'scaled_b_max':valid.scaled_b.max(),'min_remaining_df':valid.residual_df.min(),'inference':'DESCRIPTIVE_DIRECTION_AND_INFLUENCE_NOT_NEW_REPLICATION'})
    csvout(RESULTS/'NEW_RID_LEAVE_ONE_BATCH_SUMMARY.csv',summ)
    csvout(RESULTS/'INDEPENDENT_NUMERIC_CHECKS.csv',checks)
    assert all(x['passed'] for x in checks),[x for x in checks if not x['passed']][:5]
    log(f'Independent SVD/SUM3 reproduction: {len(out)} three-assay designs; five raw2019 peptides; {len(influence)} batch-deletion diagnostics; {len(checks)} checks passed')
    print(assoc[['gene','peptide','n','b_P_per_residual_SD','q_BY5']].to_string(index=False))

if __name__=='__main__':
    assert_unsealed()
    run()
