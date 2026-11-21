from config import *
from model import *
from simulation import setup,generate
import pickle
from scipy.linalg import lstsq
checks=[]
def check(name,ok,actual=None):checks.append({'check':name,'passed':bool(ok),'actual':actual})
ds=pickle.load((RESTRICTED/'DESIGNS_MAIN_RESTRICTED.pkl').open('rb'));rng=np.random.default_rng(2026102799)
d=next(iter(ds.values()));n=d['n'];Y=generate(d,setup(d,'C05'),3,rng,delta=.35)
out=batch_components(Y,d);X=np.c_[d['Z'],d['P']]
for j in range(3):
 coef,_,_,_=lstsq(X,Y[j]);res=Y[j]-X@coef
 check('independent_OLS_b_'+str(j),np.allclose(coef[-1],out['b'][j],rtol=1e-8,atol=1e-9))
 check('independent_OLS_residual_'+str(j),np.allclose(res,out['E'][j],atol=1e-8))
 c=res.T@res/d['df'];b=coef[-1];h=np.array([b[1]*c[0,2]-b[0]*c[1,2],b[2]*c[0,1]-b[0]*c[1,2]])
 check('independent_constraints_'+str(j),np.allclose(h,out['h'][j],atol=1e-9))
 V=np.zeros((2,2))
 for axis in range(3):
  labs=d['groups'][:,axis];s=np.array([out['U'][j,labs==v].sum(axis=0) for v in np.unique(labs)]);s-=s.mean(axis=0);V+=s.T@s*len(s)/(len(s)-1)
 check('independent_cluster_cov_'+str(j),np.allclose(V,out['V0'][j],atol=1e-10))
 check('independent_Q_'+str(j),np.allclose(h@np.linalg.inv(V)@h,out['Q0'][j],rtol=1e-8))
# Perturb empirical weights, jointly refit all nuisance parameters and residual covariance.
Y0=Y[0];eps=1e-4
def weighted_h(weights):
 sw=np.sqrt(weights);coef,_,_,_=lstsq(X*sw[:,None],Y0*sw[:,None]);r=Y0-X@coef
 c=(r*weights[:,None]).T@r/weights.sum()*n/d['df'];b=coef[-1]
 return np.array([b[1]*c[0,2]-b[0]*c[1,2],b[2]*c[0,1]-b[0]*c[1,2]])
for i in [0,n//2,n-1]:
 plus=np.ones(n);minus=np.ones(n);plus[i]+=eps;minus[i]-=eps
 numeric=(weighted_h(plus)-weighted_h(minus))/(2*eps)
 check('joint_IF_weight_derivative_'+str(i),np.allclose(numeric,out['U'][0,i],rtol=2e-4,atol=1e-8),float(np.max(abs(numeric-out['U'][0,i]))))
scaled=Y*np.array([.2,4.,2.])+np.array([10.,-7.,3.]);check('affine_scale_Q_invariance',np.allclose(batch_components(scaled,d)['Q0'],out['Q0'],rtol=1e-8))
one=batch_components(Y[:1],d)['Q0'][0];check('batched_matches_single',np.allclose(one,out['Q0'][0]))
try:
 from statsmodels.stats.multitest import multipletests
 pp=np.array([.001,.03,.2,.8,.0001]);check('BH_reference',np.allclose(bh(pp),multipletests(pp,method='fdr_bh')[1]));check('BY_reference',np.allclose(bh(pp,True),multipletests(pp,method='fdr_by')[1]));check('Holm_reference',np.allclose(holm(pp),multipletests(pp,method='holm')[1]))
except ImportError:check('multiple_testing_reference_available',False)
for ident,des in ds.items():
 check('design_orthonormal_'+ident,np.allclose(des['Q'].T@des['Q'],np.eye(des['rank']),atol=1e-8))
 check('P_estimating_equation_'+ident,np.allclose(des['a']@des['P'],1) and np.max(abs(des['a']@des['Z']))<1e-8)
 check('positive_df_'+ident,des['df']>=100)
df=pd.DataFrame(checks);csvout(RESULTS/'PREFLIGHT_QA.csv',df)
save_json(RESULTS/'PREFLIGHT_QA_SUMMARY.json',{'checks':len(df),'passed':int(df.passed.sum()),'failed':int((~df.passed).sum()),'clinical_protein_Y_used':False})
print(df.loc[~df.passed].to_string(index=False) if not df.passed.all() else f'PREFLIGHT {len(df)}/{len(df)} passed')
assert df.passed.all()
