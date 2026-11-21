from config import *
from engine import setup,generate,fit_noise,joint_from_z
from diagnostics import parameters
from scipy.special import ndtr
from freeze_engine import verify
import pickle,itertools

def run():
 verify();templates=pickle.load((PRIVATE/'SIMULATION_TEMPLATES_RESTRICTED.pkl').open('rb'));rng=np.random.default_rng(202610280899);rows=[]
 def check(name,ok,error=0.):rows.append({'check':name,'passed':bool(ok),'error':error})
 for i,(key,t) in enumerate(templates.items()):
  d=t['d']
  for scenario in SCENARIOS:
   spec=setup(t,scenario);Y=generate(spec,2,rng);bn,se,valid=fit_noise(Y,spec)
   maxerr=0.
   for j in range(2):
    theta,_,V,_=parameters(Y[j],d);maxerr=max(maxerr,np.max(abs(bn[j]*t['sigma']-theta[:3])),np.max(abs(se[j]*t['sigma']-np.sqrt(V.diagonal()[:3]))))
   check('reference_coefficient_and_SUM3_'+key+'_'+scenario,maxerr<1e-9 and valid.all(),maxerr)
   effect=np.array([.2,-.4,.7])*t['sigma'];eta=np.linspace(-.2,.3,d['Z'].shape[1]*3).reshape(-1,3)
   shifted=Y+d['P'][None,:,None]*effect[None,None,:]+(d['Z']@eta)[None,:,:];b2,s2,v2=fit_noise(shifted,spec)
   err=max(np.max(abs((b2-bn)*t['sigma']-effect)),np.max(abs(s2-se)))
   check('fixed_X_location_equivariance_'+key+'_'+scenario,err<1e-8 and v2.all(),err)
   check('operator_nuisance_orthogonality_'+key+'_'+scenario,np.max(abs(spec['O']@d['Z']))<1e-8)
   if i<3:
    M=np.eye(d['n'])-d['Q']@d['Q'].T;R=np.diag(spec['hetero']**2)
    for lab,G,q in spec['slopes']:
     A=np.eye(G)[lab]*q[:,None];R+=spec['slope_sd']**2*(A@A.T)
    expected=np.sum(M*R.T)/(d['df']*spec['normalizer']**2);check('independent_expected_residual_covariance_'+key+'_'+scenario,abs(expected-1)<1e-9,abs(expected-1))
 z=rng.normal(size=(100,3))*4;pp,pm,p,sgn=joint_from_z(z)
 alternate=np.where((z>0).all(axis=1)|(z<0).all(axis=1),np.max(2*ndtr(-abs(z)),axis=1),1.)
 check('independent_bidirectional_formula',np.allclose(p,alternate,atol=1e-15))
 for perm in itertools.permutations(range(3)):check('assay_permutation_'+str(perm),np.allclose(joint_from_z(z[:,perm])[2],p))
 check('global_sign_invariance',np.allclose(joint_from_z(-z)[2],p) and np.array_equal(joint_from_z(-z)[3],-sgn))
 for r in PATTERNS:
  v=np.array(r['effect']);isalt=(v>0).all() or (v<0).all();check('pattern_truth_'+r['pattern'],isalt==(r['kind']=='POWER'))
 result=pd.DataFrame(rows);csvout(RESULTS/'ENGINE_INDEPENDENT_CHECKS.csv',result)
 summary={'checks':len(result),'passed':int(result.passed.sum()),'failed':int((~result.passed).sum()),'all105_designs_all5_scenarios_verified':True}
 save_json(RESULTS/'ENGINE_QA_SUMMARY.json',summary);print(summary)
 if not result.passed.all():print(result.loc[~result.passed].to_string(index=False));raise SystemExit(1)
if __name__=='__main__':run()
