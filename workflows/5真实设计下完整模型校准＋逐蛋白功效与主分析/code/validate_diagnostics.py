from config import *
from model import batch_components
from diagnostics import *
from simulation import setup,generate
import pickle
ds=pickle.load((RESTRICTED/'DESIGNS_MAIN_RESTRICTED.pkl').open('rb'));d=next(iter(ds.values()));rng=np.random.default_rng(1729)
Y=generate(d,setup(d,'C05'),1,rng,.35)[0];result,theta,IF,V,E=summarize(Y,d,0);fast=batch_components(Y,d)
rows=[{'check':'diagnostic_vs_calibration_Q','passed':bool(np.allclose(result['Q'],fast['Q0'][0],rtol=1e-9))},
 {'check':'diagnostic_vs_calibration_h','passed':bool(np.allclose([result['H1'],result['H2']],fast['h'][0]))}]
for rho in [(0,0,0),(.2,-.1,.1),(-.2,-.2,-.2)]:
 h,J=constraints_jacobian(theta,rho);jn=np.zeros((2,9))
 for j in range(9):
  eps=1e-5*max(1.,abs(theta[j]));p=theta.copy();m=theta.copy();p[j]+=eps;m[j]-=eps
  jn[:,j]=(constraints_jacobian(p,rho)[0]-constraints_jacobian(m,rho)[0])/(2*eps)
 rows.append({'check':'rho_gradient_'+str(rho),'passed':bool(np.allclose(J,jn,rtol=1e-5,atol=1e-8))})
profile=covariance_profile(theta,V);rows.append({'check':'grid_has125_and_zero_matches','passed':len(profile)==125 and np.allclose(next(r['Q'] for r in profile if r['rho_NS']==r['rho_NT']==r['rho_ST']==0),result['Q'])})
csvout(RESULTS/'DIAGNOSTIC_PREFLIGHT_QA.csv',rows);print(pd.DataFrame(rows).to_string(index=False));assert all(r['passed'] for r in rows)
