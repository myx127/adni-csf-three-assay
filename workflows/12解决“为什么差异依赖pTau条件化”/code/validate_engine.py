from config import *
from engine import basis,stats,rcs,D
from scipy.linalg import lstsq
from scipy.interpolate import CubicSpline
from scipy.stats import chi2
import pickle
def run():
 T=pickle.load((PRIVATE/'TEST_TEMPLATES_RESTRICTED.pkl').open('rb'));c=pd.read_csv(PRIVATE/'COMMON_CLINICAL_RESTRICTED.csv');scales=pickle.load((PRIVATE/'FROZEN_SCALES_RESTRICTED.pkl').open('rb'));rng=np.random.default_rng(120090);checks=[]
 def close(name,a,b,atol=1e-9):
  ok=np.allclose(a,b,rtol=1e-7,atol=atol,equal_nan=True);checks.append({'check':name,'passed':bool(ok),'max_abs_error':float(np.max(abs(np.asarray(a,float)-np.asarray(b,float))))})
 for key,t in T.items():
  Y=rng.normal(size=(len(t['y']),5));b=[];scores=[]
  scale=scales[t['accession']]
  for j in range(t['A'].shape[0]):
   X=np.c_[t['ZZ_list'][j],t['targetx_list'][j]];fit=lstsq(X,Y,cond=max(X.shape)*np.finfo(float).eps,lapack_driver='gelsy')[0];e=Y-X@fit;b.append(fit[-1]);row=[]
   for axis,col in enumerate(['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex']):
    levels=sorted(c.iloc[scale['ids']][col].astype(str).unique());lookup={v:i for i,v in enumerate(levels)};labels=np.array([lookup[str(v)] for v in c.iloc[t['ids']][col]])
    sg=np.zeros((len(levels),Y.shape[1]));np.add.at(sg,labels,t['A'][j,:,None]*e);sg-=sg.mean(axis=0);sg*=np.sqrt(len(levels)/(len(levels)-1));row.append(sg)
   scores.append(np.concatenate(row))
  bb=np.array(b).T;U=np.stack(scores,axis=2).transpose(1,0,2);r=stats(Y,t)
  close(key+'/independent_QR_beta',r['b'],bb);close(key+'/direct_residual_cluster_scores',r['scores'],U)
  gs=U*t['sigma'][None,None,:]+bb[:,None,:]*t['Uv'][None,:,:]/(2*t['sigma'][None,None,:]);h=(bb*t['sigma'])@t['L'].T;hh=np.einsum('ap,bgp->bga',t['L'],gs);V=np.einsum('bga,bgc->bac',hh,hh)
  close(key+'/paired_crossmodel_covariance',r['V'],V);close(key+'/contrasts',r['h'],h)
  pp=[]
  for i in range(5):pp.append(chi2.sf(h[i]@np.linalg.solve(V[i],h[i]),2))
  close(key+'/independent_chi2',r['p'],pp)
  p=len(t['sigma']);v=np.r_[bb[0],t['sigma']**2];G=np.c_[t['L']*t['sigma'][None,:],t['L']*(bb[0]/(2*t['sigma']))[None,:]];num=np.empty_like(G)
  def f(x):return t['L']@(x[:p]*np.sqrt(x[p:]))
  for j in range(2*p):
   step=1e-6*max(1.,abs(v[j]));a=v.copy();z=v.copy();a[j]+=step;z[j]-=step;num[:,j]=(f(a)-f(z))/(2*step)
  close(key+'/shared_scale_gradient',G,num,atol=1e-7)
  Ah=t['Ah'];mu=Ah.T@np.linalg.solve(Ah@Ah.T,np.array([.02,-.01]));close(key+'/2d_mean_constraint',Ah@mu,np.array([.02,-.01]))
 knots=np.array(json.loads((RESULTS/'RCS_KNOTS_FROZEN.json').read_text(encoding='utf-8'))['P_LOG_knots'])
 grid=np.linspace(knots[0],knots[-1],101);vv=rng.normal(size=4);ref=CubicSpline(knots,vv,bc_type='natural')(grid);M=np.c_[np.ones(4),rcs(knots,knots)];coef=np.linalg.solve(M,vv)
 close('RCS/natural_cubic_span',np.c_[np.ones(len(grid)),rcs(grid,knots)]@coef,ref,atol=1e-9)
 for v in [knots[0]-2,knots[-1]+2]:
  h=.001;d2=(rcs([v+h],knots)-2*rcs([v],knots)+rcs([v-h],knots))/h**2;close('RCS/linear_tail/'+str(v),d2,0,atol=2e-7)
 csvout(RESULTS/'ENGINE_INDEPENDENT_CHECKS.csv',checks);out={'time':now(),'checks':len(checks),'passed':sum(r['passed'] for r in checks),'failed':sum(not r['passed'] for r in checks),'new_real_extension_targets_fitted':False};save_json(RESULTS/'ENGINE_VALIDATION_SUMMARY.json',out)
 assert out['failed']==0,[r for r in checks if not r['passed']]
 log(f'Independent synthetic validation {len(checks)}/{len(checks)} passed; no new real extension target fitted.')
if __name__=='__main__':assert_unsealed();run()

