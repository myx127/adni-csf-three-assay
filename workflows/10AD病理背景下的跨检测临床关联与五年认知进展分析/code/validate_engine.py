from config import *
from engine import joint_stats,D
import pickle
from scipy.linalg import lstsq

def run():
 templates=pickle.load((PRIVATE/'ANALYSIS_TEMPLATES_RESTRICTED.pkl').open('rb'));scales=pickle.load((PRIVATE/'EXPOSURE_SCALING_RESTRICTED.pkl').open('rb'))
 clinical=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv');rng=np.random.default_rng(202610290099)
 checks=[]
 def add(label,a,b,rtol=1e-7,atol=2e-9):
  checks.append({'check':label,'passed':bool(np.allclose(a,b,rtol=rtol,atol=atol,equal_nan=True)),'max_abs_error':float(np.max(abs(np.asarray(a)-np.asarray(b))))})
 for key,t in templates.items():
  y=rng.normal(size=t['A'].shape[1]);direct_b=[];direct_score=[]
  for j in range(3):
   T=np.c_[t['ZZ'],t['targetx'][:,j]] if t['modifier'] is None else np.c_[t['ZZ'],t['targetx'][:,j],t['targetx'][:,j]*t['modifier'][t['case_local']]]
   coef=lstsq(T,y,cond=max(T.shape)*np.finfo(float).eps,lapack_driver='gelsy')[0]
   e=y-T@coef;direct_b.append(coef[-1])
   scoreperson=t['A'][j]*e;parts=[]
   scale=scales[t['accession']]
   for axis,c in enumerate(['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex']):
    unique=sorted(clinical.loc[scale['ids'],c].astype(str).unique());lookup={s:i for i,s in enumerate(unique)}
    lab=np.array([lookup[str(clinical.loc[i,c])] for i in t['people']]);sg=np.bincount(lab,weights=scoreperson,minlength=len(unique));sg-=sg.mean();sg*=np.sqrt(len(unique)/(len(unique)-1));parts.append(sg)
   direct_score.append(np.concatenate(parts))
  direct_b=np.array(direct_b);direct_score=np.array(direct_score).T
  add(key+'_coefficients_independent_QR',t['A']@y,direct_b)
  add(key+'_joint_cluster_score_direct_residual',np.einsum('gjn,n->gj',t['K'],y),direct_score)
  stats=joint_stats(y,t)
  add(key+'_raw_joint_covariance',stats['Vraw'][0],direct_score.T@direct_score)
  scaledscore=direct_score*t['sigma'][None,:]+direct_b[None,:]*t['Uv']/(2*t['sigma'][None,:])
  add(key+'_scaled_covariance',stats['Vscaled'][0],scaledscore.T@scaledscore)
  add(key+'_diff_joint_covariance',stats['Vdiff'][0],D@(scaledscore.T@scaledscore)@D.T)
  target=np.array([0.,.03,.05]);As=t['sigma'][:,None]*t['A'];mu=As.T@np.linalg.solve(As@As.T,target)
  add(key+'_arbitrary_marginal_mean_constraint',As@mu,target)
  b0=np.r_[direct_b,t['sigma']**2];G=np.zeros((3,6))
  for j in range(3):G[j,j]=np.sqrt(b0[3+j]);G[j,3+j]=.5*b0[j]/np.sqrt(b0[3+j])
  def f(v):return v[:3]*np.sqrt(v[3:])
  numeric=np.empty((3,6))
  for j in range(6):
   step=1e-6*max(1.,abs(b0[j]));plus=b0.copy();minus=b0.copy();plus[j]+=step;minus[j]-=step;numeric[:,j]=(f(plus)-f(minus))/(2*step)
  add(key+'_SD_joint_gradient',G,numeric,atol=1e-7)
 # Exact longitudinal within-person model equivalence, independent synthetic visit outcomes.
 data=pickle.load((PRIVATE/'ASSEMBLED_CLINICAL_RESTRICTED.pkl').open('rb'))
 for acc in CANDIDATES:
  t=templates['MEM5Y__'+acc];pseudo=[];Y=[];designs=[[] for _ in range(3)]
  for local,i in enumerate(t['ids']):
   tr=data['trajectories'][(i,'ADNI_MEM',5)];tt=tr['t'];tc=tt-tt.mean();yy=rng.normal(size=len(tt));yy-=yy.mean();pseudo.append(tc@yy/np.sqrt(tr['w']));Y.extend(yy)
   Z=t['ZZ'][local]/np.sqrt(tr['w']);x=t['targetx'][local]/np.sqrt(tr['w'])
   for j in range(3):designs[j].append(tc[:,None]*np.r_[Z,x[j]][None,:])
  for j in range(3):
   T=np.concatenate(designs[j]);coef=lstsq(T,np.array(Y),cond=max(T.shape)*np.finfo(float).eps,lapack_driver='gelsy')[0][-1]
   add('VISIT_LEVEL_FIXED_INTERCEPT_EQUIVALENCE_'+acc+'_'+str(j),t['A'][j]@pseudo,coef)
 csvout(RESULTS/'ENGINE_INDEPENDENT_CHECKS.csv',checks)
 assert all(r['passed'] for r in checks),[r for r in checks if not r['passed']][:5]
 save_json(RESULTS/'ENGINE_VALIDATION_SUMMARY.json',{'time':now(),'checks':len(checks),'failed':0,'all_outcomes_synthetic':True,'new_real_target_coefficients_not_fitted':True,'designs':len(templates)})
 log(f'Independent engine validation:{len(checks)}/{len(checks)} synthetic checks passed')
if __name__=='__main__':
 assert_unsealed();run()
