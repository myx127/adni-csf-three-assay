from config import *
import itertools
from scipy.linalg import lstsq
def theory(theta):
 a=.8;bP=-.06;bU=.10;vF=1+a*a;theta=np.array(theta,float);vX=vF+theta
 g0=(a*bP+bU)/np.sqrt(vX);g1=bU*np.sqrt(vX)/(1+theta)
 return g0,g1,vF/vX
def run():
 rows=[]
 for theta in itertools.product([.1,.5,1.,2.,4.],repeat=3):
  g0,g1,r=theory(theta);row={'theta_N':theta[0],'theta_S':theta[1],'theta_T':theta[2],'one_true_signal_same_loading':True}
  for j,a in enumerate(['N','S','T']):row['gamma_NO_P_'+a]=g0[j];row['gamma_P_'+a]=g1[j];row['population_common_fraction_'+a]=r[j]
  row.update({'NO_P_S_MINUS_N':g0[1]-g0[0],'P_S_MINUS_N':g1[1]-g1[0],'CHANGE_S_MINUS_N':g1[1]-g1[0]-g0[1]+g0[0],'NO_P_T_MINUS_N':g0[2]-g0[0],'P_T_MINUS_N':g1[2]-g1[0],'CHANGE_T_MINUS_N':g1[2]-g1[0]-g0[2]+g0[0]});rows.append(row)
 csvout(RESULTS/'SINGLE_TRUE_SIGNAL_ERROR_VARIANCE_GRID.csv',rows)
 rng=np.random.default_rng(202610291299);n=200000;P=rng.normal(size=n);U=rng.normal(size=n);F=.8*P+U;Y=-.06*P+.1*U+rng.normal(scale=.2,size=n);theta=np.array([3.,.1,4.]);g0,g1,reliability=theory(theta);X=F[:,None]+rng.normal(size=(n,3))*np.sqrt(theta);example=[]
 for j,a in enumerate(['N','S','T']):
  for adjusted,expected in [(False,g0[j]),(True,g1[j])]:
   Z=np.column_stack([np.ones(n),X[:,j]]+([P] if adjusted else []));b=lstsq(Z,Y,lapack_driver='gelsy')[0];res=Y-Z@b;sigma=np.sqrt(1+.8*.8+theta[j]);est=b[1]*sigma;se=np.sqrt((res@res)/(n-Z.shape[1])*np.linalg.inv(Z.T@Z)[1,1])*sigma
   example.append({'assay':a,'pTau_adjusted':adjusted,'theta':theta[j],'analytical_standardized_coefficient':expected,'synthetic_OLS_coefficient':est,'synthetic_standard_error':se,'absolute_error':abs(est-expected),'within_6_MC_SE':abs(est-expected)<6*se+1e-4,'synthetic_people':n,'not_fitted_to_ADNI':True})
 assert all(r['within_6_MC_SE'] for r in example);csvout(RESULTS/'SINGLE_TRUE_SIGNAL_SYNTHETIC_REPRODUCTION.csv',example)
 full=pd.read_csv(RESULTS/'CALIBRATION_AND_POWER_FULL.csv',keep_default_na=False);out=[]
 assert full.mean_kind.eq('NULL').sum()==42*5*7*4
 for keys,d in full.loc[full.mean_kind.eq('NULL')].groupby(['test_id','gene','model','kind','alpha_name'],sort=False):
  r=d.loc[d.Wilson95_high.idxmax()];out.append(dict(zip(['test_id','gene','model','kind','alpha_name'],keys))|{'alpha':r.alpha,'max_null_rate':r.rate,'max_Wilson95_upper':r.Wilson95_high,'worst_scenario':r.scenario,'worst_pattern':r.pattern})
 csvout(RESULTS/'NULL_CALIBRATION_WORST_SUMMARY.csv',out)
 out=[]
 for keys,d in full.loc[full.mean_kind.eq('POWER')].groupby(['test_id','gene','model','pattern','alpha_name'],sort=False):
  r=d.loc[d.rate.idxmin()];out.append(dict(zip(['test_id','gene','model','pattern','alpha_name'],keys))|{'worst_working_rate':r.rate,'Wilson95_low':r.Wilson95_low,'worst_DGP':r.scenario,'BY_sufficient_threshold_not_rank_power':keys[-1]!='NOMINAL','not_clinical_MCID':True})
 csvout(RESULTS/'CONDITIONAL_WORKING_EFFECT_POWER.csv',out)
 save_json(RESULTS/'SYNTHETIC_COUNTEREXAMPLE_RECORD.json',{'time':now(),'F':'0.8*P+U','Y':'-0.06*P+0.10*U+eta','loading_N_S_T':[1,1,1],'eta_SD':.2,'example_error_variances':theta,'grid_values':[.1,.5,1,2,4],'grid_cells':125,'synthetic_people':200000,'mechanism_illustration_not_estimated_reliability_or_validation':True})
 log('Single-true-signal counterexample verified analytically and with synthetic data; finite-calibration summaries produced.')
if __name__=='__main__':assert_unsealed();run()

