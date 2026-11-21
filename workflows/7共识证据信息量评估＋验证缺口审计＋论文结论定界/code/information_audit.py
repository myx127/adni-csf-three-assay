from config import *
from diagnostics import parameters,constraints_jacobian
from equivalence import actual_equivalence
from scipy.stats import norm
import pickle

def load():
 ds=pickle.load((S5/'00_restricted/DESIGNS_MAIN_RESTRICTED.pkl').open('rb'));raw=np.load(S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz');meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
 return ds,raw,meta

def y_for(raw,meta,acc,d):
 j=list(raw['accessions']).index(acc);sj=list(raw['soma_features']).index(meta.loc[acc,'soma_features'])
 return np.c_[raw['N'][:,j],raw['S'][:,sj],raw['T_GIS'][:,j]][d['indices']]

def h_and_den(theta):
 A,JA=constraints_jacobian(theta);B=np.array([theta[7]*np.sqrt(theta[4]),theta[6]*np.sqrt(theta[5])]);JB=np.zeros((2,9))
 JB[0,7]=np.sqrt(theta[4]);JB[0,4]=.5*theta[7]/np.sqrt(theta[4]);JB[1,6]=np.sqrt(theta[5]);JB[1,5]=.5*theta[6]/np.sqrt(theta[5])
 return A,JA,B,JB

def fieller(A,B,va,vb,cab,alpha):
 z=norm.isf(alpha);aa=B*B-z*z*vb;bb=-2*(A*B-z*z*cab);cc=A*A-z*z*va;disc=bb*bb-4*aa*cc
 bounded=bool(aa>0 and disc>=0 and B>0)
 roots=sorted([(-bb-np.sqrt(disc))/(2*aa),(-bb+np.sqrt(disc))/(2*aa)]) if bounded else [np.nan,np.nan]
 return bounded,*roots

def p_iut(A,B,JA,JB,V,margin):
 values=[]
 for j in range(2):
  for sign in [-1,1]:
   g=margin*JB[j]+sign*JA[j];mean=margin*B[j]+sign*A[j]
   values.append(norm.sf(mean/np.sqrt(g@V@g)))
 return float(max(values)) if min(B)>0 else 1.

def run():
 ds,raw,meta=load();s6=pd.read_csv(S6/'results/OPERATIONAL_EQUIVALENCE_ALL35.csv').set_index(['accession','margin_working_SD']);rows=[];intervals=[];precision=[];templates={};moments={};reproduction=[]
 for acc,d in ds.items():
  Y=y_for(raw,meta,acc,d);theta,IF,V,E=parameters(Y,d);A,JA,B,JB=h_and_den(theta);gene=meta.loc[acc,'nulisa_target'];ref=s6.loc[(acc,.2)]
  c=theta[6:9];pos=bool((c>0).all());lam=np.sqrt([c[0]*c[1]/c[2],c[0]*c[2]/c[1],c[1]*c[2]/c[0]]) if pos else np.full(3,np.nan);psi=theta[3:6]-lam*lam
  compatible=bool(pos and (psi>=-1e-12).all());assert compatible==bool(ref.simple_cov_structure_point_compatible)
  if compatible:
   assert (psi>=0).all();templates[acc]={'d':d,'lambda':lam,'psi':psi,'C':E.T@E/d['df'],'gene':gene}
  moments[acc]={'theta':theta,'IF':IF,'V':V,'A':A,'JA':JA,'B':B,'JB':JB,'d':d}
  cov_required=max((1.959964*np.sqrt(V[6+j,6+j])/c[j])**2 for j in range(3)) if pos else np.inf
  fieller_req=max((norm.isf(.05)*np.sqrt(JB[j]@V@JB[j])/B[j])**2 for j in range(2)) if min(B)>0 else np.inf
  reasons=[]
  if not compatible:reasons.append('POINT_SINGLE_FACTOR_STRUCTURE_INCOMPATIBLE')
  if not ref.all_three_cov_CI95_positive:reasons.append('CROSS_COVARIANCE_CI_WEAK_OR_NONPOSITIVE')
  if not ref.Fieller90_both_bounded:reasons.append('RELATIVE_DENOMINATOR_FIELLER90_UNBOUNDED')
  if not ref.boundary_calibration_passed:reasons.append('FINITE_BOUNDARY_CALIBRATION_FAILED')
  for j,name in enumerate(['S_relative_to_N','T_relative_to_N']):
   va=JA[j]@V@JA[j];vb=JB[j]@V@JB[j];cab=JA[j]@V@JB[j]
   for label,alpha in [('NOMINAL',.05),('BY_RANK1_SUFFICIENT',ALPHA_BY_SUFFICIENT)]:
    bounded,lo,hi=fieller(A[j],B[j],va,vb,cab,alpha)
    intervals.append({'accession':acc,'gene':gene,'in_family34':acc in FAMILY,'contrast':name,'threshold':label,'one_sided_alpha':alpha,
     'interval_confidence':1-2*alpha,'delta_point':A[j]/B[j] if B[j]!=0 else np.nan,'bounded':bounded,'lower':lo,'upper':hi,
     'width':hi-lo if bounded else np.nan,'denominator_signal_to_SE':B[j]/np.sqrt(vb),'min_symmetric_margin_enclosing_interval':max(abs(lo),abs(hi)) if bounded else np.nan,
     'not_new_biological_margin':True})
  for margin in MARGINS:
   r=s6.loc[(acc,margin)];p=p_iut(A,B,JA,JB,V,margin);old=actual_equivalence(Y,d,margin)['p_joint_equivalence']
   reproduction.append({'accession':acc,'margin':margin,'abs_difference_Stage6':abs(p-r.p_equiv_computed),'abs_difference_reference_function':abs(p-old),'passed':np.allclose(p,r.p_equiv_computed,atol=1e-10) and np.allclose(p,old,atol=1e-10)})
   delta=A/B;inside=bool((B>0).all() and (abs(delta)<margin).all());nominal=bool(p<.05);by=bool(r.operational_equivalence_supported_BY)
   state='EXACT_ZERO_DEVIATION_AND_WITHIN_MARGIN_EQUIVALENCE_SUPPORTED' if r.S5_BY_discovery and by else 'DEVIATION_DETECTED_FROZEN_STAGE5' if r.S5_BY_discovery else 'EQUIVALENCE_SUPPORTED_FROZEN_STAGE6' if by else 'UNRESOLVED'
   bottleneck='GATE_INELIGIBLE' if not r.eligible_for_operational_equivalence else 'NOMINAL_EQUIV_NOT_BY' if nominal and not by else 'POINT_OUTSIDE_WORKING_MARGIN' if not inside else 'INTERVAL_PRECISION_INSUFFICIENT'
   rows.append({'accession':acc,'gene':gene,'margin':margin,'in_family34':acc in FAMILY,'n':d['n'],'residual_df':d['df'],
    'N_G':len(d['onehots'][0]),'S_G':len(d['onehots'][1]),'T_G':len(d['onehots'][2]),'point_structure_compatible':compatible,
    'positive_cross_cov_CI_all3':bool(r.all_three_cov_CI95_positive),'Fieller90_both_bounded':bool(r.Fieller90_both_bounded),
    'equivalence_eligible_Stage6':bool(r.eligible_for_operational_equivalence),'gate_reasons':';'.join(reasons),'delta_S':delta[0],'delta_T':delta[1],
    'point_inside_both_margins':inside,'p_equiv_computed':p,'nominal_computed_p05':nominal,'nominal_equiv_eligible':nominal and bool(r.eligible_for_operational_equivalence),
    'q_BY_Stage6':r.q_BY_equiv_exploratory,'BY_equivalence_Stage6':by,'frozen_exact_zero_deviation_BY':bool(r.S5_BY_discovery),'evidence_state':state,'main_bottleneck':bottleneck})
   for label,alpha in [('NOMINAL',.05),('BY_RANK1_SUFFICIENT',ALPHA_BY_SUFFICIENT)]:
    z=norm.isf(alpha);kr=[]
    for j in range(2):
     for sign in [-1,1]:
      g=margin*JB[j]+sign*JA[j];distance=margin*B[j]+sign*A[j]
      kr.append((z*np.sqrt(g@V@g)/distance)**2 if distance>0 else np.inf)
    numerical=max(kr);required=max(1.,numerical,cov_required,fieller_req) if compatible else np.inf
    for floor in [0.,.25,.5]:
     desired_factor=1/required if np.isfinite(required) else 0.
     k=(1-floor)/(desired_factor-floor) if np.isfinite(required) and desired_factor>floor else np.inf
     if required<=1:k=1.
     precision.append({'accession':acc,'gene':gene,'margin':margin,'in_family34':acc in FAMILY,'threshold':label,'alpha':alpha,
      'point_structure_compatible':compatible,'point_inside_both_margins':inside,'required_full_independent_information_multiplier':required,
      'fixed_cluster_variance_fraction_scenario_not_estimated':floor,'required_subject_multiplier_in_floor_scenario':k,
      'n_times_multiplier_if_new_independent_designs':int(np.ceil(d['n']*required)) if np.isfinite(required) and required<1e9 else np.nan,
      'interpretation':'CURRENT_POINT_HELD_FIXED_PRECISION_SCENARIO_NOT_POWER_OR_SAMPLE_SIZE_PROMISE',
      'limitation':'STRUCTURE_INCOMPATIBLE' if not compatible else 'POINT_OUTSIDE_MARGIN' if not inside else 'FIXED_CLUSTER_FLOOR_UNATTAINABLE' if not np.isfinite(k) else 'CONDITIONAL_PRECISION_EXTRAPOLATION'})
 csvout(RESULTS/'CONSENSUS_INFORMATION_BOTTLENECKS.csv',rows);csvout(RESULTS/'FIELLER_PRECISION_AND_DENOMINATOR_INFORMATION.csv',intervals)
 csvout(RESULTS/'CONDITIONAL_PRECISION_EXPANSION_SCENARIOS.csv',precision);csvout(RESULTS/'STAGE6_EQUIVALENCE_REPRODUCTION.csv',reproduction)
 assert all(r['passed'] for r in reproduction)
 for name,obj in [('POWER_TEMPLATES_RESTRICTED.pkl',templates),('REAL_MOMENTS_RESTRICTED.pkl',moments)]:
  with (PRIVATE/name).open('wb') as f:pickle.dump(obj,f,pickle.HIGHEST_PROTOCOL)
 table=pd.DataFrame(rows);family=table.loc[table.in_family34]
 csvout(RESULTS/'INFORMATION_BOTTLENECK_COUNTS.csv',family.groupby(['margin','main_bottleneck']).size().rename('n_proteins').reset_index())
 csvout(RESULTS/'EVIDENCE_STATE_COUNTS.csv',family.groupby(['margin','evidence_state']).size().rename('n_proteins').reset_index())
 log(f'All105 Stage6 equivalence results reproduced; {len(templates)} admissible real-covariance power templates; {M} family unchanged')
 print(pd.DataFrame(rows).query('in_family34').groupby(['margin','main_bottleneck']).size().to_string())

if __name__=='__main__':run()
