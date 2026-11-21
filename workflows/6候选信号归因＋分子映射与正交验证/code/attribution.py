from config import *
from designs import make
from diagnostics import parameters,constraints_jacobian,summarize
from model import bh,holm
import pickle

def y_for(cache,meta,acc):
 j=list(cache['accessions']).index(acc);k=list(cache['soma_features']).index(meta.loc[acc,'soma_features'])
 return np.c_[cache['N'][:,j],cache['S'][:,k],cache['T_GIS'][:,j]]
def H(b,c):return np.array([b[1]*c[1]-b[0]*c[2],b[2]*c[0]-b[0]*c[2]])
def q_of(h,V):return float(h@np.linalg.solve(V,h))

def run():
 clinical=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv');cache=np.load(S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz')
 meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession');old=pickle.load((S5/'00_restricted/DESIGNS_MAIN_RESTRICTED.pkl').open('rb'))
 sap=json.loads((S5/'results/MAIN_ANALYSIS_SAP_FROZEN.json').read_text(encoding='utf-8'));family=sap['primary_family_accessions']
 rows=[];models={};state={};calibration_designs={}
 for acc,d0 in old.items():
  Yall=y_for(cache,meta,acc);mask=np.zeros(len(clinical),bool);mask[d0['indices']]=True;state[acc]={}
  assert np.isfinite(clinical.loc[mask,['education','apoe4_carrier']]).all().all()
  for name,extra in [('M0',()),('EDU',('education',)),('APOE',('apoe4_carrier',)),('BOTH',('education','apoe4_carrier'))]:
   d=make(clinical,mask,extra=extra);assert np.array_equal(d['indices'],d0['indices']);Y=Yall[d['indices']]
   stats,theta,IF,V,E=summarize(Y,d,0);h,J=constraints_jacobian(theta);scale=d['df']/d['n'];h_fixed=h*scale;Vh_fixed=(J@V@J.T)*scale**2
   state[acc][name]={'b':theta[:3],'c':theta[6:9]*scale,'h':h_fixed,'Vh':Vh_fixed,'Q':stats['Q']}
   rows.append({'accession':acc,'gene':meta.loc[acc,'nulisa_target'],'model':name,'n':d['n'],'rank_X':d['rank'],'residual_df':d['df'],
     'P_residual_SS':float(1/(d['a']@d['a'])),'H1_on_common_1overN_scale':h_fixed[0],'H2_on_common_1overN_scale':h_fixed[1],
     'H1_SE_common_scale':float(np.sqrt(Vh_fixed[0,0])),'H2_SE_common_scale':float(np.sqrt(Vh_fixed[1,1])),**stats})
   models[(acc,name)]=d
   if name in ['EDU','APOE']:calibration_designs[name+'__'+acc]=d
 out=pd.DataFrame(rows);out['q_BY_exploratory']=np.nan
 for name,part in out.groupby('model'):
  ids=part.index[part.accession.isin(family)];out.loc[ids,'q_BY_exploratory']=bh(out.loc[ids,'p_value'].to_numpy(),True)
 csvout(RESULTS/'COVARIATE_FOUR_MODEL_RESULTS.csv',out)
 decomposition=[];cov_shap=[]
 for acc,states in state.items():
  zero=states['M0']
  for name in ['EDU','APOE','BOTH']:
   one=states[name];q00=q_of(zero['h'],zero['Vh']);q10=q_of(one['h'],zero['Vh']);q01=q_of(zero['h'],one['Vh']);q11=q_of(one['h'],one['Vh'])
   h_part=.5*((q10-q00)+(q11-q01));v_part=.5*((q01-q00)+(q11-q10))
   b_part=.5*((H(one['b'],zero['c'])-H(zero['b'],zero['c']))+(H(one['b'],one['c'])-H(zero['b'],one['c'])))
   c_part=.5*((H(zero['b'],one['c'])-H(zero['b'],zero['c']))+(H(one['b'],one['c'])-H(one['b'],zero['c'])))
   decomposition.append({'accession':acc,'gene':meta.loc[acc,'nulisa_target'],'comparison':'M0_to_'+name,'Q_M0':q00,'Q_new':q11,'delta_Q':q11-q00,
     'delta_Q_H_point_contribution':h_part,'delta_Q_variance_contribution':v_part,'Q_identity_abs_error':abs(q11-q00-h_part-v_part),
     'delta_H1_b_contribution':b_part[0],'delta_H1_c_contribution':c_part[0],'delta_H2_b_contribution':b_part[1],'delta_H2_c_contribution':c_part[1],
     'H_identity_max_abs_error':float(abs(one['h']-zero['h']-b_part-c_part).max()),'not_causal_mediation':True})
  edu=.5*((states['EDU']['Q']-zero['Q'])+(states['BOTH']['Q']-states['APOE']['Q']))
  apoe=.5*((states['APOE']['Q']-zero['Q'])+(states['BOTH']['Q']-states['EDU']['Q']))
  cov_shap.append({'accession':acc,'gene':meta.loc[acc,'nulisa_target'],'delta_Q_both':states['BOTH']['Q']-zero['Q'],
   'education_Q_contribution':edu,'APOE4_Q_contribution':apoe,'identity_abs_error':abs(states['BOTH']['Q']-zero['Q']-edu-apoe),
   'descriptive_model_contribution_not_biological_percent':True})
 csvout(RESULTS/'H_V_AND_B_C_DECOMPOSITION.csv',decomposition);csvout(RESULTS/'EDUCATION_APOE_SHAPLEY_Q.csv',cov_shap)
 with (RESTRICTED/'ATTRIBUTION_MODELS_RESTRICTED.pkl').open('wb') as f:pickle.dump(models,f,pickle.HIGHEST_PROTOCOL)
 with (RESTRICTED/'ATTRIBUTION_CALIBRATION_DESIGNS_RESTRICTED.pkl').open('wb') as f:pickle.dump(calibration_designs,f,pickle.HIGHEST_PROTOCOL)
 reference=[]
 main=pd.read_csv(S5/'results/PRIMARY_ALL43_RESULTS.csv').set_index('accession');variants=pd.read_csv(S5/'results/ALL_FIXED_VARIANT_RESULTS.csv')
 for r in out.itertuples():
  if r.model in ['M0','BOTH']:
   ref=main.loc[r.accession] if r.model=='M0' else variants.loc[variants.accession.eq(r.accession)&variants.variant.eq('FULLCOV')].iloc[0]
   reference.append({'accession':r.accession,'model':r.model,'Q_abs_difference':abs(r.Q-ref.Q),'p_abs_difference':abs(r.p_value-ref.p_value),
    'passed':bool(np.allclose(r.Q,ref.Q,rtol=1e-8) and np.allclose(r.p_value,ref.p_value,rtol=1e-8))})
 csvout(RESULTS/'S5_M0_BOTH_REPRODUCTION_QA.csv',reference);assert all(r['passed'] for r in reference)
 save_json(RESULTS/'ATTRIBUTION_DESIGN_FREEZE.json',{'before_new_EDU_APOE_simulation':True,'already_known_S5_results':True,'new_designs':70,
  'variance_version_fixed':0,'reps_per_cell':4000,'scenarios':['C01','C02','C03','C04','C05','C06'],'root_seed':SEEDS['SECONDARY'],
  'code_sha256':{p.name:sha(p) for p in [Path(__file__),ROOT/'code/designs.py',ROOT/'code/model.py',ROOT/'code/diagnostics.py',ROOT/'code/simulation_reference.py']}})
 log('Four-model attribution fitted on identical cases; 70 M0/BOTH comparisons reproduce Stage5; no primary discovery redefinition')
 print(out.loc[out.accession.isin(CANDIDATES),['gene','model','Q','p_value','q_BY_exploratory']].to_string(index=False))

if __name__=='__main__':run()
