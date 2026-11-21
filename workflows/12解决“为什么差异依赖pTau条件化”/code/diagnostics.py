from config import *
from engine import basis,nuisance,rcs,D,stats
from scipy.linalg import lstsq
import pickle
def solve(X,y):return lstsq(X,y,cond=max(X.shape)*np.finfo(float).eps,lapack_driver='gelsy')[0]
def run():
 models=pickle.load((PRIVATE/'MODELS_RESTRICTED.pkl').open('rb'));c=pd.read_csv(PRIVATE/'COMMON_CLINICAL_RESTRICTED.csv');knots=np.array(json.loads((RESULTS/'RCS_KNOTS_FROZEN.json').read_text(encoding='utf-8'))['P_LOG_knots'])
 effects=pd.read_csv(RESULTS/'PLATFORM_EFFECTS_ALL27.csv');decomp=[];components=[];pair=[];moments=[];info=[];curves=[];factor=[];loco=[]
 for acc,gene in GENES.items():
  for new,ref in [('LOG_P','NO_P'),('RAW_P','NO_P'),('RCS_P','NO_P'),('BG2_LOG_P','BG2_NO_P'),('BG2_RCS_P','BG2_NO_P'),('BG3_LOG_P','BG3_NO_P')]:
   if gene+'__'+new not in models:continue
   t=models[gene+'__'+new];r=models[gene+'__'+ref];F=t['pfun']*np.sqrt(t['w'])[:,None];ovb=[]
   for j,a in enumerate(['N','S','T']):
    X0=np.c_[r['ZZ_list'][j],r['targetx_list'][j]];X1=np.c_[X0,F];b0=solve(X0,t['y']);b1=solve(X1,t['y']);pi=solve(X0,F)[-1];delta=b1[-F.shape[1]:];before=b0[-1]*t['sigma'][j];after=b1[-F.shape[1]-1]*t['sigma'][j];contribution=pi*delta*t['sigma'][j];ovb.append(contribution.sum())
    decomp.append({'gene':gene,'reference_model':ref,'new_model':new,'assay':a,'gamma_before':before,'gamma_after':after,'after_minus_before':after-before,'OVB_before_minus_after':contribution.sum(),'identity_abs_error':abs(before-after-contribution.sum()),'not_causal_mediation':True})
    for k,name in enumerate(t['pnames']):components.append({'gene':gene,'new_model':new,'assay':a,'pTau_basis':name,'delta_outcome_coefficient':delta[k],'pi_pTau_per_protein_SD':pi[k]*t['sigma'][j],'product_before_minus_after':contribution[k]})
   predicted=-D@np.array(ovb);actual=stats(t['y'],t)['h'][0]-stats(r['y'],r)['h'][0]
   for j,label in enumerate(['S_MINUS_N','T_MINUS_N']):pair.append({'gene':gene,'reference_model':ref,'new_model':new,'contrast':label,'observed_change':actual[j],'OVB_reconstructed_change':predicted[j],'identity_abs_error':abs(actual[j]-predicted[j]),'not_a_mediation_fraction':True})
  nop=models[gene+'__NO_P'];logt=models[gene+'__LOG_P'];Q=nop['Qcommon'];xx=np.column_stack(nop['targetx_list']);E=xx-Q@(Q.T@xx);P=c.iloc[nop['ids']].P_LOG.to_numpy()*np.sqrt(nop['w']);pp=P-Q@(Q.T@P);yy=nop['y']-Q@(Q.T@nop['y']);sPP=pp@pp;sPY=pp@yy;sYY=yy@yy
  for j,a in enumerate(['N','S','T']):
   v=E[:,j];sXX=v@v;sXP=v@pp;sXY=v@yy
   b0=sXY/sXX;b1=(sXY*sPP-sXP*sPY)/(sXX*sPP-sXP*sXP);retain=1-sXP*sXP/(sXX*sPP)
   moments.append({'gene':gene,'assay':a,'weighted_partial_corr_X_pTau':sXP/np.sqrt(sXX*sPP),'weighted_partial_corr_X_Y_without_pTau':sXY/np.sqrt(sXX*sYY),'weighted_partial_corr_pTau_Y_without_X':sPY/np.sqrt(sPP*sYY),'exposure_information_retained_after_logpTau':retain,'VIF_from_logpTau':1/retain,'gamma_without_pTau_FWL':b0*nop['sigma'][j],'gamma_with_pTau_FWL':b1*nop['sigma'][j],'joint_formula_matches_full_model':bool(np.isclose(b1,stats(logt['y'],logt)['b'][0,j],rtol=1e-8,atol=1e-10))})
  for m in MODELS:
   t=models.get(gene+'__'+m)
   if t is None:continue
   for j,a in enumerate(['N','S','T']):info.append({'gene':gene,'model':m,'assay':a,'target_residual_SS':t['den'][j],'fraction_of_original_NO_P_information':t['den'][j]/nop['den'][j],'n_people':t['n_people'],'max_hat_all_assays':t['max_hat']})
   # A point single-factor compatibility calculation; no inferential or biological reliability claim.
   E=np.column_stack(t['targetx_list'])-t['Qcommon']@(t['Qcommon'].T@np.column_stack(t['targetx_list']));y=t['y']-t['Qcommon']@(t['Qcommon'].T@t['y']);df=len(y)-t['Qcommon'].shape[1];C=E.T@E/df;cy=E.T@y/df;v=np.diag(C);cross=np.array([C[0,1],C[0,2],C[1,2]])
   compatible=bool((cross>0).all());lam=np.full(3,np.nan);theta=np.full(3,np.nan);pred=np.full(3,np.nan);frac=np.full(3,np.nan);why='NONPOSITIVE_CROSS_COVARIANCE'
   if compatible:
    lam=np.sqrt([cross[0]*cross[1]/cross[2],cross[0]*cross[2]/cross[1],cross[1]*cross[2]/cross[0]]);theta=v-lam*lam;compatible=bool((theta>=0).all());why='POINT_COMPATIBLE_NOT_VALIDATED' if compatible else 'HEYWOOD_NEGATIVE_WORKING_ERROR_VARIANCE'
    if compatible:
     k=np.sum(lam*cy/v)/np.sum(lam*lam/v);pred=t['sigma']*lam*k/v;frac=lam*lam/v
   for j,a in enumerate(['N','S','T']):
    actual=effects.loc[effects.gene.eq(gene)&effects.model.eq(m)&effects.assay.eq(a),'gamma'].iloc[0]
    factor.append({'gene':gene,'model':m,'assay':a,'point_positive_single_factor_compatible':compatible,'status':why,'residual_variance':v[j],'working_loading':lam[j],'working_error_variance':theta[j],'working_common_fraction_NOT_technical_reliability':frac[j],'actual_gamma':actual,'working_factor_predicted_gamma':pred[j],'actual_minus_working_prediction':actual-pred[j],'not_mechanism_or_reliability_established':True})
  # Prespecified weighted adjusted curves, expressed relative to median pTau.
  valid=c.base_primary.astype(bool);grid=np.linspace(c.loc[valid,'P_LOG'].quantile(.01),c.loc[valid,'P_LOG'].quantile(.99),101);pmed=float(c.loc[valid,'P_LOG'].median());trans=np.polyfit(np.log2(c.loc[valid,'ptau']),c.loc[valid,'P_LOG'],1);rawgrid=2**((grid-trans[1])/trans[0])
  for new,ref in [('LOG_P','NO_P'),('RCS_P','NO_P'),('BG2_LOG_P','BG2_NO_P'),('BG2_RCS_P','BG2_NO_P')]:
   t=models[gene+'__'+new];r=models[gene+'__'+ref];F=t['pfun']*np.sqrt(t['w'])[:,None];XX=np.c_[r['ZZ_list'][0],F];fit=solve(XX,np.column_stack([t['y']]+t['targetx_list']));coef=fit[-F.shape[1]:]
   profile=rcs(grid,knots)-rcs([pmed],knots) if 'RCS' in new else (grid-pmed)[:,None];values=profile@coef
   for i in range(len(grid)):
    curves.append({'gene':gene,'model':new,'P_LOG':grid[i],'pTau_pg_per_ml':rawgrid[i],'memory_change_component_relative_median':values[i,0],'N_protein_component_in_fixedSD':values[i,1]/t['sigma'][0],'S_protein_component_in_fixedSD':values[i,2]/t['sigma'][1],'T_protein_component_in_fixedSD':values[i,3]/t['sigma'][2],'time_information_weighted_adjusted_curve_not_causal':True})
  # Delete one technical cluster: descriptive point refits only.
  f=c.iloc[nop['ids']].copy();fac=np.sqrt(nop['w']);Y=nop['y'];X=nop['X'];sigma=nop['sigma']
  for axis,col in enumerate(['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex']):
   for label in sorted(f[col].astype(str).unique()):
    keep=f[col].astype(str).ne(label).to_numpy();ff=f.iloc[np.flatnonzero(keep)];out={'gene':gene,'technical_axis':col,'omitted_cluster':label,'n_retained':int(keep.sum()),'no_new_P_values':True}
    for m in ['NO_P','LOG_P']:
     Z,_,_,_=nuisance(ff,m,knots);Zw=Z*fac[keep,None];Qz=basis(Zw);xw=X[keep]*fac[keep,None];rx=xw-Qz@(Qz.T@xw);ry=Y[keep]-Qz@(Qz.T@Y[keep]);b=(rx.T@ry)/np.sum(rx*rx,axis=0);gamma=b*sigma;contrast=D@gamma
     for j,a in enumerate(['N','S','T']):out[m+'_gamma_'+a]=gamma[j]
     out[m+'_S_MINUS_N']=contrast[0];out[m+'_T_MINUS_N']=contrast[1]
    out['CHANGE_S_MINUS_N']=out['LOG_P_S_MINUS_N']-out['NO_P_S_MINUS_N'];out['CHANGE_T_MINUS_N']=out['LOG_P_T_MINUS_N']-out['NO_P_T_MINUS_N'];loco.append(out)
 csvout(RESULTS/'OVB_ASSAY_DECOMPOSITION.csv',decomp);csvout(RESULTS/'OVB_PTAU_BASIS_COMPONENTS.csv',components);csvout(RESULTS/'OVB_PAIRED_CHANGE_RECONSTRUCTION.csv',pair);csvout(RESULTS/'FWL_WEIGHTED_MOMENT_EXPLANATION.csv',moments);csvout(RESULTS/'EXPOSURE_INFORMATION_AND_LEVERAGE.csv',info);csvout(RESULTS/'ADJUSTED_PTAU_CURVES.csv',curves);csvout(RESULTS/'WORKING_SINGLE_FACTOR_POINT_AUDIT.csv',factor);csvout(RESULTS/'LEAVE_ONE_TECHNICAL_CLUSTER_POINTS.csv',loco)
 assert max(r['identity_abs_error'] for r in decomp)<1e-8 and max(r['identity_abs_error'] for r in pair)<1e-8 and all(r['joint_formula_matches_full_model'] for r in moments)
 corr=c[['P_LOG','BG_NULISA','BG_SOMA','BG_TMT']].corr(method='spearman');csvout(RESULTS/'BACKGROUND_PTAU_SPEARMAN_DESCRIPTIVE.csv',corr.rename_axis('variable').reset_index())
 log('Exact OVB/FWL, nonlinear curves, background information, working-factor audit and leave-cluster point diagnostics complete.')
if __name__=='__main__':assert_unsealed();run()

