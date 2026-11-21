from config import *
from model import bh

models=pd.read_csv(RESULTS/'COVARIATE_FOUR_MODEL_RESULTS.csv');sap=json.loads((S5/'results/MAIN_ANALYSIS_SAP_FROZEN.json').read_text(encoding='utf-8'));family=sap['primary_family_accessions']
base=models.loc[models.model.eq('M0')].set_index('accession').reindex(family);p0=base.p_value.to_numpy();q0=bh(p0,True);rows=[]
for name in ['EDU','APOE','BOTH']:
 frame=models.loc[models.model.eq(name)].set_index('accession').reindex(family);p1=frame.p_value.to_numpy();q1=bh(p1,True)
 for i,acc in enumerate(family):
  own=p0.copy();own[i]=p1[i];other=p1.copy();other[i]=p0[i]
  q10=bh(own,True)[i];q01=bh(other,True)[i];ownpart=.5*((q10-q0[i])+(q1[i]-q01));otherpart=.5*((q01-q0[i])+(q1[i]-q10))
  rows.append({'accession':acc,'gene':base.loc[acc,'gene'],'model':name,'p_M0':p0[i],'p_new':p1[i],'q_M0':q0[i],'q_new':q1[i],
   'q_change_own_p_contribution':ownpart,'q_change_other_proteins_p_contribution':otherpart,
   'identity_abs_error':abs(q1[i]-q0[i]-ownpart-otherpart),'same_fixed_family_m':len(family),'not_causal_decomposition':True})
csvout(RESULTS/'BY_Q_OWN_VS_FAMILY_DECOMPOSITION.csv',rows)
log('BY q change separated into own-p and other-family-p contributions; fixed m=34')
