from config import *
from designs import make
from scipy.linalg import lstsq
import pickle
frame=pd.read_csv(RESTRICTED/'LEGACY_L04_CLINICAL_DESIGN_RESTRICTED.csv');cache=np.load(RESTRICTED/'LEGACY_PEPTIDE_CACHE_RESTRICTED.npz');old=np.load(S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz');meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
masks=pickle.load((RESTRICTED/'LEGACY_CASE_MASKS_RESTRICTED.pkl').open('rb'));r=pd.read_csv(RESULTS/'LEGACY_SAME_DATE_RELATIVE_DESCRIPTIVE.csv').iloc[0]
mask=masks[r.design_id];d=make(frame,mask,tech=('NULISA_PlateID','SOMASCAN_PlateId','MRM_PLATE','MRM_RUN'),additional_factors=('ptau_batch','MRM_INSTRUMENT'))
idx=frame.loc[mask,'S5_index'].to_numpy(int);j=list(old['accessions']).index(r.accession);sj=list(old['soma_features']).index(meta.loc[r.accession,'soma_features']);k=list(cache['peptides']).index(r.peptide);Y=np.c_[old['N'][idx,j],old['S'][idx,sj],cache['Y'][mask,k]]
X=np.c_[d['Z'],d['P']]
for cond in [None,max(X.shape)*np.finfo(float).eps,1e-10]:
 coef,_,rank,s=lstsq(X,Y,cond=cond);E=Y-X@coef;print('SVD cond/rank/last singulars:',cond,rank,s[-8:]);print('R:',np.corrcoef(E,rowvar=False))
print('QR rank:',d['rank']);E=Y-d['Q']@(d['Q'].T@Y);print('QR R:',np.corrcoef(E,rowvar=False));print('produced R:',r[['corr_residual_N_S','corr_residual_N_targetMS','corr_residual_S_targetMS']].to_dict())
