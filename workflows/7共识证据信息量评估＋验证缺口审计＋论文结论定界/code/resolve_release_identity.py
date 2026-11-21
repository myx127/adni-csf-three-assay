from config import *
from profile_release_precision import significant

def run():
 source=pd.read_csv(S6/'results/LEGACY_INPUT_MANIFEST.csv').set_index('input_id');a=pd.read_csv(source.loc['L01','path']);b=pd.read_csv(PRIVATE/'MRM2014_CSV_EXPORT_RESTRICTED.csv');link=pd.read_csv(PRIVATE/'MRM2014_EXPORT_SAMPLEINFO_LINK_RESTRICTED.csv');rows=[]
 cols=[c for c in b if c.replace('-','_',1) in a];A=a[[c.replace('-','_',1) for c in cols]].apply(pd.to_numeric,errors='coerce').to_numpy(float);B=b[cols].apply(pd.to_numeric,errors='coerce').to_numpy(float)
 for i,r in a.iterrows():
  cand=link.loc[link.RID.eq(r.RID)&link.VISCODE.eq(str(r.VISCODE)),'export_row'].to_numpy(int);out=[]
  for j in cand:
   valid=np.isfinite(A[i])&np.isfinite(B[j]);x=B[j,valid];y=A[i,valid]
   if len(x)<50:continue
   slope,intercept=np.linalg.lstsq(np.c_[x,np.ones(len(x))],y,rcond=None)[0];err=y-(slope*x+intercept);corr=np.corrcoef(x,y)[0,1]
   steps=10.**(np.floor(np.log10(np.maximum(abs(y),1e-300)))-5)
   out.append({'source_row':i+2,'RID':r.RID,'export_row':j,'n_common_features':len(x),'corr_all_features':corr,'slope':slope,'intercept':intercept,
    'affine_RMSE':np.sqrt(np.mean(err*err)),'affine_max_abs_error':np.max(abs(err)),
    'same_missing_pattern':np.array_equal(np.isfinite(A[i]),np.isfinite(B[j])),'log2_from_ln_max_abs_error':np.max(abs(y-x/np.log(2))),
    'ln_from_log2_max_abs_error':np.max(abs(y-x*np.log(2))),
    'max_raw_error_in_6digit_resolution_units':np.max(abs(y-x)/steps),
    'max_rounded_error_in_6digit_resolution_units':np.max(abs(y-significant(x,6))/steps),
    'n_features_differ_after_6digit_rounding':int((abs(y-significant(x,6))>1e-10).sum())})
  rows.extend(out)
 f=pd.DataFrame(rows);csvout(PRIVATE/'MRM2014_FULL_VECTOR_IDENTITY_DIAGNOSTICS_RESTRICTED.csv',f)
 if len(f):
  best=f.sort_values(['source_row','affine_RMSE']).drop_duplicates('source_row');summary=best[['n_common_features','corr_all_features','slope','intercept','affine_RMSE','affine_max_abs_error','max_raw_error_in_6digit_resolution_units','max_rounded_error_in_6digit_resolution_units','n_features_differ_after_6digit_rounding']].describe(percentiles=[.05,.5,.95]).reset_index().rename(columns={'index':'statistic'})
  csvout(RESULTS/'MRM2014_RELEASE_SCALE_IDENTITY_DIAGNOSTICS.csv',summary);print(summary.to_string(index=False))
  print('Rows with any identity candidates:',best.source_row.nunique(),'of',len(a),'maximum candidates per source row:',int(f.groupby('source_row').size().max()))
 else:print('No common identity candidates')
 log('2014 full-vector affine/log-base identity audit complete without any pathology association fit')

if __name__=='__main__':run()
