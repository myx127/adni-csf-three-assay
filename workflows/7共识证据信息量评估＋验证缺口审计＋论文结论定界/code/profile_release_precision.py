from config import *

def significant(x,p):
 x=np.asarray(x,float);scale=np.ones_like(x);valid=np.isfinite(x)&(x!=0);scale[valid]=10.**(p-1-np.floor(np.log10(abs(x[valid]))));return np.round(x*scale)/scale

def run():
 source=pd.read_csv(S6/'results/LEGACY_INPUT_MANIFEST.csv').set_index('input_id');a=pd.read_csv(source.loc['L01','path']);b=pd.read_csv(PRIVATE/'MRM2014_CSV_EXPORT_RESTRICTED.csv');cols=[c for c in b if c.replace('-','_',1) in a]
 rows=[]
 for name,frame in [('L01',a[[c.replace('-','_',1) for c in cols]]),('L02_EXPORT',b[cols])]:
  x=frame.apply(pd.to_numeric,errors='coerce').to_numpy(float);x=x[np.isfinite(x)]
  for p in range(4,11):
   for kind,rounded in [('DECIMAL_PLACES',np.round(x,p)),('SIGNIFICANT_DIGITS',significant(x,p))]:
    error=abs(x-rounded);rows.append({'source':name,'grid_kind':kind,'digits':p,'max_grid_error':error.max(),'fraction_within_1e_minus10':np.mean(error<1e-10)})
 out=pd.DataFrame(rows);csvout(RESULTS/'MRM2014_RELEASE_PRECISION_PROFILE.csv',out);print(out.to_string(index=False))

if __name__=='__main__':run()
