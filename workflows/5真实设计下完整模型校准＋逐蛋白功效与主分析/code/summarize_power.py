from config import *
from model import wilson

sap=json.loads((RESULTS/'MAIN_ANALYSIS_SAP_FROZEN.json').read_text(encoding='utf-8'));m=sap['m_primary'];version=sap['version']
data=pd.read_csv(RESULTS/'POWER_SIMULATION_SUMMARY.csv');data=data.loc[data.version.eq(version)&data.alpha.eq(.05)]
screen=pd.read_csv(RESULTS/'PRE_EFFECT_PROTEIN_SCREEN.csv').set_index('accession')
rows=[]
for r in data.itertuples():
 ddf=float(screen.loc[r.design_id,'residual_df']);n=int(screen.loc[r.design_id,'n']);mult=[1.,(n-1)/ddf,2*(n-1)/ddf][version]
 q=pd.read_parquet(SIM/r.file).Q0.to_numpy();p=np.exp(-q/(2*mult));valid=np.isfinite(p)
 for label,alpha in [('NOMINAL05',.05),('BONFERRONI',.05/m),('BY_FIRST_DISCOVERY',.05/(m*np.sum(1/np.arange(1,m+1))))]:
  k=int((valid&(p<=alpha)).sum());lo,hi=wilson(k,len(p))
  # Exact Wilson endpoints; avoid floating point 1-ulp exclusions of 0 or 1.
  if k==0:lo=0.
  if k==len(p):hi=1.
  rows.append({'accession':r.design_id,'gene':screen.loc[r.design_id,'gene'],'role':screen.loc[r.design_id,'role'],'n':n,'target':r.target,'delta_SD':r.delta,
   'threshold':label,'alpha':alpha,'n_replicates':len(p),'n_failed':int((~valid).sum()),'power':k/len(p),'wilson95_low':lo,'wilson95_high':hi,
   'noise_parameters':'sigma=.6/.8/1.2; common T=.55P+U; three-way slopes=.35; expected residual SD normalized1','observed_effect_not_used':True})
out=pd.DataFrame(rows);csvout(RESULTS/'PER_PROTEIN_CONDITIONAL_POWER.csv',out)
summary=[]
for (acc,target,threshold),part in out.groupby(['accession','target','threshold']):
 part=part.sort_values('delta_SD');met=part.loc[part.power.ge(.8)]
 first=float(met.delta_SD.iloc[0]) if len(met) else np.nan
 prev=part.loc[part.delta_SD.lt(first),'delta_SD'].max() if np.isfinite(first) else part.delta_SD.max()
 summary.append({'accession':acc,'gene':part.gene.iloc[0],'target':target,'threshold':threshold,'first_tested_delta_power_ge80':first,
  'previous_tested_delta':prev,'no_interpolated_MDE':True,'not_reached_on_grid':not len(met),'crosses80_within_MC_CI':bool(len(met) and met.wilson95_low.iloc[0]<.8)})
csvout(RESULTS/'PER_PROTEIN_DELTA80_GRID.csv',summary)
core=out.loc[out.role.eq('PRIMARY_DISCOVERY_CANDIDATE')]
agg=core.groupby(['target','delta_SD','threshold']).agg(n_proteins=('accession','size'),min_power=('power','min'),median_power=('power','median'),max_power=('power','max')).reset_index()
csvout(RESULTS/'POWER_SUMMARY_BY_EFFECT.csv',agg)
print(agg.query("delta_SD in [.35,.5,.75] and threshold in ['NOMINAL05','BY_FIRST_DISCOVERY']").to_string(index=False))
