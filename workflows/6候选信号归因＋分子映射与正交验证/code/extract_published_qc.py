from config import *
import re

f=pd.read_csv(RESULTS/'TARGETED_FEATURE_VALIDITY_COVERAGE.csv');wanted=f.loc[f.accession.isin(['P00441','P36222','P40925','P47972'])].drop_duplicates(['source','feature_key'])
allpeps=sorted(set(f.peptide.dropna()),key=len,reverse=True);rows=[]
for source,method,pools in [('V01','VM01',['ATneg','ATpos']),('V02','VM02',['QCpool','ATneg','ATpos'])]:
 text=re.sub(r'\s+',' ',(EVIDENCE/(method+'_METHODS_EXTRACTED.txt')).read_text(encoding='utf-8'))
 for r in wanted.loc[wanted.source.eq(source)].itertuples():
  found=[]
  for m in re.finditer(r'\b'+re.escape(r.peptide)+r'\b',text):
   tail=text[m.end():m.end()+500];nextpep=re.search(r'\b(?:'+'|'.join(map(re.escape,allpeps))+r')\b',tail)
   if nextpep:tail=tail[:nextpep.start()]
   pct=re.findall(r'(\d+(?:\.\d+)?)\s*%',tail)
   if len(pct)==len(pools):found.append((pct,tail[:300]))
  if found:
   values,context=found[0]
   for pool,v in zip(pools,values):rows.append({'source':source,'method':method,'accession':r.accession,'peptide':r.peptide,'published_QC_pool':pool,
    'published_CV_percent':float(v),'values_from_method_not_recomputed_from_raw':True,'method_context_after_peptide':context})
  else:rows.append({'source':source,'method':method,'accession':r.accession,'peptide':r.peptide,'published_QC_pool':'NOT_MACHINE_EXTRACTED',
    'published_CV_percent':np.nan,'values_from_method_not_recomputed_from_raw':True,'method_context_after_peptide':'See original method table; no inferred pool labels'})
csvout(RESULTS/'TARGETED_PUBLISHED_QC_EVIDENCE.csv',rows)
