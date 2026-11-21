from config import *
import re

def run():
 f=pd.read_csv(PRIVATE/'MRM2014_CSV_EXPORT_RESTRICTED.csv');x=f.iloc[:,0].dropna().astype(str)
 # Only format patterns with all digits removed are printed, not individual IDs.
 print('CSV export subject-key formats:',x.map(lambda s:re.sub(r'\d+','{N}',s)).value_counts().head(10).to_dict())
 m=pd.read_csv(PRIVATE/'MRM2014_SAMPLEINFO_RESTRICTED.csv');print('Study categories:',m.Study.value_counts(dropna=False).to_dict());print('Visit categories:',m.Viscode.value_counts(dropna=False).to_dict())
 legacy=pd.read_csv(S6/'results/LEGACY_INPUT_MANIFEST.csv').set_index('input_id');l01=pd.read_csv(legacy.loc['L01','path']);print('L01 visit categories:',l01.VISCODE.value_counts(dropna=False).to_dict())
 print('Stage5 design metadata fields:',[c for c in pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv',nrows=0).columns if re.search('date|visit|viscode|RID|key|plate|plex',c,re.I)])
 # Header/schema only; collection dates are not replaced by scheduled visit dates.
 for p in [S1/'00_restricted',S4/'00_restricted']:
  print('Potential connecting metadata:',[x.name for x in p.glob('*.csv') if re.search('JOIN|REGIST|LINK|INDEX|SUBJECT|CLIN|BIO|SAMPLE',x.name,re.I)])

if __name__=='__main__':run()
