from collections import Counter
import re
from audit_common import *

def shape(s):
 return re.sub(r'[0-9]','d',text(s))

counts=Counter()
for _,r in rows('I13',['PTDOB','PTDOBYY','PTGENDER']):
 d=parse_date(r.get('PTDOB'));y=numeric(r.get('PTDOBYY'))
 counts['PTDOB_shape_'+shape(r.get('PTDOB'))]+=1
 counts['PTDOBYY_kind_'+y.kind]+=1
 counts['PTDOBYY_shape_'+shape(r.get('PTDOBYY'))]+=1
 if d:counts['parse_full_date']+=1
 if d and y.exact:counts['year_consistent_'+str(d.year==y.low)]+=1
write_csv(RESULTS/'DEMOGRAPHIC_FORMAT_AUDIT.csv',[{'pattern':k,'n_records':v} for k,v in counts.items()])
print(dict(counts))
