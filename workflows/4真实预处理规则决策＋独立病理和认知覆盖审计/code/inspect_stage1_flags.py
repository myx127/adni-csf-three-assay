from audit_common import *
import pyarrow.parquet as pq
p=STAGE1/'00_restricted/TRIPLE_CANDIDATE_AVAILABILITY_RESTRICTED.parquet'
print('MASK_SCHEMA',pq.read_schema(p))
for name in ['PROTEIN_EVALUABILITY_70.csv','STRICT_CORE_MAPPING_DRAFT.csv','EXTENDED_CORE_MAPPING_DRAFT.csv']:
 d=pd.read_csv(STAGE1/'results'/name)
 print(name,'shape',d.shape,'columns',list(d.columns))
 print(d.iloc[:2].to_dict('records'))
