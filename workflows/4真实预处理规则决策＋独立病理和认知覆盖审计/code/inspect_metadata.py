import csv,json,collections
from pathlib import Path
from inspect_inputs import ROOT,RESULTS,MANIFEST

def main():
 inputs={r['input_id']:r for r in csv.DictReader(MANIFEST.open(encoding='utf-8-sig'))}
 extras={r['input_id']:r for r in csv.DictReader((RESULTS/'ADDITIONAL_INPUT_MANIFEST.csv').open(encoding='utf-8-sig'))}
 plans={'I10':['BATCH','KIT','STDS'],'E01':['BATCH'],'I11':['qc_flag','TRACER'],'I12':['qc_flag','TRACER'],
         'I13':['PHASE','HAS_QC_ERROR'],'I14':['APUSABLE'],'I15':['PHASE','HAS_QC_ERROR']}
 rows=[]
 for key,fields in plans.items():
  entry=inputs.get(key,extras.get(key));counts={f:collections.Counter() for f in fields};total=0
  with Path(entry['path']).open(encoding='utf-8-sig',newline='') as f:
   reader=csv.DictReader(f)
   for row in reader:
    total+=1
    for field in fields:counts[field][row.get(field,'').strip()]+=1
  for field,count in counts.items():
   for value,n in count.items():rows.append({'input_id':key,'field':field,'metadata_value':value,'n_records':n})
  print(key,'rows',total,{f:dict(c) for f,c in counts.items()})
 with (RESULTS/'METADATA_CATEGORY_COUNTS.csv').open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

if __name__=='__main__':main()
