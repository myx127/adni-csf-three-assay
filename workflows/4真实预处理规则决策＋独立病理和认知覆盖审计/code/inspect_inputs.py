from pathlib import Path
import csv,hashlib,json

ROOT=Path(__file__).resolve().parents[1]
PROJECT=ROOT.parent
RESULTS=ROOT/'results';RESULTS.mkdir(exist_ok=True)
MANIFEST=PROJECT/'规划和详细数据情况以及分析方案和研究计划/01_数据清单/CORE_INPUT_MANIFEST.csv'

def main():
 rows=list(csv.DictReader(MANIFEST.open(encoding='utf-8-sig')));checks=[];schemas=[]
 for row in rows:
  path=Path(row['path']);h=hashlib.sha256()
  with path.open('rb') as f:
   for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
  checks.append({'input_id':row['input_id'],'path':str(path),'bytes':path.stat().st_size,'sha256':h.hexdigest(),
    'expected_sha256':row['sha256'],'passed':path.stat().st_size==int(row['bytes']) and h.hexdigest()==row['sha256']})
 (RESULTS/'INPUT_HASH_CHECKS.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf-8')
 assert all(x['passed'] for x in checks)
 print('15/15 input identities match')
 for row in rows:
  if int(row['input_id'][1:])>=9:
   with Path(row['path']).open(encoding='utf-8-sig',newline='') as f:header=next(csv.reader(f))
   schemas.append({'input_id':row['input_id'],'columns':header})
   print(row['input_id'],'header_columns',len(header))
 (RESULTS/'CLINICAL_FILE_SCHEMAS.json').write_text(json.dumps(schemas,ensure_ascii=False,indent=2),encoding='utf-8')
 boundary=RESULTS/'ACCESS_BOUNDARY_INITIAL.json'
 protocol_sha=hashlib.sha256((ROOT/'AUDIT_PROTOCOL_v1.0.md').read_bytes()).hexdigest()
 if boundary.exists():
  assert json.loads(boundary.read_text(encoding='utf-8'))['protocol_sha256_before_clinical_content']==protocol_sha
  return
 boundary.write_text(json.dumps({'protocol_sha256_before_clinical_content':protocol_sha,
   'clinical_effect_estimation':False,'protein_quantitative_files_opened_for_new_analysis':False,'first_pass':'hashes and headers only'},indent=2),encoding='utf-8')

if __name__=='__main__':main()
