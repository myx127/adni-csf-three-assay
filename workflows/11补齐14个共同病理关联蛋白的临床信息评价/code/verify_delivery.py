from pathlib import Path
import csv,json,hashlib,sys
ROOT=Path(__file__).resolve().parents[1]
if hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8',errors='backslashreplace')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def verify():
 mf=ROOT/'SHA256_MANIFEST.csv';errors=[]
 if not mf.exists():print('Delivery is not sealed.');return 1
 rows=list(csv.DictReader(mf.open(encoding='utf-8-sig')))
 for i,r in enumerate(rows,1):
  p=(ROOT/r['relative_path']).resolve()
  if not p.is_relative_to(ROOT.resolve()) or not p.is_file():errors.append('Missing/escaping '+r['relative_path']);continue
  if p.stat().st_size!=int(r['bytes']) or sha(p)!=r['sha256']:errors.append('Changed '+r['relative_path'])
  if i%400==0:print(f'Verified {i}/{len(rows)} files',flush=True)
 s=json.loads((ROOT/'STAGE11_COMPLETION_STATUS.json').read_text(encoding='utf-8'))
 if s.get('status')!='COMPLETE' or s.get('fixed_family')!=34 or s.get('fixed_focus_subset')!=14:errors.append('Invalid completion record')
 if s['scientific_QA']['failed'] or s['delivery_QA']['failed']:errors.append('QA not passed')
 print(json.dumps({'status':s.get('status'),'files_checked':len(rows),'errors':len(errors),'primary_focus14':s['results']['focus14_primary_supported'],'secondary_focus14':s['results']['focus14_secondary_supported']},ensure_ascii=False))
 for e in errors[:20]:print(e)
 return int(bool(errors))
if __name__=='__main__':sys.exit(verify())

