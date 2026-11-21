from pathlib import Path
import hashlib,csv,json,sys
from datetime import datetime
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4194304),b''):h.update(b)
 return h.hexdigest()
def run():
 manifest=ROOT/'SHA256_MANIFEST.csv';receipt=ROOT/'SEAL_RECEIPT.json'
 if manifest.exists():raise SystemExit('Already sealed. Use verify_delivery.py.')
 for n in ['SCIENTIFIC_QA_SUMMARY.json','DELIVERY_QA_SUMMARY.json']:
  d=json.loads((ROOT/'results'/n).read_text(encoding='utf-8'));assert d['failed']==0
  if n.startswith('DELIVERY'):assert d['visual_review_completed'] and not d['allow_unreviewed_copy']
 files=sorted(p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p not in [manifest,receipt]);rows=[{'relative_path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files]
 with manifest.open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['relative_path','bytes','sha256']);w.writeheader();w.writerows(rows)
 out={'sealed_at':datetime.now().astimezone().isoformat(),'files':len(rows),'bytes':sum(r['bytes'] for r in rows),'manifest_sha256':sha(manifest),'excluded':['SHA256_MANIFEST.csv','SEAL_RECEIPT.json','__pycache__'],'default_entry':'python -B code/verify_delivery.py','scope':'All deliverables including restricted matrices and saved simulation P values; prior stages and raw inputs checked separately.'};receipt.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out,ensure_ascii=True))
if __name__=='__main__':run()
