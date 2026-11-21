"""Read-only standard-library verifier for sealed Stage8."""
from pathlib import Path
import csv,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1];MANIFEST=ROOT/'SHA256_MANIFEST.csv'
def excluded(path):return path==MANIFEST or '__pycache__' in path.parts or path.suffix=='.pyc'
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def verify():
 if not MANIFEST.is_file():print('No seal manifest found');return 1
 rows=list(csv.DictReader(MANIFEST.open(encoding='utf-8-sig')));listed={r['relative_path'] for r in rows};errors=[]
 if len(rows)!=len(listed):errors.append('duplicate paths')
 actual={str(p.relative_to(ROOT)).replace('\\','/') for p in ROOT.rglob('*') if p.is_file() and not excluded(p)}
 errors.extend('unlisted: '+n for n in sorted(actual-listed))
 for r in rows:
  p=(ROOT/r['relative_path']).resolve()
  if not p.is_relative_to(ROOT.resolve()):errors.append('outside root: '+r['relative_path']);continue
  if not p.is_file():errors.append('missing: '+r['relative_path']);continue
  if p.stat().st_size!=int(r['bytes']) or digest(p)!=r['sha256']:errors.append('mismatch: '+r['relative_path'])
 print(json.dumps({'entries':len(rows),'errors':len(errors),'read_only':True},ensure_ascii=False))
 for e in errors[:20]:print(e)
 return int(bool(errors))
if __name__=='__main__':sys.exit(verify())
