"""Read-only, standard-library-only verifier for the sealed Stage7 directory."""
from pathlib import Path
import csv, hashlib, json, sys

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/'SHA256_MANIFEST.csv'
def excluded(path):
 return path==MANIFEST or '__pycache__' in path.parts or path.suffix=='.pyc'
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as handle:
  for chunk in iter(lambda:handle.read(4*1024*1024),b''):h.update(chunk)
 return h.hexdigest()
def verify():
 if not MANIFEST.is_file():print('No seal manifest found.');return 1
 records=list(csv.DictReader(MANIFEST.open(encoding='utf-8-sig')))
 listed={r['relative_path'] for r in records};errors=[]
 if len(listed)!=len(records):errors.append('duplicate manifest paths')
 actual={str(p.relative_to(ROOT)).replace('\\','/') for p in ROOT.rglob('*') if p.is_file() and not excluded(p)}
 for name in sorted(actual-listed):errors.append('unlisted: '+name)
 for r in records:
  p=(ROOT/r['relative_path']).resolve()
  if not p.is_relative_to(ROOT.resolve()):errors.append('path outside root: '+r['relative_path']);continue
  if not p.is_file():errors.append('missing: '+r['relative_path']);continue
  if p.stat().st_size!=int(r['bytes']) or digest(p)!=r['sha256']:errors.append('mismatch: '+r['relative_path'])
 print(json.dumps({'sealed_entries':len(records),'verified':len(records)-sum(e.startswith(('missing:','mismatch:')) for e in errors),'errors':len(errors),'read_only':True},ensure_ascii=False))
 for e in errors[:30]:print(e)
 return int(bool(errors))
if __name__=='__main__':sys.exit(verify())
