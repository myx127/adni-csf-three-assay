"""Read-only verification; intentionally does not import config or write caches/logs."""
from pathlib import Path
import csv,hashlib,sys,json
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1];MANIFEST=ROOT/'STAGE6_INPUT_OUTPUT_SHA256.csv'

def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()

def excluded(path):return path==MANIFEST or '__pycache__' in path.parts or path.suffix in ['.pyc','.tmp','.part']

def run():
 rows=list(csv.DictReader(MANIFEST.open(encoding='utf-8-sig')));errors=[]
 for r in rows:
  p=ROOT/r['path']
  if not p.is_file() or p.stat().st_size!=int(r['bytes']) or digest(p)!=r['sha256']:errors.append(r['path'])
 actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and not excluded(p)};expected={r['path'] for r in rows}
 extra=sorted(actual-expected);missing=sorted(expected-actual)
 print(json.dumps({'manifest_entries':len(rows),'hash_or_size_failures':len(errors),'unregistered_files':len(extra),'missing_files':len(missing),'read_only':True},ensure_ascii=False))
 if errors or extra or missing:print({'failed':errors,'extra':extra,'missing':missing});raise SystemExit(1)

if __name__=='__main__':run()
