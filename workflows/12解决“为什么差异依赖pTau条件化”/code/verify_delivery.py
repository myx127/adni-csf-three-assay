"""Read-only default verification; imports only the Python standard library."""
from pathlib import Path
import csv,json,hashlib,argparse,sys
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4194304),b''):h.update(b)
 return h.hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--current-only',action='store_true');args=ap.parse_args();errors=[];manifest=ROOT/'SHA256_MANIFEST.csv';receipt=ROOT/'SEAL_RECEIPT.json';cache={}
 def verify(p,expected,label):
  try:
   if str(p) not in cache:cache[str(p)]=sha(p)
   if cache[str(p)]!=expected:errors.append(label+': HASH_MISMATCH')
  except Exception as e:errors.append(label+': '+str(e))
 try:rows=list(csv.DictReader(manifest.open(encoding='utf-8-sig')));rec=json.loads(receipt.read_text(encoding='utf-8'))
 except Exception as e:print(json.dumps({'errors':1,'detail':str(e)},ensure_ascii=True));return 1
 verify(manifest,rec['manifest_sha256'],'MANIFEST');expected=set()
 for r in rows:
  p=ROOT/r['relative_path'];expected.add(p.resolve());verify(p,r['sha256'],r['relative_path'])
  if p.exists() and p.stat().st_size!=int(r['bytes']):errors.append(r['relative_path']+': SIZE_MISMATCH')
 actual={p.resolve() for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p not in [manifest,receipt]}
 for p in actual-expected:errors.append('EXTRA_FILE: '+str(p))
 if len(rows)!=rec['files'] or sum(int(r['bytes']) for r in rows)!=rec['bytes']:errors.append('RECEIPT_COUNTS')
 previous=0;inputs=0
 if not args.current_only:
  for r in csv.DictReader((ROOT/'results'/'PRIOR_ELEVEN_STAGE_SEAL_CHECKS_FINAL.csv').open(encoding='utf-8-sig')):verify(Path(r['path']),r['expected_sha256'],'PRIOR/'+r['path']);previous+=1
  for r in csv.DictReader((ROOT/'results'/'INPUT_SOURCE_MANIFEST.csv').open(encoding='utf-8-sig')):verify(Path(r['path']),r['sha256'],'INPUT/'+r['path']);inputs+=1
 out={'mode':'READ_ONLY_VERIFY','current_files':len(rows),'previous_stage_entries':previous,'raw_inputs':inputs,'unique_hashed_paths':len(cache),'errors':len(errors),'details':errors[:30]};print(json.dumps(out,ensure_ascii=True,indent=2));return int(bool(errors))
if __name__=='__main__':sys.exit(main())
