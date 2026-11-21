from pathlib import Path
import csv,json,hashlib,sys
ROOT=Path(__file__).resolve().parents[1]
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def verify(verbose=True):
 errors=[];manifest=ROOT/'SHA256_MANIFEST.csv'
 if not manifest.exists():return ['Missing SHA256_MANIFEST.csv']
 with manifest.open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
 for row in rows:
  p=(ROOT/row['relative_path']).resolve()
  if not p.is_relative_to(ROOT.resolve()):errors.append('Escaping path '+row['relative_path']);continue
  if not p.is_file():errors.append('Missing '+row['relative_path']);continue
  if p.stat().st_size!=int(row['bytes']) or digest(p)!=row['sha256']:errors.append('Changed '+row['relative_path'])
 status=json.loads((ROOT/'STAGE10_COMPLETION_STATUS.json').read_text(encoding='utf-8'))
 if status.get('status')!='COMPLETE':errors.append('Execution status not COMPLETE')
 if status.get('new_BY_clinical_discoveries')!=0:errors.append('Incorrect clinical discovery claim')
 if status.get('clinical_prognostic_validation_established') is not False:errors.append('Incorrect validation claim')
 if verbose:print(json.dumps({'files_checked':len(rows),'errors':len(errors),'status':status.get('status'),'clinical_BY_discoveries':status.get('new_BY_clinical_discoveries'),'unqualified_amyloid_interaction_designs':3},ensure_ascii=False))
 for e in errors[:20]:print(e)
 return errors
if __name__=='__main__':sys.exit(1 if verify() else 0)
