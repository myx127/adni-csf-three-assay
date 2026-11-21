from pathlib import Path
import json,hashlib,sys

ROOT=Path(__file__).resolve().parents[1]
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b''):h.update(chunk)
    return h.hexdigest()
def verify(verbose=True):
    import csv
    errors=[];manifest=ROOT/'SHA256_MANIFEST.csv'
    if not manifest.exists():return ['Missing SHA256_MANIFEST.csv']
    with manifest.open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    for row in rows:
        p=(ROOT/row['relative_path']).resolve()
        if not p.is_relative_to(ROOT.resolve()):errors.append('Path escapes delivery root: '+row['relative_path']);continue
        if not p.is_file():errors.append('Missing '+row['relative_path']);continue
        if p.stat().st_size!=int(row['bytes']) or digest(p)!=row['sha256']:errors.append('Changed '+row['relative_path'])
    status=json.loads((ROOT/'STAGE9_COMPLETION_STATUS.json').read_text(encoding='utf-8'))
    if status.get('status')!='COMPLETE':errors.append('Execution status is not COMPLETE')
    if status.get('strict_relative_pattern_validation_completed') is not False:errors.append('Incorrect strict validation claim')
    if status.get('external_individual_replication_completed') is not False:errors.append('Incorrect external replication claim')
    if verbose:print(json.dumps({'files_checked':len(rows),'errors':len(errors),'status':status.get('status'),'strict_validation':'NOT_ESTABLISHED','external_replication':'NOT_COMPLETED'},ensure_ascii=False))
    for e in errors[:20]:print(e)
    return errors
if __name__=='__main__':sys.exit(1 if verify() else 0)
