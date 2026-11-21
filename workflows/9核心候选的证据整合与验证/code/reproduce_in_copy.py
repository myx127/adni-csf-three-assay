from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import argparse,shutil,subprocess,sys

def main():
    source=Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser(description='Recompute Stage9 in a new sibling directory; sealed source remains unchanged.')
    parser.add_argument('--destination',type=Path)
    args=parser.parse_args()
    if (source/'SHA256_MANIFEST.csv').exists():
        result=subprocess.run([sys.executable,'-X','utf8',str(source/'code/verify_delivery.py')])
        if result.returncode:raise SystemExit('Source seal did not verify')
    destination=(args.destination or source.parent/('9_reproduction_'+datetime.now(ZoneInfo('Asia/Shanghai')).strftime('%Y%m%d_%H%M%S'))).resolve()
    if destination.parent!=source.parent.resolve() or destination==source.resolve():raise SystemExit('Choose a new sibling directory within the project root')
    if destination.exists():raise SystemExit('Destination already exists; choose a new directory')
    destination.mkdir()
    shutil.copytree(source/'code',destination/'code',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    shutil.copytree(source/'evidence',destination/'evidence')
    for name in ['STAGE9_PROTOCOL_v1.0.md','RAW_QC_EXTENSION_v1.1.md']:shutil.copy2(source/name,destination/name)
    steps=['bootstrap.py','independent_analysis.py','audit_metadata.py','audit_replicates_qc.py','audit_molecular.py','triangulate_dates.py','audit_public_sources.py','integrate_evidence.py','make_figures.py','validate_science.py','build_delivery.py','validate_delivery.py','seal_delivery.py','verify_delivery.py']
    for step in steps:
        print('RUN',step,flush=True)
        result=subprocess.run([sys.executable,'-X','utf8',str(destination/'code'/step)],cwd=destination)
        if result.returncode:raise SystemExit(f'Reproduction stopped at {step}; inspect new directory {destination}')
    print('Reproduction complete:',destination)

if __name__=='__main__':main()
