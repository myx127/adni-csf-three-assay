from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import argparse,subprocess,shutil,json,sys
def run():
 source=Path(__file__).resolve().parents[1];parser=argparse.ArgumentParser(description='Recompute Stage10 in a new sibling directory; preserve all sealed results.')
 parser.add_argument('--destination',type=Path);args=parser.parse_args()
 if (source/'SHA256_MANIFEST.csv').exists():
  if subprocess.run([sys.executable,'-X','utf8',str(source/'code/verify_delivery.py')]).returncode:raise SystemExit('Source verification failed')
 destination=(args.destination or source.parent/('10_reproduction_'+datetime.now(ZoneInfo('Asia/Shanghai')).strftime('%Y%m%d_%H%M%S'))).resolve()
 if destination.parent!=source.parent.resolve() or destination==source.resolve():raise SystemExit('Use a new sibling directory under the project root')
 if destination.exists():raise SystemExit('Destination already exists')
 destination.mkdir();shutil.copytree(source/'code',destination/'code',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 for name in ['STAGE10_PROTOCOL_v1.0.md','CALIBRATION_EXECUTION_DETAIL_v1.1.md','DESCRIPTIVE_CONDITIONING_EXTENSION_v1.2.md']:shutil.copy2(source/name,destination/name)
 (destination/'REPRODUCTION_CONTEXT.json').write_text(json.dumps({'source':str(source),'original_Stage10_results_already_known':True,'purpose':'Computational reproduction, not new prospective analysis'},ensure_ascii=False,indent=2),encoding='utf-8')
 steps=['bootstrap.py','prepare_clinical.py','prepare_templates.py','validate_engine.py','run_calibration.py','describe_followup.py','analyze_real.py','describe_conditioning.py','summarize_calibration.py','make_figures.py','validate_science.py','build_delivery.py','validate_delivery.py','seal_delivery.py','verify_delivery.py']
 for name in steps:
  print('RUN',name,flush=True);result=subprocess.run([sys.executable,'-X','utf8',str(destination/'code'/name)],cwd=destination)
  if result.returncode:raise SystemExit('Stopped at '+name+'; inspect '+str(destination))
 print('Reproduction completed:',destination)
if __name__=='__main__':run()
