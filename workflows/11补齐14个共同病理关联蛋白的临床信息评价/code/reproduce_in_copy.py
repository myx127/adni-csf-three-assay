from pathlib import Path
from datetime import datetime
import shutil,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
def run():
 dest=Path(sys.argv[1]) if len(sys.argv)>1 else ROOT.with_name(ROOT.name+'_reproduce_'+datetime.now().strftime('%Y%m%d_%H%M%S'))
 dest=dest.resolve()
 if not dest.is_relative_to(ROOT.parent.resolve()) or dest==ROOT.resolve() or dest.exists():raise SystemExit('Choose a fresh sibling directory within the project.')
 dest.mkdir()
 shutil.copytree(ROOT/'code',dest/'code',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 shutil.copyfile(ROOT/'STAGE11_PROTOCOL_v1.0.md',dest/'STAGE11_PROTOCOL_v1.0.md')
 for script in ['bootstrap','prepare_templates','validate_engine','run_calibration','analyze_real','summarize_and_describe','validate_science','make_figures','build_delivery']:
  subprocess.run([sys.executable,'-B',str(dest/'code'/(script+'.py'))],cwd=dest,check=True)
 print('Fresh numerical reproduction complete. Review eight figures, create FIGURE_REVIEW_RECORD.csv, then validate_delivery and seal_delivery.')
if __name__=='__main__':run()

