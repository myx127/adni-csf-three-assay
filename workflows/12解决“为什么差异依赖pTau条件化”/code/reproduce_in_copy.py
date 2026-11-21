"""Explicit fresh-copy reproduction; never changes the sealed source or stages 1-11."""
from pathlib import Path
import shutil,subprocess,sys,argparse,re
ROOT=Path(__file__).resolve().parents[1]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--destination',required=True);args=ap.parse_args();dest=Path(args.destination).resolve();project=ROOT.parent.resolve()
 if dest.parent!=project or dest==ROOT.resolve() or not re.match(r'^12(?!\d)',dest.name):raise SystemExit('Destination must be a new sibling folder starting with 12, within this project.')
 if dest.exists() and any(dest.iterdir()):raise SystemExit('Destination must be empty; no files will be removed.')
 dest.mkdir(parents=True,exist_ok=True);(dest/'code').mkdir(exist_ok=True)
 for p in (ROOT/'code').glob('*.py'):shutil.copy2(p,dest/'code'/p.name)
 for n in ['STAGE12_PROTOCOL_v1.0.md','EXECUTION_AND_INFERENCE_BOUNDARY_v1.1.md','run_all.ps1']:shutil.copy2(ROOT/n,dest/n)
 pipeline=['bootstrap','background','prepare','validate_engine','calibrate','analyze','diagnostics','toy_and_summary','validate_science','make_figures','build_delivery','validate_delivery']
 for n in pipeline:
  command=[sys.executable,'-B',str(dest/'code'/(n+'.py'))]+(['--allow-unreviewed'] if n=='validate_delivery' else []);print('RUN '+n,flush=True);subprocess.run(command,cwd=dest,check=True)
 print('Analysis reproduced in a fresh copy. Review new figures and record their hashes before default delivery QA and sealing. No visual review was copied or invented; copy remains unsealed.')
if __name__=='__main__':main()
