"""Reproduction from fixed designs in a differently named sibling copy only."""
from config import *
import subprocess,argparse

def run():
 args=argparse.ArgumentParser();args.add_argument('--simulate',action='store_true');opt=args.parse_args()
 if ROOT.name=='6候选信号归因＋分子映射与正交验证':raise SystemExit('Use a differently named sibling copy under the same project directory; sealed original protected.')
 if not S5.exists():raise SystemExit('The reproduction copy must be a sibling of the original Stage6 directory so Stage1/4/5 remain available.')
 frozen=[p for p in RESULTS.glob('*.json') if any(w in p.name for w in ['FROZEN','FREEZE','FIRST_','ACCESS_BOUNDARY'])]
 initial={p:p.read_bytes() for p in frozen}
 save_json(RESULTS/'REPRODUCTION_CONTEXT.json',{'postdiscovery_reproduction':True,'not_first_quantitative_access':True,
  'original_freeze_bytes_preserved':True,'rerun_simulations':opt.simulate,'time':datetime.now().astimezone().isoformat()})
 def invoke(script):
  subprocess.run([sys.executable,'-X','utf8',str(ROOT/'code'/script)],cwd=ROOT,check=True)
 invoke('attribution.py')
 if opt.simulate:
  for script in ['calibrate_attribution.py','calibrate_equivalence.py','calibrate_support.py','four_assay.py','four_direction_boundary.py','calibrate_legacy.py']:invoke(script)
 for script in ['multiplicity_attribution.py','analyze_equivalence.py','analyze_support.py','analyze_legacy.py','finalize_evidence.py']:invoke(script)
 # No frozen artifact is redefined by a reproduction.
 changed=[str(p.relative_to(ROOT)) for p,b in initial.items() if p.read_bytes()!=b]
 if changed:raise RuntimeError('Frozen artifacts changed in reproduction copy: '+str(changed))
 for script in ['validate_stage6.py','make_figures.py','build_report.py','build_workbook.py','delivery_qa.py','seal_delivery.py','verify_delivery.py']:invoke(script)

if __name__=='__main__':run()
