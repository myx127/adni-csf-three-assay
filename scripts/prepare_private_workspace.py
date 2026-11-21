"""Prepare source layout and input manifests; never execute participant analyses."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re
import shutil

import yaml

REPO = Path(__file__).resolve().parents[1]


def run(config, output, check_only=False):
    paths = yaml.safe_load(Path(config).read_text(encoding='utf-8'))['inputs']
    rows = json.loads((REPO / 'metadata/input_manifest_templates.json').read_text(encoding='utf-8'))
    records = []
    for item in rows:
        row = dict(item['row']); key = row['input_id']
        path = Path(paths[key]).expanduser().resolve()
        if not path.is_file(): raise FileNotFoundError('Missing authorized input ' + key + ': provide its local path in your private config.')
        h = hashlib.sha256()
        with path.open('rb') as f:
            for block in iter(lambda: f.read(4*1024*1024), b''): h.update(block)
        if h.hexdigest() != row['sha256']: raise ValueError('Input ' + key + ' differs from the study release. Review release compatibility before reanalysis.')
        row.update(path=str(path), bytes=path.stat().st_size)
        records.append((item['origin'], row))
    if check_only: return {'inputs_verified': len(records), 'workspace_created': False}
    output = Path(output).resolve()
    if output == REPO or REPO in output.parents: raise ValueError('The private workspace must be outside the release repository.')
    if output.exists(): raise FileExistsError('Use a new private workspace directory.')
    shutil.copytree(REPO / 'workflows', output)
    # Stage 9 originally ran before stages 10-12 existed. Its prefix-only lookup
    # would otherwise confuse stage 1 and stage 10 in a newly complete workspace.
    config9 = next(p for p in output.iterdir() if p.is_dir() and re.match(r'^9(?!\d)', p.name)) / 'code/config.py'
    text9 = config9.read_text(encoding='utf-8-sig')
    old = "p.name.startswith(str(i))"
    assert text9.count(old) == 1
    text9 = text9.replace(old, "__import__('re').match(r'^'+str(i)+r'(?!\\d)',p.name)", 1)
    config9.write_text(text9, encoding='utf-8')
    for origin in ['core', 'additional']:
        rr = [r for o, r in records if o == origin]
        path = output / ('规划和详细数据情况以及分析方案和研究计划/01_数据清单/CORE_INPUT_MANIFEST.csv' if origin == 'core' else next(p.name for p in output.iterdir() if p.is_dir() and p.name.startswith('4')) + '/results/ADDITIONAL_INPUT_MANIFEST.csv')
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('w', encoding='utf-8-sig', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rr[0])); writer.writeheader(); writer.writerows(rr)
    return {'inputs_verified': len(records), 'workspace_created': True, 'analyses_executed': False}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--check-only', action='store_true')
    a = p.parse_args()
    print(json.dumps(run(a.config, a.output, a.check_only), indent=2))
