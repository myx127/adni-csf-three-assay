"""Recreate figures from published aggregates in a separate output directory."""
from pathlib import Path
import argparse
import csv
import importlib.util
import json
import shutil
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties as OriginalFontProperties, findfont

REPO = Path(__file__).resolve().parents[1]


def portable_font(*args, **kwargs):
    filename = kwargs.get('fname')
    if filename and not Path(filename).exists():
        kwargs.pop('fname', None)
        kwargs['family'] = 'DejaVu Sans'
    return OriginalFontProperties(*args, **kwargs)


def render(keys, output, lang='EN'):
    output = Path(output).resolve()
    if output == REPO or REPO / 'visualization' == output or REPO / 'visualization' in output.parents:
        raise ValueError('Choose an output directory separate from visualization sources.')
    index = {r['figure']: r for r in csv.DictReader((REPO / 'metadata/figure_code_crosswalk.csv').open(encoding='utf-8-sig'))}
    reports = []
    for key in keys:
        row = index[key]
        dest = output / key
        if dest.exists(): raise FileExistsError('Choose a fresh output directory: ' + str(dest))
        shutil.copytree(REPO / 'visualization' / key, dest)
        for name in ['editable', 'high_resolution', 'build', 'source_data']:
            (dest / name).mkdir(exist_ok=True)
        script = dest / 'source' / Path(row['draw_script']).name
        spec = importlib.util.spec_from_file_location('display_' + key, script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        # Preserve available original fonts; provide a portable English fallback.
        module.FontProperties = portable_font
        klass = module.Diagram if key == 'supp1' else module.Figure
        fig = klass(lang)
        if row['requires_draw_call'] == 'True': fig = fig.draw()
        report = fig.save()
        plt.close(fig.fig)
        reports.append({'figure': key, 'language': lang, 'report': report})
        print('Rendered ' + key + ' ' + lang, flush=True)
    (output / 'RENDER_REPORT.json').write_text(json.dumps(reports, ensure_ascii=False, indent=2), encoding='utf-8')
    return reports


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--figures', nargs='+', default=['main' + str(i) for i in range(1,6)])
    parser.add_argument('--output', type=Path, default=REPO / 'outputs/figures')
    parser.add_argument('--lang', choices=['EN', 'ZH'], default='EN')
    a = parser.parse_args()
    render(a.figures, a.output, a.lang)
