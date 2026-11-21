"""Check released aggregate P/q values and the article's principal result counts."""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'src'))
from adni_csf.multiplicity import by


def check():
    root = REPO / 'aggregate_results'
    pathology = pd.read_csv(root / 'D013.csv')
    for _, f in pathology.groupby('model', sort=False):
        assert len(f) == 34
        np.testing.assert_allclose(by(f.p_used), f.q, rtol=1e-11, atol=1e-14)
        assert (f.q <= .05).sum() == 14
        assert set(f.loc[f.q <= .05, 'direction']) == {'Positive'}
    memory = pd.read_csv(root / 'D023.csv')
    assert len(memory) == 68
    np.testing.assert_allclose(by(memory.p_used), memory.q, rtol=1e-11, atol=1e-14)
    failed = memory.loc[memory.gene == 'CST3']
    assert len(failed) == 2 and (failed.p_used == 1).all()
    same = memory.loc[memory.test == 'Shared direction']
    assert len(same) == 34 and not (same.q <= .05).any()
    difference = memory.loc[(memory.test == 'Assay differences') & (memory.q <= .05)]
    assert set(difference.gene) == {'NPTX2', 'TAFA5', 'TIMP3'}
    return {'pathology_models': 3, 'shared_positive_proteins_per_model': 14, 'memory_family': 68, 'failed_calibration_slots_retained': 2, 'shared_memory_min_q': float(same.q.min()), 'memory_difference_candidates': sorted(difference.gene), 'status': 'PASS'}


if __name__ == '__main__':
    print(json.dumps(check(), indent=2))
