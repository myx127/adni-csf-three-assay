from pathlib import Path
import ast
import csv
import numpy as np
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

from adni_csf import pathology, clinical, conditioning
from adni_csf.multiplicity import by
from run_synthetic_demo import run
from check_published_results import check

ROOT = Path(__file__).resolve().parents[1]


def test_shared_direction_requires_all_assays_and_accounts_for_direction_choice():
    plus, minus, p, direction = pathology.joint_from_z(np.array([[3.,3.,3.],[-3.,-3.,-3.],[3.,3.,-3.],[3.,3.,0.]]))
    assert p[0] < .01 and p[1] < .01
    np.testing.assert_allclose(p[0], p[1])
    assert p[2] == 1 and p[3] == 1
    assert direction[:2].tolist() == [1,-1]


def test_by_matches_independent_statsmodels_implementation_with_p1_slots():
    p = np.r_[np.geomspace(1e-7, .9, 66), 1., 1.]
    np.testing.assert_allclose(by(p), multipletests(p, method='fdr_by')[1], atol=1e-14)


def test_wald_matches_direct_quadratic_form():
    b = np.array([[.1,.2], [.3,-.2]])
    V = np.array([[[.04,.01],[.01,.09]], [[.09,-.02],[-.02,.04]]])
    q, p, valid = clinical.wald_batch(b, V)
    expected = np.array([one @ np.linalg.inv(v) @ one for one,v in zip(b,V)])
    np.testing.assert_allclose(q, expected)
    np.testing.assert_allclose(p, chi2.sf(expected, 2))
    assert valid.all()


def test_coefficients_match_independent_weighted_least_squares():
    example = run()
    sw = np.sqrt(example['weights'])
    for model, Z in enumerate(example['Z']):
        expected = []
        for j in range(3):
            design = np.c_[Z, example['X'][:,j]] * sw[:,None]
            expected.append(np.linalg.lstsq(design, example['y'] * sw, rcond=None)[0][-1])
        np.testing.assert_allclose(example['fits'][model]['b'][0], expected, rtol=1e-10, atol=1e-12)


def test_adjustment_change_preserves_between_model_covariance():
    e = run()
    U = np.c_[e['scores'][0], e['scores'][1]]
    L = np.c_[-clinical.D, clinical.D]
    expected = L @ (U.T @ U) @ L.T
    np.testing.assert_allclose(e['change']['V'][0], expected, rtol=1e-10, atol=1e-12)
    # Treating the two models as independent must not reproduce this paired variance.
    independent = clinical.D @ (e['scores'][0].T @ e['scores'][0] + e['scores'][1].T @ e['scores'][1]) @ clinical.D.T
    assert not np.allclose(independent, expected, rtol=.01, atol=1e-9)


def test_figure_table_aggregate_values_preserve_fixed_families_and_results():
    assert check()['status'] == 'PASS'


def test_reusable_kernel_functions_are_unchanged_from_frozen_workflows():
    crosswalk = list(csv.DictReader((ROOT / 'metadata/source_code_crosswalk.csv').open(encoding='utf-8-sig')))
    rows = [r for r in crosswalk if r['release_path'].startswith('src/adni_csf/')]
    assert len(rows) == 3
    for r in rows:
        old = ast.parse((ROOT / 'workflows' / r['original_relative_path']).read_text(encoding='utf-8-sig'))
        new = ast.parse((ROOT / r['release_path']).read_text(encoding='utf-8'))
        f = lambda tree: [ast.dump(x, include_attributes=False) for x in tree.body if isinstance(x, ast.FunctionDef)]
        assert f(old) == f(new)
