"""Illustrate paired assays using simulated records and the final study kernels."""
from pathlib import Path
import argparse
import json
import sys

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'src'))
from adni_csf import clinical, conditioning


def template(X, Z, groups, weights):
    """Build fixed-scale coefficient/cluster operators for the synthetic example."""
    sw = np.sqrt(weights)
    Z = sw[:, None] * Z
    A, axes = [], [[] for _ in range(3)]
    Qz = clinical.basis(Z)
    for j in range(3):
        x = sw * X[:, j]
        r = x - Qz @ (Qz.T @ x)
        a = r / (r @ r)
        Q = np.c_[Qz, r / np.sqrt(r @ r)]
        A.append(a)
        for axis in range(3):
            lab = groups[:, axis]
            G = lab.max() + 1
            B = np.eye(G)[lab].T * a[None, :]
            K = B - (B @ Q) @ Q.T
            K = (K - K.mean(axis=0, keepdims=True)) * np.sqrt(G / (G - 1))
            axes[axis].append(K)
    K = np.concatenate([np.stack(one, axis=1) for one in axes], axis=0)
    # Scales are known from this synthetic generator and fixed across models.
    # Production templates retain nonzero scale-estimation influence scores.
    return {'A': np.asarray(A), 'K': K, 'sigma': np.ones(3), 'Uv': np.zeros((len(K), 3))}


def run(output=None, seed=2026112101, n=360):
    rng = np.random.default_rng(seed)
    p = rng.normal(size=n)
    u = rng.normal(size=n)
    common = .75 * p + np.sqrt(1 - .75**2) * u
    errors = np.asarray([.45, .15, .8])
    X = (common[:, None] + rng.normal(size=(n, 3)) * errors) / np.sqrt(1 + errors**2)
    y = .20 * common - .50 * p + rng.normal(scale=.25, size=n)
    age = rng.normal(size=n)
    sex = rng.integers(0, 2, size=n)
    groups = np.column_stack([rng.integers(0, g, size=n) for g in [12, 13, 15]])
    weights = rng.uniform(.5, 2, size=n)
    Z0 = np.c_[np.ones(n), age, sex]
    for axis in range(3):
        Z0 = np.c_[Z0, pd.get_dummies(groups[:, axis], drop_first=True, dtype=float)]
    Z1 = np.c_[Z0, p]
    ts = [template(X, Z, groups, weights) for Z in [Z0, Z1]]
    yw = y * np.sqrt(weights)
    fits = [clinical.joint_stats(yw, t) for t in ts]
    scores = [np.einsum('gjn,n->gj', t['K'], yw) for t in ts]
    L = np.c_[-clinical.D, clinical.D]
    change_t = {'L': L, 'sigma': np.ones(6), 'Uv': np.zeros((scores[0].shape[0], 6))}
    change = conditioning.moment_stats(np.c_[fits[0]['b'], fits[1]['b']], np.c_[scores[0], scores[1]][None, :, :], change_t, full=True)
    rows = []
    for label, fit in zip(['without_p_tau', 'with_p_tau'], fits):
        for j, assay in enumerate(['Assay_A', 'Assay_B', 'Assay_C']):
            rows.append({'data': 'SYNTHETIC', 'model': label, 'assay': assay, 'coefficient': float(fit['b'][0, j]), 'shared_direction_p': float(fit['p_same'][0]), 'assay_difference_p': float(fit['p_diff'][0])})
    summary = {'data': 'SYNTHETIC', 'seed': seed, 'n': n, 'scope': 'Illustration using simulated annual memory changes and known fixed exposure scales; not a reanalysis of ADNI or an article estimate', 'paired_adjustment_p': float(change['p'][0]), 'paired_changes': change['h'][0].tolist()}
    if output is not None:
        output = Path(output); output.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(output / 'synthetic_assay_associations.csv', index=False)
        (output / 'synthetic_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return {'X': X, 'y': y, 'Z': [Z0, Z1], 'weights': weights, 'templates': ts, 'fits': fits, 'scores': scores, 'change': change, 'summary': summary}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=REPO / 'outputs/synthetic_demo')
    args = parser.parse_args()
    print(json.dumps(run(args.output)['summary'], indent=2))
