from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
HANDOFF = PROJECT / "规划和详细数据情况以及分析方案和研究计划"
RESULTS = ROOT / "results"
RESTRICTED = ROOT / "00_restricted"
FIGURES = ROOT / "figures"
LOGS = ROOT / "logs"
for directory in [RESULTS, RESTRICTED, FIGURES, LOGS]:
    directory.mkdir(exist_ok=True, parents=True)

MISSING_TEXT = {"", "NA", "N/A", "NAN", "NONE", "NULL", "NAT"}
UNIPROT = re.compile(r"(?:(?:[OPQ][0-9][A-Z0-9]{3}[0-9])|(?:[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}))(?:-\d+)?$")
ALIASES = {f"A{variant}{length}": f"Aβ{length}" for variant in ["?", "Î²"] for length in ["38", "40", "42"]}
SPECIAL_TARGETS = {"APOE4", "Aβ38", "Aβ40", "Aβ42", "BD-MAPT", "BD-pTau-181", "BD-pTau-217", "BD-pTau-231",
                   "pTau-181", "pTau-217", "pTau-231", "Oligo-SNCA", "pSNCA-129", "pTDP43-409", "IL12p70"}
CSF_UNAVAILABLE = {"APOE", "CRP"}


def log(message):
    line = datetime.now().strftime("%H:%M:%S") + " " + str(message)
    print(line, flush=True)
    with (LOGS / "AUDIT_RUN.log").open("a", encoding="utf-8") as f:
        f.write(line + "\n")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def clean_text(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    value = str(value).strip()
    return "" if value.upper() in MISSING_TEXT else value


def join_unique(values):
    return "|".join(sorted({clean_text(v) for v in values if clean_text(v)}))


def accession(value):
    return re.sub(r"-\d+$", "", clean_text(value))


def tokens(value):
    return {accession(x) for x in re.split(r"[|;,\s]+", clean_text(value)) if x}


def canonical_target(value):
    value = clean_text(value)
    return ALIASES.get(value, value)


def write_csv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(exist_ok=True, parents=True)
    if isinstance(rows, pd.DataFrame):
        rows.to_csv(path, index=False, encoding="utf-8-sig")
        return
    rows = list(rows)
    if fields is None:
        fields = list(dict.fromkeys(k for r in rows for k in r)) if rows else ["status"]
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def safe_json(value):
    if isinstance(value, dict):
        return {str(k): safe_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [safe_json(v) for v in value]
    if isinstance(value, np.ndarray):
        return safe_json(value.tolist())
    if isinstance(value, np.generic):
        return safe_json(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(path, value):
    Path(path).write_text(json.dumps(safe_json(value), ensure_ascii=False, indent=2), encoding="utf-8")


def make_keys(frame, label):
    rid = pd.to_numeric(frame["RID"], errors="coerce")
    valid_rid = rid.notna() & (rid > 0) & rid.eq(np.floor(rid))
    dates = pd.to_datetime(frame["EXAMDATE"].replace({"NA": "", "": None}), errors="coerce")
    valid = valid_rid & dates.notna()
    key = pd.Series("", index=frame.index, dtype=object)
    key.loc[valid] = rid.loc[valid].astype(int).astype(str) + "|" + dates.loc[valid].dt.strftime("%Y-%m-%d")
    return key, {"platform": label, "rows": len(frame), "invalid_rid_rows": int((~valid_rid).sum()),
                 "invalid_date_with_valid_rid_rows": int((valid_rid & dates.isna()).sum()),
                 "unique_valid_rid": int(rid.loc[valid].nunique()), "unique_valid_keys": int(key.loc[valid].nunique())}


def numerical_stats(values):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    out = {"n_valid": n, "n_unique_values": len(np.unique(x)), "mean": np.nan, "sd": np.nan,
           "variance": np.nan, "mad_unscaled": np.nan, "iqr": np.nan, "minimum": np.nan, "maximum": np.nan,
           "p01": np.nan, "p05": np.nan, "p25": np.nan, "p50": np.nan, "p75": np.nan, "p95": np.nan, "p99": np.nan}
    if n:
        ps = np.quantile(x, [.01, .05, .25, .5, .75, .95, .99])
        out.update({"mean": float(x.mean()), "minimum": float(x.min()), "maximum": float(x.max()),
                    "mad_unscaled": float(np.median(np.abs(x - np.median(x)))), "iqr": float(ps[4] - ps[2])})
        out.update(dict(zip(["p01", "p05", "p25", "p50", "p75", "p95", "p99"], ps)))
    if n >= 2:
        out["variance"] = float(np.var(x, ddof=1))
        out["sd"] = float(np.std(x, ddof=1))
    out["zero_variance"] = bool(n >= 2 and out["sd"] == 0)
    out["near_zero_variance_sd_le_1e_6"] = bool(n >= 2 and out["sd"] <= 1e-6)
    out["low_distinct_values_le_5"] = bool(n and out["n_unique_values"] <= 5)
    return out


def batch_effect(values, batches):
    values = np.asarray(values, dtype=float)
    batches = np.asarray([clean_text(x) for x in batches], dtype=object)
    keep = np.isfinite(values) & (batches != "")
    x, b = values[keep], batches[keep]
    labels = np.unique(b)
    n, k = len(x), len(labels)
    result = {"n_batch_effect": n, "n_batches_with_values": k, "batch_r2": np.nan,
              "batch_adjusted_r2": np.nan, "smallest_batch_n_valid": 0, "largest_batch_n_valid": 0}
    if n and k:
        sizes = [int((b == v).sum()) for v in labels]
        result.update({"smallest_batch_n_valid": min(sizes), "largest_batch_n_valid": max(sizes)})
    if n >= 2:
        total = float(np.sum((x - x.mean())**2))
        if total > 0:
            between = sum(float(np.sum(b == label)) * float((x[b == label].mean() - x.mean())**2) for label in labels)
            result["batch_r2"] = min(1.0, max(0.0, between / total))
            if n > k:
                result["batch_adjusted_r2"] = 1 - (1 - result["batch_r2"]) * (n - 1) / (n - k)
    return result


def per_batch_rows(values, batches, base):
    values = np.asarray(values, dtype=float)
    batches = np.asarray([clean_text(x) for x in batches], dtype=object)
    rows = []
    for label in sorted(set(batches) - {""}):
        mask = batches == label
        st = numerical_stats(values[mask])
        rows.append({**base, "batch": label, "n_records": int(mask.sum()), "n_valid": st["n_valid"],
                     "missing_fraction": 1 - st["n_valid"] / int(mask.sum()), "mean_log_scale": st["mean"],
                     "sd_log_scale": st["sd"], "median_log_scale": st["p50"]})
    return rows
