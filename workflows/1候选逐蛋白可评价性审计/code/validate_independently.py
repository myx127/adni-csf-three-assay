"""Independent raw-CSV verification, without importing the main audit code.

Checks set membership, detection, missing-value codes and feature-level moments
using csv/math/statistics rather than pandas transformations and groupby.
Only platform measurements and annotation are opened; no clinical files.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import statistics
import zipfile
from collections import defaultdict
from datetime import date, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
RESULTS = ROOT / "results"


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def number(value, positive=False):
    try:
        x = float(value)
        return x if math.isfinite(x) and (not positive or x > 0) else None
    except (ValueError, TypeError):
        return None


def key_for(rid, examdate):
    n = number(rid)
    if n is None or n <= 0 or n != int(n):
        return None
    value = str(examdate).strip()
    try:
        d = date.fromisoformat(value[:10])
    except ValueError:
        try:
            # TMT release uses US month/day/year; 12/30 and 1/17 are unambiguous.
            d = datetime.strptime(value, "%m/%d/%Y").date()
        except ValueError:
            return None
    return str(int(n)) + "|" + d.isoformat()


def moments(values):
    values = [v for v in values if v is not None and math.isfinite(v)]
    return {"n_valid": len(values), "mean": statistics.mean(values) if values else None,
            "sd": statistics.stdev(values) if len(values) >= 2 else None,
            "variance": statistics.variance(values) if len(values) >= 2 else None}


def close(expected, observed):
    e = number(expected)
    if e is None or observed is None:
        return e is None and observed is None
    return math.isclose(e, observed, rel_tol=1e-10, abs_tol=1e-10)


def main():
    manifest = {r["input_id"]: r for r in read_csv(PROJECT / "规划和详细数据情况以及分析方案和研究计划/01_数据清单/CORE_INPUT_MANIFEST.csv")}
    mapping = read_csv(RESULTS / "EXTENDED_CORE_MAPPING_DRAFT.csv")
    feature_qc = read_csv(RESULTS / "FEATURE_LEVEL_QC.csv")
    main_summary = {r["accession"]: r for r in read_csv(RESULTS / "PROTEIN_EVALUABILITY_70.csv")}
    targets = {r["nulisa_target"] for r in mapping}
    soma_features = {f for r in mapping for f in r["soma_features"].split("|")}
    tmt_features = {r["tmt_columns"] for r in mapping}
    nkeys = set()
    nvalues = defaultdict(list)
    n_all_detect = {}
    n_any_detect = defaultdict(bool)
    n_known = defaultdict(bool)
    with Path(manifest["I01"]["path"]).open(encoding="utf-8-sig", newline="") as handle:
        for r in csv.DictReader(handle):
            if r["SampleMatrixType"] != "CSF" or r["SampleType"] != "Sample" or r["SampleQC"] != "passed":
                continue
            key = key_for(r["RID"], r["EXAMDATE"])
            if key is None:
                continue
            nkeys.add(key)
            target = r["Target"]
            if target not in targets:
                continue
            pair = key, target
            npq, lod = number(r["NPQ"]), number(r["TargetLOD_NPQ"])
            if npq is not None:
                nvalues[pair].append(npq)
            known = npq is not None and lod is not None
            n_known[pair] |= known
            n_any_detect[pair] |= known and npq >= lod
            n_all_detect[pair] = n_all_detect.get(pair, True) and known and npq >= lod
    nmeans = {pair: statistics.mean(v) for pair, v in nvalues.items() if v}
    srows = {}
    with Path(manifest["I04"]["path"]).open(encoding="utf-8-sig", newline="") as handle:
        for r in csv.DictReader(handle):
            key = key_for(r["RID"], r["EXAMDATE"])
            if key is None:
                continue
            if key in srows:
                raise RuntimeError("Independent validator found a SomaScan RID/date duplicate")
            srows[key] = {f: math.log2(v) if (v := number(r[f], positive=True)) is not None else None for f in soma_features}
    trows, tplex, gis = {}, {}, {}
    with zipfile.ZipFile(manifest["I07"]["path"]) as archive:
        with archive.open("EMORY_CSF_TMT_MS.csv") as binary:
            for r in csv.DictReader(io.TextIOWrapper(binary, encoding="utf-8-sig")):
                # Parsing independently from the primary regular expression.
                batch_label, channel = r["ProteomicID"].rsplit(".", 1)
                plex = int(batch_label.split("_b", 1)[1])
                key = key_for(r["RID"], r["EXAMDATE"])
                vals = {f: math.log2(v) if (v := number(r[f], positive=True)) is not None else None for f in tmt_features}
                if key is not None:
                    if key in trows:
                        raise RuntimeError("Independent validator found TMT RID/date duplicate")
                    trows[key], tplex[key] = vals, plex
                else:
                    if channel != "126" or plex in gis:
                        raise RuntimeError("Independent GIS standard/channel check failed")
                    gis[plex] = vals
    triple_all = sorted(nkeys & set(srows) & set(trows), key=lambda k: (int(k.split("|")[0]), k.split("|")[1]))
    seen, triple = set(), []
    for key in triple_all:
        rid = key.split("|")[0]
        if rid not in seen:
            seen.add(rid)
            triple.append(key)
    check_rows = [{"check": "independent_triple_registry_count", "accession": "ALL", "metric": "n_triple",
                   "expected": 571, "observed": len(triple), "passed": len(triple) == 571}]
    raw_pairings = read_csv(RESULTS / "PAIRWISE_SAMPLE_COUNTS.csv")
    sets = {"NULISA": nkeys, "SOMASCAN": set(srows), "TMT_MS": set(trows)}
    for row in raw_pairings:
        names = row["comparison"].split(" + ")
        shared = set.intersection(*(sets[name] for name in names))
        check_rows.append({"check": "independent_pairing_keys", "accession": "ALL", "metric": row["comparison"],
                           "expected": int(row["same_date_keys"]), "observed": len(shared), "passed": int(row["same_date_keys"]) == len(shared)})
    per_route = {}
    for r in mapping:
        protein, target, tf = r["accession"], r["nulisa_target"], r["tmt_columns"]
        sf = r["soma_features"].split("|")
        unambiguous = [f for f in r["soma_single_accession_features"].split("|") if f]
        observed = defaultdict(int)
        for key in triple:
            n = (key, target) in nmeans
            nd_any, nd_all = n_any_detect[(key, target)], n_all_detect.get((key, target), False)
            sv = [srows[key][f] is not None for f in sf]
            uv = [srows[key][f] is not None for f in unambiguous]
            sa, sall = any(sv), bool(sv) and all(sv)
            ua, uall = any(uv), bool(uv) and all(uv)
            tv = trows[key][tf] is not None
            gv = gis[tplex[key]][tf] is not None
            vals = {"N_numeric": n, "N_detect_any_repeat": nd_any, "N_detect_all_repeats": nd_all,
                    "S_any": sa, "S_all": sall, "S_unambiguous_any": ua, "S_unambiguous_all": uall,
                    "T_raw": tv, "GIS_positive": gv, "T_GIS_usable": tv and gv,
                    "complete_raw_S_any": n and sa and tv, "complete_raw_S_all": n and sall and tv,
                    "complete_GIS_S_any": n and sa and tv and gv, "complete_GIS_S_all": n and sall and tv and gv,
                    "complete_GIS_unambiguous_S_any": n and ua and tv and gv,
                    "complete_GIS_unambiguous_S_all": n and uall and tv and gv,
                    "complete_GIS_detect_all_N_S_any": n and nd_all and sa and tv and gv,
                    "complete_GIS_detect_all_N_S_all": n and nd_all and sall and tv and gv}
            for metric, val in vals.items():
                observed["n_" + metric] += int(val)
        for metric, count in observed.items():
            expected = int(main_summary[protein][metric])
            check_rows.append({"check": "independent_candidate_availability", "accession": protein, "metric": metric,
                               "expected": expected, "observed": count, "passed": expected == count})
        per_route[("NULISA", target)] = {key: nmeans.get((key, target)) for key in nkeys}
        for feature in sf:
            per_route[("SOMASCAN", feature)] = {key: row[feature] for key, row in srows.items()}
        per_route[("TMT_MS_RAW", tf)] = {key: row[tf] for key, row in trows.items()}
        per_route[("TMT_MS_GIS_RATIO", tf)] = {
            key: row[tf] - gis[tplex[key]][tf] if row[tf] is not None and gis[tplex[key]][tf] is not None else None
            for key, row in trows.items()}
    moment_rows = []
    for row in feature_qc:
        vals = per_route[(row["platform"], row["feature"])]
        selected = list(vals.values()) if row["scope"] == "NATIVE_RECORDS" else [vals.get(k) for k in triple]
        st = moments(selected)
        for metric, value in st.items():
            expected = row[metric]
            passed = close(expected, value)
            moment_rows.append({"accession": row["accession"], "platform": row["platform"], "feature": row["feature"],
                                "scope": row["scope"], "metric": metric, "expected": expected, "observed": value,
                                "passed": passed})
    write_csv(RESULTS / "INDEPENDENT_COUNT_VALIDATION.csv", check_rows)
    write_csv(RESULTS / "INDEPENDENT_NUMERIC_VALIDATION.csv", moment_rows)
    failed_counts = [x for x in check_rows if not x["passed"]]
    failed_moments = [x for x in moment_rows if not x["passed"]]
    summary = {"method": "Independent csv/math/statistics raw-file reading; no primary-audit helper imported",
               "clinical_files_opened": 0, "count_checks": len(check_rows), "numeric_moment_checks": len(moment_rows),
               "failed_count_checks": len(failed_counts), "failed_numeric_checks": len(failed_moments),
               "triple_participants": len(triple), "all_passed": not failed_counts and not failed_moments}
    (RESULTS / "INDEPENDENT_VALIDATION_SUMMARY.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    if not summary["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
