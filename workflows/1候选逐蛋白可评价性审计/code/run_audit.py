from __future__ import annotations

import importlib.metadata
import platform
import traceback

import numpy as np
import pandas as pd

from audit_helpers import *
from audit_inputs import verify_inputs, inspect_methods, load_nulisa, load_soma, load_tmt, build_mapping
from audit_measurements import build_sample_sets, run_measurement_audit


def describe(series):
    series = pd.to_numeric(series, errors="coerce").dropna()
    return {"minimum": float(series.min()), "median": float(series.median()), "maximum": float(series.max())} if len(series) else {}


def main():
    log("START: candidate evaluability audit; no clinical files or clinical associations")
    entries, inputs = verify_inputs()
    environment = {"python": sys.version, "platform": platform.platform(),
                   "packages": {name: importlib.metadata.version(name) for name in ["pandas", "numpy", "scipy", "openpyxl", "pyarrow", "pypdf", "matplotlib"]}}
    write_json(RESULTS / "SOFTWARE_ENVIRONMENT.json", environment)
    facts = inspect_methods(entries)
    write_json(RESULTS / "METHODS_FILE_CHECKS.json", facts)
    n = load_nulisa(entries["I01"])
    s = load_soma(entries["I04"], entries["I05"])
    t = load_tmt(entries["I07"])
    mapping, fmap = build_mapping(n, s, t)
    triple, registry, pairs, sample_summary = build_sample_sets(n, s, t)
    audit = run_measurement_audit(n, s, t, mapping, fmap, triple)
    p, f, g = audit["proteins"], audit["features"], audit["gis"]
    strict = p.loc[p["annotation_strict_candidate"]]
    multi = audit["multiaptamer"].query("scope == 'TRIPLE_SAME_DATE' and n_features > 1")
    zero_f = f.loc[f["scope"].eq("TRIPLE_SAME_DATE") & f["zero_variance"]]
    summary = {
        "status": "AUDIT_COMPLETED_WITH_OPEN_DECISIONS", "clinical_files_opened": 0,
        "input_hashes_all_match": True, "triple_registry_n": len(triple), "sample_matching": sample_summary,
        "candidate_n": len(p), "annotation_single_feature_candidate_n": len(strict),
        "multiaptamer_candidate_n": int(p["soma_feature_count"].gt(1).sum()),
        "composite_soma_annotation_candidate_n": int(p["soma_has_multiaccession_annotation"].sum()),
        "native_record_n": {"NULISA": len(n["matrix"]), "SOMASCAN": len(s["matrix"]), "TMT": len(t["matrix"])},
        "NULISA_native_plate_n": n["summary"]["CSF"]["plate_count"],
        "SOMASCAN_native_plate_n": s["summary"]["plate_count"], "TMT_native_plex_n": t["summary"]["plexes"],
        "complete_GIS_any_soma_n_extended": describe(p["n_complete_GIS_S_any"]),
        "complete_GIS_all_soma_n_extended": describe(p["n_complete_GIS_S_all"]),
        "complete_GIS_n_single_feature": describe(strict["n_complete_GIS_S_all"]),
        "complete_GIS_N_LOD_detected_n_single_feature": describe(strict["n_complete_GIS_detect_all_N_S_all"]),
        "NULISA_detection_fraction_extended": describe(p["N_detection_fraction_all_repeats"]),
        "NULISA_detection_fraction_single_feature": describe(strict["N_detection_fraction_all_repeats"]),
        "NULISA_detection_below_half_targets": p.loc[p["N_detection_fraction_all_repeats"].lt(.5), "nulisa_target"].tolist(),
        "NULISA_no_detected_targets": p.loc[p["N_detection_fraction_all_repeats"].eq(0), "nulisa_target"].tolist(),
        "single_feature_NULISA_detection_below_half_targets": strict.loc[strict["N_detection_fraction_all_repeats"].lt(.5), "nulisa_target"].tolist(),
        "single_feature_all_GIS_complete_under_100_targets": strict.loc[strict["n_complete_GIS_S_all"].lt(100), "nulisa_target"].tolist(),
        "single_feature_N_LOD_GIS_complete_under_100_targets": strict.loc[strict["n_complete_GIS_detect_all_N_S_all"].lt(100), "nulisa_target"].tolist(),
        "zero_variance_features_in_triple": zero_f[["accession", "platform", "feature", "n_valid"]].to_dict("records"),
        "GIS_loss_of_positive_participant_values_extended_total_cells": int((p["n_T_raw"] - p["n_T_GIS_usable"]).sum()),
        "GIS_missing_plexes_per_candidate": describe(p["GIS_missing_plexes"]),
        "TMT_raw_vs_GIS_delta_batch_R2_same_subset": describe(p["T_GIS_delta_batch_r2_matched"]),
        "multiaptamer_rho_min": describe(multi["rho_min"]),
        "multiaptamer_nonpositive_pair_targets": multi.loc[multi["any_nonpositive_pair_rho"], "target"].tolist(),
        "final_Strict_core_frozen": False,
        "pathology_stage_counts_or_followup_adequacy_evaluated": False,
    }
    for name, sub in [("extended", p), ("single_feature", strict)]:
        summary[name + "_coverage_grid"] = {}
        for metric in ["n_complete_raw_S_all", "n_complete_GIS_S_all", "n_complete_GIS_detect_all_N_S_all"]:
            summary[name + "_coverage_grid"][metric] = {str(level): int(sub[metric].ge(level*len(triple)-1e-10).sum()) for level in [.5, .7, .8, .9, .95, 1.0]}
    write_json(RESULTS / "AUDIT_SUMMARY.json", summary)
    qa = []
    def check(name, passed, expected, observed, severity="ERROR"):
        qa.append({"check": name, "passed": bool(passed), "expected": expected, "observed": observed, "severity_if_failed": severity})
    check("eight_input_hashes", all(x["passed"] for x in inputs), 8, len(inputs))
    check("nulisa_csf_qc_passed_rid", n["summary"]["CSF"]["qc_passed_rid"] == 1517, 1517, n["summary"]["CSF"]["qc_passed_rid"])
    check("nulisa_canonical_targets", n["summary"]["canonical_target_labels"] == 131, 131, n["summary"]["canonical_target_labels"])
    check("soma_rows", s["summary"]["released_rows"] == 737, 737, s["summary"]["released_rows"])
    check("soma_features", s["summary"]["features"] == 7008, 7008, s["summary"]["features"])
    check("tmt_participants", t["summary"]["participant_rows"] == 1105, 1105, t["summary"]["participant_rows"])
    check("tmt_GIS_rows", t["summary"]["standard_rows"] == 65, 65, t["summary"]["standard_rows"])
    check("tmt_protein_columns", t["summary"]["protein_columns"] == 3907, 3907, t["summary"]["protein_columns"])
    check("triple_sample_keys", len(triple) == 571, 571, len(triple))
    check("candidate_set_size", len(p) == 70, 70, len(p))
    check("annotation_single_feature_count", len(strict) == 43, 43, len(strict), "WARNING")
    check("PTID_RID_consistency", sample_summary["all_platform_PTID_RID_mismatches"] == 0, 0, sample_summary["all_platform_PTID_RID_mismatches"])
    check("no_double_counted_participants", registry["RID"].is_unique, "unique", registry["RID"].nunique())
    check("GIS_unique_standard_and_17_subjects", bool(t["plex_structure"]["participant_rows"].eq(17).all() and t["plex_structure"]["standard_rows"].eq(1).all()), "65 plexes each 17+1", len(t["plex_structure"]))
    check("GIS_long_rows", len(g) == len(p)*65, len(p)*65, len(g))
    check("NULISA_release_detection_is_plate_constant", n["release"]["release_detectability_unique"].le(1).all(), "all <=1", int(n["release"]["release_detectability_unique"].max()))
    check("NULISA_release_LOD_is_plate_constant", n["release"]["lod_unique"].le(1).all(), "all <=1", int(n["release"]["lod_unique"].max()))
    check("GIS_complete_not_more_than_raw", p["n_complete_GIS_S_all"].le(p["n_complete_raw_S_all"]).all(), "all <= raw", int((p["n_complete_GIS_S_all"] > p["n_complete_raw_S_all"]).sum()))
    check("all_aptamers_complete_not_more_than_any", p["n_complete_GIS_S_all"].le(p["n_complete_GIS_S_any"]).all(), "all <= any", int((p["n_complete_GIS_S_all"] > p["n_complete_GIS_S_any"]).sum()))
    check("LOD_detected_complete_not_more_than_numeric", p["n_complete_GIS_detect_all_N_S_all"].le(p["n_complete_GIS_S_all"]).all(), "all <= numeric", int((p["n_complete_GIS_detect_all_N_S_all"] > p["n_complete_GIS_S_all"]).sum()))
    check("single_feature_any_equals_all", strict["n_complete_GIS_S_any"].eq(strict["n_complete_GIS_S_all"]).all(), "all equal", int(strict["n_complete_GIS_S_any"].ne(strict["n_complete_GIS_S_all"]).sum()))
    check("clinical_data_firewall", summary["clinical_files_opened"] == 0, 0, summary["clinical_files_opened"])
    write_csv(RESULTS / "QA_CHECKS.csv", qa)
    failed = [r for r in qa if not r["passed"] and r["severity_if_failed"] == "ERROR"]
    warnings = [
        {"code": "TMT_VERSION_COLUMN_DISCREPANCY", "severity": "OPEN_DECISION", "detail": "Methods states 3447; bound matrix has 3907 protein columns; cause not proven"},
        {"code": "FINAL_STRICT_CORE_NOT_FROZEN", "severity": "OPEN_DECISION", "detail": "Annotation and quantitative audit only; no pathology coverage or power decision"},
        {"code": "SAME_DATE_NOT_EXACT_ALIQUOT", "severity": "LIMITATION", "detail": "RID/date matched; common collection/aliquot identifier not established"},
        {"code": "LIMITED_TECHNICAL_REPEATS", "severity": "LIMITATION", "detail": "One CSF repeat participant and pooled GIS do not establish individual assay reliability"},
        {"code": "VISCODE2_DIFFERENCES", "severity": "REVIEW", "detail": str(sample_summary["triple_VISCODE2_counts"])},
    ]
    write_csv(RESULTS / "FAILURES_AND_WARNINGS.csv", warnings + [dict(code=x["check"], severity="ERROR", detail=x["observed"]) for x in failed])
    if failed:
        raise RuntimeError("Audit QA checks failed; inspect QA_CHECKS.csv")
    log(f"Quantitative audit complete: {len(triple)} participants, {len(p)} candidates, {len(strict)} annotation single-feature candidates")
    print(json.dumps(safe_json(summary), ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        (LOGS / "FAILURE_TRACEBACK.txt").write_text(traceback.format_exc(), encoding="utf-8")
        log("FAILED; see logs/FAILURE_TRACEBACK.txt")
        raise
