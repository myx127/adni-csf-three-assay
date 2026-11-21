from __future__ import annotations

import csv
import json
import math
import re
from datetime import datetime
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook

from audit_helpers import ROOT, RESULTS, RESTRICTED, digest, write_csv, write_json


def main():
    checks = []
    def check(name, passed, detail):
        checks.append({"check": name, "passed": bool(passed), "detail": detail})

    q = pd.read_csv(RESULTS / "QA_CHECKS.csv")
    check("primary_QA", q.passed.all(), {"passed": int(q.passed.sum()), "total": len(q)})
    independent = json.loads((RESULTS / "INDEPENDENT_VALIDATION_SUMMARY.json").read_text(encoding="utf-8"))
    check("independent_raw_validation", independent["all_passed"] and not independent["failed_count_checks"] and not independent["failed_numeric_checks"], independent)
    access = json.loads((RESULTS / "ACCESS_BOUNDARY_AUDIT.json").read_text(encoding="utf-8"))
    check("protocol_unchanged_since_analysis", digest(ROOT / "AUDIT_PROTOCOL_v1.0.md") == access["protocol_sha256_before_analysis"], access["protocol_sha256_before_analysis"])
    check("clinical_input_boundary", access["clinical_files_opened"] == 0 and not access["opened_clinical_input_ids"] and not access["clinical_outcome_linkage"], "Code-boundary record; not an OS file-access audit")
    inputs = pd.read_csv(RESULTS / "INPUT_HASH_VERIFICATION.csv")
    input_checks = []
    for row in inputs.itertuples(index=False):
        path = Path(row.path)
        observed = digest(path)
        passed = observed == row.expected_sha256 and path.stat().st_size == row.expected_bytes
        input_checks.append({"input_id": row.input_id, "still_matches": passed, "bytes": path.stat().st_size,
                             "sha256": observed, "path": str(path)})
    check("eight_inputs_still_unchanged", len(input_checks) == 8 and all(r["still_matches"] for r in input_checks), input_checks)
    p = pd.read_csv(RESULTS / "PROTEIN_EVALUABILITY_70.csv")
    check("70_candidates_unique_accessions", len(p) == 70 and p.accession.is_unique, len(p))
    check("43_annotation_single_feature_only", int(p.annotation_strict_candidate.sum()) == 43, int(p.annotation_strict_candidate.sum()))
    check("no_SAP_or_molecular_identity_claim", not p.final_SAP_membership_frozen.any() and not p.molecular_identity_validated.any(), "All frozen/identity-validated flags false")
    registry = pd.read_csv(RESTRICTED / "TRIPLE_SAMPLE_REGISTRY_RESTRICTED.csv")
    check("restricted_registry_571_unique", len(registry) == 571 and registry.RID.is_unique and registry.key.is_unique, len(registry))
    masks = pd.read_parquet(RESTRICTED / "TRIPLE_CANDIDATE_AVAILABILITY_RESTRICTED.parquet")
    check("restricted_boolean_matrix_dimensions", len(masks) == 70*571 and not masks.duplicated(["key", "accession"]).any(), {"rows": len(masks), "columns": len(masks.columns)})
    feature = pd.read_csv(RESULTS / "FEATURE_LEVEL_QC.csv")
    check("feature_scope_and_dimensions", len(feature) == 630 and set(feature.scope) == {"NATIVE_RECORDS", "TRIPLE_SAME_DATE"}, "315 feature/processing routes x 2 scopes")
    check("GIS_70_by65", len(pd.read_csv(RESULTS / "TMT_GIS_BATCH_AUDIT.csv")) == 4550, "4,550 rows")
    example = pd.read_csv(RESULTS / "DESCRIPTIVE_80PCT_PRIORITY_CANDIDATES.csv")
    expected_example = set(p.loc[p.annotation_strict_candidate & p.n_complete_GIS_detect_all_N_S_all.ge(math.ceil(571*.8)), "accession"])
    check("illustrative_30_subset", len(example) == 30 and set(example.accession) == expected_example, "Descriptive >=457/571 only; not frozen")

    schema = pd.read_csv(RESULTS / "WORKBOOK_SHEET_COUNTS.csv")
    wb = load_workbook(ROOT / "候选逐蛋白可评价性审计总表.xlsx", read_only=True, data_only=False)
    workbook_counts = []
    all_workbook_ok = True
    forbidden_headers = {"RID", "PTID", "EXAMDATE", "key", "NULISA_PTID", "TMT_PTID", "GUSPECID", "ProteomicID"}
    aggregate_only = True
    for row in schema.itertuples(index=False):
        ws = wb[row.sheet]
        passed = ws.max_row == row.data_rows+1 and ws.max_column == row.columns
        all_workbook_ok &= passed
        headers = list(next(ws.iter_rows(min_row=1, max_row=1, values_only=True)))
        aggregate_only &= not bool(set(headers) & forbidden_headers)
        workbook_counts.append({"sheet": row.sheet, "rows": ws.max_row-1, "columns": ws.max_column, "passed": passed})
    check("workbook_sheet_dimensions", all_workbook_ok and len(wb.sheetnames) == len(schema), workbook_counts)
    check("workbook_no_participant_headers", aggregate_only, "Aggregate outputs only; source paths and batch labels retained")
    ws = wb["70候选总览"]
    rows = list(ws.values)
    header = list(rows[0])
    accession_index = header.index("UniProt")
    complete_index = header.index("三平台完整人数_GIS_Soma全部")
    lod_index = header.index("完整人数_GIS_Soma全部_NULISA达到LOD")
    actual = {r[accession_index]: (r[complete_index], r[lod_index]) for r in rows[1:]}
    expected = {r.accession: (r.n_complete_GIS_S_all, r.n_complete_GIS_detect_all_N_S_all) for r in p.itertuples()}
    check("workbook_counts_match_csv_for_all70", actual == expected, "Checked GIS and GIS+LOD counts for all70")
    wb.close()
    required = ["README_审计总览.md", "STAGE1_QA_REPORT_ZH.md", "TMT_COLUMN_COUNT_RECONCILIATION.md", "候选逐蛋白可评价性审计总表.xlsx",
                "figures/01_CANDIDATE_AVAILABILITY_70.png", "figures/01_CANDIDATE_AVAILABILITY_70.pdf",
                "figures/02_TMT_BATCH_R2_BEFORE_AFTER_GIS.png", "figures/02_TMT_BATCH_R2_BEFORE_AFTER_GIS.pdf",
                "code/run_all.ps1", "code/validate_independently.py"]
    check("required_deliverables_exist", all((ROOT / name).is_file() and (ROOT / name).stat().st_size for name in required), required)
    missing_links = []
    for filename in ["README_审计总览.md", "STAGE1_QA_REPORT_ZH.md", "TMT_COLUMN_COUNT_RECONCILIATION.md"]:
        for link in re.findall(r"\]\(([^)]+)\)", (ROOT / filename).read_text(encoding="utf-8")):
            if link in {"STAGE1_INPUT_OUTPUT_SHA256.csv", "results/AUDIT_COMPLETION_QA.json"}:
                continue  # Generated by this script immediately below.
            if not (ROOT / link).is_file():
                missing_links.append({"file": filename, "link": link})
    check("report_local_links", not missing_links, missing_links)
    check("batch_allocations_have_full571", pd.read_csv(RESULTS / "CROSS_PLATFORM_BATCH_ALLOCATION_SUMMARY.csv").n.eq(571).all(), "All three metadata-only pairwise allocations contain571")
    passed = all(r["passed"] for r in checks)
    qa = {"status": "STAGE1_PASS_WITH_OPEN_DECISIONS" if passed else "DELIVERABLE_QA_FAILED",
          "completed_at_local": datetime.now().isoformat(), "checks_passed": sum(r["passed"] for r in checks), "checks_total": len(checks),
          "clinical_files_opened": 0, "final_Strict_core_frozen": False,
          "open_items": ["TMT3447/3907 source-version explanation", "LOD/left-censoring and final eligibility", "GIS/batch rule",
                         "multiaptamer rule", "power/pathology coverage", "measurement model identification"],
          "manifest_scope": "8 source files plus all audit deliverables, restricted files included; manifest itself and __pycache__ excluded",
          "checks": checks}
    write_json(RESULTS / "AUDIT_COMPLETION_QA.json", qa)
    if not passed:
        raise RuntimeError("Deliverable QA failed; inspect AUDIT_COMPLETION_QA.json")
    (ROOT / "logs/FINALIZATION.log").write_text(datetime.now().isoformat()+f" Completion QA {len(checks)}/{len(checks)} passed. Manifest generated after this log.\n", encoding="utf-8")
    manifest_rows = []
    for row in input_checks:
        manifest_rows.append({"type": "SOURCE_I01_I08", "input_id": row["input_id"], "path": row["path"],
                              "bytes": row["bytes"], "sha256": row["sha256"], "restricted_participant_detail": False})
    manifest_path = ROOT / "STAGE1_INPUT_OUTPUT_SHA256.csv"
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or path == manifest_path or "__pycache__" in path.parts:
            continue
        relative = path.relative_to(ROOT)
        manifest_rows.append({"type": "AUDIT_OUTPUT", "input_id": "", "path": relative.as_posix(),
                              "bytes": path.stat().st_size, "sha256": digest(path),
                              "restricted_participant_detail": relative.parts[0] == "00_restricted"})
    write_csv(manifest_path, manifest_rows)
    # Verify serialized manifest hashes after writing; no audited artifact is mutated afterwards.
    with manifest_path.open(encoding="utf-8-sig", newline="") as f:
        serialized = list(csv.DictReader(f))
    for row in serialized:
        path = Path(row["path"]) if row["type"] == "SOURCE_I01_I08" else ROOT / row["path"]
        if path.stat().st_size != int(row["bytes"]) or digest(path) != row["sha256"]:
            raise RuntimeError(f"Serialized manifest verification failed: {path.name}")
    print(json.dumps({"status": qa["status"], "completion_checks": f"{len(checks)}/{len(checks)}",
                      "manifest_rows_verified": len(serialized), "workbook_sheets": len(schema)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
