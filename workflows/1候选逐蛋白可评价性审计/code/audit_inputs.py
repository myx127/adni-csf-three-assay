from __future__ import annotations

import csv
import gc
import io
import logging
import re
import zipfile

import numpy as np
import pandas as pd
from pypdf import PdfReader

from audit_helpers import *


def verify_inputs():
    manifest_path = HANDOFF / "01_数据清单/CORE_INPUT_MANIFEST.csv"
    with manifest_path.open(encoding="utf-8-sig", newline="") as f:
        all_entries = list(csv.DictReader(f))
    entries = {r["input_id"]: r for r in all_entries if int(r["input_id"][1:]) <= 8}
    checks = []
    for key, entry in entries.items():
        path = Path(entry["path"])
        present = path.is_file()
        observed = digest(path) if present else ""
        checks.append({"input_id": key, "role": entry["role"], "path": str(path),
                       "expected_bytes": int(entry["bytes"]), "observed_bytes": path.stat().st_size if present else None,
                       "expected_sha256": entry["sha256"], "observed_sha256": observed,
                       "passed": bool(present and path.stat().st_size == int(entry["bytes"]) and observed == entry["sha256"])})
    write_csv(RESULTS / "INPUT_HASH_VERIFICATION.csv", checks)
    if not all(row["passed"] for row in checks):
        raise RuntimeError("I01-I08 identity verification failed; no quantitative audit run")
    write_json(RESULTS / "ACCESS_BOUNDARY_AUDIT.json", {
        "allowed_input_ids": list(entries), "opened_clinical_input_ids": [],
        "excluded_input_ids": [r["input_id"] for r in all_entries if int(r["input_id"][1:]) > 8],
        "clinical_outcome_linkage": False, "cross_platform_consensus_or_residual_estimation": False,
        "protocol_sha256_before_analysis": digest(ROOT / "AUDIT_PROTOCOL_v1.0.md"),
        "clinical_files_opened": 0})
    log("I01-I08 input hashes and byte counts all match")
    return entries, checks


def inspect_methods(entries):
    logging.getLogger("pypdf").setLevel(logging.ERROR)
    output = ["# 方法文件核对记录", "", "只核对检测尺度、LOD、GIS、缺失编码与列数说明。", ""]
    facts = {}
    for key, limit, terms in [("I02", 5, ["NPQ", "limit of detection", "detectability"]),
                              ("I08", 6, ["3,447", "not batch corrected", "missing values", "global internal"])]:
        reader = PdfReader(entries[key]["path"])
        text_pages = []
        for page in reader.pages[:limit]:
            txt = " ".join((page.extract_text() or "").split())
            text_pages.append(txt)
            for term in terms:
                match = re.search(re.escape(term), txt, re.I)
                if match:
                    excerpt = txt[max(0, match.start()-90):match.end()+430]
                    output += [f"## {key} PDF第{page.page_number + 1}页：{term}", "", excerpt, ""]
        facts[key] = {"pages": len(reader.pages), "terms_found": {term: any(term.lower() in x.lower() for x in text_pages) for term in terms}}
    (RESULTS / "METHODS_EVIDENCE_EXCERPTS.md").write_text("\n".join(output), encoding="utf-8")
    dictionary = pd.read_csv(entries["I03"]["path"], dtype=str, keep_default_na=False)
    subset = dictionary.loc[dictionary["TBLNAME"].eq("BSHRI_PLA_CSF_NULISA_CNS")]
    write_csv(RESULTS / "NULISA_DICTIONARY_FIELDS.csv", subset)
    facts["I03"] = {"dictionary_rows": len(dictionary), "nulisa_dictionary_rows": len(subset)}
    return facts


def load_nulisa(entry):
    log("Read NULISA released table and audit exact duplicates")
    raw = pd.read_csv(entry["path"], dtype=str, keep_default_na=False)
    original_columns = list(raw.columns)
    raw["source_csv_row"] = raw.index + 2
    duplicated = raw.duplicated(subset=original_columns)
    counts = {"raw_rows": len(raw), "raw_columns": len(original_columns),
              "matrix_type_rows": raw["SampleMatrixType"].value_counts().to_dict(),
              "raw_target_labels": raw["Target"].nunique(),
              "canonical_target_labels": raw["Target"].map(canonical_target).nunique(),
              "update_stamps": sorted(raw["update_stamp"].unique().tolist())}
    for fluid in ["PLASMA", "CSF"]:
        in_fluid = raw["SampleMatrixType"].eq(fluid)
        has_rid = pd.to_numeric(raw["RID"], errors="coerce").notna()
        frame = raw.loc[in_fluid & has_rid]
        dupe_participants = frame.loc[frame.duplicated(subset=original_columns, keep=False), "RID"].nunique()
        counts[fluid] = {"unique_rid": frame["RID"].nunique(),
                         "qc_passed_rid": frame.loc[frame["SampleQC"].eq("passed"), "RID"].nunique(),
                         "qc_warning_rid": frame.loc[frame["SampleQC"].eq("warning"), "RID"].nunique(),
                         "exact_rows_removed": int((in_fluid & has_rid & duplicated).sum()),
                         "participants_with_exact_duplicate_rows": dupe_participants}
    annotation = raw[["Target", "UniProtID", "ProteinName"]].drop_duplicates().copy()
    annotation["canonical_target"] = annotation["Target"].map(canonical_target)
    annotation["accession_tokens"] = annotation["UniProtID"].map(tokens)
    annotation["accession"] = annotation["accession_tokens"].map(lambda x: next(iter(x)) if len(x) == 1 else "")
    annotation["target_class"] = annotation["canonical_target"].map(
        lambda x: "SPECIAL_FORM_OR_COMPOSITE_EXCLUDED" if x in SPECIAL_TARGETS else
        "CSF_NO_NPQ_EXCLUDED" if x in CSF_UNAVAILABLE else "TOTAL_LIKE_TARGET_ANNOTATION_ONLY")
    annotation_out = annotation.drop(columns=["accession_tokens"]).copy()
    write_csv(RESULTS / "NULISA_TARGET_CLASSIFICATION.csv", annotation_out)

    csf_all = raw.loc[raw["SampleMatrixType"].eq("CSF") & raw["SampleType"].eq("Sample") & ~duplicated].copy()
    no_rid = pd.to_numeric(csf_all["RID"], errors="coerce").isna()
    counts["CSF"]["nonrid_sample_name_plate_blocks"] = len(csf_all.loc[no_rid, ["SampleName", "PlateID"]].drop_duplicates())
    counts["CSF"]["nonrid_sample_names"] = csf_all.loc[no_rid, "SampleName"].nunique()
    csf = csf_all.loc[~no_rid & csf_all["SampleQC"].eq("passed")].copy()
    csf["key"], key_quality = make_keys(csf, "NULISA_CSF_QC_PASSED_LONG")
    csf = csf.loc[csf["key"].ne("")].copy()
    csf["target"] = csf["Target"].map(canonical_target)
    csf["npq"] = pd.to_numeric(csf["NPQ"], errors="coerce").replace([np.inf, -np.inf], np.nan)
    csf["lod"] = pd.to_numeric(csf["TargetLOD_NPQ"], errors="coerce").replace([np.inf, -np.inf], np.nan)
    csf["detectability_fraction"] = pd.to_numeric(csf["TargetDetectability"].str.replace("%", "", regex=False), errors="coerce") / 100
    known = csf["npq"].notna() & csf["lod"].notna()
    csf["detection"] = np.where(known, csf["npq"].ge(csf["lod"]).astype(float), np.nan)
    csf["n_known_lod_and_npq"] = known.astype(int)
    csf["npq_equals_lod"] = (known & csf["npq"].eq(csf["lod"])).astype(int)
    duplicate_rows = csf.loc[csf.duplicated(["key", "target"], keep=False)].copy()
    write_csv(RESTRICTED / "NULISA_NONIDENTICAL_TECHNICAL_REPEATS_RESTRICTED.csv", duplicate_rows)
    counts["CSF"]["nonidentical_repeat_participants"] = duplicate_rows["RID"].nunique()
    counts["CSF"]["nonidentical_repeat_rows"] = len(duplicate_rows)
    counts["CSF"]["nonidentical_repeat_key_target_groups"] = duplicate_rows.groupby(["key", "target"]).ngroups
    counts["CSF"]["plate_count"] = csf["PlateID"].nunique()
    counts["CSF"]["key_quality"] = key_quality

    grouped = csf.groupby(["key", "target"], sort=False)
    agg = grouped.agg(npq=("npq", "mean"), npq_first=("npq", "first"), npq_last=("npq", "last"),
                      npq_min=("npq", "min"), npq_max=("npq", "max"), n_replicate_records=("npq", "size"),
                      n_numeric_npq=("npq", "count"), n_known=("n_known_lod_and_npq", "sum"),
                      detect_any=("detection", "max"), detect_min=("detection", "min"),
                      n_equal_lod=("npq_equals_lod", "sum"), lod_mean=("lod", "mean"))
    agg["detect_all"] = np.where(agg["n_known"].eq(agg["n_replicate_records"]), agg["detect_min"], np.nan)
    agg["replicate_npq_range"] = agg["npq_max"] - agg["npq_min"]
    matrix = agg["npq"].unstack("target")
    any_detect = agg["detect_any"].unstack("target").reindex(matrix.index)
    all_detect = agg["detect_all"].unstack("target").reindex(matrix.index)
    lod_known = agg["n_known"].gt(0).unstack("target").reindex(matrix.index)
    meta_cols = ["key", "PTID", "RID", "EXAMDATE", "VISCODE2", "VISCODE", "SampleName", "PlateID", "PHASE"]
    meta = csf[meta_cols].drop_duplicates().groupby("key", sort=False).agg({c: join_unique for c in meta_cols if c != "key"})
    meta = meta.reindex(matrix.index)
    # Never assign a cross-plate technical repeat arbitrarily to one plate.
    batch = meta["PlateID"].where(~meta["PlateID"].str.contains("|", regex=False), "")
    release = csf.groupby(["target", "PlateID"], sort=False).agg(
        n_records=("npq", "size"), n_numeric_npq=("npq", "count"), n_known=("n_known_lod_and_npq", "sum"),
        n_above_or_equal_lod=("detection", "sum"), release_detectability=("detectability_fraction", "first"),
        release_detectability_unique=("TargetDetectability", "nunique"), lod_npq=("lod", "first"),
        lod_unique=("TargetLOD_NPQ", "nunique")).reset_index()
    release["recomputed_fraction_in_QC_passed_RID_records"] = release["n_above_or_equal_lod"] / release["n_known"].replace(0, np.nan)
    release["note"] = "Release percentage denominator may include other sample/QC categories; not an individual detection flag"
    write_csv(RESULTS / "NULISA_PLATE_TARGET_DETECTION.csv", release)
    write_csv(RESTRICTED / "NULISA_TECHNICAL_REPEAT_AGGREGATION_RESTRICTED.csv", agg.loc[agg["n_replicate_records"].gt(1)].reset_index())
    write_json(RESULTS / "NULISA_STRUCTURAL_AUDIT.json", counts)
    del raw, csf_all
    gc.collect()
    return {"summary": counts, "annotation": annotation, "long": csf, "aggregated": agg,
            "matrix": matrix, "any_detect": any_detect, "all_detect": all_detect,
            "lod_known": lod_known, "meta": meta, "batch": batch, "release": release}


def load_soma(entry, annotation_entry):
    log("Read SomaScan release, metadata and analyte annotation")
    frame = pd.read_csv(entry["path"], low_memory=False)
    columns = list(frame.columns)
    features = [c for c in columns if re.fullmatch(r"X\d+\.\d+", c)]
    annotations = pd.read_csv(annotation_entry["path"], dtype=str, keep_default_na=False)
    annotations["tokens"] = annotations["UniProt"].map(tokens)
    annotations["n_accessions"] = annotations["tokens"].map(len)
    keys, quality = make_keys(frame, "SOMASCAN_CSF_RELEASE")
    frame["key"] = keys
    before = len(frame)
    frame = frame.drop_duplicates(subset=columns)
    bad = frame["key"].eq("")
    duplicates = frame.loc[frame.duplicated("key", keep=False) & ~bad]
    if len(duplicates):
        write_csv(RESTRICTED / "SOMASCAN_UNRESOLVED_DUPLICATES_RESTRICTED.csv", duplicates[["key", "RID", "EXAMDATE", "VISCODE2", "PlateId"]])
        raise RuntimeError("Unresolved SomaScan same RID/date duplicate; quantitative audit stopped")
    frame = frame.loc[~bad].set_index("key")
    meta_columns = [c for c in columns if c not in features]
    meta = frame[meta_columns].copy()
    for col in meta_columns:
        meta[col] = meta[col].map(clean_text)
    inferred_ptid = []
    ptid_conflicts = []
    for _, row in meta.iterrows():
        found = set()
        for value in row:
            found.update(re.findall(r"\d{3}_S_\d{4,6}", clean_text(value)))
        inferred_ptid.append(next(iter(found)) if len(found) == 1 else "")
        ptid_conflicts.append(len(found) > 1)
    meta["PTID_from_available_metadata"] = inferred_ptid
    meta["conflicting_PTID_patterns"] = ptid_conflicts
    raw_values = frame[features].apply(pd.to_numeric, errors="coerce")
    invalid = ~np.isfinite(raw_values)
    summary = {"released_rows": before, "rows_after_exact_dedup": len(frame), "features": len(features),
               "unique_rid": int(pd.to_numeric(meta["RID"], errors="coerce").nunique()),
               "plate_count": meta["PlateId"].nunique(), "metadata_columns": meta_columns,
               "visit_record_counts": meta["VISCODE2"].value_counts(dropna=False).to_dict(),
               "metadata_PTID_available_rows": sum(bool(x) for x in inferred_ptid),
               "metadata_PTID_conflict_rows": sum(ptid_conflicts),
               "global_missing_nonnumeric_nonfinite_cells": int(invalid.sum().sum()),
               "global_nonpositive_cells": int(raw_values.le(0).sum().sum()),
               "annotation_rows": len(annotations), "key_quality": quality}
    log_values = np.log2(raw_values.where(raw_values.gt(0)))
    write_json(RESULTS / "SOMASCAN_STRUCTURAL_AUDIT.json", summary)
    return {"summary": summary, "matrix": log_values, "raw_matrix": raw_values,
            "meta": meta, "batch": meta["PlateId"], "annotation": annotations, "features": features}


def load_tmt(entry):
    log("Read TMT archive; distinguish participant channels and pooled standards")
    with zipfile.ZipFile(entry["path"]) as z:
        member = next(n for n in z.namelist() if Path(n).name == "EMORY_CSF_TMT_MS.csv")
        with z.open(member) as f:
            frame = pd.read_csv(f, low_memory=False)
        info = z.getinfo(member)
        with z.open(member) as f:
            inner_hash = hashlib.sha256(f.read()).hexdigest()
    metadata = ["PHASE", "PTID", "ProteomicID", "RID", "VISCODE", "VISCODE2", "EXAMDATE"]
    features = [c for c in frame.columns if c not in metadata]
    parsed = frame["ProteomicID"].astype(str).str.extract(r"^adni_b(?P<plex>\d+)\.(?P<channel>\d+[NC]?)$")
    if parsed.isna().any().any():
        raise RuntimeError("Unparsed TMT ProteomicID; GIS linkage must be resolved")
    parsed["plex"] = parsed["plex"].map(lambda v: f"P{int(v):02d}")
    frame["plex"] = parsed["plex"]
    frame["channel"] = parsed["channel"]
    rid = pd.to_numeric(frame["RID"], errors="coerce")
    real = rid.notna() & rid.gt(0) & rid.eq(np.floor(rid))
    standard = ~real
    if not frame.loc[standard, "channel"].eq("126").all():
        raise RuntimeError("A no-RID TMT row is not reference channel 126")
    if frame.loc[standard, "plex"].duplicated().any():
        raise RuntimeError("Multiple TMT pooled standards per plex")
    values = frame[features].apply(pd.to_numeric, errors="coerce")
    participants = frame.loc[real].copy()
    keys, quality = make_keys(participants, "TMT_MS_CSF_PARTICIPANTS")
    participants["key"] = keys
    if participants["key"].eq("").any() or participants["key"].duplicated().any():
        raise RuntimeError("TMT participant RID/date invalid or nonunique")
    meta = participants[["key"] + metadata + ["plex", "channel"]].set_index("key")
    for col in meta.columns:
        meta[col] = meta[col].map(clean_text)
    raw = values.loc[real].copy()
    raw.index = meta.index
    raw = raw.where(np.isfinite(raw) & raw.gt(0))
    gis = values.loc[standard].copy()
    gis.index = frame.loc[standard, "plex"]
    gis = gis.where(np.isfinite(gis) & gis.gt(0))
    protein_annotations = []
    for feature in features:
        suffix = feature.rsplit("_", 1)[-1]
        protein_annotations.append({"feature": feature, "gene_label": feature.rsplit("_", 1)[0], "raw_suffix": suffix,
                                    "normalized_accession": accession(suffix), "valid_uniprot_suffix": bool(UNIPROT.fullmatch(suffix)),
                                    "isoform_suffix": bool(re.search(r"-\d+$", suffix))})
    protein_annotations = pd.DataFrame(protein_annotations)
    summary = {"inner_csv": member, "inner_bytes": info.file_size, "inner_sha256": inner_hash,
               "total_rows": len(frame), "participant_rows": int(real.sum()), "unique_rid": int(rid.loc[real].nunique()),
               "standard_rows": int(standard.sum()), "plexes": frame["plex"].nunique(), "protein_columns": len(features),
               "metadata_columns": metadata, "participant_visit_counts": meta["VISCODE2"].value_counts().to_dict(),
               "sentinel_minus4_cells_all": int(values.eq(-4).sum().sum()),
               "sentinel_minus4_cells_participants": int(values.loc[real].eq(-4).sum().sum()),
               "other_nonpositive_cells": int(((values <= 0) & values.ne(-4)).sum().sum()),
               "missing_nonnumeric_nonfinite_cells": int((~np.isfinite(values)).sum().sum()),
               "valid_uniprot_suffix_columns": int(protein_annotations["valid_uniprot_suffix"].sum()),
               "unique_valid_uniprot_accessions": protein_annotations.loc[protein_annotations["valid_uniprot_suffix"], "normalized_accession"].nunique(),
               "special_columns": protein_annotations.loc[~protein_annotations["valid_uniprot_suffix"], "feature"].tolist(),
               "key_quality": quality}
    plex_structure = frame.groupby("plex").agg(rows=("ProteomicID", "size"))
    plex_structure["participant_rows"] = participants.groupby("plex").size()
    plex_structure["standard_rows"] = frame.loc[standard].groupby("plex").size()
    write_csv(RESULTS / "TMT_PLEX_STRUCTURE.csv", plex_structure.reset_index())
    write_csv(RESULTS / "TMT_PROTEIN_COLUMN_ANNOTATION.csv", protein_annotations)
    write_json(RESULTS / "TMT_STRUCTURAL_AUDIT.json", summary)
    return {"summary": summary, "raw_matrix": raw, "matrix": np.log2(raw), "gis_raw": gis, "gis_log": np.log2(gis),
            "meta": meta, "batch": meta["plex"], "annotation": protein_annotations, "features": features,
            "plex_structure": plex_structure}


def build_mapping(nulisa, soma, tmt):
    log("Rebuild target/analyte/protein mapping from released annotations")
    nann = nulisa["annotation"].copy()
    eligible = nann.loc[nann["target_class"].eq("TOTAL_LIKE_TARGET_ANNOTATION_ONLY") & nann["accession"].ne("")]
    sanni, tanni = soma["annotation"], tmt["annotation"]
    sset = set().union(*sanni["tokens"])
    tset = set(tanni.loc[tanni["valid_uniprot_suffix"], "normalized_accession"])
    nset = set(eligible["accession"])
    candidates = sorted(nset & sset & tset)
    rows, feature_map = [], {}
    for protein in candidates:
        nmatch = eligible.loc[eligible["accession"].eq(protein)]
        smatch = sanni.loc[sanni["tokens"].map(lambda s: protein in s)]
        tmatch = tanni.loc[tanni["normalized_accession"].eq(protein) & tanni["valid_uniprot_suffix"]]
        ntargets = sorted(set(nmatch["canonical_target"]))
        sfeatures = list(dict.fromkeys(smatch["Analytes"]))
        tfeatures = list(dict.fromkeys(tmatch["feature"]))
        simple_s = smatch.loc[smatch["n_accessions"].eq(1), "Analytes"].tolist()
        unresolved_isoform = bool(nmatch["UniProtID"].str.contains(r"-\d+", regex=True).any()
                                  or smatch["UniProt"].str.contains(r"-\d+", regex=True).any()
                                  or tmatch["isoform_suffix"].any())
        strict = (len(ntargets) == 1 and len(sfeatures) == 1 and len(simple_s) == 1 and len(tfeatures) == 1 and not unresolved_isoform)
        row = {"accession": protein, "nulisa_target": "|".join(ntargets), "nulisa_raw_target_labels": join_unique(nmatch["Target"]),
               "nulisa_raw_accession_fields": join_unique(nmatch["UniProtID"]), "nulisa_target_count": len(ntargets),
               "soma_features": "|".join(sfeatures), "soma_feature_count": len(sfeatures),
               "soma_single_accession_features": "|".join(simple_s), "soma_single_accession_feature_count": len(simple_s),
               "soma_has_multiaccession_annotation": bool(smatch["n_accessions"].gt(1).any()),
               "soma_raw_accession_fields": "; ".join(sorted(set(smatch["UniProt"]))),
               "soma_target_labels": join_unique(smatch["Target"]), "soma_units": join_unique(smatch["Units"]),
               "tmt_columns": "|".join(tfeatures), "tmt_column_count": len(tfeatures),
               "isoform_mapping_requires_review": unresolved_isoform, "annotation_strict_candidate": strict,
               "mapping_status": "ANNOTATION_ONLY_SINGLE_FEATURE_DRAFT" if strict else "EXTENDED_MAPPING_REQUIRES_REVIEW",
               "molecular_identity_validated": False, "final_SAP_membership_frozen": False}
        rows.append(row)
        feature_map[protein] = {"NULISA": ntargets, "SOMASCAN": sfeatures, "SOMASCAN_SINGLE_ACCESSION": simple_s, "TMT": tfeatures}
    mapping = pd.DataFrame(rows)
    write_csv(RESULTS / "EXTENDED_CORE_MAPPING_DRAFT.csv", mapping)
    write_csv(RESULTS / "STRICT_CORE_MAPPING_DRAFT.csv", mapping.loc[mapping["annotation_strict_candidate"]])
    prior = PROJECT / "数据熟悉记录_20261026/CANDIDATE_OVERLAP_NAME_CHECK.csv"
    prior_frame = pd.read_csv(prior, keep_default_na=False)
    write_csv(RESULTS / "CANDIDATE_MAPPING_CROSSCHECK.csv", [
        {"check": "candidate_accession_set_matches_prior", "passed": set(mapping["accession"]) == set(prior_frame["accession"]),
         "current_count": len(mapping), "prior_count": len(prior_frame)},
        {"check": "annotation_single_feature_unambiguous", "passed": int(mapping["annotation_strict_candidate"].sum()) == 43,
         "current_count": int(mapping["annotation_strict_candidate"].sum()), "prior_count": 43}])
    return mapping, feature_map
