from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from audit_helpers import *


def visit(value):
    value = clean_text(value).lower()
    return {"baseline": "bl", "m00": "bl"}.get(value, value)


def check_ptid(ptid, rid):
    ptid = clean_text(ptid)
    if not ptid:
        return "NOT_AVAILABLE"
    match = re.fullmatch(r"\d{3}_S_(\d{4,6})", ptid)
    if not match:
        return "UNPARSED"
    return "MATCH" if int(match.group(1)) == int(rid) else "MISMATCH"


def build_sample_sets(n, s, t):
    log("Build exact RID/date pairings and cross-check PTID/VISCODE2")
    platforms = {"NULISA": n, "SOMASCAN": s, "TMT_MS": t}
    keys = {name: set(v["meta"].index) for name, v in platforms.items()}
    pairs = []
    for size in [2, 3]:
        for names in itertools.combinations(keys, size):
            common = set.intersection(*(keys[name] for name in names))
            common_rid = set.intersection(*({k.split("|")[0] for k in keys[name]} for name in names))
            matches, mismatches, unknown = 0, 0, 0
            for key in common:
                codes = []
                for name in names:
                    col = "VISCODE2"
                    codes.append(visit(platforms[name]["meta"].loc[key, col]))
                if any(not x for x in codes):
                    unknown += 1
                elif len(set(codes)) == 1:
                    matches += 1
                else:
                    mismatches += 1
            pairs.append({"comparison": " + ".join(names), "common_RID_any_date": len(common_rid),
                          "same_date_keys": len(common), "same_date_participants": len({k.split('|')[0] for k in common}),
                          "same_date_VISCODE2_match": matches, "same_date_VISCODE2_mismatch": mismatches,
                          "same_date_VISCODE2_unknown": unknown})
    triple_all = sorted(set.intersection(*keys.values()), key=lambda k: (int(k.split("|")[0]), k.split("|")[1]))
    triple = []
    seen_rids = set()
    for key in triple_all:
        rid = key.split("|")[0]
        if rid not in seen_rids:
            triple.append(key)
            seen_rids.add(rid)
    rows = []
    for key in triple:
        rid, date = key.split("|")
        nr, sr, tr = n["meta"].loc[key], s["meta"].loc[key], t["meta"].loc[key]
        nv, sv, tv = [visit(x) for x in [nr["VISCODE2"], sr["VISCODE2"], tr["VISCODE2"]]]
        codes = [nv, sv, tv]
        pstatus = "UNKNOWN" if not nr["PTID"] or not tr["PTID"] else "MATCH" if nr["PTID"] == tr["PTID"] else "MISMATCH"
        row = {"key": key, "RID": int(rid), "EXAMDATE": date, "NULISA_PTID": nr["PTID"], "TMT_PTID": tr["PTID"],
               "SOMASCAN_PTID_from_metadata": sr["PTID_from_available_metadata"], "NULISA_TMT_PTID_check": pstatus,
               "NULISA_PTID_RID_check": check_ptid(nr["PTID"], rid), "TMT_PTID_RID_check": check_ptid(tr["PTID"], rid),
               "SOMASCAN_PTID_RID_check": check_ptid(sr["PTID_from_available_metadata"], rid),
               "NULISA_VISCODE2": nv, "SOMASCAN_VISCODE2": sv, "TMT_VISCODE2": tv,
               "VISCODE2_status": "UNKNOWN" if any(not x for x in codes) else "MATCH" if len(set(codes)) == 1 else "MISMATCH",
               "NULISA_SampleName": nr["SampleName"], "NULISA_PlateID": nr["PlateID"], "SOMASCAN_PlateId": sr["PlateId"],
               "SOMASCAN_GUSPECID": sr.get("GUSPECID", ""), "SOMASCAN_ExtIdentifier": sr.get("ExtIdentifier", ""),
               "TMT_ProteomicID": tr["ProteomicID"], "TMT_plex": tr["plex"], "TMT_channel": tr["channel"],
               "exact_aliquot_identity_confirmed": False}
        rows.append(row)
    registry = pd.DataFrame(rows)
    write_csv(RESTRICTED / "TRIPLE_SAMPLE_REGISTRY_RESTRICTED.csv", registry)
    write_csv(RESULTS / "PAIRWISE_SAMPLE_COUNTS.csv", pairs)
    identity_rows = []
    for name, data in platforms.items():
        for key, row in data["meta"].iterrows():
            rid, date = key.split("|")
            ptid = row.get("PTID", row.get("PTID_from_available_metadata", ""))
            identity_rows.append({"platform": name, "key": key, "RID": int(rid), "EXAMDATE": date,
                                  "PTID_available_or_inferred": ptid, "PTID_RID_check": check_ptid(ptid, rid)})
    identity = pd.DataFrame(identity_rows)
    write_csv(RESTRICTED / "ALL_PLATFORM_SAMPLE_ID_CHECKS_RESTRICTED.csv", identity)
    identity_summary = identity.groupby(["platform", "PTID_RID_check"]).size().reset_index(name="n_records")
    write_csv(RESULTS / "SAMPLE_ID_CHECK_SUMMARY.csv", identity_summary)
    summary = {"triple_all_same_date_keys": len(triple_all), "triple_earliest_per_RID": len(triple),
               "multiple_common_dates_removed_by_prespecified_rule": len(triple_all) - len(triple),
               "triple_VISCODE2_counts": registry["VISCODE2_status"].value_counts().to_dict(),
               "triple_NULISA_TMT_PTID_counts": registry["NULISA_TMT_PTID_check"].value_counts().to_dict(),
               "all_platform_PTID_RID_mismatches": int(identity["PTID_RID_check"].eq("MISMATCH").sum()),
               "exact_aliquot_identity_confirmed": False,
               "NULISA_plates_in_triple": len(set(v for x in registry["NULISA_PlateID"] for v in x.split('|') if v)),
               "SOMASCAN_plates_in_triple": registry["SOMASCAN_PlateId"].nunique(),
               "TMT_plexes_in_triple": registry["TMT_plex"].nunique()}
    write_json(RESULTS / "SAMPLE_MATCHING_SUMMARY.json", summary)
    if summary["all_platform_PTID_RID_mismatches"] or registry["NULISA_TMT_PTID_check"].eq("MISMATCH").any():
        raise RuntimeError("Hard sample ID conflict; see restricted identity checks")
    return triple, registry, pairs, summary


def availability(matrix, features, index, how):
    if not features:
        return pd.Series(False, index=index)
    good = np.isfinite(matrix.reindex(index)[features])
    return good.any(axis=1) if how == "any" else good.all(axis=1)


def run_measurement_audit(n, s, t, mapping, feature_map, triple):
    log("Audit candidate distributions, LOD, per-batch missingness and same-plex GIS support")
    scopes = {"NATIVE_RECORDS": None, "TRIPLE_SAME_DATE": triple}
    feature_rows, batch_rows, protein_rows, gis_rows = [], [], [], []
    pair_rows, multia_rows, restricted_masks = [], [], []
    tmt_features = sorted({f for d in feature_map.values() for f in d["TMT"]})
    sample_gis = t["gis_log"].reindex(t["meta"]["plex"])[tmt_features]
    sample_gis.index = t["meta"].index
    ratio_matrix = t["matrix"][tmt_features] - sample_gis
    data_routes = {"NULISA": n, "SOMASCAN": s, "TMT_MS_RAW": t,
                   "TMT_MS_GIS_RATIO": {"matrix": ratio_matrix, "batch": t["batch"], "meta": t["meta"]}}

    for m in mapping.to_dict("records"):
        protein = m["accession"]
        fm = feature_map[protein]
        if len(fm["NULISA"]) != 1 or len(fm["TMT"]) != 1:
            raise RuntimeError("Candidate has multiple NULISA/TMT routes: explicit aggregation required")
        nt, tf = fm["NULISA"][0], fm["TMT"][0]
        routes = {"NULISA": [nt], "SOMASCAN": fm["SOMASCAN"], "TMT_MS_RAW": [tf], "TMT_MS_GIS_RATIO": [tf]}
        all_stats = {}
        for platform, features in routes.items():
            data = data_routes[platform]
            for scope, wanted in scopes.items():
                index = data["matrix"].index if wanted is None else pd.Index(wanted)
                for feature in features:
                    v = data["matrix"].reindex(index)[feature]
                    b = data["batch"].reindex(index)
                    st = numerical_stats(v)
                    effect = batch_effect(v, b)
                    row = {"accession": protein, "target": nt, "annotation_strict_candidate": m["annotation_strict_candidate"],
                           "platform": platform, "feature": feature, "scope": scope, "n_records": len(index),
                           **st, "n_missing": len(index) - st["n_valid"], "missing_fraction": 1 - st["n_valid"] / len(index),
                           "scale": "released_NPQ_log2_like" if platform == "NULISA" else "log2_RFU" if platform == "SOMASCAN" else "log2_sample_over_GIS" if platform.endswith("RATIO") else "log2_positive_intensity",
                           **effect}
                    if platform == "NULISA":
                        nd = n["any_detect"].reindex(index)[nt]
                        ad = n["all_detect"].reindex(index)[nt]
                        row.update({"n_individual_detection_known": int(nd.notna().sum()), "n_detect_any_repeat": int(nd.eq(1).sum()),
                                    "n_detect_all_repeats": int(ad.eq(1).sum()), "detect_any_fraction_all_records": float(nd.eq(1).mean()),
                                    "detect_all_fraction_all_records": float(ad.eq(1).mean()),
                                    "n_with_npq_equal_lod_records": int(n["aggregated"]["n_equal_lod"].unstack("target").reindex(index)[nt].sum()),
                                    "detection_basis": "NPQ>=matching_record_plate_LOD; release TargetDetectability is plate aggregate"})
                    else:
                        row["detection_basis"] = "Individual LOD not provided; numerical availability only"
                    if platform == "SOMASCAN":
                        annotation = s["annotation"].loc[s["annotation"]["Analytes"].eq(feature)].iloc[0]
                        row.update({"feature_raw_accession_field": annotation["UniProt"],
                                    "feature_is_multiaccession": annotation["n_accessions"] > 1})
                    if platform == "TMT_MS_GIS_RATIO":
                        raw_paired = t["matrix"].reindex(index)[tf].where(v.notna())
                        before = batch_effect(raw_paired, b)
                        row.update({"batch_r2_before_GIS_on_identical_available_subset": before["batch_r2"],
                                    "batch_r2_after_GIS": effect["batch_r2"],
                                    "delta_batch_r2_GIS_minus_raw_same_subset": effect["batch_r2"] - before["batch_r2"],
                                    "GIS_ratio_is_audit_sensitivity_not_frozen_processing": True})
                    feature_rows.append(row)
                    all_stats[(platform, scope, feature)] = row
                    batch_rows.extend(per_batch_rows(v, b, {"accession": protein, "target": nt, "platform": platform,
                                                           "feature": feature, "scope": scope}))

        ix = pd.Index(triple)
        nv = n["matrix"].reindex(ix)[nt].notna()
        nd_any = n["any_detect"].reindex(ix)[nt].eq(1)
        nd_all = n["all_detect"].reindex(ix)[nt].eq(1)
        sa = availability(s["matrix"], fm["SOMASCAN"], ix, "any")
        sall = availability(s["matrix"], fm["SOMASCAN"], ix, "all")
        su_any = availability(s["matrix"], fm["SOMASCAN_SINGLE_ACCESSION"], ix, "any")
        su_all = availability(s["matrix"], fm["SOMASCAN_SINGLE_ACCESSION"], ix, "all")
        tv = t["matrix"].reindex(ix)[tf].notna()
        gv = sample_gis.reindex(ix)[tf].notna()
        tg = tv & gv
        masks = {"N_numeric": nv, "N_detect_any_repeat": nd_any, "N_detect_all_repeats": nd_all,
                 "S_any": sa, "S_all": sall, "S_unambiguous_any": su_any, "S_unambiguous_all": su_all,
                 "T_raw": tv, "GIS_positive": gv, "T_GIS_usable": tg,
                 "complete_raw_S_any": nv & sa & tv, "complete_raw_S_all": nv & sall & tv,
                 "complete_GIS_S_any": nv & sa & tg, "complete_GIS_S_all": nv & sall & tg,
                 "complete_GIS_unambiguous_S_any": nv & su_any & tg,
                 "complete_GIS_unambiguous_S_all": nv & su_all & tg,
                 "complete_GIS_detect_all_N_S_any": nv & nd_all & sa & tg,
                 "complete_GIS_detect_all_N_S_all": nv & nd_all & sall & tg}
        pr = {**m, "n_triple_registry": len(ix)}
        for name, mask in masks.items():
            pr["n_" + name] = int(mask.sum())
            pr["fraction_" + name] = float(mask.mean())
        ns = all_stats[("NULISA", "TRIPLE_SAME_DATE", nt)]
        ts = all_stats[("TMT_MS_RAW", "TRIPLE_SAME_DATE", tf)]
        gs = all_stats[("TMT_MS_GIS_RATIO", "TRIPLE_SAME_DATE", tf)]
        pr.update({"N_NPQ_sd": ns["sd"], "N_NPQ_iqr": ns["iqr"], "N_batch_r2": ns["batch_r2"],
                   "T_log2_sd": ts["sd"], "T_GIS_log2_sd": gs["sd"], "T_raw_batch_r2": ts["batch_r2"],
                   "T_raw_batch_r2_on_GIS_available": gs["batch_r2_before_GIS_on_identical_available_subset"],
                   "T_GIS_batch_r2": gs["batch_r2"], "T_GIS_delta_batch_r2_matched": gs["delta_batch_r2_GIS_minus_raw_same_subset"],
                   "N_detection_fraction_all_repeats": ns["detect_all_fraction_all_records"],
                   "N_n_detection_known": ns["n_individual_detection_known"]})
        soma_stats = [all_stats[("SOMASCAN", "TRIPLE_SAME_DATE", f)] for f in fm["SOMASCAN"]]
        pr["S_log2_sd_min"] = min((x["sd"] for x in soma_stats if np.isfinite(x["sd"])), default=np.nan)
        pr["S_log2_sd_max"] = max((x["sd"] for x in soma_stats if np.isfinite(x["sd"])), default=np.nan)
        pr["S_batch_r2_min"] = min((x["batch_r2"] for x in soma_stats if np.isfinite(x["batch_r2"])), default=np.nan)
        pr["S_batch_r2_max"] = max((x["batch_r2"] for x in soma_stats if np.isfinite(x["batch_r2"])), default=np.nan)
        flags = []
        if not m["annotation_strict_candidate"]:
            flags.append("MAPPING_OR_MULTIAPTAMER_DECISION_REQUIRED")
        if m["soma_has_multiaccession_annotation"]:
            flags.append("COMPOSITE_SOMA_REAGENT_PRESENT")
        if ns["zero_variance"] or ts["zero_variance"] or any(x["zero_variance"] for x in soma_stats):
            flags.append("ZERO_VARIANCE_IN_AT_LEAST_ONE_FEATURE")
        if ns["n_detect_all_repeats"] == 0:
            flags.append("NO_NULISA_LOD_DETECTED_PARTICIPANTS")
        elif ns["detect_all_fraction_all_records"] < .5:
            flags.append("NULISA_DETECTION_BELOW_HALF_DESCRIPTIVE")
        if pr["n_complete_GIS_S_all"] == 0:
            flags.append("NO_COMPLETE_ALL_FEATURE_GIS_OBSERVATIONS")
        if pr["n_complete_GIS_S_all"] < 100:
            flags.append("COMPLETE_ALL_FEATURE_GIS_N_LT_100_DESCRIPTIVE")
        pr["audit_flags"] = "|".join(flags) if flags else "NO_HARD_STRUCTURE_FLAG_IDENTIFIED"
        pr["final_evaluability_status"] = "FROZEN_SAP_POWER_PATHOLOGY_COVERAGE_PENDING"
        protein_rows.append(pr)
        part = pd.DataFrame({"key": ix, "accession": protein})
        for name, mask in masks.items():
            part[name] = mask.to_numpy(dtype=bool)
        restricted_masks.append(part)

        gis_raw = t["gis_raw"][tf]
        gst = numerical_stats(t["gis_log"][tf])
        positive_gis = gis_raw.dropna()
        pr["GIS_positive_plexes"] = len(positive_gis)
        pr["GIS_missing_plexes"] = len(gis_raw) - len(positive_gis)
        pr["GIS_across_plex_raw_CV"] = float(positive_gis.std(ddof=1) / positive_gis.mean()) if len(positive_gis) > 1 and positive_gis.mean() > 0 else np.nan
        pr["GIS_across_plex_log2_sd"] = gst["sd"]
        for plex in t["gis_raw"].index:
            native_keys = t["meta"].index[t["meta"]["plex"].eq(plex)]
            triple_keys = [k for k in native_keys if k in set(triple)]
            valid_native = t["matrix"].reindex(native_keys)[tf].notna()
            valid_triple = t["matrix"].reindex(triple_keys)[tf].notna()
            gis_positive = bool(pd.notna(gis_raw.loc[plex]))
            gis_rows.append({"accession": protein, "target": nt, "tmt_feature": tf, "plex": plex,
                             "n_native_participants": len(native_keys), "n_triple_participants": len(triple_keys),
                             "GIS_positive": gis_positive, "GIS_scaled_intensity": gis_raw.loc[plex],
                             "GIS_log2_intensity": t["gis_log"].loc[plex, tf],
                             "n_native_raw_positive": int(valid_native.sum()), "n_triple_raw_positive": int(valid_triple.sum()),
                             "n_native_GIS_usable": int(valid_native.sum()) if gis_positive else 0,
                             "n_triple_GIS_usable": int(valid_triple.sum()) if gis_positive else 0,
                             "all_native_participants_missing": not bool(valid_native.any()),
                             "triple_mean_raw_log2": numerical_stats(t["matrix"].reindex(triple_keys)[tf])["mean"],
                             "triple_sd_raw_log2": numerical_stats(t["matrix"].reindex(triple_keys)[tf])["sd"],
                             "triple_mean_GIS_log2_ratio": numerical_stats(ratio_matrix.reindex(triple_keys)[tf])["mean"]})

        protein_pair_rows = []
        for f1, f2 in itertools.combinations(fm["SOMASCAN"], 2):
            for scope, wanted in scopes.items():
                sf = s["matrix"] if wanted is None else s["matrix"].reindex(wanted)
                a, b = sf[f1], sf[f2]
                valid = a.notna() & b.notna()
                nvalid = int(valid.sum())
                rho = float(spearmanr(a[valid], b[valid]).statistic) if nvalid >= 3 and a[valid].nunique() > 1 and b[valid].nunique() > 1 else np.nan
                pa = {"accession": protein, "target": nt, "feature_1": f1, "feature_2": f2,
                      "scope": scope, "n_pairwise_complete": nvalid, "spearman_rho": rho,
                      "selection_used_for_primary_aptamer": False}
                pair_rows.append(pa)
                protein_pair_rows.append(pa)
        for scope in scopes:
            own = [x for x in protein_pair_rows if x["scope"] == scope]
            rhos = [x["spearman_rho"] for x in own if np.isfinite(x["spearman_rho"])]
            ns_pair = [x["n_pairwise_complete"] for x in own]
            multia_rows.append({"accession": protein, "target": nt, "scope": scope, "n_features": len(fm["SOMASCAN"]),
                               "n_single_accession_features": len(fm["SOMASCAN_SINGLE_ACCESSION"]), "n_pairs": len(own),
                               "n_evaluable_pairs": len(rhos), "rho_min": min(rhos, default=np.nan),
                               "rho_median": float(np.median(rhos)) if rhos else np.nan, "rho_max": max(rhos, default=np.nan),
                               "n_pairwise_complete_min": min(ns_pair, default=np.nan),
                               "n_pairwise_complete_max": max(ns_pair, default=np.nan),
                               "any_nonpositive_pair_rho": any(x <= 0 for x in rhos),
                               "aggregation_not_required": len(fm["SOMASCAN"]) == 1,
                               "multiaptamer_aggregation_rule_frozen": False,
                               "note": "within-platform descriptive QC only; no cross-platform correlation selection"})

    feature_df, protein_df = pd.DataFrame(feature_rows), pd.DataFrame(protein_rows)
    batch_df, gis_df = pd.DataFrame(batch_rows), pd.DataFrame(gis_rows)
    multi_df, pair_df = pd.DataFrame(multia_rows), pd.DataFrame(pair_rows)
    write_csv(RESULTS / "FEATURE_LEVEL_QC.csv", feature_df)
    write_csv(RESULTS / "PROTEIN_PLATFORM_EVALUABILITY.csv", feature_df)
    write_csv(RESULTS / "PROTEIN_EVALUABILITY_70.csv", protein_df)
    write_csv(RESULTS / "PROTEIN_EVALUABILITY_SINGLE_FEATURE_DRAFT.csv", protein_df.loc[protein_df["annotation_strict_candidate"]])
    write_csv(RESULTS / "PLATFORM_BATCH_QC_LONG.csv", batch_df)
    write_csv(RESULTS / "TMT_GIS_BATCH_AUDIT.csv", gis_df)
    write_csv(RESULTS / "MULTIAPTAMER_AUDIT.csv", multi_df)
    write_csv(RESULTS / "MULTIAPTAMER_PAIRWISE_QC.csv", pair_df)
    pd.concat(restricted_masks, ignore_index=True).to_parquet(RESTRICTED / "TRIPLE_CANDIDATE_AVAILABILITY_RESTRICTED.parquet", index=False)
    write_csv(RESULTS / "TMT_GIS_SUMMARY.csv", protein_df[["accession", "nulisa_target", "GIS_positive_plexes", "GIS_missing_plexes",
                                                       "GIS_across_plex_raw_CV", "GIS_across_plex_log2_sd", "n_T_raw", "n_T_GIS_usable",
                                                       "T_raw_batch_r2_on_GIS_available", "T_GIS_batch_r2", "T_GIS_delta_batch_r2_matched"]])
    threshold_rows = []
    metrics = ["n_complete_raw_S_any", "n_complete_raw_S_all", "n_complete_GIS_S_any", "n_complete_GIS_S_all",
               "n_complete_GIS_detect_all_N_S_any", "n_complete_GIS_detect_all_N_S_all"]
    for group, sub in [("EXTENDED_ANNOTATION_CANDIDATES", protein_df),
                       ("SINGLE_FEATURE_ANNOTATION_CANDIDATES", protein_df.loc[protein_df["annotation_strict_candidate"]])]:
        for metric in metrics:
            for fraction in [.5, .7, .8, .9, .95, 1.0]:
                threshold_rows.append({"candidate_set": group, "n_candidates": len(sub), "metric": metric,
                                       "threshold_type": "fraction_of_triple_registry", "threshold": fraction,
                                       "n_proteins_at_or_above": int(sub[metric].ge(len(triple)*fraction - 1e-10).sum()),
                                       "threshold_is_frozen_inclusion_rule": False})
            for number in [100, 200, 300, 400, 500]:
                threshold_rows.append({"candidate_set": group, "n_candidates": len(sub), "metric": metric,
                                       "threshold_type": "participant_count", "threshold": number,
                                       "n_proteins_at_or_above": int(sub[metric].ge(number).sum()),
                                       "threshold_is_frozen_inclusion_rule": False})
    write_csv(RESULTS / "DESCRIPTIVE_COVERAGE_THRESHOLDS.csv", threshold_rows)
    return {"proteins": protein_df, "features": feature_df, "batches": batch_df, "gis": gis_df,
            "multiaptamer": multi_df, "multiaptamer_pairs": pair_df, "thresholds": pd.DataFrame(threshold_rows)}
