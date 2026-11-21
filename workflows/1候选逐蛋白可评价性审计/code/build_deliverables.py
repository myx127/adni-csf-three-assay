from __future__ import annotations

import itertools
import math
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.formatting.rule import ColorScaleRule, DataBarRule
from openpyxl.utils import get_column_letter

from audit_helpers import ROOT, RESULTS, RESTRICTED, FIGURES, write_csv, write_json


def load_json(name):
    import json
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def read(name):
    return pd.read_csv(RESULTS / name)


def md_table(frame):
    def cell(x):
        if pd.isna(x):
            return "—"
        return str(x).replace("|", "/").replace("\n", " ")
    rows = ["| " + " | ".join(map(str, frame.columns)) + " |",
            "| " + " | ".join(["---"] * len(frame.columns)) + " |"]
    rows += ["| " + " | ".join(cell(x) for x in row) + " |" for row in frame.itertuples(index=False, name=None)]
    return "\n".join(rows)


def metadata_batch_audit():
    # Only batch labels are used. No quantitative cross-platform associations.
    registry = pd.read_csv(RESTRICTED / "TRIPLE_SAMPLE_REGISTRY_RESTRICTED.csv", dtype=str, keep_default_na=False)
    cols = {"NULISA": "NULISA_PlateID", "SOMASCAN": "SOMASCAN_PlateId", "TMT_MS": "TMT_plex"}
    labels = {}
    sizes = []
    for platform, column in cols.items():
        label = registry[column].where(~registry[column].str.contains("|", regex=False), "")
        labels[platform] = label
        counts = label.loc[label.ne("")].value_counts().sort_index()
        for batch, n in counts.items():
            sizes.append({"platform": platform, "batch": batch, "n_triple_participants": int(n),
                          "excluded_multibatch_or_unknown_records": int(label.eq("").sum())})
    long_rows, summaries = [], []
    for a, b in itertools.combinations(cols, 2):
        valid = labels[a].ne("") & labels[b].ne("")
        table = pd.crosstab(labels[a][valid], labels[b][valid])
        for aa in table.index:
            for bb in table.columns:
                if table.loc[aa, bb]:
                    long_rows.append({"platform_A": a, "platform_B": b, "batch_A": aa, "batch_B": bb,
                                      "n_participants": int(table.loc[aa, bb])})
        # Directional allocation descriptions; neither is a technical effect estimate.
        summaries.append({"platform_A": a, "platform_B": b, "n": int(table.to_numpy().sum()),
                          "n_A_batches": len(table), "n_B_batches": len(table.columns),
                          "nonempty_batch_pairs": int((table > 0).to_numpy().sum()),
                          "A_to_dominant_B_weighted_share": float(table.max(axis=1).sum()/table.to_numpy().sum()),
                          "B_to_dominant_A_weighted_share": float(table.max(axis=0).sum()/table.to_numpy().sum()),
                          "interpretation": "Batch allocation counts only; no clinical or protein association testing"})
    write_csv(RESULTS / "TRIPLE_BATCH_SAMPLE_SIZES.csv", sizes)
    write_csv(RESULTS / "CROSS_PLATFORM_BATCH_ALLOCATION.csv", long_rows)
    write_csv(RESULTS / "CROSS_PLATFORM_BATCH_ALLOCATION_SUMMARY.csv", summaries)
    return pd.DataFrame(sizes), pd.DataFrame(summaries)


def decisions():
    rows = [
        ("D01", "输入身份", "AUDIT_CLOSED", "I01–I08大小与SHA-256均匹配", "后续阶段继续绑定同一输入，换版本须重审计"),
        ("D02", "同日配对", "AUDIT_CLOSED_WITH_LIMITATION", "571人同RID同日；570访视一致、1缺访视；无已观察到的ID冲突", "同日未证明同一aliquot；SomaScan无可用PTID，不能声称三平台PTID均已验证"),
        ("D03", "候选映射", "PREFREEZE_OPEN", "70扩展；43注释单特征；26多适配体；3涉及复合注释", "43不等于最终Strict core；须锁定分子解释范围与正式集合"),
        ("D04", "多适配体", "PREFREEZE_OPEN", "IL6R、UBB、SQSTM1、PARK7含非正的适配体对相关", "保留逐特征结果，不能按与其他平台/病理的相关性择优；任何聚合均须另定规则"),
        ("D05", "NULISA检出及左删失", "PREFREEZE_OPEN", "TargetDetectability为板汇总；9蛋白检出不足50%；43中4个", "确定LOD以下保留、删失/敏感性策略，不能将NPQ存在视为可靠检出"),
        ("D06", "TMT版本3,447/3,907", "PREFREEZE_OPEN", "方法3,447，绑定矩阵3,907；其中3,904普通UniProt、3特殊列", "来源差异原因未证实。冻结前取得与此矩阵匹配的官方说明/版本依据；勿自行删460列"),
        ("D07", "GIS及批次处理", "PREFREEZE_OPEN", "68蛋白批次R²下降，2上升；可用性总损失25个参与者×蛋白单元", "GIS比值目前是审计敏感性口径；比较加性plex调整等并用模拟评估，不能按结果逐蛋白挑法"),
        ("D08", "技术重复", "AUDIT_CLOSED_WITH_LIMITATION", "CSF仅1位跨板重复参与者；GIS为池标准", "不足以估计个体技术可靠性或识别真实丰度、表位/片段和测量误差"),
        ("D09", "病理/认知可用性和功效", "NOT_IN_STAGE1", "未打开I09–I15；未计A/T组、随访和协变量缺失", "下一阶段按冻结前边界审计表型可用性并模拟功效；覆盖人数不是功效证明"),
        ("D10", "共识/不一致模型识别", "PREFREEZE_OPEN", "本轮没有运行PCA/共识/不一致或临床关联", "延续研究指导：PC1与正交坐标不能自动分离真丰度和病理特异测量效应；先通过人工零效应模拟"),
        ("D11", "TMT基线描述", "PREFREEZE_OPEN", "方法称基线；实际VISCODE2为1,104 bl + 1 m06", "此次按同日主键保留。正式访视范围需预先决定，不擅自把全部1,105标记为基线"),
    ]
    d = pd.DataFrame(rows, columns=["decision_id", "问题", "状态", "本轮证据", "对冻结或解释的影响"])
    write_csv(RESULTS / "OPEN_DECISION_IMPACTS.csv", d)
    return d


def dictionary(proteins, features):
    definitions = {
        "accession": ("用于候选连接的规范化UniProt accession", "文本", "注释"),
        "nulisa_target": ("NULISA规范化靶标名；Aβ编码别名已合并", "文本", "注释"),
        "target": ("对应NULISA规范化靶标名", "文本", "注释"),
        "nulisa_raw_target_labels": ("未规范化的原始靶标标签，竖线分隔", "文本", "注释"),
        "nulisa_raw_accession_fields": ("原始NULISA accession字段", "文本", "注释"),
        "nulisa_target_count": ("此accession匹配的NULISA总量型候选靶标数", "个", "注释"),
        "soma_features": ("所有匹配SomaScan特征；竖线分隔，未择优", "文本", "注释"),
        "soma_feature_count": ("所有匹配适配体特征数", "个", "注释"),
        "soma_single_accession_features": ("仅含一个accession的SomaScan匹配特征", "文本", "注释"),
        "soma_single_accession_feature_count": ("仅含一个accession的匹配特征数", "个", "注释"),
        "soma_has_multiaccession_annotation": ("是否至少一个适配体注释包含多个accession", "布尔", "注释"),
        "soma_raw_accession_fields": ("原始SomaScan accession字段；保留复合注释", "文本", "注释"),
        "soma_target_labels": ("原始SomaScan靶标描述", "文本", "注释"),
        "soma_units": ("SomaScan原始注释单位；审计统计另取log2", "文本", "注释"),
        "tmt_columns": ("匹配TMT原始列名", "文本", "注释"),
        "tmt_column_count": ("匹配TMT普通UniProt列数", "个", "注释"),
        "isoform_mapping_requires_review": ("去除isoform后缀连接是否带来需复核的映射", "布尔", "注释"),
        "annotation_strict_candidate": ("唯一NULISA、一个单accession适配体、唯一普通TMT列；仅注释候选", "布尔", "70候选；True共43"),
        "mapping_status": ("注释映射状态，不是最终可评价性结论", "文本", "注释"),
        "molecular_identity_validated": ("是否独立证明三平台相同分子物种；本轮均False", "布尔", "注释"),
        "final_SAP_membership_frozen": ("是否冻结最终SAP集合；本轮均False", "布尔", "注释"),
        "n_triple_registry": ("三平台精确同RID同日注册表人数", "人", "固定571"),
        "N_NPQ_sd": ("三平台集有效NULISA NPQ样本标准差", "NPQ", "有效NPQ；ddof=1"),
        "N_NPQ_iqr": ("三平台集NULISA NPQ的75%-25%分位差", "NPQ", "有效NPQ"),
        "N_batch_r2": ("NULISA板间均值差异平方和/总平方和", "0–1", "三平台集有效值及单一板标签"),
        "T_log2_sd": ("TMT正强度log2的样本标准差", "log2强度", "三平台集原始正值"),
        "T_GIS_log2_sd": ("TMT同plex log2(样本/GIS)的样本标准差", "log2比值", "三平台集样本及GIS均正"),
        "T_raw_batch_r2": ("全部原始可用子集的TMT plex描述性R²", "0–1", "三平台集原始正值"),
        "T_raw_batch_r2_on_GIS_available": ("限制到GIS同时可用子集后的原始TMT plex R²", "0–1", "与GIS后R²完全相同子集"),
        "T_GIS_batch_r2": ("同plex log2(样本/GIS)的plex R²", "0–1", "三平台集样本及GIS均正"),
        "T_GIS_delta_batch_r2_matched": ("同子集GIS后R²减GIS前R²；负值表示描述性下降", "R²差", "每蛋白匹配子集"),
        "N_detection_fraction_all_repeats": ("所有技术记录均NPQ≥其板LOD的人数/571；LOD未知计不满足", "比例", "固定571；不是仅以LOD已知人数作分母"),
        "N_n_detection_known": ("至少一条技术记录NPQ与LOD同时已知的人数", "人", "三平台集571"),
        "S_log2_sd_min": ("同蛋白所有匹配适配体log2 RFU标准差的最小值", "log2 RFU", "三平台集各适配体有效值"),
        "S_log2_sd_max": ("同蛋白所有匹配适配体log2 RFU标准差的最大值", "log2 RFU", "三平台集各适配体有效值"),
        "S_batch_r2_min": ("同蛋白适配体板描述性R²最小值", "0–1", "三平台集"),
        "S_batch_r2_max": ("同蛋白适配体板描述性R²最大值", "0–1", "三平台集"),
        "audit_flags": ("描述性审计标记；不是自动删除/纳入规则", "文本", "逐蛋白"),
        "final_evaluability_status": ("最终可评价性仍待SAP、功效和病理覆盖决定", "文本", "逐蛋白"),
        "GIS_positive_plexes": ("该蛋白GIS原始强度为正的plex数量", "plex", "全部65 plex"),
        "GIS_missing_plexes": ("65减去GIS正值plex数量", "plex", "全部65 plex"),
        "GIS_across_plex_raw_CV": ("GIS正强度样本SD/均值；不是个体技术可靠性", "无量纲", "GIS正值plex"),
        "GIS_across_plex_log2_sd": ("GIS正强度log2的样本标准差", "log2强度", "GIS正值plex"),
        "platform": ("NULISA/SOMASCAN/TMT_MS_RAW/TMT_MS_GIS_RATIO", "类别", "平台或处理路径"),
        "feature": ("原始特征列名；NULISA用规范化靶标名", "文本", "逐特征"),
        "scope": ("NATIVE_RECORDS=平台原生审计记录；TRIPLE_SAME_DATE=571人", "类别", "Soma原生737记录来自735人"),
        "n_records": ("当前scope总记录数", "记录", "N原生1517/S原生737/T原生1105/三平台571"),
        "n_valid": ("当前特征当前scope有限且有效的数值数", "记录", "当前scope"),
        "n_unique_values": ("有效值中的不同数值数", "个", "有效数值"),
        "mean": ("有效值算术均值；按scale指定尺度", "对应scale", "有效数值"),
        "sd": ("样本标准差；ddof=1", "对应scale", "有效数值；至少2值"),
        "variance": ("样本方差；ddof=1", "对应scale平方", "有效数值；至少2值"),
        "mad_unscaled": ("median(abs(x-median(x)))；未乘1.4826", "对应scale", "有效数值"),
        "iqr": ("第75分位减第25分位", "对应scale", "有效数值；线性插值分位数"),
        "minimum": ("有效值最小值", "对应scale", "有效数值"),
        "maximum": ("有效值最大值", "对应scale", "有效数值"),
        "zero_variance": ("至少2有效值且SD=0", "布尔", "逐特征"),
        "near_zero_variance_sd_le_1e_6": ("至少2有效值且log尺度SD≤1e-6；仅数值稳定性标记", "布尔", "逐特征"),
        "low_distinct_values_le_5": ("有值且不同数值≤5；不自动删除", "布尔", "逐特征"),
        "n_missing": ("n_records-n_valid", "记录", "当前scope"),
        "missing_fraction": ("1-n_valid/n_records", "比例", "当前scope总记录数"),
        "scale": ("发布NPQ、log2 RFU、log2正强度或log2样本/GIS", "文本", "逐特征"),
        "n_batch_effect": ("同时有有效定量值及单一批次标签的记录数", "记录", "批次R²实际分母"),
        "n_batches_with_values": ("至少一个有效值的批次数k", "批次", "批次R²子集"),
        "batch_r2": ("批次组间平方和/总平方和；含未测样本组成差异", "0–1", "n_batch_effect"),
        "batch_adjusted_r2": ("1-(1-R²)*(n-1)/(n-k)；n>k才定义", "可为负", "n=n_batch_effect，k=有值批次数"),
        "smallest_batch_n_valid": ("有值批次中最少的有效值数", "记录", "当前特征与scope"),
        "largest_batch_n_valid": ("有值批次中最多的有效值数", "记录", "当前特征与scope"),
        "n_individual_detection_known": ("NULISA至少一条记录NPQ和板LOD均已知的记录/人数", "记录/人", "当前scope"),
        "n_detect_any_repeat": ("至少一条技术记录达到其板LOD", "记录/人", "当前scope"),
        "n_detect_all_repeats": ("所有技术记录的NPQ/LOD均已知且均NPQ≥LOD", "记录/人", "当前scope"),
        "detect_any_fraction_all_records": ("任一重复达到LOD人数/n_records", "比例", "全部当前scope记录，不排除未知"),
        "detect_all_fraction_all_records": ("全部重复达到LOD人数/n_records", "比例", "全部当前scope记录，不排除未知"),
        "n_with_npq_equal_lod_records": ("NPQ恰等于其板LOD的原始技术记录条数；可能非独立人数", "条", "当前scope原始技术记录"),
        "detection_basis": ("检出口径；Soma/TMT无个体LOD，只审计数值可用", "文本", "逐特征"),
        "feature_raw_accession_field": ("SomaScan此适配体的原始accession字段", "文本", "Soma特征"),
        "feature_is_multiaccession": ("Soma适配体注释包含多个accession", "布尔", "Soma特征"),
        "batch_r2_before_GIS_on_identical_available_subset": ("GIS后同一子集上原始TMT的plex R²", "0–1", "当前scope且样本/GIS均正"),
        "batch_r2_after_GIS": ("同子集GIS比值的plex R²", "0–1", "当前scope且样本/GIS均正"),
        "delta_batch_r2_GIS_minus_raw_same_subset": ("同子集后减前的plex R²差", "R²差", "当前scope且样本/GIS均正"),
        "GIS_ratio_is_audit_sensitivity_not_frozen_processing": ("GIS比值是审计敏感性处理、尚未冻结", "布尔", "TMT GIS路径"),
    }
    masks = {
        "N_numeric": "NULISA有有效NPQ",
        "N_detect_any_repeat": "NULISA任一技术记录达到对应板LOD",
        "N_detect_all_repeats": "NULISA所有技术记录均已知且达到对应板LOD",
        "S_any": "Soma至少一个匹配适配体有正RFU",
        "S_all": "Soma全部匹配适配体均有正RFU",
        "S_unambiguous_any": "Soma仅单accession特征中至少一个可用；没有此类特征则False",
        "S_unambiguous_all": "Soma仅单accession特征全部可用；没有此类特征则False",
        "T_raw": "TMT样本原始强度为正",
        "GIS_positive": "该参与者所在plex的此蛋白GIS强度为正",
        "T_GIS_usable": "TMT样本及同plex GIS强度均为正",
        "complete_raw_S_any": "N数值+S任一+T原始正值同时可用",
        "complete_raw_S_all": "N数值+S全部+T原始正值同时可用",
        "complete_GIS_S_any": "N数值+S任一+T样本及GIS均正同时可用",
        "complete_GIS_S_all": "N数值+S全部+T样本及GIS均正同时可用",
        "complete_GIS_unambiguous_S_any": "N数值+S单accession特征任一+T及GIS均正",
        "complete_GIS_unambiguous_S_all": "N数值+S单accession特征全部+T及GIS均正",
        "complete_GIS_detect_all_N_S_any": "N数值且全部技术重复≥LOD+S任一+T及GIS均正",
        "complete_GIS_detect_all_N_S_all": "N数值且全部技术重复≥LOD+S全部+T及GIS均正",
    }
    for mask, description in masks.items():
        definitions["n_"+mask] = (description+"的人数", "人", "三平台集固定571")
        definitions["fraction_"+mask] = (description+"的人数/571", "比例", "固定571")
    for quantile in ["01", "05", "25", "50", "75", "95", "99"]:
        definitions["p"+quantile] = (f"第{int(quantile)}百分位；线性插值", "对应scale", "有效数值")
    rows = []
    for filename, frame in [("PROTEIN_EVALUABILITY_70.csv", proteins), ("FEATURE_LEVEL_QC.csv", features)]:
        for field in frame.columns:
            if field not in definitions:
                raise ValueError(f"Undocumented field: {field}")
            description, unit, denominator = definitions[field]
            rows.append({"数据表": filename, "字段": field, "中文定义": description, "单位": unit, "分母或范围": denominator})
    d = pd.DataFrame(rows)
    write_csv(RESULTS / "METRIC_DICTIONARY_ZH.csv", d)
    return d


def friendly_table(p, n):
    columns = {
        "nulisa_target": "蛋白靶标", "accession": "UniProt",
        "annotation_strict_candidate": "注释单特征候选_非最终Strict",
        "soma_feature_count": "Soma适配体数", "soma_has_multiaccession_annotation": "含复合注释",
        "n_N_numeric": "NULISA数值人数", "n_N_detect_all_repeats": "NULISA全部重复达到LOD人数",
        "N_detection_fraction_all_repeats": "NULISA检出比例_分母571",
        "n_S_any": "Soma任一可用人数", "n_S_all": "Soma全部可用人数",
        "n_T_raw": "TMT原始正值人数", "n_T_GIS_usable": "TMT及GIS均正人数",
        "n_complete_raw_S_all": "三平台完整人数_原始TMT_Soma全部",
        "n_complete_GIS_S_any": "三平台完整人数_GIS_Soma任一",
        "n_complete_GIS_S_all": "三平台完整人数_GIS_Soma全部",
        "n_complete_GIS_detect_all_N_S_all": "完整人数_GIS_Soma全部_NULISA达到LOD",
        "fraction_complete_GIS_detect_all_N_S_all": "完整比例_GIS及LOD_分母571",
        "GIS_positive_plexes": "GIS正值plex数_分母65",
        "T_raw_batch_r2_on_GIS_available": "TMT批次R2_GIS前_匹配子集",
        "T_GIS_batch_r2": "TMT批次R2_GIS后",
        "T_GIS_delta_batch_r2_matched": "批次R2变化_后减前",
        "N_NPQ_sd": "NULISA_NPQ_SD", "S_log2_sd_min": "Soma_log2_SD最小",
        "S_log2_sd_max": "Soma_log2_SD最大", "T_GIS_log2_sd": "TMT_GIS_log2_SD",
        "audit_flags": "审计标记_非删除规则",
    }
    d = p.sort_values(["annotation_strict_candidate", "nulisa_target"], ascending=[False, True])[list(columns)].rename(columns=columns)
    d["80%示例覆盖_非纳入规则"] = np.where(p.loc[d.index, "n_complete_GIS_detect_all_N_S_all"].ge(math.ceil(n*.8)),
                                            f"达到示例线≥{math.ceil(n*.8)}人", "未达到示例线")
    d["最终状态"] = "待SAP/功效/病理覆盖冻结"
    return d.reset_index(drop=True)


def figures(p, n):
    plt.rcParams.update({"font.size": 10, "pdf.fonttype": 42})
    d = p.sort_values(["annotation_strict_candidate", "nulisa_target"], ascending=[False, True])
    fields = ["n_N_numeric", "n_N_detect_all_repeats", "n_S_all", "n_T_raw", "n_T_GIS_usable",
              "n_complete_GIS_S_all", "n_complete_GIS_detect_all_N_S_all"]
    names = ["NULISA\nnumeric", "NULISA\n>= LOD", "SomaScan\nall features", "TMT\npositive", "TMT + GIS\npositive", "Complete\nGIS", "Complete\nGIS + N LOD"]
    values = d[fields].to_numpy()/n
    fig, ax = plt.subplots(figsize=(11.4, 20.8))
    im = ax.imshow(values, aspect="auto", cmap="YlGnBu", norm=Normalize(0, 1), interpolation="none")
    ax.set_xticks(range(len(fields)), names)
    ax.set_yticks(range(len(d)), [str(x.nulisa_target)+(" *" if x.annotation_strict_candidate else "") for x in d.itertuples()])
    ax.tick_params(axis="y", labelsize=9, length=0)
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    for i in range(len(d)):
        for j in range(len(fields)):
            ax.text(j, i, f"{int(d.iloc[i][fields[j]])}", ha="center", va="center", fontsize=8,
                    color="white" if values[i, j] > .62 else "black")
    split = int(d.annotation_strict_candidate.sum()) - .5
    ax.axhline(split, color="#dd8452", lw=2)
    ax.set_title(f"Candidate availability audit | {n} same-date participants\n* Annotation single-feature candidate (43); remaining 27 require mapping / feature decisions", pad=58, fontsize=12)
    fig.colorbar(im, ax=ax, fraction=.025, pad=.025, label=f"Available fraction of {n}")
    fig.text(.12, .015, "Cells show participant counts. All SomaScan features must be available in 'Complete'.\nLOD rule: every NULISA replicate NPQ >= its plate LOD. Availability does not establish assay reliability.\nNo candidate set or preprocessing method was frozen by this audit.", fontsize=9)
    fig.subplots_adjust(left=.15, right=.92, top=.91, bottom=.055)
    for extension in ["png", "pdf"]:
        fig.savefig(FIGURES / f"01_CANDIDATE_AVAILABILITY_70.{extension}", dpi=190)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.8, 7.6))
    for status, marker, color, label in [(True, "o", "#286b9c", "43 annotation single-feature candidates"),
                                        (False, "s", "#cf8449", "27 other extended candidates")]:
        q = p.loc[p.annotation_strict_candidate.eq(status)]
        ax.scatter(q.T_raw_batch_r2_on_GIS_available, q.T_GIS_batch_r2, s=42, marker=marker,
                   color=color, alpha=.8, label=label)
    ax.plot([0, 1], [0, 1], "--", color=".5", lw=1)
    for row in p.loc[p.T_GIS_delta_batch_r2_matched.gt(0)].itertuples():
        ax.annotate(row.nulisa_target, (row.T_raw_batch_r2_on_GIS_available, row.T_GIS_batch_r2),
                    xytext=(10, 8), textcoords="offset points", fontsize=10)
    ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="Plex R² before GIS (identical available subset)",
           ylabel="Plex R² after log2(sample / same-plex GIS)",
           title="TMT batch description: 68 decrease, 2 increase\nOne dot per candidate; available subset varies by protein")
    ax.set_aspect("equal")
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    ax.grid(alpha=.15)
    fig.text(.1, .035, "Plex R² also reflects sample composition; it is not a pure technical variance estimate.\nGIS ratios are an audit sensitivity analysis, not a frozen processing rule.", fontsize=9)
    fig.subplots_adjust(left=.13, bottom=.17, top=.87, right=.96)
    for extension in ["png", "pdf"]:
        fig.savefig(FIGURES / f"02_TMT_BATCH_R2_BEFORE_AFTER_GIS.{extension}", dpi=190)
    plt.close(fig)


def write_workbook(p, feature, dec, metric_dictionary, batch_sizes, batch_summary):
    summary = load_json("AUDIT_SUMMARY.json")
    n = summary["triple_registry_n"]
    overview = friendly_table(p, n)
    guide = pd.DataFrame([
        ("状态", "Stage 1审计完成；最终Strict core及SAP尚未冻结。"),
        ("核心分母", "三平台精确同RID同日571人；Soma原生737条记录来自735人，不应混作737个独立参与者。"),
        ("70与43", "70为扩展注释候选；43为注释上单特征匹配候选，均不是已经完成效能验证的蛋白集合。"),
        ("主要完整人数", "NULISA有效NPQ、全部匹配Soma适配体正RFU、TMT样本及同plex GIS正值的共同人数。"),
        ("LOD严格人数", "在主要完整人数基础上要求NULISA所有技术重复NPQ≥各自板LOD。"),
        ("80%示例", "80%×571向上取整为457；该列仅描述覆盖，未用于冻结纳入，也不是功效证明。"),
        ("原始尺度", "NULISA发布NPQ已是对数尺度；Soma正RFU取log2；TMT -4视为缺失，正值取log2。"),
        ("GIS处理", "log2样本-log2同plex GIS；原始可用与GIS可用分别报告。仅审计敏感性，未冻结。"),
        ("TargetDetectability", "NULISA靶标×板汇总百分比，不能作为逐个体检出标志。"),
        ("批次R2", "板/plex组间平方和占比，含样本组成差异。GIS前后用完全相同子集比较。"),
        ("Soma多适配体", "不聚合、不按跨平台/病理相关择优。任一/全部可用并列；内部相关仅作描述性QC。"),
        ("未运行", "未打开I09–I15临床/病理文件；未运行跨平台共识、不一致或任何病理/认知关联。"),
        ("个体明细", "Excel不含RID/PTID/标本号/精确日期。参与者映射和缺失布尔矩阵在00_restricted。"),
        ("复合注释", "IL6R、UBB、HBA1涉及复合Soma注释；同accession不证明三平台检测同一种分子物种。"),
        ("版本开放项", "TMT方法3,447与矩阵3,907的差异原因未证实；详见TMT_COLUMN_COUNT_RECONCILIATION.md。"),
        ("独立验证", "1,265项计数及2,520项数值矩核对全部通过；详细CSV在results。"),
        ("科学解释", "有数值、达到LOD、技术可靠性、病理相关性分别需要证据，不能互相替代。"),
    ], columns=["项目", "说明"])
    source_sheets = [
        ("逐蛋白完整字段70", p),
        ("映射70", read("EXTENDED_CORE_MAPPING_DRAFT.csv")),
        ("映射43草案", read("STRICT_CORE_MAPPING_DRAFT.csv")),
        ("逐特征QC", feature),
        ("适配体汇总", read("MULTIAPTAMER_AUDIT.csv")),
        ("适配体配对QC", read("MULTIAPTAMER_PAIRWISE_QC.csv")),
        ("GIS逐蛋白", read("TMT_GIS_SUMMARY.csv")),
        ("GIS逐plex", read("TMT_GIS_BATCH_AUDIT.csv")),
        ("NULISA板检出", read("NULISA_PLATE_TARGET_DETECTION.csv")),
        ("原生批次QC", read("PLATFORM_BATCH_QC_LONG.csv").query("scope=='NATIVE_RECORDS'")),
        ("同日批次QC", read("PLATFORM_BATCH_QC_LONG.csv").query("scope=='TRIPLE_SAME_DATE'")),
        ("同日批次样本数", batch_sizes),
        ("跨平台批次分配", batch_summary),
        ("精确日期配对人数", read("PAIRWISE_SAMPLE_COUNTS.csv")),
        ("ID核对汇总", read("SAMPLE_ID_CHECK_SUMMARY.csv")),
        ("描述阈值网格", read("DESCRIPTIVE_COVERAGE_THRESHOLDS.csv")),
        ("QA检查", read("QA_CHECKS.csv")),
        ("输入校验", read("INPUT_HASH_VERIFICATION.csv")),
        ("开放决定", dec),
        ("字段字典", metric_dictionary),
        ("TMT全部列注释", read("TMT_PROTEIN_COLUMN_ANNOTATION.csv")),
        ("NULISA靶标分类", read("NULISA_TARGET_CLASSIFICATION.csv")),
    ]
    wb = Workbook()
    wb.remove(wb.active)
    schemas = []
    all_sheets = [("阅读说明", guide), ("70候选总览", overview),
                  ("43单特征候选", overview.loc[overview["注释单特征候选_非最终Strict"]].reset_index(drop=True))] + source_sheets
    for name, frame in all_sheets:
        ws = wb.create_sheet(name)
        ws.append(list(frame.columns))
        for row in frame.itertuples(index=False, name=None):
            values = []
            for value in row:
                if pd.isna(value):
                    values.append(None)
                elif isinstance(value, np.generic):
                    values.append(value.item())
                elif isinstance(value, str) and value.startswith("="):
                    values.append("'"+value)
                else:
                    values.append(value)
            ws.append(values)
        ws.freeze_panes = "C2" if name in {"70候选总览", "43单特征候选"} else "A2"
        ws.auto_filter.ref = ws.dimensions
        ws.sheet_view.zoomScale = 85
        ws.row_dimensions[1].height = 45
        for c in ws[1]:
            c.fill = PatternFill("solid", fgColor="214D6B")
            c.font = Font(color="FFFFFF", bold=True)
            c.alignment = Alignment(wrap_text=True, vertical="center")
        for j, column in enumerate(frame.columns, 1):
            sampled = frame[column].dropna().astype(str).head(100)
            longest = max([len(str(column)), *[len(x) for x in sampled]], default=12)
            width = min(42, max(14, longest * (1.2 if any(ord(c)>127 for c in str(column)) else .95)))
            if name == "阅读说明":
                width = 26 if j == 1 else 112
            if name == "字段字典":
                width = {1:38, 2:48, 3:85, 4:20, 5:50}.get(j, 22)
            if name == "开放决定":
                width = {1:16, 2:28, 3:35, 4:85, 5:95}.get(j, 22)
            ws.column_dimensions[get_column_letter(j)].width = width
            percent = (str(column).startswith("fraction_") or "比例" in str(column) or
                       str(column) in {"missing_fraction", "N_detection_fraction_all_repeats"})
            if pd.api.types.is_float_dtype(frame[column]):
                fmt = "0.0%" if percent else "0.0000"
                for cell in list(ws.iter_cols(min_col=j, max_col=j, min_row=2))[0]:
                    cell.number_format = fmt
            if name in {"70候选总览", "43单特征候选"} and ("人数" in str(column) or "比例" in str(column)):
                ws.conditional_formatting.add(f"{get_column_letter(j)}2:{get_column_letter(j)}{len(frame)+1}",
                                              DataBarRule(start_type="num", start_value=0, end_type="num",
                                                          end_value=1 if percent else n, color="77AED4", showValue=True))
        if name in {"阅读说明", "字段字典", "开放决定"}:
            for row in ws.iter_rows(min_row=2):
                for cell in row:
                    cell.alignment = Alignment(wrap_text=True, vertical="top")
                ws.row_dimensions[row[0].row].height = 45 if name != "开放决定" else 65
        schemas.append({"sheet": name, "data_rows": len(frame), "columns": len(frame.columns)})
    wb.properties.title = "候选逐蛋白可评价性审计（Stage 1）"
    wb.properties.subject = "无临床结局审计；候选集合未冻结"
    wb.properties.creator = "Research audit"
    wb.save(ROOT / "候选逐蛋白可评价性审计总表.xlsx")
    write_csv(RESULTS / "WORKBOOK_SHEET_COUNTS.csv", schemas)
    return overview


def write_reports(p, features, dec, batch_sizes):
    summary = load_json("AUDIT_SUMMARY.json")
    independent = load_json("INDEPENDENT_VALIDATION_SUMMARY.json")
    n, nstruct = load_json("NULISA_STRUCTURAL_AUDIT.json"), load_json("TMT_STRUCTURAL_AUDIT.json")
    s = load_json("SOMASCAN_STRUCTURAL_AUDIT.json")
    strict = p.loc[p.annotation_strict_candidate]
    population = summary["triple_registry_n"]
    qacount = len(read("QA_CHECKS.csv"))
    low = p.sort_values("n_complete_GIS_detect_all_N_S_all").head(13)[[
        "nulisa_target", "annotation_strict_candidate", "n_N_detect_all_repeats", "n_complete_GIS_S_all", "n_complete_GIS_detect_all_N_S_all"]].copy()
    low.columns = ["蛋白", "注释单特征", "NULISA达到LOD人数", "GIS三平台完整人数", "再要求NULISA达到LOD人数"]
    low["注释单特征"] = low["注释单特征"].map({True: "是", False: "否"})
    grid_rows = []
    for fraction in [.5, .7, .8, .9, .95, 1.0]:
        grid_rows.append({"描述覆盖比例": f"{fraction:.0%}", "对应至少人数": math.ceil(population*fraction),
                          "70候选_GIS完整": int(p.n_complete_GIS_S_all.ge(population*fraction-1e-10).sum()),
                          "70候选_GIS及LOD完整": int(p.n_complete_GIS_detect_all_N_S_all.ge(population*fraction-1e-10).sum()),
                          "43候选_GIS完整": int(strict.n_complete_GIS_S_all.ge(population*fraction-1e-10).sum()),
                          "43候选_GIS及LOD完整": int(strict.n_complete_GIS_detect_all_N_S_all.ge(population*fraction-1e-10).sum())})
    grid = pd.DataFrame(grid_rows)
    write_csv(RESULTS / "COVERAGE_GRID_READING_ZH.csv", grid)
    example = strict.loc[strict.n_complete_GIS_detect_all_N_S_all.ge(math.ceil(population*.8))].copy()
    example["descriptive_rule_only_not_frozen"] = "43 annotation single-feature AND GIS/all-Soma/N-LOD complete >=457 of571; illustrative only"
    write_csv(RESULTS / "DESCRIPTIVE_80PCT_PRIORITY_CANDIDATES.csv", example)
    pairs = read("PAIRWISE_SAMPLE_COUNTS.csv")[["comparison", "common_RID_any_date", "same_date_participants", "same_date_VISCODE2_match", "same_date_VISCODE2_unknown"]]
    pairs.columns = ["平台组合", "共同RID不论日期", "精确同日人数", "访视一致", "访视未知"]
    exceptions = p.loc[p.T_GIS_delta_batch_r2_matched.gt(0), ["nulisa_target", "T_raw_batch_r2_on_GIS_available", "T_GIS_batch_r2", "T_GIS_delta_batch_r2_matched"]].copy()
    exceptions.iloc[:, 1:] = exceptions.iloc[:, 1:].round(4)
    exceptions.columns = ["蛋白", "同子集GIS前R²", "GIS后R²", "后减前"]
    mult = read("MULTIAPTAMER_AUDIT.csv").query("scope=='TRIPLE_SAME_DATE' and any_nonpositive_pair_rho")[["target", "n_features", "n_evaluable_pairs", "rho_min", "rho_max"]].copy()
    mult[["rho_min", "rho_max"]] = mult[["rho_min", "rho_max"]].round(4)
    mult.columns = ["蛋白", "适配体数", "可评价对数", "最小Spearman相关", "最大Spearman相关"]
    batches = batch_sizes.groupby("platform").n_triple_participants.agg(["count", "min", "median", "max"]).reset_index()
    batches.columns = ["平台", "批次数", "每批最少人数", "每批中位人数", "每批最多人数"]

    report = f"""# 候选逐蛋白可评价性审计报告

审计日期：2026-10-26。状态：**Stage 1 PASS（审计交付通过，冻结前开放项保留）**。

本轮确认了**571名精确同RID同日的三平台参与者、70个扩展注释候选、43个注释上单特征匹配候选**。43个候选的GIS三平台完整人数中位数为559；再要求NULISA全部技术记录达到对应板LOD后，中位数为542。以80%覆盖（至少457/571人）仅作描述，分别有35个和30个候选达到该示例线。**30与43均不是最终Strict core；本轮没有冻结纳入规则，也没有证明功效或病理覆盖充足。**

## 1. 审计边界和验证

- 仅打开移交清单I01–I08的检测数据、注释、方法文件和数据字典；8个文件的大小、SHA-256均与原清单一致。
- 未打开I09–I15临床、病理、诊断、人口学和APOE数据。未计算跨平台共识/不一致、平台间蛋白相关、A/T或认知关联。
- 定量读取仅用于缺失、分布、方差、LOD、批次与可评价性。多适配体相关限SomaScan平台内部，未据此择优。
- 审计前方案为[AUDIT_PROTOCOL_v1.0.md](AUDIT_PROTOCOL_v1.0.md)，其事前SHA保存在[ACCESS_BOUNDARY_AUDIT.json](results/ACCESS_BOUNDARY_AUDIT.json)。原始文件、原移交包、原SAP和既有论文未改写。
- 主流程{qacount}项QA全部通过。独立程序用Python标准库csv/math/statistics直接重读原始数据，未导入主审计的加载或计算函数；**{independent['count_checks']:,}项人数、{independent['numeric_moment_checks']:,}项有效数/均值/样本SD/方差核对均通过**。独立验证读取主流程候选映射，不构成独立的分子靶标特异性认证。
- 本轮“未打开临床文件”依据代码限定输入与执行记录，并非操作系统级取证日志。最终交付核验见[AUDIT_COMPLETION_QA.json](results/AUDIT_COMPLETION_QA.json)。

## 2. 数据结构、QC与重复

| 项目 | 本轮核实情况 | 审计处理 |
| --- | --- | --- |
| NULISA全表 | 574,435长表行：PLASMA 295,012、CSF 217,853、CONTROL 61,570 | 本研究候选审计使用CSF；血浆仅复现结构与重复情况 |
| NULISA CSF | 1,530个RID；1,517通过QC、13 warning；19板 | 限定CSF、SampleType=Sample、QC=passed、有效RID与日期 |
| NULISA重复 | 血浆824人有完全重复行，冗余105,472行；CSF无完全重复，1人131靶标跨两板，涉及262条非相同记录 | 完全重复先去除；CSF非相同重复NPQ取均值仅用于审计，保留每条原值及任一/全部检出状态 |
| NULISA无RID样本 | CSF有132个SampleName×PlateID块，102个SampleName | 不虚构参与者身份、不计入1,517人 |
| NULISA靶标 | 134原始字符串归一为131靶标；Aβ38/40/42存在编码别名 | 保留原标签；APOE/CRP在CSF无NPQ；本文件未见KNG1靶标 |
| SomaScan | 737条发布后QC记录、735人、7,008适配体、36板；无同RID同日未解决重复 | 原生分布分母为737条记录，不写成737名独立参与者 |
| SomaScan访视 | 708 bl、22其他访视、7访视缺失 | 精确日期配对；不默认所有记录是基线 |
| TMT | 1,170行=1,105参与者+65 GIS；65 plex，每plex 17参与者+1个126通道GIS | GIS不作参与者；每人绑定其同plex GIS |
| TMT列与编码 | 3,907定量列+7元数据；参与者`-4`有1,492,044单元，全表1,579,629；未见其他非正/不可数值化强度 | `-4`转NA，仅正值进入log2；缺失不插补 |
| TMT访视 | 1,104 bl、1 m06；方法称基线样本 | 不将全部1,105行默认为基线；正式访视范围留待冻结 |

NULISA跨板重复的参与者不被强行指定到某一板，板效应汇总排除其复合板标签；定量与检出审计仍保留其聚合记录。该参与者不在571人的三平台同日集合内，因此三平台批次描述未因此少人。

结构证据见[NULISA](results/NULISA_STRUCTURAL_AUDIT.json)、[SomaScan](results/SOMASCAN_STRUCTURAL_AUDIT.json)、[TMT](results/TMT_STRUCTURAL_AUDIT.json)。技术重复的个体明细单独保存在`00_restricted`。

## 3. 身份、日期和同日配对

主键是有效整数RID＋精确EXAMDATE。若同一人存在多个三平台共同日期，方案事前规定取最早共同日期；实际没有这种重复共同日期。

{md_table(pairs)}

三平台共同RID为583，其中571能精确同日配对。571人的NULISA–TMT PTID完全一致；可用PTID的RID后缀核对无冲突。**SomaScan现有元数据没有可用PTID，无法独立完成SomaScan的PTID交叉核验。**VISCODE2有570人一致、1人未知，未观察到已知访视冲突。日期才是配对主键，访视不能替代日期。

TMT原始日期是美国月/日/年（例如12/30/2010），NULISA/SomaScan主要为ISO日期。独立验证最初只接受ISO，补入明确的`%m/%d/%Y`解析后完成复核；这属于验证器兼容性修正，未改动源日期，也未改变主流程的571人结果。

**同RID同日不等于已证明同一aliquot。**未找到足以跨三平台证明同一分装标本的共同标识。不能据此将跨平台差异全归为检测技术差异。

配对计数见[PAIRWISE_SAMPLE_COUNTS.csv](results/PAIRWISE_SAMPLE_COUNTS.csv)；含RID/PTID/精确日期的注册表仅在`00_restricted/TRIPLE_SAMPLE_REGISTRY_RESTRICTED.csv`。

## 4. 70个候选与43个单特征候选的含义

映射由原始NULISA/SomaScan注释及TMT列重新生成，并与既有名称核对表比较。排除Aβ片段、APOE4、位点/修饰型tau、BD-MAPT、Oligo-SNCA、pSNCA、pTDP43、IL12p70等特异形式/复合靶标，以及CSF无NPQ的APOE/CRP。“总量型”只是注释类别，不能证明检测完整全长蛋白。

- 扩展注释候选70个。44个只有一个匹配Soma特征，26个有多个适配体。
- 其中43个同时满足NULISA唯一、Soma仅一个且单accession、TMT唯一普通UniProt列及无未决isoform映射。
- HBA1虽只有一个适配体，但注释含P69905/P68871，因此不在43个中。IL6R和UBB也涉及多accession特征；本轮另列“只保留单accession特征”的任一/全部可用人数。
- **26个多适配体加1个HBA1形成其余27个候选。没有通过与其他平台或病理的相关性选适配体。**

映射草案：[70候选](results/EXTENDED_CORE_MAPPING_DRAFT.csv)、[43单特征候选](results/STRICT_CORE_MAPPING_DRAFT.csv)。分子身份并未独立验证，所有最终SAP集合冻结标志均为False。

## 5. 有数值、达到LOD和完整观测的区别

NULISA使用发布NPQ，不再取log2。方法文件明确其已经过对数变换。`TargetDetectability`在靶标×板内恒定，本轮观察到2,489个靶标×板组合；它是汇总百分比，**不是某个参与者的检出标志**。逐记录审计检出定义为NPQ≥本条记录对应的板LOD，等于LOD的条数单列。技术重复分别报告任一与全部达到LOD；“全部”还要求每条重复的NPQ/LOD均已知。

SomaScan正RFU取log2。没有提供可作此审计的逐个体LOD，因此仅称数值可用。TMT样本正值取log2；GIS口径要求该样本及同plex的该蛋白GIS都为正，使用log2(样本/GIS)。未缩尾、未插补、未按板残差化。

主要完整人数是**NULISA数值可用＋全部匹配Soma适配体可用＋TMT样本和GIS均正**的交集；LOD严格人数再要求NULISA所有技术重复达到各自板LOD。多适配体“任一可用”也完整保留，但不把“任一可用”当成已经完成一种蛋白聚合。

| 候选范围与口径 | 最少人数 | 中位人数 | 最多人数 |
| --- | --- | --- | --- |
| 70个扩展：GIS＋任一Soma特征 | 53 | 565 | 571 |
| 70个扩展：GIS＋全部Soma特征 | 53 | 557 | 571 |
| 43个单特征：GIS完整 | 53 | 559 | 571 |
| 43个单特征：GIS完整且NULISA全部重复达到LOD | 4 | 542 | 571 |

以下阈值均为预先列出的描述网格，**不是正式纳入标准或功效合格线**。每格是达到该覆盖线的蛋白数；“GIS完整”使用全部Soma特征。

{md_table(grid)}

例如80%对应至少457人：43候选中35个达到GIS完整示例线，30个再满足NULISA LOD要求。30个示例候选见[描述性80%名单](results/DESCRIPTIVE_80PCT_PRIORITY_CANDIDATES.csv)，保留它是为了便于后续模拟，不可当作已经冻结的主分析名单。

## 6. 低覆盖与低检出蛋白

{md_table(low)}

在43个单特征候选中，MME、S100A12、TEK的GIS三平台完整人数不足100。YWHAZ和BASP1的数值完整人数很高，但分别只有4/571和29/571达到NULISA LOD，不能因数值完整而认定检测充分。

70个候选的NULISA检出不足50%者为MME、UBB、SAA1、ARSA、BDNF、YWHAZ、S100A12、BASP1、PARK7；43候选中的对应4个是MME、YWHAZ、S100A12、BASP1。未出现完全零检出靶标。三平台集的逐特征检查未见零方差、SD≤1e-6或不同数值≤5的特征；这不抵消低于LOD的限制。

低检出蛋白仍保留在审计全集，后续需预先决定左删失/敏感性及范围，不能看病理效应后决定去留。

## 7. TMT GIS与批次

70个候选共形成4,550条蛋白×plex记录。GIS缺失plex数中位0，最多58/65；稀疏蛋白的可评价性集中在少数plex。要求GIS为正较仅要求样本为正，总共少25个“参与者×蛋白”可用单元：SLIT2少6、S100A12少12、PDLIM5少7。这不是少25名独立参与者。

GIS前后的one-way plex R²严格在各蛋白相同的样本/GIS正值子集比较。70个中68个下降，2个上升；后减前的R²差中位数−0.4602，范围−0.7359至+0.2891。例外为：

{md_table(exceptions)}

R²反映plex间均值差异，包含样本组成差异，**不能直接称为纯技术方差**。GIS比值和跨plex CV也不能替代个体测量可靠性估计。不能仅据68个下降便宣布GIS已被证明为所有蛋白的最优正式处理。

三平台571人覆盖NULISA 19板、SomaScan 36板及TMT 65 plex。各批次参与者人数如下；完整批次分配计数另存CSV。

{md_table(batches)}

目前只审计批次分配和分布，未检查其与病理/诊断的混杂。TMT的3,447/3,907差异尚未解决，详见[TMT_COLUMN_COUNT_RECONCILIATION.md](TMT_COLUMN_COUNT_RECONCILIATION.md)。当前矩阵可用于本轮描述审计，不能据方法文档任意裁到3,447列。

## 8. 多适配体与技术可靠性边界

26个多适配体蛋白按原生737记录及三平台571人分别计算平台内部两两Spearman相关。以下蛋白在三平台集至少有一对相关≤0：

{md_table(mult)}

存在相关强的适配体对并不抵消同蛋白另一些适配体的弱/非正相关，也不证明任意平均合理。上述相关未用于选“最佳”适配体。43单特征候选避免了聚合问题，但仍不证明相同分子物种或测量可靠性。

NULISA CSF仅1位跨板技术重复参与者，无法在逐蛋白层面稳定估计ICC；TMT GIS来自池标准，不是个体技术复测。因而本审计无法把真实丰度、蛋白片段/表位与测量误差充分分开。

## 9. 下一阶段应解决的具体问题

{md_table(dec[["decision_id", "问题", "状态", "对冻结或解释的影响"]])}

建议先以43单特征候选的全部覆盖结构作功效和稳定性模拟，再比较35个GIS覆盖≥80%与30个GIS/LOD覆盖≥80%的描述子集。模拟必须包含左删失、蛋白依赖缺失、19/36/65批次及平台误差不等；只有覆盖阈值远远不够。

延续既有研究指导的模型问题：PC1和与其正交的不一致坐标不能自动识别“真实共同丰度”与“病理相关检测差异”。本轮没有运行这些分解，正式分析前仍须完成无病理特异检测效应的人工零模型检验，避免把平台误差不等造成的几何结果作机制解释。

## 10. 阅读与复现

先看[Excel总表](候选逐蛋白可评价性审计总表.xlsx)的“阅读说明”“70候选总览”“43单特征候选”“开放决定”；所有原始字段仍在CSV与Excel完整字段页。

- [逐蛋白70候选](results/PROTEIN_EVALUABILITY_70.csv)：每蛋白完整人数、映射、检出、SD、R²及审计标记。
- [逐特征QC](results/FEATURE_LEVEL_QC.csv)：630行=315个特征/处理路径×2个scope；其中TMT原始和GIS路径分开，不能称为315个不同检测。实际候选映射为70个NULISA、105个Soma和70个TMT列。
- [逐批次QC](results/PLATFORM_BATCH_QC_LONG.csv)、[GIS逐plex](results/TMT_GIS_BATCH_AUDIT.csv)、[适配体配对QC](results/MULTIAPTAMER_PAIRWISE_QC.csv)、[字段字典](results/METRIC_DICTIONARY_ZH.csv)。
- [70候选可用性图](figures/01_CANDIDATE_AVAILABILITY_70.png)、[GIS前后批次R²图](figures/02_TMT_BATCH_R2_BEFORE_AFTER_GIS.png)，均附PDF。
- [主流程QA](results/QA_CHECKS.csv)、[独立验证摘要](results/INDEPENDENT_VALIDATION_SUMMARY.json)、[输入输出SHA-256](STAGE1_INPUT_OUTPUT_SHA256.csv)。
- 运行`code/run_all.ps1`依次重建主审计、独立验证、交付与封存核验。固定输入可复现数值；Excel/PDF创建时间、执行日志及其文件哈希会随重跑变化。
- `00_restricted`保存参与者映射、技术重复和布尔可用性矩阵，仅作本地授权科研使用。目录名不是额外操作系统权限控制；未上传至外部服务。

本报告完成了Stage 1的无结局审计。列出的冻结前开放项需要在后续方法决定中处理；本报告没有替代那些决定。
"""
    (ROOT / "STAGE1_QA_REPORT_ZH.md").write_text(report, encoding="utf-8")
    memo = f"""# TMT 3,447 / 3,907列差异核对备忘录

2026-10-26；状态：**本地证据已核对，来源差异原因仍未解决**。

## 绑定证据

- I07为原清单绑定的`ADNI_EMORY_CSF_TMT_MS.zip`，其大小及SHA-256均与移交清单一致。zip内`EMORY_CSF_TMT_MS.csv`为{nstruct['inner_bytes']:,}字节，SHA-256：`{nstruct['inner_sha256']}`。
- I08方法文件`ADNI_EMORY_CSF_TMT_MS_Methods_20251007.pdf`第1页称3,447 proteins及未作batch correction；第3页说明0强度按缺失编码为−4；第4页解释HGNC_UniProtID列名与GIS。
- 实际CSV共有3,914列：7元数据列（PHASE/PTID/ProteomicID/RID/VISCODE/VISCODE2/EXAMDATE）加**3,907定量列**；1,170行包含1,105参与者和65 GIS。
- 3,907定量列中3,904列具有普通UniProt后缀，对应3,904个不同规范化accession；其余为`ABETA40_ABETA40`、`APOE2_APOE2`、`APOE4_APOE4`。

## 能确定和不能确定的内容

数值差为460列，但本地文件未提供3,447列清单或对应版本映射，不能确定这460列的身份或来源。将3个特殊列排除后仍有3,904列，差异不是由3个特殊列解释。亦没有证据足以把差异归因于GIS、元数据列、isoform合并或某种QC过滤。

不得将任意缺失率阈值、选前3,447列或主观去列当作“复现方法文档”。也不能因为3907个列名均存在就声称它们已全通过正式生物学可评价性。

## 本轮采用的口径及影响

审计以已核实SHA的3,907列矩阵作为**实际输入**，将−4视为缺失，GIS独立于参与者处理。三平台70候选及43单特征候选均从这份矩阵的注释连接得到；本轮结果的可复现性不依赖猜测3,447列版本。

在冻结正式研究数据版本前，应取得ADNI/数据提交方针对该矩阵的版本说明或与其一致的方法文档/发布记录，确认3,447是否为另一批筛选或另一版矩阵。若输入版本改变，须重新做哈希、映射、配对和可评价性审计。此次未对外发送询问信息。

证据表：[TMT_STRUCTURAL_AUDIT.json](results/TMT_STRUCTURAL_AUDIT.json)、[全部列注释](results/TMT_PROTEIN_COLUMN_ANNOTATION.csv)、[方法摘录](results/METHODS_EVIDENCE_EXCERPTS.md)、[输入身份](results/INPUT_HASH_VERIFICATION.csv)。
"""
    (ROOT / "TMT_COLUMN_COUNT_RECONCILIATION.md").write_text(memo, encoding="utf-8")

    checklist = []
    evidence = [
        ("S01", "验证I01–I08身份", "results/INPUT_HASH_VERIFICATION.csv", "PASS"),
        ("S02", "NULISA结构、重复、QC和特殊靶标", "results/NULISA_STRUCTURAL_AUDIT.json", "PASS"),
        ("S03", "SomaScan结构和多适配体", "results/SOMASCAN_STRUCTURAL_AUDIT.json", "PASS"),
        ("S04", "TMT参与者、GIS、plex及-4", "results/TMT_STRUCTURAL_AUDIT.json", "PASS"),
        ("S05", "精确同日注册表", "00_restricted/TRIPLE_SAMPLE_REGISTRY_RESTRICTED.csv", "PASS_WITH_ALIQUOT_LIMITATION"),
        ("S06", "两两与三平台配对计数", "results/PAIRWISE_SAMPLE_COUNTS.csv", "PASS"),
        ("S07", "单特征映射草案", "results/STRICT_CORE_MAPPING_DRAFT.csv", "PASS_DRAFT_NOT_FROZEN"),
        ("S08", "扩展映射草案", "results/EXTENDED_CORE_MAPPING_DRAFT.csv", "PASS_DRAFT_NOT_FROZEN"),
        ("S09", "逐蛋白可评价性", "results/PROTEIN_PLATFORM_EVALUABILITY.csv", "PASS"),
        ("S10", "3,447/3,907备忘录", "TMT_COLUMN_COUNT_RECONCILIATION.md", "PASS_WITH_UNRESOLVED_SOURCE_VERSION"),
        ("S11", "QA与SHA封存", "STAGE1_QA_REPORT_ZH.md;STAGE1_INPUT_OUTPUT_SHA256.csv", "PASS_WITH_OPEN_DECISIONS"),
    ]
    for sid, task, evidence_path, status in evidence:
        checklist.append({"step_id": sid, "任务": task, "本目录证据": evidence_path, "Stage1状态": status})
    write_csv(RESULTS / "STAGE1_COMPLETION_CHECKLIST.csv", checklist)
    readme = f"""# 候选逐蛋白可评价性审计

2026-10-26。**Stage 1已完成；保留冻结前开放项。**

已核实571名三平台同日参与者、70个扩展候选、43个注释单特征候选。43候选中，35个GIS完整覆盖≥80%，再要求NULISA达到LOD后为30个；80%仅作描述，不是冻结纳入标准。

建议按以下顺序阅读：

1. [中文审计报告](STAGE1_QA_REPORT_ZH.md)：结果、限制与下一阶段决定。
2. [Excel总表](候选逐蛋白可评价性审计总表.xlsx)：70/43候选、逐特征/批次、适配体、GIS、QA及字段解释。
3. [TMT版本差异备忘录](TMT_COLUMN_COUNT_RECONCILIATION.md)：方法3,447与实际3,907的差异尚未解释。
4. [可用性图](figures/01_CANDIDATE_AVAILABILITY_70.png)、[GIS批次图](figures/02_TMT_BATCH_R2_BEFORE_AFTER_GIS.png)。
5. [交付完成核验](results/AUDIT_COMPLETION_QA.json)、[输入输出SHA-256](STAGE1_INPUT_OUTPUT_SHA256.csv)。

`results`包含机器可读聚合结果；`code/run_all.ps1`可复现全流程；`logs`保留执行记录。`00_restricted`为参与者映射及技术重复明细，未设额外系统ACL，请按已有数据授权管理。

本轮未打开临床/病理文件，未计算跨平台共识、不一致或病理/认知关联。原始数据、旧稿和原SAP未修改。
"""
    (ROOT / "README_审计总览.md").write_text(readme, encoding="utf-8")


def main():
    p, feature = read("PROTEIN_EVALUABILITY_70.csv"), read("FEATURE_LEVEL_QC.csv")
    def flattened(value, prefix=""):
        for key, item in value.items():
            field = f"{prefix}.{key}" if prefix else key
            if isinstance(item, dict):
                yield from flattened(item, field)
            else:
                yield {"metric": field, "value": "|".join(map(str, item)) if isinstance(item, list) else item}
    for platform in ["NULISA", "SOMASCAN", "TMT"]:
        write_csv(RESULTS / f"{platform}_STRUCTURAL_AUDIT.csv", flattened(load_json(f"{platform}_STRUCTURAL_AUDIT.json")))
    batch_sizes, batch_summary = metadata_batch_audit()
    dec = decisions()
    metric_dictionary = dictionary(p, feature)
    figures(p, load_json("AUDIT_SUMMARY.json")["triple_registry_n"])
    write_workbook(p, feature, dec, metric_dictionary, batch_sizes, batch_summary)
    write_reports(p, feature, dec, batch_sizes)
    (ROOT / "logs/DELIVERABLE_BUILD.log").write_text(datetime.now().isoformat()+" Report, workbook, dictionaries, figures and local completion checklist generated.\n", encoding="utf-8")
    print("Reports, workbook, figures and supporting CSVs generated.")


if __name__ == "__main__":
    main()
