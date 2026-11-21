from __future__ import annotations

import json
import math
import numpy as np
import pandas as pd
from scipy.stats import norm
from sim_config import *


def wilson(success,total):
    if total==0: return np.nan,np.nan
    z=1.959963984540054;p=success/total
    midpoint=(p+z*z/(2*total))/(1+z*z/total)
    width=z*math.sqrt(p*(1-p)/total+z*z/(4*total*total))/(1+z*z/total)
    return midpoint-width,midpoint+width


def method_registry():
    labels=[
        ('PCA_IN_SAMPLE','同样本PCA＋不一致坐标','描述基线','坐标关联不等于额外病理检测效应'),
        ('PCA_CROSSFIT_2FOLD','两折交叉拟合PCA','描述基线','交叉拟合不修复总体方向错配；此HC3亦不包含所有学习步骤的不确定性'),
        ('COMMON_LOADING_MEAN_WALD','共同方向均值约束','条件候选','线性共同响应，平台误差条件不相关，无选择性缺失；2df'),
        ('INVARIANCE_JOINT_IID','均值＋协方差联合约束','条件候选','允许共同均值/方差及误差方差变化；依赖线性、条件不相关误差及个体独立；8df'),
        ('INVARIANCE_MEAN_IID','均值约束诊断','联合检验的组成诊断','4df；不另建发现家族'),
        ('INVARIANCE_COV_IID','协方差约束诊断','联合检验的组成诊断','4df；不另建发现家族'),
        ('INVARIANCE_JOINT_TMT_CLUSTER','联合约束＋第三平台聚类','有限簇候选','仅处理第三平台批次依赖；渐近chi-square需校准'),
        ('INVARIANCE_JOINT_3WAY_CLUSTER','联合约束＋三向聚类','有限簇候选','容斥协方差可能非PSD，失败明确保留，未修补'),
        ('ORACLE_COMPLETE_TRUTH','已知T和真实响应的完整数据基准','非现实信息基准','使用现实没有的T、f(T)及未删失值；不能作为可部署方法'),
    ]
    d=pd.DataFrame(labels,columns=['method','中文名称','角色','假设或限制']);csv(RESULTS/'METHOD_REGISTRY.csv',d);return d


def main():
    grid=design_grid(scenarios());lookup={x['scenario_id']:x for x in scenarios()}
    rows=[];statusrows=[];paired=[]
    for cell in grid:
        d=pd.read_parquet(RESULTS/'replicates'/(cell['cell_id']+'.parquet'))
        for method,q in d.groupby('method',observed=True):
            total=len(q);success=q.status.eq('OK');rejected=success&q.p_value.lt(.05)
            rejected_n=int(rejected.sum());lo,hi=wilson(rejected_n,total)
            lower_cond,upper_cond=wilson(rejected_n,int(success.sum()))
            scenario=lookup[cell['scenario_id']]
            row={**cell,'method':str(method),'scenario_label':scenario['label_zh'],'truth_class':scenario['truth_class'],
                 'compatible_linear_independent_null':scenario['compatible_linear_independent_null'],
                 'n_replicates':total,'n_successful':int(success.sum()),'n_failed':int((~success).sum()),
                 'failure_rate':float((~success).mean()),'n_rejected':rejected_n,
                 'rejection_rate_all_replicates':rejected_n/total,'wilson95_lower':lo,'wilson95_upper':hi,
                 'monte_carlo_se':math.sqrt((rejected_n/total)*(1-rejected_n/total)/total),
                 'rejection_rate_successful_only':rejected_n/int(success.sum()) if success.any() else np.nan,
                 'conditional_wilson95_lower':lower_cond,'conditional_wilson95_upper':upper_cond,
                 'complete_n_min':int(q.n_complete.min()),'complete_n_median':float(q.n_complete.median()),
                 'complete_n_max':int(q.n_complete.max()),'smallest_group_n_median':float(q.smallest_group_n.median()),
                 'nonpositive_cross_covariance_rate':float(q.nonpositive_pooled_cross_covariance.mean()),
                 'interpretation_warning':'PCA: coordinate association, not scientific type-I error; pressure cases may violate constraint assumptions; failures retained in denominator'}
            rows.append(row)
            for status,count in q.status.value_counts().items():
                if count:statusrows.append({'cell_id':cell['cell_id'],'method':str(method),'status':str(status),'n':int(count)})
        # Paired comparisons use identical synthetic samples, not independent rate SEs.
        wide=d.pivot(index='replicate',columns='method',values='p_value')
        for baseline in ['PCA_IN_SAMPLE','PCA_CROSSFIT_2FOLD']:
            both=wide[[baseline,'INVARIANCE_JOINT_IID']].notna().all(axis=1)
            a=wide.loc[both,baseline].lt(.05);b=wide.loc[both,'INVARIANCE_JOINT_IID'].lt(.05)
            diff=a.astype(float)-b.astype(float)
            paired.append({'cell_id':cell['cell_id'],'baseline':baseline,'comparison':'INVARIANCE_JOINT_IID','n_pairs':int(both.sum()),
                           'baseline_only_reject':int((a&~b).sum()),'comparison_only_reject':int((~a&b).sum()),
                           'both_reject':int((a&b).sum()),'neither_reject':int((~a&~b).sum()),
                           'paired_rate_difference':float(diff.mean()),'paired_mc_se':float(diff.std(ddof=1)/math.sqrt(len(diff)))})
    summary=pd.DataFrame(rows)
    csv(RESULTS/'MONTE_CARLO_SUMMARY.csv',summary)
    csv(RESULTS/'METHOD_STATUS_COUNTS.csv',statusrows)
    csv(RESULTS/'PAIRED_METHOD_COMPARISONS.csv',paired)
    gates=summary.loc[summary.phase.eq('MAIN')&summary.scenario_id.isin(COMPATIBLE)].copy()
    gates['upper_bound_le_0_075']=gates.wilson95_upper.le(.075)
    gates['failure_le_0_01']=gates.failure_rate.le(.01)
    gates['conditional_size_upper_gate_pass']=gates.upper_bound_le_0_075&gates.failure_le_0_01
    csv(RESULTS/'PREDECLARED_CALIBRATION_GATES.csv',gates)
    aggregate=gates.groupby('method',observed=True).agg(compatible_cells=('scenario_id','size'),
        passed_cells=('conditional_size_upper_gate_pass','sum'),maximum_null_rejection=('rejection_rate_all_replicates','max'),
        maximum_null_upper95=('wilson95_upper','max'),maximum_failure=('failure_rate','max')).reset_index()
    csv(RESULTS/'CALIBRATION_GATE_BY_METHOD.csv',aggregate)
    method_registry()
    mainrows=summary.loc[summary.phase.eq('MAIN')]
    csv(RESULTS/'MAIN_REJECTION_RATE_MATRIX.csv',mainrows.pivot(index='scenario_id',columns='method',values='rejection_rate_all_replicates').reset_index())
    csv(RESULTS/'MAIN_FAILURE_RATE_MATRIX.csv',mainrows.pivot(index='scenario_id',columns='method',values='failure_rate').reset_index())
    save_json(RESULTS/'SUMMARY_COUNTS.json',{'cells':len(grid),'independent_replicates':sum(c['repetitions'] for c in grid),
        'summary_rows':len(summary),'method_result_rows':int(summary.n_replicates.sum()),'clinical_files_opened':0,
        'alpha':.05,'rate_denominator':'All planned replicates, including failures; successful-only rates separately labeled'})
    print(aggregate.to_string(index=False))


if __name__=='__main__':main()
