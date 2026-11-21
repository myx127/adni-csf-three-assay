from __future__ import annotations

import json
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font,PatternFill,Alignment
from openpyxl.utils import get_column_letter
from sim_config import *

SHORT={'PCA_IN_SAMPLE':'PCA same sample','PCA_CROSSFIT_2FOLD':'PCA cross-fit',
       'COMMON_LOADING_MEAN_WALD':'Mean constraint','INVARIANCE_JOINT_IID':'Joint constraint',
       'INVARIANCE_JOINT_TMT_CLUSTER':'Joint / 65 clusters','INVARIANCE_JOINT_3WAY_CLUSTER':'Joint / 3-way clusters',
       'ORACLE_COMPLETE_TRUTH':'Oracle (not available)'}

def read(name):return pd.read_csv(RESULTS/name)

def md(frame):
    def cell(v):return str(v).replace('|','/').replace('\n',' ')
    result=['| '+' | '.join(map(str,frame.columns))+' |','| '+' | '.join(['---']*len(frame.columns))+' |']
    result+=['| '+' | '.join(map(cell,row))+' |' for row in frame.itertuples(index=False,name=None)]
    return '\n'.join(result)

def pct(x):return f'{x:.1%}'

def rate_text(row):return f"{row.rejection_rate_all_replicates:.1%} [{row.wilson95_lower:.1%}, {row.wilson95_upper:.1%}]"

def make_figures(summary):
    main=summary.loc[summary.phase.eq('MAIN')]
    order=[r['scenario_id'] for r in scenarios()];methods=list(SHORT)
    matrix=main.pivot(index='scenario_id',columns='method',values='rejection_rate_all_replicates').reindex(index=order,columns=methods)
    failed=main.pivot(index='scenario_id',columns='method',values='failure_rate').reindex(index=order,columns=methods)
    fig,ax=plt.subplots(figsize=(12.8,10.5))
    data=matrix.to_numpy().copy();data[failed.to_numpy()>.01]=np.nan
    cmap=plt.colormaps['YlOrRd'].copy();cmap.set_bad('#d5d9de')
    im=ax.imshow(data,aspect='auto',vmin=0,vmax=1,cmap=cmap)
    ax.set_xticks(range(len(methods)),[SHORT[x] for x in methods],rotation=30,ha='left')
    ax.xaxis.tick_top()
    ax.set_yticks(range(24),order)
    for i in range(24):
        for j in range(len(methods)):
            text=f'F {failed.iloc[i,j]:.0%}' if failed.iloc[i,j]>.01 else f'{matrix.iloc[i,j]:.1%}'
            ax.text(j,i,text,ha='center',va='center',fontsize=9,
                    color='white' if np.isfinite(data[i,j]) and data[i,j]>.65 else 'black')
    ax.axhline(19.5,color='#276c95',lw=2);ax.axhline(22.5,color='#276c95',lw=2)
    ax.set_title('Synthetic pathology association / constraint rejection rates | n=571, R=2,000 per scenario',pad=95,fontsize=12)
    fig.colorbar(im,ax=ax,fraction=.025,pad=.025,label='Rejections / all planned replicates')
    fig.text(.08,.025,'N01-N20: no extra pathology-dependent response; some violate the linear independent-error assumptions.\nA01-A03: relative response effects. U01: direct effect parallel to the shared direction (not identifiable).\nGrey cells: >1% method failure; F is failure rate, not power. PCA rejection is coordinate association, not proof of a new mechanism.\nOracle uses true latent values and full unselected data; those inputs do not exist in the current study.',fontsize=9)
    fig.subplots_adjust(left=.09,right=.94,bottom=.13,top=.78)
    for ext in ['png','pdf']:fig.savefig(FIGURES/f'01_MAIN_REJECTION_AND_FAILURE.{ext}',dpi=190)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,5.8))
    colors=['#216b9c','#c76237','#777777']
    for method,color in zip(['COMMON_LOADING_MEAN_WALD','INVARIANCE_JOINT_IID','ORACLE_COMPLETE_TRUTH'],colors):
        q=summary.loc[summary.phase.eq('EFFECT_CURVE')&summary.method.eq(method)].sort_values('offset_effect_override')
        x=q.offset_effect_override.to_numpy();y=q.rejection_rate_all_replicates.to_numpy()
        ax.plot(x,y,'o-' if method!='ORACLE_COMPLETE_TRUTH' else 'o--',color=color,label=SHORT[method])
        ax.fill_between(x,q.wilson95_lower.to_numpy(),q.wilson95_upper.to_numpy(),color=color,alpha=.12)
    ax.axhline(.05,color='.6',ls=':',label='Nominal 0.05')
    ax.set(xlabel='Additional M2 offset per unit P (synthetic log units)',ylabel='Rejection proportion',ylim=(0,1.04),
           title='Relative offset effect curve | n=571, R=1,000\nUnequal error SDs; conditional independent-error model')
    ax.grid(alpha=.15);ax.legend(frameon=False,fontsize=9)
    fig.text(.12,.025,'Bands: Wilson 95% Monte Carlo intervals. This is not clinical effect-size or power calibration.',fontsize=8.5)
    fig.subplots_adjust(bottom=.17,top=.84)
    for ext in ['png','pdf']:fig.savefig(FIGURES/f'02_OFFSET_EFFECT_CURVE.{ext}',dpi=190)
    plt.close(fig)
    family=read('OBSERVATIONALLY_EQUIVALENT_PARAMETER_FAMILY.csv').sort_values('error_correlation12')
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.6))
    for ax,field,label in zip(axes,['relative_loading3','direct_pathology_effect3'],['Relative M3 loading','Direct P-to-M3 coefficient']):
        ax.plot(family.error_correlation12,family[field],color='#286f8d',lw=2)
        ax.set(xlabel='Assumed M1/M2 residual error correlation',ylabel=label)
        ax.grid(alpha=.2)
    axes[0].axhline(1,color='.5',ls='--');axes[1].axhline(0,color='.5',ls='--')
    fig.suptitle('Identical observed conditional means and covariance for every point',fontsize=12)
    fig.text(.08,.025,'A constructed equivalence family, not an ADNI estimate or confidence interval.\nChanging the latent/error decomposition changes the inferred direct effect while leaving the observable distribution unchanged.',fontsize=9)
    fig.subplots_adjust(bottom=.25,top=.85,wspace=.35)
    for ext in ['png','pdf']:fig.savefig(FIGURES/f'03_OBSERVATIONAL_EQUIVALENCE.{ext}',dpi=190)
    plt.close(fig)
    recovery=read('PARAMETER_RECOVERY_REPLAY.csv')
    fig,axes=plt.subplots(1,2,figsize=(10.5,5.3))
    for ax,cases,title in [(axes[0],['N03','N05','N07'],'Relative M3 loading: strong common signal'),
                           (axes[1],['N19'],'Weak common signal; all values retained')]:
        arrays=[recovery.loc[recovery.scenario_id.eq(s),'relative_loading3_estimate'].to_numpy() for s in cases]
        ax.boxplot(arrays,labels=cases,showfliers=True,flierprops={'markersize':2,'alpha':.4})
        ax.axhline(1,color='#b65e35',ls='--');ax.set(title=title,ylabel='Estimated loading ratio (true = 1)')
        ax.grid(axis='y',alpha=.15)
    axes[1].set_yscale('symlog',linthresh=1)
    fig.text(.08,.025,'Replays of original main-experiment seeds; no new independent replicates.\nN07 violates the estimator\'s uncorrelated-error assumption. Weak-signal ratios are not trimmed or winsorized.',fontsize=9)
    fig.subplots_adjust(bottom=.20,top=.88,wspace=.3)
    for ext in ['png','pdf']:fig.savefig(FIGURES/f'04_PARAMETER_RECOVERY.{ext}',dpi=190)
    plt.close(fig)


def key_table(summary):
    main=summary.loc[summary.phase.eq('MAIN')]
    rows=[]
    for case in scenarios():
        q=main.loc[main.scenario_id.eq(case['scenario_id'])].set_index('method')
        row={'情景':case['scenario_id'],'情景说明':case['label_zh']}
        for method,label in [('PCA_IN_SAMPLE','PCA坐标关联'),('PCA_CROSSFIT_2FOLD','交叉拟合PCA坐标关联'),
                             ('COMMON_LOADING_MEAN_WALD','共同方向均值约束'),('INVARIANCE_JOINT_IID','联合稳定性约束')]:
            row[label]=rate_text(q.loc[method])
        row['第三平台聚类拒绝率']=pct(q.loc['INVARIANCE_JOINT_TMT_CLUSTER','rejection_rate_all_replicates'])
        row['三向聚类失败率']=pct(q.loc['INVARIANCE_JOINT_3WAY_CLUSTER','failure_rate'])
        row['完整人数中位数']=int(q.loc['INVARIANCE_JOINT_IID','complete_n_median'])
        rows.append(row)
    table=pd.DataFrame(rows);csv(RESULTS/'KEY_RESULTS_ZH.csv',table);return table


def workbook(summary,keys):
    notes=pd.DataFrame([
        ('状态','67,000次人工重复、603,000行方法结果已完成。模拟完成不等于正式模型已获准。'),
        ('真实数据边界','仅使用Stage1聚合人数和批次规模；没有读取个体蛋白、临床/病理或restricted数据。'),
        ('平台含义','M1/M2/M3与人工SD组合不代表三种真实平台的精度排序。'),
        ('主要拒绝比例','拒绝数/所有计划重复；失败保留在分母。另列成功条件下比例，不能混作功效。'),
        ('PCA口径','PCA的D-P关联是真实构造坐标关联，不能称为对统计零假设的已证实错误检验；它不证明额外分子效应。'),
        ('约束方法','依赖固定线性共同方向、条件不相关误差等。压力情景拒绝不必然意味着病理改变检测关系。'),
        ('连续P','均值约束直接使用连续P；三组约束使用事前固定三分位切点，仅为分段诊断，未实现全自由MNLFA。'),
        ('Oracle','使用已知T和真实固定响应及未删失值，仅说明额外信息的作用；不可部署。HC3有限样本亦非完美校准。'),
        ('三向聚类','负/奇异方差矩阵明确记失败，未做PSD修补。低的无条件拒绝比例不表示该方法可用。'),
        ('参数恢复','8,000次为主实验种子重放，不能加算到67,000个独立重复中。'),
        ('等价参数族','同一观测分布可对应不同直接项；所画曲线是构造族，不是ADNI估计、置信区间或完整识别集合。'),
        ('后续','当前可保留条件性的测量关系检验，需补识别信息/约束并解决有限簇校准，尚不冻结机制推断SAP。'),
    ],columns=['项目','说明'])
    qa=json.loads((RESULTS/'SIMULATION_VALIDATION_QA.json').read_text(encoding='utf-8'))
    qa_table=pd.DataFrame(qa['checks']);qa_table['detail']=qa_table.detail.map(lambda v:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else str(v))
    sheets=[('阅读说明',notes),('关键结果中文',keys),('全部MonteCarlo汇总',summary),
            ('校准上限检查',read('CALIBRATION_GATE_BY_METHOD.csv')),('校准逐情景',read('PREDECLARED_CALIBRATION_GATES.csv')),
            ('截距效应曲线',summary.loc[summary.phase.eq('EFFECT_CURVE')]),('小样本',summary.loc[summary.phase.eq('N_GRID')]),
            ('连续病理变量',summary.loc[summary.phase.eq('CONTINUOUS')]),('失败状态计数',read('METHOD_STATUS_COUNTS.csv')),
            ('配对方法比较',read('PAIRED_METHOD_COMPARISONS.csv')),('参数恢复汇总',read('PARAMETER_RECOVERY_SUMMARY.csv')),
            ('参数恢复种子重放',read('PARAMETER_RECOVERY_REPLAY.csv')),('观测等价参数族',read('OBSERVATIONALLY_EQUIVALENT_PARAMETER_FAMILY.csv')),
            ('情景定义',read('SCENARIO_REGISTRY.csv')),('候选方法',read('METHOD_REGISTRY.csv')),('实验与种子',read('EXPERIMENT_GRID.csv')),
            ('输入聚合哈希',read('INPUT_AGGREGATE_HASH_CHECKS.csv')),('实现验证',qa_table)]
    wb=Workbook();wb.remove(wb.active);schema=[]
    for name,frame in sheets:
        if 'seed' in frame.columns:
            frame=frame.copy();frame['seed']=frame.seed.map(str)
        ws=wb.create_sheet(name);ws.append(list(frame.columns))
        for row in frame.itertuples(index=False,name=None):
            values=[]
            for value in row:
                if pd.isna(value):values.append(None)
                elif isinstance(value,np.generic):values.append(value.item())
                else:values.append(value)
            ws.append(values)
        ws.freeze_panes='C2' if name=='关键结果中文' else 'A2';ws.auto_filter.ref=ws.dimensions;ws.sheet_view.zoomScale=85
        ws.row_dimensions[1].height=42
        for cell in ws[1]:
            cell.font=Font(bold=True,color='FFFFFF');cell.fill=PatternFill('solid',fgColor='23536A');cell.alignment=Alignment(wrap_text=True,vertical='center')
        for j,col in enumerate(frame.columns,1):
            width=18 if pd.api.types.is_numeric_dtype(frame[col]) else min(70,max(20,len(str(col))*1.1))
            if name=='阅读说明':width=25 if j==1 else 110
            if name=='关键结果中文':width=27 if j==2 else 31
            ws.column_dimensions[get_column_letter(j)].width=width
            if pd.api.types.is_float_dtype(frame[col]):
                fmt='0.0%' if any(x in str(col) for x in ['rate','coverage','wilson','upper95','maximum_null_rejection','maximum_failure']) else '0.0000'
                for cells in ws.iter_cols(min_col=j,max_col=j,min_row=2):
                    for cell in cells:cell.number_format=fmt
        if name=='阅读说明':
            for row in ws.iter_rows(min_row=2):
                ws.row_dimensions[row[0].row].height=45
                for cell in row:cell.alignment=Alignment(wrap_text=True,vertical='top')
        schema.append({'sheet':name,'rows':len(frame),'columns':len(frame.columns)})
    wb.properties.title='候选模型识别模拟';wb.properties.subject='全部个体为人工生成，未冻结正式测量模型'
    wb.save(ROOT/'候选模型识别模拟总表.xlsx');csv(RESULTS/'WORKBOOK_SCHEMA.csv',schema)


def write_reports(summary,keys):
    main=summary.loc[summary.phase.eq('MAIN')].set_index(['scenario_id','method'])
    def r(s,m):return rate_text(main.loc[(s,m)])
    proof=json.loads((RESULTS/'EXACT_IDENTIFICATION_AND_EQUIVALENCE.json').read_text(encoding='utf-8'))
    gates=read('CALIBRATION_GATE_BY_METHOD.csv')
    gates['校准通过格数']=gates.passed_cells.astype(str)+'/'+gates.compatible_cells.astype(str)
    small=gates[['method','校准通过格数','maximum_null_rejection','maximum_failure']].copy()
    small.maximum_null_rejection=small.maximum_null_rejection.map(pct);small.maximum_failure=small.maximum_failure.map(pct)
    recovery=read('PARAMETER_RECOVERY_SUMMARY.csv')
    table=keys.loc[keys['情景'].isin(['N02','N03','N05','N09','N14','N16','N17','N18','N19','A01','A02','U01']),
                   ['情景','情景说明','PCA坐标关联','交叉拟合PCA坐标关联','共同方向均值约束','联合稳定性约束']]
    power=summary.loc[summary.phase.eq('EFFECT_CURVE')&summary.method.isin(['COMMON_LOADING_MEAN_WALD','INVARIANCE_JOINT_IID'])].pivot(index='offset_effect_override',columns='method',values='rejection_rate_all_replicates').reset_index()
    for col in power.columns[1:]:power[col]=power[col].map(pct)
    power.columns=['人工M2截距偏移','共同方向均值约束','联合稳定性约束']
    report=f'''# 候选模型的识别模拟报告

状态：**预定模拟及交付完成；尚无模型获得不依赖识别假设的正式机制推断许可。** 执行日期：{datetime.now().date()}。

完成24个主情景、49个实验格、**67,000个独立人工重复、9条方法/诊断路径、603,000行逐重复结果**。另对原有种子重放8,000次检查参数恢复，不计为新增独立重复。所有个体、病理P、共同信号T及蛋白读数均为人工生成。

最有价值的进展是：均值及协方差约束方法避开了本轮“不等精度共同信号”情景中的PCA几何问题，并能检出特定相对偏移。最关键的边界是：**三个单次读出不能无条件区分真实共同信号与额外平台效应；相关误差、固定非线性和技术混杂仍可改变解释。**

## 1. 范围、预先规定及真实数据使用

方案[SIMULATION_PROTOCOL_v1.0.md](SIMULATION_PROTOCOL_v1.0.md)在本目录定量模拟前写定，事前SHA见[访问与方案记录](results/ACCESS_AND_PROTOCOL_RECORD.json)。主情景n=571来自Stage1聚合人数；人工批次规模分别为19/36/65，沿用聚合人数分布，但参与者对应关系全部重新生成。只读取聚合AUDIT_SUMMARY、批次人数及封存清单用于核验，未打开个体蛋白、restricted或临床/病理文件。

主实验每情景2,000次；小样本n=200/400每格500次；连续P和截距效应曲线每格1,000次。根种子20261027，每格及每次重复可独立重放。平台SD=[0.25,0.8,2.0]是人工压力参数，不能当作NULISA/SomaScan/TMT真实精度排名；三组比例不是实际ADNI A/T覆盖。

显著性0.05，比例的分母为全部计划重复，失败未被静默删除；方括号为Wilson 95% Monte Carlo区间。R=2,000且概率0.05时MC SE约0.49个百分点。检验比较的是各自定义的统计对象，不能把所有N情景的拒绝一律称为统计一类错误。

## 2. 实际比较的模型

- **PCA及两折交叉拟合PCA**：将两个正交不一致坐标对P与PC1回归，HC3联合Wald。交叉拟合诊断训练泄漏，不能改变总体共同方向错配；此HC3也不涵盖所有训练步骤的不确定性。
- **共同方向均值约束（2df）**：从Y~P残差的非对角协方差估计相对共同方向，检验三平台P斜率是否与之平行；方向估计的不确定性进入sandwich。
- **联合稳定性约束（8df）**：三组均值差应与共同方向平行，组内非对角协方差应只按共同方差比例变化。允许共同均值/方差及误差方差随组改变；不要求潜在信号严格正态。均值与协方差的4df检验只作组成诊断，不另增发现家族。
- **一向/三向聚类版本**：在同一8df约束下改变影响函数协方差估计；负或奇异矩阵明确标失败，没有PSD修补。
- **Oracle基准**：使用当前真实数据没有的T、真实固定响应、完整未删失读数及已知技术扰动。它不是可部署方法；其有限样本HC3检验亦不是完美校准标准。

这里的约束方法是候选可观测矩检验，**没有把它称为已经验证的全自由MNLFA，也没有估出个体的真实蛋白丰度**。将所有直接效应、相关误差及潜在分布参数同时放自由的模型，先由结构识别证明判定，不用“优化收敛”替代识别。

## 3. 主结果

{md(table)}

完整24情景及失败率见[关键中文表](results/KEY_RESULTS_ZH.csv)和[全部汇总](results/MONTE_CARLO_SUMMARY.csv)。

**共同信号＋不等精度（N03）**：没有任何病理到平台的直接项，PCA坐标关联比例为{r('N03','PCA_IN_SAMPLE')}，交叉拟合为{r('N03','PCA_CROSSFIT_2FOLD')}。均值约束为{r('N03','COMMON_LOADING_MEAN_WALD')}，联合约束为{r('N03','INVARIANCE_JOINT_IID')}。PCA结果针对构造坐标的关联可以是真的，但不能据此解释为新的分子过程。

**真实相对偏移**：A01的M2额外0.35P，均值约束检出率{r('A01','COMMON_LOADING_MEAN_WALD')}，联合约束{r('A01','INVARIANCE_JOINT_IID')}。A02的M3额外0.40P*T，均值约束只有{r('A02','COMMON_LOADING_MEAN_WALD')}，联合约束为{r('A02','INVARIANCE_JOINT_IID')}。所以均值检验不能覆盖载荷变化，8df方法亦未在这一人工载荷效应下达到80%检出率。

**固定响应与处理规则**：N09的二次响应不随P改变，但联合线性约束拒绝率为{r('N09','INVARIANCE_JOINT_IID')}；N14低值地板替代为{r('N14','INVARIANCE_JOINT_IID')}。N10固定饱和响应在本组参数下未明显增加拒绝，不能据此宣称该方法对任意非线性稳健。

**相关误差**：N07/N08在当前人工SD及相关强度下，两种约束方法没有出现很高的拒绝率；这只是有限情景表现，不能验证误差不相关假设。下面的精确等价构造表明，即便无限样本也存在解释不可区分的问题。

**缺失和批次**：30%/70%单平台LOD删除在本组参数下主要降低完整人数和检出能力，并非必然增加拒绝；N18选择性缺失使均值/联合约束分别达到{r('N18','COMMON_LOADING_MEAN_WALD')}和{r('N18','INVARIANCE_JOINT_IID')}。未校正plex混杂N16的联合约束为{r('N16','INVARIANCE_JOINT_IID')}，它发现的是技术混杂下的观测关系差异，不能被解释成病理特异分子效应。

## 4. 校准、聚类及失败

事前规定对7个兼容线性、独立误差零情景逐格检查：拒绝比例的Wilson上限≤7.5%、失败≤1%。这是上限检查，不是证明精确5%校准，也不能奖励没有功效的保守方法。

{md(small)}

共同方向均值约束7/7通过，联合IID约束7/7通过，但后者总体偏保守。一向聚类版本7/7未通过预定上限检查；N03即使误差独立，其拒绝比例仍为{r('N03','INVARIANCE_JOINT_TMT_CLUSTER')}，说明简单使用渐近chi-square加65簇sandwich还不足以完成校准。

三向聚类在24个主情景的失败率为86.6%–91.0%，N03为{pct(main.loc[('N03','INVARIANCE_JOINT_3WAY_CLUSTER'),'failure_rate'])}。主要原因是容斥方差矩阵非正定，另有奇异矩阵。**这表示本轮“19/36/65批次＋8df容斥sandwich＋渐近Wald”实现不适用，不能外推为所有多向聚类方法均不适用。**其约5%–12%的无条件拒绝率被大量失败压低，不可叫校准良好或功效不足。

Oracle也有部分兼容情景未通过上限检查，反映有限样本HC3联合Wald的近似限制；没有将它当作统计检验已完美校准的真值算法。

## 5. 严格的结构识别反例

构造A：载荷(1,1,1)，共同残差方差1，独立误差方差(0.4,0.6,1)，共同P效应0.55，没有直接平台项。

构造B：载荷(1,1,1.25)，共同残差方差0.8，误差方差(0.6,0.8,0.75)，M1/M2误差协方差0.2（相关约0.289），M3直接P项−0.1375。

**两者对每个P产生完全相同的Y条件均值和协方差；正态构造下就是完全相同的条件分布。**协方差最大差{proof['observational_equivalence']['maximum_covariance_difference']:.3g}，P斜率最大差{proof['observational_equivalence']['maximum_slope_difference']:.3g}。两组误差协方差均正定。一个“无额外直接项”，另一个“有额外直接项”，但没有任何只用这些观测数据的检验能将它们区分。

即使载荷已知，beta=lambda*alpha+delta也有alpha→alpha+c、delta→delta−c*lambda的等价变换；3个观测斜率对应4个未知量，Jacobian秩为3。U01的三个平台平行直接项因此可被共同均值变化吸收，不能要求观测约束方法把它检出为相对偏离。Oracle能检出是因为额外知道T，不代表现实可做到。

可复核证据：[精确识别构造](results/EXACT_IDENTIFICATION_AND_EQUIVALENCE.json)、[等价参数族](results/OBSERVATIONALLY_EQUIVALENT_PARAMETER_FAMILY.csv)、[图3](figures/03_OBSERVATIONAL_EQUIVALENCE.png)。参数族不是ADNI估计、置信区间或完整识别集合。

## 6. 约束检验可用，不等于丰度/可靠性参数可稳定恢复

主实验种子重放显示，N03虽有良好的约束检验校准，直接矩分解仍有20.95%的重复出现至少一个负残差方差；N07为55.0%。这不等于约束检验本身失败，而是不能把该矩分解的个体分数或可靠性估计直接当真。

弱信号N19有40.55%的重复出现非正的汇总跨平台协方差，载荷比的RMSE超过45；联合约束拒绝率为0也不能被称为测量关系稳定。必须保留“信号太弱、无法有力评价”的类别。

参数诊断未截尾、未删异常估计，全部8,000次重放见[参数恢复原表](results/PARAMETER_RECOVERY_REPLAY.csv)；这些不是新增独立模拟次数。详见[补充诊断范围](PARAMETER_DIAGNOSTIC_SCOPE.md)。

## 7. 效应强度、小样本和连续P

以下仅为人工M2截距效应曲线，平台噪声、共同水平和n固定；它不是临床效应大小或ADNI实际功效估计。

{md(power)}

连续P补充实验使用事前固定三分位切点作联合约束诊断，没有按结果优化分组，也没有实现全自由连续MNLFA。小样本n=200/400、连续P的全部结果单独列入Excel及汇总CSV，不以它们选择一个最有利结果代替主实验。

## 8. 对论文设计的具体建议

1. **保留“病理相关相对测量关系”作为可检验目标。**将共同方向均值约束与联合稳定性约束作为条件候选；PCA及D保留描述作用，不能独自支撑“已经分离真丰度与独立病理机制”的主张。
2. **先修正有限簇推断。**评估在预先限定零模型下的bootstrap/有限簇校准，并用新种子独立验证；不能把本轮失败的聚类Wald直接用于真实论文推断。
3. **明确误差相关、非线性、左删失的识别限制。**合理范围内结论不稳定者保留为不确定。不能把误差相关全部放自由后宣称模型已识别；强先验/正则化也不是新增观测证据。
4. **增加可核验的信息。**优先取得匹配肽段、同一样本技术复测或有依据的检测响应信息；多适配体须保留同平台依赖，不能当独立重复。对正式候选的分子机制解释需要这些额外证据。
5. **独立病理锚与功效仍待后续。**本轮没有真实A/T、PET、随访及协变量覆盖，因此不冻结蛋白集，也不声称真实研究有80%功效；稳定关系的阳性证据还需要预设等效界值和足够精度，不以P>0.05替代。

下一步方法开发应围绕“有限簇校准＋误差相关/非线性敏感性＋可补充识别信息”推进，而不是继续将旋转后的残差自动生物学化。若额外信息无法取得，论文仍可报告透明的相对测量关系及其病理相关性，但须限定解释。

## 9. 交付与验证

- [Excel总表](候选模型识别模拟总表.xlsx)：18个工作表，含全部汇总、失败、参数恢复、种子与验证。
- [逐重复Parquet目录](results/replicates)：49文件共603,000行；包含每方法每重复P值、统计量、状态和完整人数，无真实参与者信息。
- [独立实现核验](results/SIMULATION_VALIDATION_QA.json)：约束Jacobian、总体零约束、OLS/HC3、聚类容斥公式、等价构造及种子重放共99项通过。种子重放共享生成器，只证明可复现；不冒充独立科学验证。
- [图1](figures/01_MAIN_REJECTION_AND_FAILURE.png)、[图2](figures/02_OFFSET_EFFECT_CURVE.png)、[图3](figures/03_OBSERVATIONAL_EQUIVALENCE.png)、[图4](figures/04_PARAMETER_RECOVERY.png)，均附PDF。
- [输入输出SHA-256](SIMULATION_INPUT_OUTPUT_SHA256.csv)、[交付核验](results/DELIVERY_COMPLETION_QA.json)、[一键重跑](code/run_all.ps1)。重跑数值应相同，日志、Excel/PDF创建时间及其哈希会改变。

方法背景：[MNLFA及识别限制](https://doi.org/10.1037/met0000077)、[固定非线性](https://pubmed.ncbi.nlm.nih.gov/16221030/)、[相关误差模拟](https://pmc.ncbi.nlm.nih.gov/articles/PMC9979199/)、[多向聚类](https://doi.org/10.1198/jbes.2010.07136)。本课题矩约束程序及其实际表现是本轮模拟的结果，不是声称这些文献已经验证其在ADNI中的适用性。
'''
    (ROOT/'MODEL_IDENTIFICATION_SIMULATION_REPORT_ZH.md').write_text(report,encoding='utf-8')
    decision='''# 模型选择与冻结前决定

本轮识别模拟已完成。正式机制推断SAP：**不签发**。

| 方法或主张 | 本轮决定 | 下一步条件 |
| --- | --- | --- |
| PCA/正交D作为新分子机制证据 | 不支持 | 保留描述，不能以调整PC1后显著为充分证据 |
| 共同方向均值约束 | 条件候选 | 线性/误差相关/选择性缺失及批次敏感性；不能覆盖所有载荷变化 |
| 均值＋协方差联合约束 | 条件候选 | 保守性与载荷变化功效；固定非线性及相关误差解释边界 |
| 本轮65簇渐近Wald | 不通过校准上限检查 | 有限簇或零模型bootstrap，并用独立新种子验证 |
| 本轮三向容斥sandwich Wald | 不可部署 | 解决高比例非正定/奇异协方差，不能删失败后声称成功 |
| 全自由直接项＋误差相关＋共同信号分解 | 结构不可识别 | 增加独立信息或有依据的限制；强先验不能变成数据识别 |
| 真实丰度或proteoform机制已被分离 | 本轮不支持 | 需匹配的肽段/复测/检测响应及独立验证证据 |
| 最终蛋白集、独立病理终点和临床功效 | 尚未评估/冻结 | 表型覆盖与随访审计、协变量及多重检验层级后再模拟 |

无需等待额外分子信息也可继续开发相对测量关系检验；但结论必须受识别假设约束。本轮未打开临床或个体蛋白数据。
'''
    (ROOT/'MODEL_RECOMMENDATION_AND_FREEZE_GATES.md').write_text(decision,encoding='utf-8')
    readme='''# 候选模型的识别模拟

已完成67,000次独立人工重复、9条方法/诊断路径及精确识别反例。模拟完成，正式机制推断模型尚未冻结。

1. [中文结果报告](MODEL_IDENTIFICATION_SIMULATION_REPORT_ZH.md)
2. [Excel总表](候选模型识别模拟总表.xlsx)
3. [模型选择与冻结条件](MODEL_RECOMMENDATION_AND_FREEZE_GATES.md)
4. [主结果图](figures/01_MAIN_REJECTION_AND_FAILURE.png)及[观测等价图](figures/03_OBSERVATIONAL_EQUIVALENCE.png)
5. [交付核验](results/DELIVERY_COMPLETION_QA.json)与[哈希清单](SIMULATION_INPUT_OUTPUT_SHA256.csv)

均值/协方差约束方法避开了不等精度下的PCA几何问题，但依赖识别假设；现有三个单次读出无法无条件分离共同信号和额外平台效应。聚类推断、非线性及选择性缺失仍需处理。

`results/replicates`保存49个逐重复Parquet文件；`code/run_all.ps1`可重建。没有读取真实个体蛋白或临床/病理数据，没有修改Stage1封存文件。
'''
    (ROOT/'README_阅读顺序.md').write_text(readme,encoding='utf-8')


def main():
    summary=read('MONTE_CARLO_SUMMARY.csv');keys=key_table(summary)
    make_figures(summary);workbook(summary,keys);write_reports(summary,keys)
    (LOGS/'DELIVERY_BUILD.log').write_text(datetime.now().isoformat()+' Reports, workbook and figures built.\n',encoding='utf-8')
    print('Report, 18-sheet workbook and four PNG/PDF figure pairs generated.')


if __name__=='__main__':main()
