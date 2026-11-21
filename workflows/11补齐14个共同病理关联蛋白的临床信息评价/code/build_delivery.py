from config import *
import platform,importlib.metadata
import xlsxwriter

def run():
 formal=pd.read_csv(RESULTS/'FORMAL204_CLINICAL_RESULTS.csv');focus=formal.loc[formal.focus14];audit=pd.read_csv(RESULTS/'ACTUAL_DESIGN_AUDIT.csv');cal=pd.read_csv(RESULTS/'CALIBRATION_QUALIFICATION_DECISIONS.csv');states=pd.read_csv(RESULTS/'EVIDENCE_STATE_COUNTS.csv');science=json.loads((RESULTS/'SCIENTIFIC_QA_SUMMARY.json').read_text(encoding='utf-8'));runinfo=json.loads((RESULTS/'CALIBRATION_RUN_SUMMARY.json').read_text(encoding='utf-8'));figs=pd.read_csv(RESULTS/'FIGURE_INDEX.csv')
 def chosen(branch,test,combined=False,onlyfocus=True):
  d=focus if onlyfocus else formal;d=d.loc[d.branch.eq(branch)&d.test.eq(test)&d['supported_combined204' if combined else 'supported_within_family']]
  return '、'.join(d.gene.tolist()) or '无'
 def grouptext():
  out=[]
  for br in FORMAL_BRANCHES:
   for test in ['SAME','DIFF']:
    d=focus.loc[focus.branch.eq(br)&focus.test.eq(test)]
    out.append(f'| {br} | {test} | {int(d.calibration_qualified.sum())}/14 | {chosen(br,test)} | {chosen(br,test,True)} |')
  return '\n'.join(out)
 a=audit.loc[audit.branch.eq('MEM5Y')&audit.gene.isin(FOCUS14)]
 rows=focus.loc[focus.branch.eq('MEM5Y')&focus.test.eq('SAME')].sort_values('p_computed')
 main_table=rows[['gene','n_people','direction','p_computed','q_BY_within_family','q_BY_combined204','calibration_qualified']].to_markdown(index=False,floatfmt='.6g')
 diffrows=focus.loc[focus.test.eq('DIFF')&focus.supported_within_family]
 diff_table=diffrows[['gene','branch','n_people','p_computed','q_BY_within_family','q_BY_combined204']].to_markdown(index=False,floatfmt='.6g')
 sens=pd.read_csv(RESULTS/'PRESPECIFIED_SENSITIVITY_RESULTS.csv')
 selected_genes=focus.loc[focus.branch.eq('MEM5Y')&focus.test.eq('DIFF')&focus.supported_within_family,'gene'].tolist()
 sens_table=sens.loc[sens.gene.isin(selected_genes)&sens.test.eq('DIFF'),['gene','branch','p_computed','q_BY28_diagnostic']].to_markdown(index=False,floatfmt='.6g')
 failed_designs='、'.join(cal.loc[~cal.qualified,'design_id'].tolist()) or '无'
 nowtext=now()
 report=f'''# 第11阶段：共同病理关联蛋白的临床信息评价

完成时间：{nowtext}。本轮以原34蛋白固定家族完成204项正式问题，14个既有三检测共同pTau信号作为重点解释子集。所有前十阶段和原始来源保持只读。

## 主要结果

14重点蛋白主记忆进展的共同临床关联，完整主68项BY支持：**{chosen('MEM5Y','SAME')}**；检测间标准化临床关联差异支持：**{chosen('MEM5Y','DIFF')}**。
基线记忆共同关联，完整次136项BY支持：**{chosen('MEM_BASE','SAME')}**；执行功能进展共同关联支持：**{chosen('EF5Y','SAME')}**。
全部204项合并BY下，14重点蛋白记忆进展共同关联支持：**{chosen('MEM5Y','SAME',True)}**。
本轮新增的主要正面证据是检测间临床关联系数差异：记忆进展支持{chosen('MEM5Y','DIFF')}，执行功能进展支持{chosen('EF5Y','DIFF')}。这是规定蛋白尺度、调整pTau的条件模型结果，不能解释为某检测更准确或某一分子形式已被识别。

{diff_table}

| 终点 | 问题 | 校准合格 | 完整家族BY支持蛋白 | 合并204项BY支持 |
|---|---|---:|---|---|
{grouptext()}

主记忆进展三检测共同关联逐蛋白结果如下。排序仅用于阅读，分析名单和检验家族均预先固定。

{main_table}

## 本轮补齐了什么

第10阶段只评价NPTX1、KLK6、PDGFRB。本轮扩展到原34家族，重点回答14个共同病理关联信号是否同样具有三检测共同的临床关联。其余20蛋白均保留在正式家族和交付中。原三个候选的重叠设计重新计算并与第10阶段核对；因家族大小改变，q值不能要求与第10阶段相同。
所有结果仍来自同一ADNI脑脊液人群，是发现后的临床延伸，不构成独立人群验证。原M0/APOE/BOTH病理结果仅作为既有证据整合，不在本轮改选。

## 模型、单位与问题边界

主终点是最长五年观察窗口的ADNI_MEM年变化；次终点是基线ADNI_MEM与ADNI_EF年变化。
主纵向模型含每人固定截距，时间乘年龄、性别、教育、APOE、独立Elecsys log2 pTau181、三检测技术因子和Elecsys批次，以及每一项蛋白测量。三方程使用相同病例、共同协变量，保留配对相关。
时间中心化充分统计量与逐访视固定截距OLS等价，不将单人斜率当作无误差结局。
蛋白单位固定为原完整临床病例中、去除人口学和技术因素但不含pTau/新结局的残差SD；跨终点和敏感性保持固定，并传播定标不确定度。
正向认知年变化系数表示条件下更正向的变化/相对较慢下降，不意味着实际认知改善或因果保护。

SAME检验全正或全负三边际关联的并集，使用交并检验及两个方向校正；不是平均系数显著。DIFF是三个标准化临床系数相等的2df配对检验，不是原病理H约束检验。
共同关联、检测间差异、实用等效和真实共同蛋白丰度是不同问题。本轮不建立临床等效、个体共识/不一致评分、因果中介或风险预测模型。

## 多重检验与固定子集

主家族68项BY，次家族136项BY，合并204项BY作为跨主次问题敏感性。14子集摘录的全部q均来自完整家族，不按子集另行缩小正式校正范围。
BY控制目标是所声明的完整家族，不是整个项目，也不额外保证按pTau结果筛出的14子集的条件FDR。
本轮未联合模拟34蛋白全家族的经验FDR；BY的依赖稳健性仍要求输入P有效，有限设计校准只核查预定工作情景。
不合格设计正式P置1并保留家族。单项原始P、点置信区间或单个检测的显著性不替代正式共同/差异结论。

## 样本、随访与信息量

144个实际设计均通过预设几何门槛。14重点蛋白主记忆样本范围为{int(a.n_people.min())}—{int(a.n_people.max())}人。详细逐蛋白完整交集、设计秩和批次见ACTUAL_DESIGN_AUDIT.csv。
五年指观察窗口，真实随访跨度和末次达到4.5年的比例另表交付。时间权重有效人数仅为权重集中度诊断，不是统计学真实有效样本量。
基线和执行功能使用各自实际完整病例；不能从某一个蛋白或终点的人数推断全部设计。

## 有限设计校准与功效

102正式设计每个5情景×100000次；42敏感性设计每个5情景×10000次。共{runinfo['independent_noise_draws']:,}次独立噪声重复，复用37均值模式形成{runinfo['pattern_evaluations_reusing_noise']:,}次模式评价，后者不是独立重复。
校准合格{int(cal.qualified.sum())}/144设计，其中正式{int(cal.loc[cal.formal].qualified.sum())}/102，敏感性{int(cal.loc[~cal.formal].qualified.sum())}/42。
未合格设计：{failed_designs}。当前CST3主记忆设计仅小尾部资格未通过，两项正式P均置1；未移除蛋白或缩小家族。14重点蛋白的全部正式及敏感性设计均合格。
情景覆盖独立访视投影、AR1及个体随机斜率、pTau异方差、技术簇年龄随机斜率、t5及技术斜率。工作方差沿三技术轴求和G/(G−1)簇得分矩阵，称SUM3，不称CGM或CR2。
正式设计核验BY68、BY136、BY204第一步充分阈值的小尾部；敏感性核验BY28阈值。有限情景合格不证明所有误差结构、随机暴露或信息性失访下有效。
工作效应0.02/0.05/0.10记忆分数单位/年/蛋白残差SD不是临床MCID，也不是按观测效应推得的事后功效。充分阈值通过率只提供实际BY排序功效的充分条件下界。

## 固定敏感性与条件化解释

对14重点蛋白完成不调整pTau、高信息随访、个体等权三个分支，84项诊断问题，各分支BY28仅用于诊断。
不调整pTau保持相同病例，其他分支蛋白尺度不变。pTau调整前后系数变化的遗漏变量代数恒等式另表交付；它不是因果中介、保护机制或信号解释比例。
高信息子集和等权模型评价估计对观察安排的依赖，不能自动解决非随机失访。原第10阶段其他敏感性仍限于原三候选。

新增主记忆差异信号的预设敏感性结果：

{sens_table}

三个主记忆差异在不调整pTau时均未获支持。高信息子集下三项均达到预设诊断BY28；个体等权下仅NPTX2达到该诊断阈值。敏感性与主模型回答的条件/权重目标不同，不能用较小诊断家族升级正式结论。主要结论应明确限定在pTau条件模型中；TIMP3未通过合并204项BY，TAFA5等权结果较弱，NPTX2有相对更多跨终点和权重诊断支持，均不称独立复现。

## 结论定界

本轮的结论必须与上面的完整家族q和校准资格对应。发现共同临床关联时，可报告三检测在规定条件模型中共同提供临床关联信息；不能因此宣称测量等效、疾病特异、病因作用或可直接临床预测。
未获支持的蛋白应表述为当前数据和检验范围未确认，不能一概称为无生物学或临床作用。
即使共同病理关联和共同临床关联同时成立，也不能推断病理介导临床效应。
血浆、FDG、全局脑脊液蛋白背景、分子机制及其他队列均不属于本轮执行范围；既有分子与独立验证缺口继续保留。

## 可复核性与结果入口

科学检查{science['passed']}/{science['checks']}通过；全部保存模拟检验P值逐项重数、种子回放及实测设计独立QR/直接残差协方差复核通过。
结果和图见Excel与FIGURE_INDEX.csv；逐人临床和设计矩阵仅保存在00_restricted。前十阶段{science['prior_sealed_entries']}条封存记录在阶段末再次核验。
'''
 (ROOT/'STAGE11_MAIN_REPORT_ZH.md').write_text(report,encoding='utf-8')
 key=f'''# 第11阶段关键决策
1. 固定原34家族，14个共同pTau蛋白作为解释子集；主68项、次136项、合并204项，不缩小正式校正家族。
2. 14重点蛋白记忆进展共同关联BY68支持：{chosen('MEM5Y','SAME')}；差异支持：{chosen('MEM5Y','DIFF')}。
3. 14重点蛋白基线记忆共同关联BY136支持：{chosen('MEM_BASE','SAME')}；执行功能进展共同关联支持：{chosen('EF5Y','SAME')}，检测差异支持：{chosen('EF5Y','DIFF')}。
4. 合并204项BY下重点记忆共同关联支持：{chosen('MEM5Y','SAME',True)}。
5. 14重点主样本{int(a.n_people.min())}—{int(a.n_people.max())}；{int(cal.qualified.sum())}/144设计有限校准合格；5310万独立噪声，37模式复用不另计独立重复。
6. 三个主记忆差异在不调pTau时均未确认；高信息诊断保留三项，等权诊断仅NPTX2保留。临床差异结论限定于规定的pTau条件模型与蛋白尺度，不是检测准确性优劣。
7. 同ADNI事后临床延伸，不是独立验证；共同关联不等于测量等效、因果保护或个体预测用途。
8. 科学检查{science['passed']}/{science['checks']}通过；旧阶段及原始来源只读并已核验。
'''
 (ROOT/'关键决策摘要.md').write_text(key,encoding='utf-8')
 (ROOT/'模型原理与结果解释.md').write_text('''# 模型原理与结果解释
每人先去除认知水平固定截距，再利用真实访视时间估计蛋白与年变化的关联。以时间中心化平方和的平方根加权，可精确约化为每人一行的回归；此权重来自原始OLS估计目标，不是把估计斜率当成准确观测。
三检测分别进入相同病例与相同协变量的方程，病例和技术簇得分联合计算协方差。不能把三个系数当作独立研究做普通荟萃。
共同检验要求三边际关联全部同向成立：正方向P为三单侧P的最大值，负方向同理，再为两个方向取2倍较小值并截断到1。
差异检验将三个系数乘固定蛋白残差SD，并以S−N、T−N构成2df Wald。标准化的SD也是估计量，联合影响函数保留其误差。
三轴SUM3为三套簇得分矩阵外积之和，每套乘G/(G−1)；没有CGM包含排除，也没有CR2偏差校正。资格来自预定有限情景校准，不能用BY补救不合格的单项推断。
原病理H约束检验与这里的临床系数差异检验不同。相同方向允许不同强度；没有差异证据也不能证明等效。
临床系数条件于pTau。增大或符号改变可以由线性条件化解释，不能自动认定保护、中介或分离出某一分子成分。
''',encoding='utf-8')
 (ROOT/'参考文献与来源说明.md').write_text('''# 方法与来源
- Benjamini & Yekutieli, 2001. DOI: https://doi.org/10.1214/aos/1013699998 。BY要求输入P有效，控制范围为声明家族。
- Nichols et al., 2005. DOI: https://doi.org/10.1016/j.neuroimage.2004.12.005 。共同关联的交并检验逻辑；本轮为双方向并集，另作两个方向校正。
- MacKinnon, Nielsen & Webb, 2023. DOI: https://doi.org/10.1016/j.jeconom.2022.04.001 。有限簇推断的局限；本项目SUM3并非该文某一估计器的普遍保证。
- 数据来自已封存第5阶段三检测变换矩阵、第8阶段34/14清单、第10阶段已核实的ADNI认知装配。精确路径和散列见INPUT_SOURCE_MANIFEST.csv。
- UWNPSYCHSUM使用固定2023版ADNI_MEM/EF复合分数，源文件2026年下载；数据处理继承第10阶段已验证的实际日期与冲突处置。
- 既有三个候选的临床结果和14个pTau共同信号已知，本轮明确为发现后延伸。原JTM及其他工程共享数据背景按既有发表边界披露，不复用网络结果作新发现。
''',encoding='utf-8')
 (ROOT/'可复现与运行说明.md').write_text('''# 运行说明
默认运行code/run_all.ps1：已封存则只读核验；未封存才执行完整流程。也可python -B code/verify_delivery.py只读核验。
完整流程：bootstrap → prepare_templates → validate_engine → run_calibration → analyze_real → summarize_and_describe → validate_science → make_figures → build_delivery → validate_delivery → seal_delivery。
新复现可用python -B code/reproduce_in_copy.py，先创建新的同级目录，仅复制代码和方案，再冷启动完整分析。不会删除或覆盖原目录。源1—10阶段必须保持原结构和散列。
计算为102正式设计×5×100000，加42敏感性设计×5×10000；8进程、每进程BLAS单线程。种子由固定根种子、设计ID和情景确定；按文件散列及运行签名检查恢复点。完整模拟检验P存parquet，边际组件只保存汇总计数。
复现软件依赖见results/ENVIRONMENT.json；相同版本、批大小和硬件数值条件下可重放模拟。日期、运行日志和交付文件时间戳不是数值重现目标。
逐人资料在00_restricted，不应放入公开论文附件。
''',encoding='utf-8')
 english=f'''# Suggested methods and results wording
This post-discovery extension evaluated a fixed family of 34 CSF proteins, with 14 previously supported same-direction pTau-associated proteins designated as an interpretation subset. Prior clinical results for NPTX1, KLK6 and PDGFRB were known. This was neither an independent validation cohort nor a wholly prospective study.

The primary endpoint was ADNI-MEM change within an observation window of up to five years. The 68 primary hypotheses comprised a three-assay same-direction intersection-union test and a two-degree-of-freedom test of equality of standardized marginal clinical coefficients for each protein. Baseline ADNI-MEM and longitudinal ADNI-EF formed a combined secondary family of 136 hypotheses. BY adjustment was applied to the full primary and secondary families; all 204 hypotheses were additionally adjusted together. The 14-protein subset retained these full-family q values.

Longitudinal estimation used person fixed intercepts and actual visit times, with an exact within-person sufficient-statistic reduction. Models adjusted age, sex, education, APOE, independent Elecsys pTau181 and all assay technical factors; nuisance effects acted on time slopes. Paired covariance retained dependence across the three marginal models. Exposure units were fixed residual standard deviations estimated without the new outcome or pTau, with their uncertainty propagated jointly. Three-axis SUM3 working covariance was qualified under prespecified finite scenarios rather than claimed as universally valid.

Within the focus subset, primary-family BY-supported common memory-change associations: {chosen('MEM5Y','SAME')}. Primary paired clinical-coefficient differences: {chosen('MEM5Y','DIFF')}. Secondary-family supported common baseline-memory associations: {chosen('MEM_BASE','SAME')}; common executive-change associations: {chosen('EF5Y','SAME')}.
Secondary-family supported executive-change coefficient differences: {chosen('EF5Y','DIFF')}. The primary memory-coefficient differences were not supported in the same-case models omitting pTau; diagnostic equal-person weighting retained NPTX2 alone. The findings are therefore explicitly conditional and do not establish assay superiority.

Same-direction clinical association does not establish assay equivalence, true protein abundance, mediation, causal protection or individual predictive utility. Nonsignificant clinical-coefficient differences likewise do not establish equivalence. All conclusions remain conditional on the stated observational models and missing-data patterns.
'''
 (ROOT/'论文方法与结果建议_EN.md').write_text(english,encoding='utf-8')
 # Workbook: compact tables plus explicit CSV-only large calibration/recount tables.
 preferred=['EVIDENCE_STATE_COUNTS','FOCUS14_PRIMARY28_FROM_FULL_FAMILY','FOCUS14_SECONDARY56_FROM_FULL_FAMILY','FOCUS14_INTEGRATED_EVIDENCE','FOCUS14_PLATFORM_EFFECTS','FOCUS14_PAIRED_CONTRASTS','PRIMARY68_MEMORY_PROGRESSION_RESULTS','SECONDARY136_RESULTS','FORMAL204_CLINICAL_RESULTS','INTEGRATED_PATHOLOGY_CLINICAL_FAMILY34','FROZEN_FAMILY34_AND_FOCUS14','ACTUAL_DESIGN_AUDIT','EXPOSURE_SCALING_AUDIT','FIXED_FORMAL_HYPOTHESIS_REGISTRY','PRESPECIFIED_SENSITIVITY_RESULTS','FOCUS14_PTAU_CONDITIONING_DECOMPOSITION','FOLLOWUP_INFORMATION_FAMILY34','MAIN_CASE_PROFILE_FAMILY34','CALIBRATION_QUALIFICATION_DECISIONS','NULL_CALIBRATION_WORST_SUMMARY','PRIMARY_WORKING_EFFECT_POWER','STAGE10_OVERLAPPING_DESIGN_REPRODUCTION','FIGURE_INDEX','INPUT_SOURCE_MANIFEST','SCIENTIFIC_QA_CHECKS','ENGINE_INDEPENDENT_CHECKS','REAL_CLUSTER_INFLUENCE_DIAGNOSTICS','SIMULATION_FILE_INVENTORY','COMPONENT_NULL_CALIBRATION','CONDITIONAL_POWER_WORST_DGP','SIMULATION_MEAN_PATTERN_GRID','INDEPENDENT_SEED_REPLAY_CHECKS','CALIBRATION_AND_POWER_FULL','INDEPENDENT_SIMULATION_RECOUNTS']
 entries=[];book=ROOT/'14个共同病理关联蛋白_临床信息评价_结果与审计.xlsx'
 with pd.ExcelWriter(book,engine='xlsxwriter') as writer:
  wb=writer.book;head=wb.add_format({'bold':True,'bg_color':'#DCEAF7','text_wrap':True,'border':1});num=wb.add_format({'num_format':'0.000000'});sci=wb.add_format({'num_format':'0.0000E+00'})
  for i,name in enumerate(preferred,1):
   f=RESULTS/(name+'.csv');d=pd.read_csv(f);large=d.size>180000
   sheet=f'T{i:02d}_'+name[:25] if not large else 'CSV_ONLY'
   entries.append({'table_id':i,'csv':f.name,'rows':len(d),'columns':len(d.columns),'sheet':sheet,'delivery':'CSV_ONLY' if large else 'EXCEL_AND_CSV'})
   if large:continue
   d.to_excel(writer,sheet_name=sheet,index=False);ws=writer.sheets[sheet];ws.freeze_panes(1,0);ws.autofilter(0,0,len(d),len(d.columns)-1)
   for j,c in enumerate(d.columns):
    ws.write(0,j,c,head);fmt=sci if any(k in c.lower() for k in ['p_computed','q_by','p_for','p_saved','p_replayed']) else num if pd.api.types.is_float_dtype(d[c]) else None
    ws.set_column(j,j,min(42,max(13,len(c)*.8)),fmt)
   ws.set_row(0,38)
  pd.DataFrame(entries).to_excel(writer,sheet_name='INDEX',index=False);ws=writer.sheets['INDEX'];ws.freeze_panes(1,0);ws.set_column(0,0,10);ws.set_column(1,1,68);ws.set_column(2,3,12);ws.set_column(4,5,35)
 csvout(RESULTS/'WORKBOOK_SHEET_INDEX.csv',entries)
 env={'time':now(),'python':sys.version,'platform':platform.platform(),'packages':{n:importlib.metadata.version(n) for n in ['numpy','pandas','scipy','matplotlib','openpyxl','xlsxwriter','pyarrow']}}
 save_json(RESULTS/'ENVIRONMENT.json',env)
 (ROOT/'README.md').write_text('''# 第11阶段交付入口
- [关键决策摘要](关键决策摘要.md)
- [完整中文报告](STAGE11_MAIN_REPORT_ZH.md)
- [结果与审计Excel](14个共同病理关联蛋白_临床信息评价_结果与审计.xlsx)
- [冻结方案](STAGE11_PROTOCOL_v1.0.md)
- [模型原理](模型原理与结果解释.md)
- [运行说明](可复现与运行说明.md)
结果以原34家族校正，14共同病理信号为重点子集。全部正式、敏感性及不合格结果保留。封存完成标记与文件清单见STAGE11_COMPLETION_STATUS.json和SHA256_MANIFEST.csv。
''',encoding='utf-8')
 (ROOT/'字段说明与结果导航.md').write_text('''# 字段说明
gene/accession为冻结映射；focus14为既有共同pTau子集。
N=NULISA，S=SomaScan，T=TMT-MS。
SAME是三检测全正或全负共同临床关联；DIFF是固定残差SD尺度临床系数相等的2df检验。
p_computed是工作计算P；p_for_family在资格不合格或数学无效时为1。
q_BY_within_family对应主68或次136家族；q_BY_combined204对应全部正式问题。敏感性q_BY28_diagnostic仅用于诊断。
standardized_effect为ADNI复合分数单位/年/蛋白残差SD，基线模型则无/年。不是浓度、相关系数、HR或OR。
standardized_SE_full_joint包含回归与定标的联合不确定度。CI95_pointwise为点区间，不是所有蛋白/分支的同时区间。
qualified是固定实际设计有限情景资格。nominal_pass/tail_pass/marginal_pass分别为名义、小尾部和零边际组件检查。
模拟pattern的37次复用不独立；BY充分阈值通过率不等于实际BY排序功效。全部表的行列与Excel索引见WORKBOOK_SHEET_INDEX.csv。
''',encoding='utf-8')
 log(f'Chinese report, decision memo and workbook built: {len(entries)} indexed tables.')
if __name__=='__main__':assert_unsealed();run()

