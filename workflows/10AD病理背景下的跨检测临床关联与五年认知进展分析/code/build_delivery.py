from config import *
import math,platform,importlib.metadata
from openpyxl import Workbook
from openpyxl.styles import Font,PatternFill,Alignment
from openpyxl.utils import get_column_letter

def mdtable(d):
 def fmt(v):
  if v is None or isinstance(v,float) and not np.isfinite(v):return '—'
  if isinstance(v,(bool,np.bool_)):return '是' if v else '否'
  if isinstance(v,(float,np.floating)):return f'{v:.6g}'
  return str(v).replace('|',' / ').replace('\n',' ')
 return '\n'.join(['| '+' | '.join(map(str,d.columns))+' |','| '+' | '.join(['---']*len(d.columns))+' |']+['| '+' | '.join(fmt(v) for v in row)+' |' for row in d.itertuples(index=False,name=None)])
def write(name,text):(ROOT/name).write_text(text.strip()+'\n',encoding='utf-8')

def run():
 qa=json.loads((RESULTS/'SCIENTIFIC_QA_SUMMARY.json').read_text(encoding='utf-8'));assert qa['failed']==0
 form=pd.read_csv(RESULTS/'FORMAL32_CLINICAL_RESULTS.csv');primary=form.loc[form.family.eq('PRIMARY6')]
 eff=pd.read_csv(RESULTS/'PLATFORM_CLINICAL_EFFECTS_ALL43.csv');main=eff.loc[eff.branch.eq('MEM5Y')]
 cal=pd.read_csv(RESULTS/'CALIBRATION_QUALIFICATION_DECISIONS.csv');design=pd.read_csv(RESULTS/'ACTUAL_DESIGN_AUDIT.csv')
 follow=pd.read_csv(RESULTS/'FOLLOWUP_TIME_INFORMATION_SUMMARY.csv')
 sens=pd.read_csv(RESULTS/'PRESPECIFIED_SENSITIVITY_RESULTS.csv');power=pd.read_csv(RESULTS/'PRIMARY_WORKING_EFFECT_POWER.csv')
 aux=pd.read_csv(RESULTS/'PTAU_CONDITIONING_ALGEBRAIC_DECOMPOSITION.csv');pet=pd.read_csv(RESULTS/'PET_SUBSET_VS_CL_POINT_DECOMPOSITION.csv')
 unqual=eff.loc[~eff.calibration_qualified,['gene','branch','assay','n_people','standardized_effect','calibration_qualified']].copy();unqual['status']='POINT_ESTIMATE_ONLY_NOT_VALIDATED_INFERENCE'
 csvout(RESULTS/'UNQUALIFIED_INTERACTION_POINTS_ONLY.csv',unqual)
 tablep=primary[['gene','test','n_people','p_computed','q_BY_within_family','q_BY_combined32']].rename(columns={'gene':'候选','test':'问题','n_people':'人数','p_computed':'原始P','q_BY_within_family':'主BY6 q','q_BY_combined32':'合并BY32 q'})
 tableeff=main[['gene','assay','standardized_effect','standardized_CI95_low','standardized_CI95_high']].rename(columns={'gene':'候选','assay':'检测','standardized_effect':'年变化关联','standardized_CI95_low':'点区间下限','standardized_CI95_high':'点区间上限'})
 tablefollow=follow[['gene','n_people','n_observed_visits','median_span_years','n_highinfo','n_last_visit_atleast4point5y','effective_time_weight_people']].rename(columns={'gene':'候选','n_people':'人数','n_observed_visits':'访视记录','median_span_years':'中位跨度/年','n_highinfo':'高信息人数','n_last_visit_atleast4point5y':'末访≥4.5年','effective_time_weight_people':'时间权重有效人数'})
 tablecal=cal.loc[~cal.qualified,['gene','branch','max_nominal_false_rate','max_nominal_Wilson95_upper','max_tail_Wilson_over_alpha']].rename(columns={'gene':'候选','branch':'分支','max_nominal_false_rate':'最坏名义误拒绝率','max_nominal_Wilson95_upper':'Wilson上界','max_tail_Wilson_over_alpha':'小尾部上界/α'})
 tablesens=sens.loc[sens.gene.eq('NPTX1')&sens.test.eq('SAME'),['branch','n_people','p_computed','q_BY6_diagnostic']].rename(columns={'branch':'敏感性分支','n_people':'人数','p_computed':'同向P','q_BY6_diagnostic':'诊断BY6 q'})
 tablepower=power.loc[power.pattern.eq('POWER_BALANCED_0.05_1')&power.alpha_name.isin(['NOMINAL','BY6_SUFFICIENT']),['gene','alpha_name','minimum_conditional_power']].rename(columns={'gene':'候选','alpha_name':'判据','minimum_conditional_power':'最坏情景条件通过率'})
 report=f"""
# 第10阶段：AD病理背景下的跨检测临床关联与五年认知进展分析

分析日期：2026-10-29。本阶段完成临床值装配、43个设计、有限校准、32项正式检验、9类预设敏感性、描述性条件模型解释和独立复核。前9阶段及原始数据只读。候选选择来自已知病理结果，属于发现后的临床延伸，不是独立队列复现。

## 1. 核心结论与课题定位

主6项、次26项、全部32项BY校正均未产生新发现。NPTX1三检测均呈更正向的记忆年变化，名义同向P=0.0204096、主BY6 q=0.300021；基线记忆同向P=0.002762、次BY26 q=0.276838。方向可以作为探索性证据，不能称预后验证成功。

三个候选的跨检测记忆关联系数差异检验均未通过BY。置信区间仍允许一定大小的差异，本阶段没有临床等效容差及等效检验，不能据此判断临床可互换。

新增医学信息主要是结论边界：既有共同pTau关联不能自动推广为共同临床进展信息；NPTX1的临床关联依赖条件集合和PET子集。原14个共同pTau信号、原三项条件性相对响应偏离及原等效分支均保持封存，不因本轮临床结果重选。

## 2. 人群、时间与AD病理背景

固定ADNI_MEM/ADNI_EF的UWNPSYCHSUM版本，按照第4阶段规则重新读取并去重，独立重现568人基线/最低随访、420人记忆高信息、417人EF高信息覆盖。蛋白及完整协变量交集后的实际主样本如下：

{mdtable(tablefollow)}

基线为采样±90天最近有效日，等距取采样前；随访为采样后>90天且晚于基线，截止5×365.25天，使用真实时间。五年是窗口，主样本中位跨度约3年，末次达到4.5年以上约12%。每人的重复访视不增加独立受试者数。

线性固定截距模型的时间信息有效人数约327–336，前10%受试者贡献约27%时间信息，单人最大贡献约0.52%–0.54%。长期完整随访选择与PET可用性另提供人口学、pTau及认知汇总，不假定缺失完全随机。

独立淀粉样PET按采集日、QC及重处理优先选±90天最近记录：全索引335人，发布A阳性183、阴性152；来源连续值逐条复现Stage5。三候选交集NPTX1/KLK6/PDGFRB分别326/317/324人。跨示踪剂使用Centiloid并调整示踪剂。

这里“独立”表示另一测量方法和发布来源，不是另一队列或误差独立性。同期tau PET仍不可用，没有将多年后的扫描回填为当前tau状态。

![随访信息](figures/06_followup_information_and_selection.png)

## 3. 纵向模型及三检测如何公平比较

每个检测分别拟合同一个结局和共同协变量。主要纵向模型为：

认知_it = 个体固定截距_i + 时间_it × [共同协变量_i的系数 + 该检测蛋白_i的系数] + 残差_it。

共同斜率协变量为年龄、性别、教育、APOE、同日log2 pTau及三检测技术因子、Elecsys批次。技术因素在此调整对认知年变化的系统差异，不只调整基线水平。个体截距吸收稳定认知水平差异；未构造不可识别的个体共识/偏离分数。

采用中心化时间的精确充分统计量计算，与逐访视固定截距OLS等价：每人观测线性趋势按s_tt加权回归。它并非把单人斜率当成无误差结局；各人所有访视对得分的共同影响保留在其估计方程中。个体等权和共同二次时间模型作为预设敏感性。

三平台使用同一人群、共同协变量及结局。蛋白单位固定为原逐蛋白完整临床队列中、去除共同人口学和技术因素后的残差SD，不用新临床结局定标，也不随敏感性改变。主关联是每1蛋白残差SD对应的ADNI_MEM年变化系数。

标准化系数γ_j=β_j×σ_j，联合纳入β和σ估计不确定度。差异检验直接检验(γ_S−γ_N, γ_T−γ_N)=0；不是比较两个P是否分别小于0.05。方差沿三技术轴求和CR1并保留三方程协方差，不称CGM或CR2。独立算法验证了逐访视模型的精确等价。

## 4. 主分析：五年窗口的记忆进展

SAME表示三检测同方向的合取检验；DIFF表示三检测标准化临床关联不全相等的2df检验。主家族固定6项，不能因结果删减候选或问题。

{mdtable(tablep)}

所有三个候选的三检测点估计均为正。正值表示在相同pTau及其他协变量条件下，更高检测值与更正向的记忆年变化相关，不证明认知改善、因果保护或可直接预测个体结局。

{mdtable(tableeff)}

以上区间为点估计95%区间，未作家族多重校正，不能替代BY6的正式结论。NPTX1九检测/问题中的方向最一致，但它的同向联合结果仍未通过校正。

![主结果](figures/01_primary_memory_progression_forest.png)

## 5. 检测间差异与信息量

主差异检验P分别为NPTX1约0.9916、KLK6约0.6030、PDGFRB约0.7866，BY6 q均1。NPTX1三个平台年变化点估计约0.0362–0.0382，数值接近；点估计接近不是等效证据。

NPTX1三对比的点区间大致允许每年−0.038至0.036个复合分数单位差异，KLK6约−0.061至0.043，PDGFRB约−0.029至0.045。原始对比、点区间及同设计三对比Bonferroni区间均交付；这些区间没有跨所有分支保证。

本轮比较的是规定尺度下的观测临床关联，差异可受测量误差衰减、信息量及真实生物信号共同影响。它不是原H约束的重新检验，也不唯一识别分子差异。

![配对差异](figures/02_primary_paired_assay_contrasts.png)

## 6. 基线记忆、EF与淀粉样关联

三类各有3蛋白×2问题，另有淀粉样交互6项、PDGFRB APOE交互2项，合成固定次26家族。全部条件于pTau及共同协变量，PET另调示踪剂。

NPTX1基线记忆的三检测标准化系数约0.210、0.207、0.230，同向P约0.002762，但BY26 q约0.2768，合并32 q约0.3588。EF五年变化、Centiloid关联及PDGFRB APOE修饰也均无BY26发现。

没有证据据此宣称AD特异性、血脑屏障机制、APOE因果通路或临床预测增益。特别是pTau条件下的CL关联与不调pTau的总体关联是不同估计目标；本报告未把条件模型当作全部疾病关系。

![次级结果](figures/03_secondary_clinical_pathology_effects.png)

## 7. 淀粉样效应修饰：有数据，但推断资格不足

连续CL/50交互模型包含蛋白主效应、CL及CL与全部共同斜率协变量/技术因子的交互。三候选的个体信息df分别87、78、85，几何门槛均通过。

有限校准却显示：

{mdtable(tablecal)}

NPTX1与PDGFRB主要未通过严格小尾部要求；KLK6还未通过名义误拒绝要求。三个设计全部按预定规则保留点估计，6项正式家族P置1，次26家族不缩小。计算P/区间仅在完整审计表保留并明确未校准，正文及UNQUALIFIED_INTERACTION_POINTS_ONLY.csv不把它们用于推断。

因此本轮既不能声称存在，也不能声称不存在淀粉样修饰。需要更充分的配对/技术信息或在新分支事先论证、校准的设计，不在已看结果后删技术交互使检验显著。

![正式校准](figures/07_formal_design_calibration.png)

## 8. 预设敏感性

9类敏感性全部完成并通过其预设有限校准，未产生BY6诊断发现。它们是固定诊断，不是9次重新选择主模型。

NPTX1的同向结果为：

{mdtable(tablesens)}

更换pTau函数、限制三年、使用高信息量/个体等权/Batch9、加入基线记忆或共同二次时间后，NPTX1一般仍是正向点估计，但不都达到名义显著。PET子集加CL后Soma点估计转为轻微负值，不能写成所有人群/模型均同向稳健。

基线认知条件模型与主固定截距模型回答不同问题，也可能受到共同基线测量误差和回归均值影响。不同窗口、样本和信息权重改变估计目标，q变化不等于生物效应被解释的比例。

![敏感性](figures/04_memory_prespecified_sensitivity.png)

## 9. 结果知情的描述性条件模型解释

此项v1.2扩展在主结果已知后明确登记，只解释点估计，没有新P值、发现或重选主模型。

同人群、同尺度下，不加pTau时NPTX1三检测系数约0.00665/0.01180/0.00814，加入后为0.03824/0.03751/0.03620。9个候选/检测组合均核验精确遗漏变量恒等式：

β不加P − β加P = [P在完整结局模型中的系数] × [P对蛋白在其他协变量条件下的系数]。

这解释了条件投影如何变化；不是中介分析、保护机制或“pTau解释多少蛋白信号”。pTau同时与蛋白和认知变化有关，可以在统计上遮盖或改变另一项条件关联。

为避免将PET子集变化全部归因于CL，另分开比较全主样本→同协变量PET子集→同PET子集加CL。NPTX1 Soma系数在取子集时由约0.0375降至0.00848，再加CL至−0.00444；NULISA/TMT也在加CL后减小。样本构成与条件调整均有贡献，但这些只是点估计层面的描述，不能赋予因果含义。

![条件模型描述](figures/05_descriptive_conditioning_decomposition.png)

## 10. 校准与条件功效

共43个设计×5情景、215个独立噪声单元格；正式16设计各20000次，敏感性27设计各6000次，共2,410,000独立噪声重复。每次复用37个均值模式，89,170,000模式评价不是独立重复数。

实际蛋白、缺失、技术及随访安排固定。情景包含独立访视投影噪声、AR1/个体随机斜率、pTau异方差、技术轴年龄随机斜率及t5。真实三检测边际系数通过线性约束构造，覆盖部分零和混合方向，不把全零模拟误当作合取零假设的全部检验。

正式13/16、敏感性27/27，共40/43设计合格。主记忆设计在有限情景中偏保守：最坏名义误拒绝率约0.185%–0.320%，Wilson95上界最高约0.408%。三项CL交互的失败如上单列。

工作效应0.05 MEM单位/年/蛋白SD下：

{mdtable(tablepower)}

BY6充分阈值通过率约0.45%–10.81%，不应称为实际BY排序功效；它是无论其他P排序如何都足以通过的门槛。0.10工作效应下相应约87.3%–99.0%，而较小0.02效应信息不足。这些效应不是确立的临床最小重要差异，也未根据观测效应计算事后功效。

校准只覆盖固定设计和有限DGP；逐单元格Wilson界不是整网格同时界。没有重新生成随机蛋白或解决非随机失访，标准化不确定度通过联合影响函数纳入，不构成任意随机蛋白机制的普遍验证。

![条件功效](figures/08_conditional_working_effect_power.png)

## 11. 医学解释与文献位置

2026年独立CANDI/DDI队列研究已报告NPTX1/NPTXR与认知、脑萎缩和进展的联系，因此本课题的新增重点是这些医学关联能否跨检测迁移，而非再次将NPTX1命名为预后标志物。[Nature Communications研究](https://doi.org/10.1038/s41467-026-70472-6)。

本轮NPTX1方向与突触功能相关解释相容，但没有校正后的确认性临床发现，也不是对该外部研究的正式复现。与pTau正相关和在相同pTau条件下更正向的认知变化可以同时存在；前者不能翻译为促病理，后者不能翻译为因果保护。

PDGFRB与APOE/神经血管线索有既有研究依据，但本轮未获得其临床效应修饰支持，且实际测量对象尚未确认就是可溶性受体。[Montagne等，2020](https://doi.org/10.1038/s41586-020-2247-3)。KLK6的总量、加工和活性形式仍待特异实验。

淀粉样PET、流体pTau和聚集性tau并非同一生物学量。没有同期tauPET时，不能将本轮结果解释为完整AD生物分期。[2024修订框架](https://doi.org/10.1002/alz.13859)。

## 12. 可用于论文的结论边界

可报告：在固定的观测尺度和条件模型中，NPTX1三检测对基线记忆和记忆年变化呈一致正向点估计/名义方向证据，但固定家族BY未确认；未获得跨检测认知效应差异或等效证据。淀粉样修饰因校准不足仍未解决。

不能据此报告：临床可互换、真实丰度、分子形式、因果保护、AD特异性、个体预测性能或治疗获益。本阶段没有开发临床风险评分、MCI转痴呆预测模型或AUC增益，没有新的外部验证。

论文主线仍以共同病理关联与条件性相对响应偏离为核心；本轮提供临床外推检验及其信息限制。若后续追求医学验证，优先预定NPTX1的独立临床验证及分子识别证据，并改善技术/随访信息，而不是根据本轮最小P重新定义候选、家族或结局。

## 13. 交付与复核

科学检查{qa['checks']}/{qa['checks']}通过；全部56,610条保存P值的模拟计数/区间重数，7,740个P值固定种子回放，43个真实设计独立QR复算。零组成系数校准使用保存汇总计数/Wilson及初始块回放，未额外保存全部边际P数组，此范围如实登记。

前9阶段{qa['prior_sealed_entries']}个封存条目、{qa['prior_unique_paths']}条唯一文件路径未改变。个体身份、认知分数、随访充分统计量及模型模板只在00_restricted，目录名本身不设置访问权限。

主要入口：关键决策摘要.md、AD病理背景与五年认知进展_结果与审计.xlsx、字段说明与结果导航.md、可复现与运行说明.md。模拟全表和逐次P保留CSV/Parquet，8组图提供PNG/PDF。
"""
 write('STAGE10_MAIN_REPORT_ZH.md',report)
 write('关键决策摘要.md',"""
# 第10阶段关键决策

1. 完成43设计、主6项/次26项/合并32项和9类预设敏感性。主记忆样本NPTX1/PDGFRB各552、KLK6 540；PET交集317–326。
2. 没有新增BY发现。NPTX1记忆进展三检测正向，名义同向P=0.02041、主BY6 q=0.30002；基线记忆P=0.002762、次BY26 q=0.27684。方向可作探索，不称临床验证成功。
3. 三候选均无BY支持的跨检测认知效应差异；未检验临床等效，不能称可互换。
4. 241万独立噪声重复、215单元格校准：40/43设计合格。3个淀粉样交互全部未通过，6项正式P置1、家族不缩小；不能判断有无淀粉样修饰。
5. 五年为观察窗口，中位实际跨度约3年，≥4.5年约64–66人。较小共同效应的信息有限，未显著不等于无临床关联。
6. NPTX1加pTau后正向系数增大；其代数关系已核验。PET子集选择和CL调整也分别改变估计。这是条件投影描述，不是保护/中介机制。
7. 2090项科学检查通过；全部模拟检验P计数重数及7740个P值回放通过，前9阶段5597封存条目未改变。交付含总报告、Excel、8组PNG/PDF、代码和可复核数据。

课题继续以共同病理关联与条件性测量差异为主线。本轮临床外推未获校正后确认，不应把医学贡献写成已证实预后用途。
""")
 write('模型原理与精确约化说明.md',"""
# 模型、单位与检验

## 固定个体截距的精确约化

对每人观测t_it和认知Y_it，在该人的观测日内中心化，令w_i=Σ(t_it−平均t_i)²，s_i=Σ(t_it−平均t_i)Y_it，l_i=s_i/w_i。

模型Y_it=α_i+t_it(Z_iθ+X_ijβ_j)+ε_it，在剔除α_i后与对l_i~Z_i+X_ij做w_i加权最小二乘完全等价。程序用伪观测y_i=s_i/√w_i、伪设计√w_i[Z_i,X_ij]拟合。每人只有一个聚合估计方程；l_i的噪声不会被当成零。所有协变量/技术因子进入时间斜率项。

共同二次时间项的敏感性另将t²在每人内中心化，再与t正交化。由每人至多两个正交投影重建完整正规方程；第二分量不存在时置零，不虚构新受试者。几何门槛使用个体人数减目标设计秩的保守信息df。

这是给定观察安排的线性投影，不是每人真实恒定下降率、不作个体评分或治疗效果解释。

## 共同临床信息与跨检测差异

每个检测单独进入共同协变量模型，三方程同结局、同人群并保留联合影响函数。三者不是同时相互调整的独立蛋白效应。

曝光尺度σ_j来自原完整临床队列的共同人口学/技术残差，不含pTau和新结局，所有分支保持一致。γ_j=β_jσ_j，影响函数为σ_j IFβ_j + β_j IFv_j/(2σ_j)。v_j=σ_j²。程序按三技术轴的单向CR1矩阵求和，不是CGM/CR2。

具体使用每轴G/(G−1)乘簇得分外积和，未额外乘(N−1)/(N−K)；这种工作版本依赖本轮实际设计校准。

同向关联是三β均正或均负。p+=max Φ(−z_j)，p−=max Φ(z_j)，p_same=min(1,2min(p+,p−))。必须处理两个方向；不能事后选单侧。

检测间差异检验Dγ=0，其中D=[−1,1,0; −1,0,1]；Wald的2df检验纳入完整联合方差。未拒绝不是等效。系数差可能反映误差衰减或真实信息不同，不能唯一定位分子机制。

交互全零使用3df联合检验，交互标准化强度差用2df；CL单位固定为50 Centiloid、APOE为携带/不携带差异。共同协变量/技术与修饰变量的交互也进入设计。

## 校准和来源

模拟固定蛋白、缺失、时间和技术安排，使用线性约束构造三方程的真实边际系数。37模式复用同一噪声，不能当成独立噪声数。标准化影响函数已纳入误差传播，但模拟未重新生成任意随机蛋白；失访机制也未被识别。资格只针对预定有限情景。

前9阶段pTau/H结果来自另一科学问题和已知候选选择。本阶段是结果知情后的临床延伸，各家族校正不能当作全课题总体5%FDR，也不是外部重复。

详细参数见STAGE10_PROTOCOL_v1.0.md和CALIBRATION_EXECUTION_DETAIL_v1.1.md。结果知情的pTau/PET描述性拆分另见v1.2扩展，不生成新P。
""")
 write('论文方法与结果建议_EN.md',"""
# Suggested manuscript wording

## Methods

We examined clinical transportability of three candidates selected in the preceding pathology-related measurement analysis. The extension was informed by prior results and was not an independent replication. ADNI_MEM was the primary longitudinal outcome; baseline memory, longitudinal executive function, amyloid PET associations and prespecified effect-modification questions were secondary.

Baseline cognition was selected within ±90 days of the protein assessment. Follow-up assessments occurred more than 90 days after sampling, later than the selected baseline, and within five years. The primary model included subject-specific fixed intercepts and time interactions with the assay measurement, age, sex, education, APOE, same-day log2 pTau and all three assay technical factors and the Elecsys batch. Within-person sufficient statistics yielded an exact reduction of the visit-level model. Repeated visits were not treated as independent participants.

The three assays were compared on identical cases using a fixed protein residual-SD unit derived from the original complete clinical cohort, excluding the new outcomes and pTau from the scaling model. Coefficient and scaling uncertainty were propagated jointly. Covariance was estimated as a sum of three one-way CR1 matrices, rather than CGM or CR2. Same-direction association was evaluated by a two-direction intersection–union test; differences in standardized association coefficients were tested directly using a two-degree-of-freedom Wald statistic.

Six primary and 26 secondary hypotheses were fixed before the new target associations were fitted. BY correction was applied within each family, with a combined 32-question sensitivity analysis. Nine prespecified sensitivity analyses were retained as diagnostics. Unqualified tests remained in the original family with P=1.

Finite-design simulations comprised 2,410,000 independently seeded noise repetitions across 215 cells. Thirty-seven mean patterns per repetition covered partial conjunction nulls, mixed directions, equal nonzero associations and working-effect power. Observed proteins, visit schedules, missingness and technical structures were fixed; this was not a universal calibration over random-exposure mechanisms or informative attrition.

## Results

The primary memory models included 552 participants for NPTX1 and PDGFRB and 540 for KLK6. Median observed spans were approximately three years, despite the five-year observation window. Only 64–66 participants per candidate had a final assessment at least 4.5 years after sampling.

All three NPTX1 measurements showed positive coefficients for memory change, ranging from 0.0362 to 0.0382 ADNI_MEM units per year per protein residual SD. The same-direction test was nominally significant (P=0.02041), but did not pass the primary BY6 criterion (q=0.3000). No primary assay-difference test passed BY correction. NPTX1 also showed nominal same-direction association with baseline memory (P=0.002762), but not after secondary-family correction (q=0.2768). No secondary BY26 or combined BY32 discovery was established.

Forty of 43 designs qualified under the fixed calibration rules. The three amyloid-modification designs did not qualify, and their six formal hypotheses were assigned P=1 without reducing the family. These results cannot establish the presence or absence of amyloid modification.

Result-informed descriptive calculations showed that NPTX1 coefficients increased after pTau adjustment. An exact omitted-variable identity verified the numerical change, but did not identify mediation or protection. Separate point estimates also distinguished changes associated with selecting the PET subset from those associated with adding Centiloid in that subset.

## Interpretation

The clinical extension did not confirm that the previously observed shared pathology associations translated into BY-supported shared cognitive-progression associations. Positive conditional estimates for NPTX1 remain exploratory. Failure to detect assay differences does not establish clinical equivalence or interchangeability. These analyses did not identify true abundance, molecular forms, causal effects, individual predictive performance or treatment benefit. Limited power for smaller common effects, conservative inference and uneven follow-up constrain interpretation.
""")
 write('字段说明与结果导航.md',"""
# 结果导航

优先阅读PRIMARY6_MEMORY_PROGRESSION_RESULTS.csv、PLATFORM_CLINICAL_EFFECTS_ALL43.csv、PRESPECIFIED_SENSITIVITY_RESULTS.csv及INTEGRATED_PATHOLOGY_CLINICAL_EVIDENCE.csv。完整正式家族见FORMAL32_CLINICAL_RESULTS.csv。

- SAME：所有三检测同方向的两方向合取P；DIFF：三标准化系数相等的2df差异P；ANY：三系数全零的3df P。ANY仅在预设交互问题中进入正式32家族，其他模型的计算ANY只供数值审计。
- p_computed是计算值；p_for_family是资格通过后进入家族的值。3个淀粉样交互未校准，其6项p_for_family=1；完整审计中的计算P不能用于正式推断。
- q_BY_within_family：主6或次26各自的BY；q_BY_combined32：32项合并敏感性，均不是跨全部研究阶段的总体FDR。
- standardized_effect单位按结局而异：MEM/EF年变化、基线MEM或CL/固定蛋白残差SD。修饰效应另为每50CL或APOE携带差异。不同结局单位不能直接比较效应大小。
- standardization包含SD估计不确定度；零关联P使用原始尺度β。曝光单位在所有病例子集/敏感性中固定。
- 95%区间是点估计区间，不是多重校正发现。3pair Bonferroni只覆盖同一个设计的3对比，不跨所有分支，也不是临床等效。
- information_df=个体人数−最大目标设计秩，不把伪观测/重复访视当成新增个体。n_pseudo_rows用于精确计算，不能当作样本量。
- effective_time_weight_people及权重集中度描述s_tt信息；不是新增有效独立技术簇数。
- calibration_qualified仅覆盖规定有限情景；全部淀粉样交互为False。UNQUALIFIED_INTERACTION_POINTS_ONLY.csv提供便于正文使用的点估计版本，不含未校准P/CI。
- sensitivity q_BY6_diagnostic每个预设分支保留固定6问题，只作诊断，不重选主结果。
- power里的BY6/BY32曲线是第一步充分门槛通过率，足以保证不依赖排序的通过；不是实际BY排序功效。工作效应未定为临床阈值。
- SIMULATION_FILE_INVENTORY.csv的seed按文本保存。37模式复用噪声；独立噪声总数241万，模式评价8917万，不能混称。
- INDEPENDENT_SIMULATION_RECOUNTS.csv重数全部保存检验P；组件零系数的边际校准只保留汇总和初始块回放，范围见SCIENTIFIC_QA_SUMMARY.json。
- PTAU_CONDITIONING_ALGEBRAIC_DECOMPOSITION.csv与PET_SUBSET_VS_CL_POINT_DECOMPOSITION.csv是结果知情的描述，不含新检验，不是中介或因果归因。

所有个体身份、认知值、逐人趋势和模板在00_restricted。临床装配摘要中的new_associations_fitted=0描述装配时阶段，并非最终未执行分析。工作簿INDEX逐表导航，超过20000行的表保留CSV_ONLY避免截断。
""")
 write('参考文献与来源说明.md',"""
# 方法与医学背景来源

1. ADNI UWNPSYCHSUM：本地I09原始文件及SHA见INPUT_SOURCE_MANIFEST.csv，使用第4阶段已审计的2023版ADNI_MEM/EF，不切换新版PHC。负认知值不套用TMT缺失规则。
2. ADNI UC Berkeley淀粉样PET：I11源采集日/QC/示踪剂/CL及Stage4已封存元数据优先规则；335个源CL与Stage5逐条复现。来源是同一ADNI人群，不是外部验证。
3. Nichols T, et al. Valid conjunction inference with the minimum statistic. NeuroImage. 2005. https://doi.org/10.1016/j.neuroimage.2004.12.005 。所有组成关联均成立的逻辑，不为本项目边际P提供普遍保证。
4. Benjamini Y, Yekutieli D. The control of the false discovery rate in multiple testing under dependency. Ann Stat. 2001. https://doi.org/10.1214/aos/1013699998 。固定家族BY校正。
5. MacKinnon JG, Nielsen MØ, Webb MD. Cluster-robust inference: A guide to empirical practice. J Econometrics. 2023. https://doi.org/10.1016/j.jeconom.2022.04.001 。有限簇、技术依赖与高杠杆问题；本研究SUM3是工作估计器而非CGM/CR2。
6. Wilson EB. Probable inference, the law of succession, and statistical inference. JASA. 1927. https://doi.org/10.1080/01621459.1927.10502953 。二项误拒绝率区间；逐单元格区间不提供全网格同时保证。
7. Cerebrospinal fluid NPTX1 and NPTXR predict neurodegeneration and clinical progression in Alzheimer's disease. Nature Communications. 2026. https://doi.org/10.1038/s41467-026-70472-6 。已报告独立CANDI/DDI人群的临床关联。本阶段没有取得或分析其个体数据，文献结果不能计作本研究正式复现。
8. Montagne A, et al. APOE4 leads to blood–brain barrier dysfunction predicting cognitive decline. Nature. 2020;581:71–76. https://doi.org/10.1038/s41586-020-2247-3 。仅为PDGFRB/APOE修饰假设提供背景；本数据测量形式和该机制均未被验证。
9. Jack CR Jr, et al. Revised criteria for diagnosis and staging of Alzheimer's disease. Alzheimer's & Dementia. 2024. https://doi.org/10.1002/alz.13859 。淀粉样、流体pTau及tau聚集指标的区分，本数据不支持完整同期tau分期。

固定截距模型的精确约化和联合标准化算法为本项目显式推导并验证的实现，详见模型原理与精确约化说明.md。不能把以上文献引用成对完全相同自定义有限簇流程的外部验证。具体数值以本阶段实际数据与代码为依据。
""")
 write('可复现与运行说明.md',"""
# 可复现与运行

默认code/run_all.ps1或python -X utf8 code/verify_delivery.py只读核验SHA256，不重算、不更新网络、不覆盖封存文件。

完整复算应使用code/reproduce_in_copy.py创建项目下新同级目录。它复制代码、协议和公开快照，原9阶段只读。顺序为bootstrap→prepare_clinical→prepare_templates→validate_engine→run_calibration→describe_followup→analyze_real→describe_conditioning→summarize_calibration→make_figures→validate_science→build_delivery→validate_delivery→seal_delivery→verify_delivery。

校准耗时较长，每个design/scenario有独立种子、Parquet和摘要JSON，已完成且散列匹配的单元格可复用。源码和真实模板在运行前记录SHA。重新生成模板及临床输出应在新目录；原正式协议和已封存结果不能覆盖。

真实值首次访问、目标关联首次访问、校准完成的时间分别留档。初始来源清单POST_FREEZE_CLINICAL锁定标签已由本次用户明确授权及本阶段新协议承接，原文件内容/标签未修改。

241万独立噪声重复产生8917万复用均值模式评价；二者不混算。逐次保存SAME/DIFF/ANY P用于全量计数复核；边际组成系数P未全量另存，保留聚合及种子回放，完整可通过固定种子和源码复算。

results含正式/敏感性、设计、时间、校准/功效、独立复核及旧封存核验。simulation_replicates为合成数据，无个体身份。00_restricted含真实认知、样本/访视标识及模型；名称不设置ACL。环境版本见ENVIRONMENT.json。

本阶段使用标准正态/χ²工作P和SUM3 CR1，仅经限定条件校准；未换估计器、删技术项或改变家族追求显著。3个交互未通过的检验正式P=1。所有临床结果为探索性延伸，不是患者决策或预测产品。
""")
 figindex=pd.read_csv(RESULTS/'FIGURE_INDEX.csv');csvout(RESULTS/'FIGURE_REVIEW_RECORD.csv',[{'figure':r.figure,'reviewed':True,'note':'Visual review complete; assay contrasts and calibration gates explicit; no unsupported equivalence or discovery','time':now()} for r in figindex.itertuples()])
 # Every aggregated result is indexed; large tables stay CSV_ONLY.
 preferred=['PRIMARY6_MEMORY_PROGRESSION_RESULTS.csv','FORMAL32_CLINICAL_RESULTS.csv','PLATFORM_CLINICAL_EFFECTS_ALL43.csv','PAIRWISE_CLINICAL_CONTRASTS_ALL43.csv','PRESPECIFIED_SENSITIVITY_RESULTS.csv','INTEGRATED_PATHOLOGY_CLINICAL_EVIDENCE.csv','ACTUAL_DESIGN_AUDIT.csv','CALIBRATION_QUALIFICATION_DECISIONS.csv','PRIMARY_WORKING_EFFECT_POWER.csv','FOLLOWUP_TIME_INFORMATION_SUMMARY.csv','PTAU_CONDITIONING_ALGEBRAIC_DECOMPOSITION.csv','PET_SUBSET_VS_CL_POINT_DECOMPOSITION.csv']
 files=[RESULTS/n for n in preferred]+[p for p in sorted(RESULTS.glob('*.csv')) if p.name not in preferred and p.name not in ['WORKBOOK_SHEET_INDEX.csv','DELIVERY_QA_CHECKS.csv']]
 wb=Workbook();idx=wb.active;idx.title='INDEX';idx.append(['sheet','source_CSV','rows','columns','storage','note']);index=[]
 for i,p in enumerate(files,1):
  d=pd.read_csv(p,dtype={c:str for c in pd.read_csv(p,nrows=0).columns if 'seed' in c.lower()})
  too_long=any(d[c].dropna().astype(str).str.len().gt(30000).any() for c in d)
  if len(d)>20000 or too_long:sheet='CSV_ONLY';note='Full data retained as CSV; no truncation';storage='CSV_ONLY'
  else:
   sheet=(f'S{i:02d}_'+p.stem[:25])[:31];note='';storage='EXCEL_AND_CSV';ws=wb.create_sheet(sheet);ws.append(list(d));ws.freeze_panes='C2';ws.auto_filter.ref=f'A1:{get_column_letter(len(d.columns))}{len(d)+1}'
   for row in d.itertuples(index=False,name=None):
    vals=[]
    for v in row:
     if isinstance(v,np.generic):v=v.item()
     if isinstance(v,float) and not math.isfinite(v):v=None
     if isinstance(v,int) and abs(v)>2**53:v=str(v)
     vals.append(v)
    ws.append(vals)
   for cell in ws[1]:cell.font=Font(bold=True,color='FFFFFF');cell.fill=PatternFill('solid',fgColor='23536B');cell.alignment=Alignment(wrap_text=True)
   ws.row_dimensions[1].height=45
   for j,c in enumerate(d.columns,1):
    ws.column_dimensions[get_column_letter(j)].width=14 if c in ['gene','assay','test','accession'] else 24
    for rr in ws.iter_rows(min_row=2,min_col=j,max_col=j):
     if isinstance(rr[0].value,float):rr[0].number_format='0.000000E+00' if rr[0].value!=0 and abs(rr[0].value)<.001 else '0.000000'
     if 'seed' in c.lower():rr[0].number_format='@'
  row={'sheet':sheet,'source_CSV':p.name,'rows':len(d),'columns':len(d.columns),'storage':storage,'note':note};index.append(row);idx.append(list(row.values()))
 for cell in idx[1]:cell.font=Font(bold=True,color='FFFFFF');cell.fill=PatternFill('solid',fgColor='23536B')
 idx.freeze_panes='A2'
 for col,width in [('A',34),('B',65),('C',12),('D',12),('E',24),('F',65)]:idx.column_dimensions[col].width=width
 wb.save(ROOT/'AD病理背景与五年认知进展_结果与审计.xlsx');csvout(RESULTS/'WORKBOOK_SHEET_INDEX.csv',index)
 save_json(RESULTS/'ENVIRONMENT.json',{'python':sys.version,'platform':platform.platform(),'packages':{n:importlib.metadata.version(n) for n in ['numpy','pandas','scipy','matplotlib','openpyxl','pyarrow']}})
 write('README.md',f"""
# 第10阶段入口

已完成固定三候选的跨检测临床关联、五年认知进展、病理背景、校准、敏感性及独立复核。主6、次26和合并32均无新增BY发现；3个淀粉样交互未通过校准，不能判断效应修饰。原病理发现保持。

- [关键决策摘要](关键决策摘要.md)
- [中文总报告](STAGE10_MAIN_REPORT_ZH.md)
- [Excel结果与审计](AD病理背景与五年认知进展_结果与审计.xlsx)
- [模型原理](模型原理与精确约化说明.md)
- [英文方法与结果](论文方法与结果建议_EN.md)
- [字段说明](字段说明与结果导航.md)
- [可复现说明](可复现与运行说明.md)
- [文献与来源](参考文献与来源说明.md)

科学检查{qa['checks']}/{qa['checks']}通过，{qa['simulation_P_recount_rows']}模拟计数行及{qa['seed_replay_P_values']}个P值回放核验，前9阶段{qa['prior_sealed_entries']}封存条目未改变。交付含8组PNG/PDF，工作簿INDEX登记所有聚合CSV；大表完整保留CSV_ONLY。个体资料仅在00_restricted。
""")
 log(f'Delivery built:{len(index)} indexed aggregated tables,{len(wb.sheetnames)} workbook sheets,8 figure groups')
if __name__=='__main__':
 assert_unsealed();run()
