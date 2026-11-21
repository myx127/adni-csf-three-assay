from config import *
import platform

BOOK=ROOT/'第12阶段_pTau条件化解释_结果总表.xlsx'
def tab(df):return df.to_markdown(index=False,floatfmt='.6g')
def write(name,s): (ROOT/name).write_text(s.strip()+'\n',encoding='utf-8')

def run():
 qa=json.loads((RESULTS/'SCIENTIFIC_QA_SUMMARY.json').read_text(encoding='utf-8'));assert qa['failed']==0
 r=pd.read_csv(RESULTS/'ALL42_DIRECTIONAL_EXPLANATION_RESULTS.csv');ovb=pd.read_csv(RESULTS/'OVB_ASSAY_DECOMPOSITION.csv');mom=pd.read_csv(RESULTS/'FWL_WEIGHTED_MOMENT_EXPLANATION.csv');inf=pd.read_csv(RESULTS/'EXPOSURE_INFORMATION_AND_LEVERAGE.csv');lo=pd.read_csv(RESULTS/'LEAVE_ONE_TECHNICAL_CLUSTER_POINTS.csv');cal=pd.read_csv(RESULTS/'NULL_CALIBRATION_WORST_SUMMARY.csv');fac=pd.read_csv(RESULTS/'WORKING_SINGLE_FACTOR_POINT_AUDIT.csv');co=pd.read_csv(RESULTS/'COMMON_CASE_COHORT_AUDIT.csv');bgs=json.loads((RESULTS/'BACKGROUND_CONSTRUCTION_SUMMARY.json').read_text(encoding='utf-8'));ss=json.loads((RESULTS/'CALIBRATION_RUN_SUMMARY.json').read_text(encoding='utf-8'));dec=[]
 for g in GENES.values():
  d=r.loc[r.gene.eq(g)].set_index('model');o=ovb.loc[ovb.gene.eq(g)&ovb.new_model.eq('LOG_P')].set_index('assay');ll=lo.loc[lo.gene.eq(g)]
  dec.append({'gene':g,'n':int(co.loc[co.gene.eq(g),'common_stage12_n'].iloc[0]),'gamma_increment_N':o.loc['N','after_minus_before'],'gamma_increment_S':o.loc['S','after_minus_before'],'gamma_increment_T':o.loc['T','after_minus_before'],'LOG_SHIFT_q9':d.loc['LOG_SHIFT','q_BY_family'],'LOG_SHIFT_q42':d.loc['LOG_SHIFT','q_BY_combined42'],'three_P_forms_change_supported_scope':bool(d.loc[PRIMARY_COMPARISONS,'supported_family'].all()),'NO_P_DIFF_p':d.loc['NO_P','p_computed'],'LOG_P_DIFF_q33':d.loc['LOG_P','q_BY_family'],'RAW_P_DIFF_q33':d.loc['RAW_P','q_BY_family'],'RCS_P_DIFF_q33':d.loc['RCS_P','q_BY_family'],'BG2_LOG_DIFF_q33':d.loc['BG2_LOG_P','q_BY_family'],'BG2_RCS_DIFF_q33':d.loc['BG2_RCS_P','q_BY_family'],'BG3_LOG_DIFF_q33':d.loc['BG3_LOG_P','q_BY_family'],'BG2_SHIFT_q33':d.loc['BG2_SHIFT','q_BY_family'],'BG3_SHIFT_q33':d.loc['BG3_SHIFT','q_BY_family'],'LOCO_change_SN_min':ll.CHANGE_S_MINUS_N.min(),'LOCO_change_SN_max':ll.CHANGE_S_MINUS_N.max(),'priority':'Main selected candidate; background sensitivity remains' if g=='NPTX2' else ('Context-dependent supporting candidate' if g=='TAFA5' else 'Background-sensitive supporting candidate'),'post_selection_confirmatory_FDR_claimed':False})
 dec=pd.DataFrame(dec);csvout(RESULTS/'DECISION_SUMMARY.csv',dec)
 refs=[
 {'id':'R01','reference':'MacKinnon DP, Krull JL, Lockwood CM. Equivalence of the mediation, confounding and suppression effect. Prevention Science. 2000;1:173–181.','doi':'10.1023/A:1026595011371','url':'https://pubmed.ncbi.nlm.nih.gov/11523746/','role':'Regression coefficient-change algebra does not identify causal role.'},
 {'id':'R02','reference':'Addressing inter individual variability in CSF levels of brain derived proteins across neurodegenerative diseases. Scientific Reports. 2025.','doi':'10.1038/s41598-024-83281-y','url':'https://www.nature.com/articles/s41598-024-83281-y','role':'Motivates cross-protein background sensitivity, without validating our median proxy or causal interpretation.'},
 {'id':'R03','reference':'Benjamini Y, Yekutieli D. The control of the false discovery rate in multiple testing under dependency. Annals of Statistics. 2001;29:1165–1188.','doi':'10.1214/aos/1013699998','url':'https://doi.org/10.1214/aos/1013699998','role':'BY adjustment; does not cure outcome-driven candidate selection.'},
 {'id':'R04','reference':'MacKinnon JG, Nielsen MO, Webb MD. Cluster-robust inference: A guide to empirical practice. Journal of Econometrics. 2023;232:272–299.','doi':'10.1016/j.jeconom.2022.04.001','url':'https://doi.org/10.1016/j.jeconom.2022.04.001','role':'Context for finite-cluster limitations; our SUM3 implementation is explicitly not CGM or CR2.'},
 {'id':'R05','reference':'A Unification of Mediator, Confounder, and Collider Effects. Prevention Science. 2021;22:1185–1193.','doi':'10.1007/s11121-021-01268-x','url':'https://pubmed.ncbi.nlm.nih.gov/34164779/','role':'Why conditional observational coefficients cannot choose among causal structures.'}]
 csvout(RESULTS/'METHOD_AND_BACKGROUND_REFERENCES.csv',refs)
 main=dec[['gene','n','gamma_increment_N','gamma_increment_S','gamma_increment_T','LOG_SHIFT_q9','LOG_SHIFT_q42']].rename(columns={'gamma_increment_N':'NULISA增量','gamma_increment_S':'Soma增量','gamma_increment_T':'TMT增量','LOG_SHIFT_q9':'变化q9','LOG_SHIFT_q42':'变化q42'})
 bgtable=dec[['gene','LOG_P_DIFF_q33','RAW_P_DIFF_q33','RCS_P_DIFF_q33','BG2_LOG_DIFF_q33','BG2_RCS_DIFF_q33','BG3_LOG_DIFF_q33','BG2_SHIFT_q33','BG3_SHIFT_q33']]
 nominal=cal.loc[cal.alpha_name.eq('NOMINAL')];maxnom=nominal.max_null_rate.max();maxwil=nominal.max_Wilson95_upper.max();maxtail=(cal.loc[cal.alpha_name.ne('NOMINAL'),'max_Wilson95_upper']/cal.loc[cal.alpha_name.ne('NOMINAL'),'alpha']).max()
 boundary='三个候选是第11阶段使用同一记忆结局发现后选入。本轮规则虽在新增拟合前冻结，但没有校准“第11阶段筛选→第12阶段解释”的联合流程。q9/q33/q42只表示本轮固定范围的多重性整理，不能宣称选择后或全课题确认性FDR得到控制，也不是独立复现。'
 write('00_先读_关键结论与决策.md',f'''第12阶段已完成：为什么差异依赖pTau条件化

本轮已把现象推进到可检验的统计解释：同一批人、同一认知随访、固定蛋白标准化尺度下，pTau条件化对SomaScan的临床关联系数提升更大，从而扩大了跨检测差异。原值、log和样条给出一致的基础模式；背景调整显示，候选的稳健程度不同。

{tab(main)}

单位：ADNI-MEM年变化 / 第11阶段原冻结蛋白残差SD。增量为有log pTau减无pTau。变化检验是两维配对检验，并非比较显著与不显著。

{boundary}

建议将NPTX2作为最有信息的条件依赖临床例证：BG2 log/样条、BG3 log条件下的跨检测差异q33分别为{dec.iloc[0].BG2_LOG_DIFF_q33:.4g}、{dec.iloc[0].BG2_RCS_DIFF_q33:.4g}、{dec.iloc[0].BG3_LOG_DIFF_q33:.4g}。不过BG3后的“加入pTau所致变化”未通过q33，因此保留残余差异不等于保留pTau特异触发。

TAFA5作为背景及函数形式敏感的支持例证：BG2 log仍有范围内支持，BG2样条与BG3 log不足。TIMP3作为背景敏感例证：加入背景后跨检测差异不足。不能把这些未通过结果写成差异已被完全解释。

结论边界：支持条件依赖的关联非互换性；尚未识别真实丰度、不同分子形式、技术可靠性、保护作用或pTau因果机制。前11阶段封存文件保持原样；本轮未开展血浆、外队列或新临床结局分析。

详见[完整报告](01_完整分析报告.md)、[原理](02_原理与论文表述.md)、结果总表Excel和figures内8组PNG/PDF图。
''')
 report=f'''第12阶段完整分析报告：pTau条件化依赖的统计解释

完成时间：{now()}。范围：ADNI脑脊液、NPTX2/TAFA5/TIMP3、五年ADNI-MEM年变化；源于第11阶段定向解释。

研究问题与推断等级

本轮检验的是三种检测的临床关联系数差异 Dγ=(γS−γN, γT−γN)，及加入pTau前后该差异的改变 Δ=Dγ有P−Dγ无P。它不是第5阶段病理相对反应偏离H，也不是三检测病理关联IUT、共识因子F或个人D分数。不同层面的证据不能相互替代。

{boundary}

同样本与冻结设计

{tab(co)}

九种模型分别为NO_P、LOG_P、RAW_P、RCS_P、BG2_NO_P、BG2_LOG_P、BG2_RCS_P、BG3_NO_P、BG3_LOG_P。三检测使用相同受试者和同一临床时间信息权重，原冻结蛋白残差SD贯穿全部模型，避免把改换标准化分母当作系数变化。RCS使用原base_primary固定5%、35%、65%、95%四结点；未根据本轮结局挑选结点。模型均保留年龄、性别、教育、APOE、三平台技术簇及Elecsys批次。

三平台背景指标由原571人全覆盖矩阵构建：NULISA {bgs['reference_counts']['NULISA']}、SomaScan {bgs['reference_counts']['SOMA']}、TMT {bgs['reference_counts']['TMT']}个参考蛋白。排除原34候选及MAPT/APP/APOE，要求≥90%有限值；NULISA另要求≥80%高于LOD。蛋白按中位数和IQR稳健标准化，再取个人跨蛋白中位数；具体可追溯到BACKGROUND_FEATURE_AUDIT。它们是工作背景代理，不是实测CSF总蛋白或确定的清除/稀释指标。[跨蛋白变异的背景依据]({refs[1]['url']})。

主结果：差异本身确实变化

{tab(main)}

log、原值、样条三种pTau形式×3蛋白的九个配对变化检验均通过本轮q9，并均通过合并q42敏感性范围。相同样本和尺度下，SomaScan的系数增量最大。基础模型中的跨检测差异在三种pTau形式下均保留范围内支持。因此，单纯使用log线性形式不足以解释本轮模式。此结论不证明所有可能函数形式或遗漏变量均已排除。

精确解释：两个方向的关联叠加

控制原协变量后，蛋白与pTau的加权偏相关均为正，SomaScan最大；pTau与记忆年变化的加权偏相关约−0.25。高pTau对应更差的年变化方向，而pTau耦合强的蛋白在无P回归中混合了这一负向项。加入pTau后移除该统计项，蛋白的正向系数上移，而且上移幅度因检测不同。这符合统计抑制模式，不能仅凭代数将其解释为去除了真实病理混杂或识别了直接保护效应。[统计抑制与因果含义的区别]({refs[0]['url']})。

{tab(mom[['gene','assay','weighted_partial_corr_X_pTau','weighted_partial_corr_pTau_Y_without_X','exposure_information_retained_after_logpTau','VIF_from_logpTau']])}

对每种检测，β无P−β有P=δ×π：δ是同时包含蛋白及pTau时pTau的系数，π是pTau对蛋白和原协变量回归时的蛋白系数。样条为各基函数乘积求和。全部54个逐检测分解及36个配对对比重构的最大绝对误差分别为{ovb.identity_abs_error.max():.3g}及{pd.read_csv(RESULTS/'OVB_PAIRED_CHANGE_RECONSTRUCTION.csv').identity_abs_error.max():.3g}。这证明所报告系数变化由条件化的代数结构准确重构，不是生物通路贡献百分比。

单独log pTau后的VIF为{mom.VIF_from_logpTau.min():.3f}–{mom.VIF_from_logpTau.max():.3f}，不支持把模式简单归结为极端共线性；它仍消耗了蛋白信息。加入背景后信息可进一步减少，不能只追踪q值变化。

背景敏感性：三候选需要分层解释

{tab(bgtable)}

表中DIFF为各模型跨检测系数差异的q33，SHIFT为在相同背景下加入log pTau后的配对变化q33。NPTX2在BG2的log/样条及BG3 log均保留跨检测差异；BG2下pTau变化仍有范围内支持，BG3则不足。说明NPTX2残余非互换性较稳定，但“差异由pTau特异引起”的强解释仍受背景影响。

TAFA5在BG2 log保留差异，BG2样条及BG3 log不足；同背景下pTau变化也不足。TIMP3在背景模型中的差异及变化均不足。这一层级支持把NPTX2优先作为临床解释例证，将TAFA5和TIMP3放入背景敏感支持结果。

背景加入同时改变了估计对象、点估计及可用信息。BG3+log条件下，三平台蛋白残差信息相对原无P模型：NPTX2约34%–42%，TAFA5约17%–23%，TIMP3约33%–52%。这些是固定设计的残差平方和比例，不是功效百分比或技术可靠性。背景可能携带真实生物信号，未通过不等于背景解释了全部差异；本轮没有为“背景移除了多少效应”另加配对检验，因此不把模型间变弱定性为已证明的背景因果作用。

测量误差解释的边界

建立F=0.8P+U、Y=−0.06P+0.10U+η、Xj=F+ej的合成模型，三检测载荷均为1，误差方差取125组预定组合。在同一真实信号、没有检测特异分子形式的情况下，条件化P后也可产生显著不同的标准化关联系数；解析结果与200,000合成人OLS相符。其作用是提供反例，不是ADNI误差方差估计。

真实数据中的单因子独立误差工作模型在基础条件下多出现负误差方差/非正交叉协方差，不能报告为可靠性估计。背景调整后部分模型出现点上相容，也不提供拟合优度、机制确认或准确度排序。因此“检测不同分子形式”“SomaScan更准确”“真实蛋白保护认知”都超出了数据。

技术簇、校准及验证

逐一删除NULISA板、Soma板、TMT plex后，三个蛋白log条件化所致Soma−NULISA差异增量均保持正向；范围见DECISION_SUMMARY和F08。这是固定尺度、重拟合协变量的描述性点估计诊断，没有新增P值，也不能排除未测混杂。

本轮42个固定检验×5噪声情景×60,000重复，共{ss['independent_noise_draws']:,}独立噪声；7个零均值结构+24个工作效应结构复用同一噪声，共{ss['pattern_evaluations_reusing_noise']:,}模式评价。包含独立访视、AR+人随机斜率、pTau异方差、年龄相关技术斜率和t5技术斜率情景。42/42通过预定有限情景门槛。名义α=.05下最坏经验零拒绝率{maxnom:.6g}、Wilson上限{maxwil:.6g}，显示SUM3在这些情景下明显保守。BY最小充分阈值的最坏Wilson上限/α为{maxtail:.6g}。这不是任意均值结构/任意聚类相关下的普遍校准证明。检验没有达到显著时，尤其不能解释成等效。

工作效应功效依赖对比方向；功效表中的0.02/0.05/0.10是模拟差异幅度，不是临床MCID；BY最小充分阈值功效也不是完整排序BY流程或选择后功效。背景代理和冻结尺度是在现有设计下条件化的，未传播背景构建的所有随机性。

科学QA {qa['passed']}/{qa['checks']}通过：42个真实检验独立QR与残差簇得分复算；6个已知基础模型复现第11阶段完整P/系数/SE；{qa['saved_P_recount_rows']:,}行已保存P值拒绝数重计；{qa['seed_replay_P_values']:,}个固定种子P值重放；新非线性背景模型与访视固定效应形式核对；BY独立重算；原输入及前11阶段{qa['prior_sealed_entries']:,}条封存记录复核原样。初次下游读CSV曾将字面标签NULL解析为空值，已修正读取规则并完整重跑QA，未改变冻结模型、模拟或真实结果；修复记录保存在evidence。

论文整合建议

原14个共同病理关联蛋白仍作为病理证据层；原病理相对反应偏离H保持其假设与校正敏感性边界。本轮只补充临床层：病理同向证据不保证认知关联可跨检测互换；条件化改变了观察到的关联对象，且作用幅度与检测–pTau耦合及背景结构有关。

建议正文使用NPTX2的同样本配对变化与背景稳健性作为例证，TAFA5/TIMP3为支持及局限性结果。图1–2呈现条件化变化，图4解释代数，图6呈现背景与信息权衡，图7界定机制；其余可作为补充。不要把阶段12的q33替换成阶段11的q68，两个检验范围服务不同问题。表述为“与pTau条件化相容的统计抑制和检测相关投影差异”，避免写成pTau介导或检测不同蛋白形式的证实。

本阶段已回答统计层面的为什么，并定位了无法由现有观察数据识别的因果层面。后续在CSF论文内优先整合证据层级、候选排序及限制，无需用新增血浆分析掩盖识别边界。
'''
 write('01_完整分析报告.md',report)
 write('02_原理与论文表述.md',f'''如何理解“为什么依赖pTau条件化”

先区分两个研究对象。无pTau模型问：在相同人口学和技术背景下，蛋白读数高的人，其记忆年变化如何？有pTau模型问：在相同pTau水平及上述背景下，蛋白剩余变化与记忆年变化如何？第二个问题改变了条件集合，不能自动当作更接近真实因果的问题。

可将原协变量从蛋白X、pTau P、个体记忆线性趋势Y中投影出去。时间权重来自真实随访，投影后的变量记作x、p、y。定义cXY=xᵀy、cXP=xᵀp、cPY=pᵀy，vX=xᵀx、vP=pᵀp。FWL定理给出：

β无P=cXY/vX

β有P=(cXY−cXP×cPY/vP)/(vX−cXP²/vP)

本轮cXP为正、cPY为负，分子减去负项而变大；分母因投影pTau而减少。SomaScan与pTau耦合最强，因此增量更大。这是描述数据协方差结构的精确关系。即便所有β与γ都是正值，也不能读成提高蛋白即可改善记忆。

另一种完全等价的表达是β无P−β有P=δπ。δ来自Y~X+P+原协变量，π来自P~X+原协变量。对样条则为Σkδkπk。乘原冻结σ得到γ尺度的相同恒等式。本轮逐检测和配对差异重构均达到浮点精度。该代数关系本身不能区分抑制、混杂、中介或碰撞变量机制；相关原理见[MacKinnon等]({refs[0]['url']})及[因果结构区分]({refs[4]['url']})。

不能用“无P不显著、有P显著”来证明变化。我们直接构造Δ=[(γS−γN)有P−(γS−γN)无P, (γT−γN)有P−(γT−γN)无P]，联合检验Δ=0。两个模型使用同一个人同一结局，需保留跨模型协方差。实现为堆叠6个系数，使用[-D,D]对比；加上原冻结尺度估计的共享影响函数。六系数矩阵可降秩，只要求目标两维对比可识别。

方差使用三个技术轴簇得分外积和SUM3，每轴G/(G−1)，不额外乘(N−1)/(N−K)。这是项目沿用的技术方差约定，不称作CGM或CR2。用有限情景模拟检查它在本设计下是否过度拒绝，不能从模拟中推导任意数据生成机制都有效。点区间也不是同时或选择后置信区间。

背景调整可揭示差异依赖的条件集合，但同时可能删除真实生物学信息并减少残差信息。q变大可能同时来自效应缩小和精度变差；未通过不能证明完全解释。图6并列展示点估计与信息比例。

真实数据单因子模型的负误差方差是工作假设与样本矩结构不相容的警讯，不是负可靠性。合成反例只证明：不同误差方差本身足以形成条件化差异；不能反推它就是本研究唯一原因。

建议英文结果措辞：

“In this outcome-informed explanatory analysis, pTau conditioning altered cross-assay associations with memory trajectories on an unchanged sample and fixed protein scale. Paired changes persisted across raw, logarithmic and spline specifications. Background-proxy adjustment differentiated candidate robustness, with NPTX2 retaining cross-assay differences across the examined background models. These observations support conditional association non-interchangeability, without identifying assay-specific molecular forms or causal effects.”

建议局限性措辞：

“Candidates were selected using the same memory outcome in the preceding stage. The fixed families in this analysis do not account for that selection. Background proxies may contain biology, and the finite-scenario calibration conditions on their realized construction. Single-factor working assumptions and observational conditioning do not establish true abundance, technical reliability, proteoform specificity or pTau mediation.”

{boundary}
''')
 write('03_Methods_EN.md',f'''Stage 12: explaining pTau-conditioned cross-assay clinical associations

This outcome-informed explanatory analysis evaluated NPTX2, TAFA5 and TIMP3, selected from the preceding ADNI CSF memory-trajectory analysis. The candidate set and hypotheses were not independently selected from the outcome. Rules for the new fits and finite-scenario calibration were fixed before accessing the new extension coefficients. Within-stage BY-adjusted values are sensitivity summaries and do not establish selective or project-wide confirmatory FDR control.

For each candidate, all nine models used the same complete cases (553, 557 and 552 participants) and the original frozen protein residual scales. Individual fixed-effect linear trajectory models were represented exactly by time-information-weighted individual slopes. All models retained demographic, APOE and technical covariates from Stage 11. The conditioning specifications were no pTau, log pTau, raw pTau, restricted cubic spline log pTau, and selected combinations with two or three cross-protein background proxies. Four spline knots were frozen at the 5th, 35th, 65th and 95th percentiles in the original eligible base cohort.

Background proxies were the per-person medians of robustly standardized eligible proteins from each full platform matrix, excluding the entire frozen 34-protein family and MAPT, APP and APOE targets. Coverage requirements were at least 90%, with at least 80% above-LOD detection for NULISA. The reference sets contained 49, 6001 and 1906 proteins for NULISA, SomaScan and TMT. Proxies were fixed working covariates, not measured total CSF protein or validated nuisance quantities.

The model-specific two-dimensional contrast was (gamma_Soma−gamma_NULISA, gamma_TMT−gamma_NULISA). A conditioning-change contrast compared this vector before and after pTau adjustment, preserving cross-model covariance by stacking the six coefficient and scale influence functions. Variance followed the project's SUM3 convention: the sum of one-way technical-cluster score cross-products, each corrected by G/(G−1). It was not interpreted as CGM or CR2. Nine primary paired changes comprised three pTau specifications for each candidate. Thirty-three secondary tests comprised 27 model differences and six background-conditioned changes. BY correction was applied within these fixed scopes, with a combined 42-test sensitivity summary.

Before extension inference, each of 42 tests was calibrated using five noise scenarios, 60,000 independent noise draws per scenario, seven contrast-null mean structures and 24 working effect structures, totaling 12.6 million independent noise draws and 390.6 million mean-pattern evaluations reusing noise. All 42 tests met the fixed finite-scenario gates. Null rejection was conservative within these scenarios. Calibration does not cover all possible mean or correlation structures, the preceding candidate selection, or all uncertainty in constructing background proxies.

Exact omitted-variable and FWL decompositions quantified coefficient changes; these were regression identities rather than causal mediation decompositions. Leave-one-technical-cluster refits were descriptive and introduced no new P values. Single-factor covariance calculations were point compatibility checks only; incompatible variance estimates were flagged rather than interpreted as reliability. An analytical and 200,000-person synthetic single-signal example demonstrated that unequal measurement error can produce conditional cross-assay associations without distinct molecular forms. No error variance was fitted to ADNI in this example.

Independent QA refitted all 42 target tests, reproduced the six known Stage 11 base models, recounted 26,040 saved-P summary rows and replayed 2,520 seeded P values. Stages 1–11 remained unchanged. These results address conditional association structure, not true abundance, causal pTau mediation, neuroprotection, clinical prediction utility or external validation.
''')
 write('04_参考文献与证据用途.md','本阶段方法与背景参考\n\n'+ '\n\n'.join(f"{x['id']}. {x['reference']} DOI: {x['doi']}. [原始来源]({x['url']})。用途：{x['role']}" for x in refs)+'\n\n参考文献支持方法和背景；Stage12的数据结果来自本地计算。未将背景研究的发现当作本队列机制验证。检索核对日期：2026-10-29。')
 write('05_表格与字段说明.md','''结果表与字段说明

ALL42_DIRECTIONAL_EXPLANATION_RESULTS：完整42检验注册；kind区分模型DIFF与配对CHANGE；p_computed为数学计算值；p_for_family为通过设计/校准后用于固定范围校正的P，未合格时保留P=1；q_BY_family为9或33范围内BY，q_BY_combined42为合并敏感性范围。它们未处理前轮结局驱动的候选筛选。

PLATFORM_EFFECTS_ALL27：27模型×3检测。beta为原读数单位的年变化关联；gamma=beta×原冻结sigma，单位ADNI-MEM/year/original protein residual SD；所有模型固定sigma。SE_full_joint包含共享尺度影响；point_CI95是逐项区间，不是同时或选择后区间。

PAIRED_CONTRAST_ESTIMATES_ALL42：S_MINUS_N和T_MINUS_N。MODEL_DIFF行是模型内对比，CONDITIONING_CHANGE行是“有P减无P”的对比变化，二者不能混读。

OVB_ASSAY_DECOMPOSITION：before_minus_after=Σδπ；after_minus_before为相反符号。不是中介贡献比例。OVB_PAIRED_CHANGE_RECONSTRUCTION对照配对跨检测变化。

FWL_WEIGHTED_MOMENT_EXPLANATION：先控制原人口学/技术变量后的时间加权偏相关；information retained=1−r_XP²，VIF为倒数，仅针对增加log pTau。EXPOSURE_INFORMATION_AND_LEVERAGE比较所有模型的蛋白残差平方和；比值不是样本量、功效或可靠性。

WORKING_SINGLE_FACTOR_POINT_AUDIT：独立误差单因子工作模型的矩估计；负误差方差或非正交叉协方差时无有效共同方差比例。即使点上相容，也不证明技术可靠性或生物机制。

CALIBRATION_AND_POWER_FULL：五情景、31均值结构、四阈值。字面字符串NULL是均值结构标签，读取时需keep_default_na=False，避免被解析为缺失。SIMULATION_FILE_INVENTORY列出每格种子与文件哈希；replicates保留逐次P值供重放。功效参数0.02/0.05/0.10不是临床MCID。

LEAVE_ONE_TECHNICAL_CLUSTER_POINTS：逐一删除技术簇的系数与变化点估计，没有新增P值。BACKGROUND_FEATURE_AUDIT保留全部参考蛋白纳入排除过程。00_restricted含受限人级数据，不应公开共享。

所有原始CSV优先于Excel显示格式；完整表和是否入工作簿见WORKBOOK_TABLE_INDEX。PNG用于查看，PDF为矢量图。图表的q值应与对应检验范围一同引用。
''')
 # All aggregate CSVs are indexed; large audit/replicate-summary tables remain CSV to keep Excel practical.
 skip={'WORKBOOK_TABLE_INDEX.csv','DELIVERY_QA_CHECKS.csv','SCIENTIFIC_QA_CHECKS.csv','FIGURE_REVIEW_RECORD.csv','INDEPENDENT_SIMULATION_RECOUNTS.csv','CALIBRATION_AND_POWER_FULL.csv','BACKGROUND_FEATURE_AUDIT.csv','PRIOR_ELEVEN_STAGE_SEAL_CHECKS_INITIAL.csv','PRIOR_ELEVEN_STAGE_SEAL_CHECKS_FINAL.csv'}
 files=sorted(p for p in RESULTS.glob('*.csv') if p.name!='WORKBOOK_TABLE_INDEX.csv' and not p.name.startswith('DELIVERY_'));idx=[]
 selected=[]
 for i,p in enumerate(files,1):
  d=pd.read_csv(p,keep_default_na=False);sheet=(f'{i:02d}_'+p.stem)[:31] if p.name not in skip and len(d)<15000 else ''
  idx.append({'file':p.name,'rows':len(d),'columns':len(d.columns),'sheet':sheet,'delivery':'XLSX_AND_CSV' if sheet else 'CSV_ONLY','sha256':sha(p)})
  if sheet:selected.append((p,d,sheet))
 csvout(RESULTS/'WORKBOOK_TABLE_INDEX.csv',idx)
 with pd.ExcelWriter(BOOK,engine='xlsxwriter',engine_kwargs={'options':{'strings_to_formulas':False,'strings_to_urls':False}}) as w:
  pd.DataFrame(idx).to_excel(w,sheet_name='00_INDEX',index=False)
  for p,d,sheet in selected:
   d.to_excel(w,sheet_name=sheet,index=False);ws=w.sheets[sheet];ws.freeze_panes(1,0);ws.autofilter(0,0,len(d),len(d.columns)-1);ws.set_column(0,len(d.columns)-1,18)
 write('README.md',f'''第12阶段交付入口

先读00_先读_关键结论与决策.md，再读01_完整分析报告.md和02_原理与论文表述.md。03为英文方法，04为参考来源，05解释字段。结果总表为{BOOK.name}，完整机器可读结果在results，8组图在figures。

默认入口仅验证，不重跑：在PowerShell运行 `python -B .\\code\\verify_delivery.py`，或 `.\\run_all.ps1`。默认同时核对本阶段、原输入及前11阶段封存记录，验证只在终端输出。`--current-only`只核对本阶段。

本目录封存后不再覆盖；完整重现请执行 `python -B .\\code\\reproduce_in_copy.py --destination <同一项目下新的空目录>`。它复制代码和冻结协议，在新目录依次完成分析、QA、图与文档，不触碰本目录及前11阶段。全量模拟有12.6百万独立噪声、390.6百万模式计算，需足够磁盘和约几十分钟；依赖numpy/pandas/scipy/pyarrow/matplotlib/xlsxwriter/openpyxl/pypdf。新副本的图须重新审阅，不能沿用本次人工审阅记录，因此自动流程不伪造审阅或自动封存。

本轮为既有发现的定向解释。病理关联、病理相对反应偏离、临床关联差异是不同统计对象。不能把本轮q值解释成处理前轮筛选后的确认性FDR。本轮没有因果机制、蛋白形式、技术可靠性或临床预测效能的识别结论。00_restricted为受限人级材料。
''')
 save_json(RESULTS/'DELIVERY_BUILD_RECORD.json',{'time':now(),'workbook':BOOK.name,'tables':len(idx),'sheets_including_index':len(selected)+1,'csv_only':len(idx)-len(selected),'python':sys.version,'platform':platform.platform(),'scientific_qa_passed':qa['passed'],'scientific_qa_checks':qa['checks'],'references':len(refs),'figures':8,'post_selection_confirmatory_FDR_claimed':False})
 log(f'Documents and workbook built: {len(idx)} tables, {len(selected)+1} sheets.')
if __name__=='__main__':assert_unsealed();run()
