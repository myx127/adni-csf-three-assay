from config import *

def table(df,columns=None):
 if columns:df=df[columns]
 df=df.copy().fillna('—')
 def val(v):
  if isinstance(v,(float,np.floating)):return f'{v:.6g}'
  return str(v).replace('|',' / ').replace('\n',' ')
 return '| '+' | '.join(map(str,df.columns))+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(val(v) for v in r)+' |' for r in df.itertuples(index=False,name=None))

def write(name,content):(ROOT/name).write_text(content.strip()+'\n',encoding='utf-8')

def run():
 qa=json.loads((RESULTS/'INDEPENDENT_QA_SUMMARY.json').read_text(encoding='utf-8'));assert qa['failed']==0
 cov=pd.read_csv(RESULTS/'COVARIATE_FOUR_MODEL_RESULTS.csv');cov=cov.loc[cov.gene.isin(CANDIDATES.values())]
 c=cov[['gene','model','n','Q','p_value','q_BY_exploratory']].rename(columns={'gene':'蛋白','model':'模型','p_value':'原始P','q_BY_exploratory':'BY q'})
 shap=pd.read_csv(RESULTS/'EDUCATION_APOE_SHAPLEY_Q.csv');shap=shap.loc[shap.gene.isin(CANDIDATES.values()),['gene','delta_Q_both','education_Q_contribution','APOE4_Q_contribution']]
 hv=pd.read_csv(RESULTS/'H_V_AND_B_C_DECOMPOSITION.csv');hv=hv.loc[hv.gene.isin(CANDIDATES.values())&hv.comparison.eq('M0_to_APOE'),['gene','delta_Q','delta_Q_H_point_contribution','delta_Q_variance_contribution']]
 direct=pd.read_csv(RESULTS/'THREE_DISCOVERY_FINAL_EVIDENCE_STATUS.csv')
 direct=direct[['gene','direct_2019_Caprion_peptides','2019_same_day_n_min','2019_same_day_n_max','2019_same_day_residual_df_min','2019_same_day_residual_df_max','new_RID_single_assay_n']]
 assoc=pd.read_csv(RESULTS/'LEGACY_NEW_RID_PEPTIDE_PTAU_ASSOCIATION.csv');assoc=assoc.loc[assoc.cohort.eq('NEW_RID_CLINICAL'),['gene','peptide','n','b_P_per_residual_SD','p_positive_direction','q_BY5_single_assay']]
 support=pd.read_csv(RESULTS/'ORTHOGONAL_MEASUREMENT_SUPPORT_RESULTS.csv');support=support.loc[support.cohort.eq('S5_SAME_DATE'),['gene','peptide','n','corr_residual_N_targetMS','corr_residual_S_targetMS','p_value','q_BY5_exploratory']]
 molecular=pd.read_csv(RESULTS/'LEGACY_PEPTIDE_CANONICAL_MAPPING.csv');molecular=molecular[['gene','peptide','canonical_start_1based','canonical_end_1based','EXACT_sequence_matches','I_L_EQUIVALENT_gene_specific_within_reference']]
 eq=pd.read_csv(RESULTS/'EQUIVALENCE_WORKING_MARGIN_SUMMARY.csv');four=pd.read_csv(RESULTS/'FOUR_ASSAY_CROSSBLOCK_RESULTS.csv')
 four=four[['gene','peptide','n','p_tetrad','direction_ratio','direction_ratio_CI95_low','direction_ratio_CI95_high','q_BY5_equiv_factor_1.5']]
 by=pd.read_csv(RESULTS/'BY_Q_OWN_VS_FAMILY_DECOMPOSITION.csv');kl=by.loc[by.gene.eq('KLK6')&by.model.eq('APOE')].iloc[0]
 report=f'''# 第六阶段：候选信号归因、分子映射与正交验证

执行日期：2026-10-27。本报告对应本地已登记的数据与实际完成的分析；第五阶段34蛋白主家族、三项发现及其P值均保持封存。它是已知发现后的解释与验证，不是发现前预注册，也不是外部队列复现。

## 1. 对研究决策的影响

可以继续以“正交检测下与pTau相关的相对测量约束偏离”为主要问题。PDGFRB、NPTX1、KLK6仍是冻结主分析的三项BY发现，q均约0.0291。教育调整后仍通过BY；APOE调整揭示条件模型和不确定度的依赖，不能把它们写成已定位的蛋白形式差异或APOE因果通路。

新增旧版质谱数据直接覆盖NPTX1两肽段、KLK6三肽段。与主样本同日重叠90–94人，完整共同设计残差自由度20–24，未通过预先冻结的40门槛，因此只作描述，不报告正式相对约束P值。另有95名完全不在Stage5索引571人中的ADNI受试者完成肽段–pTau正向关联；NPTX1两肽段方向一致，但五肽段BY校正后最小q=0.0563，不能写验证成功。PDGFRB尚无直接靶向覆盖。

共同方向和共同病理响应需要正面证据。34家族中11个蛋白满足本阶段信息及结构门槛，在±0.2、±0.35、±0.5 SD三个事后工作容差下均无BY等效成立。这个结果既不把未显著蛋白判为共识，也不能把它们全部判为不一致。

## 2. 证据来源与边界

原三检测共同样本为571个唯一RID/日期，主模型使用独立同日Elecsys pTau及年龄、性别和技术因素。原输入、Stage1/4/5封存件只读使用。第六阶段新增Emory两份靶向MS、两份方法、六份旧版CSF来源，以及官方UniProt条目和人类参考蛋白质组；来源字节数、SHA256、URL/版本均留档。

| 来源 | 本轮可评价内容 | 独立性及限制 |
| --- | --- | --- |
| V01/V02：Emory 2022/2023靶向MS | 70候选池覆盖9个基因，Stage5合格集覆盖SOD1/CHI3L1/MDH1/NPTX2，共5肽段 | 无NPTX1/PDGFRB/KLK6；跨发布重叠按固定来源优先，不能算独立验证次数 |
| L01 CSV＋L02 2014整合XLSX＋L03方法 | 2014 MRM同一实验/发布家族，303行、287 RID；含NPTX1两肽段、KLK6三肽段 | L01没有EXAMDATE；尚无核实的逐行运行/日期交叉表。本轮完成覆盖/来源审计，未拟合候选关联；L01/L02不计为两次实验 |
| L04：2019 Caprion纵向MRM压缩档 | 750行，749有效身份/日期，198 RID、729唯一RID/日期，20组同键重复 | 与2014是不同实验，但来自ADNI；同日仅表示报告日期一致，不能证明同一分装 |
| L05旧版多重免疫数据 | 检查PDGF标签 | PDGF-BB是配体，不是P09619受体PDGFRB，不能替代 |
| L06 BACE/sAPPβ档案 | 结构和方法审计 | 不覆盖三项重点，不做无关结果扩张 |

冻结初始方案、各扩展及首次定量访问记录均保留。早期“未覆盖”仅适用于Emory V01/V02，最终总状态以`DIRECT_DISCOVERY_VALIDATION_STATUS.json`为准；没有把后来发现的旧数据追称为最初已知。

## 3. 归因的原理

每个蛋白的三种检测均回归同一个P和共同协变量集合，得到病理系数b_N、b_S、b_T及残差交叉协方差c_NS、c_NT、c_ST。在本研究条件假设下，相对均值响应与残差共同方向相容要求：

H₁ = b_S c_NT − b_N c_ST = 0；
H₂ = b_T c_NS − b_N c_ST = 0。

联合统计量Q = H′V⁻¹H。V包含系数及协方差估计的不确定度，并沿三技术轴求和CR1。H偏离反映可观察约束的偏离；任意相关测量误差仍可产生相似模式，不能据此唯一分解真实丰度与检测特异作用，也不构造个人F/D评分。

四模型使用相同病例与同一个P标准化：M0；M0＋教育；M0＋APOE ε4；M0＋二者。原始pTau选择、缺失、LOD、TMT/GIS规则不改。比较时以1/n残差乘积统一H及V的尺度，防止单纯残差df改变造成虚假的“点估计变化”；Q与原检验相同。

对于H和V的改变，分别固定其中一个、替换另一个，并平均两个替换顺序；教育/APOE添加也平均两个顺序。这是精确代数Shapley分解，不是中介分析。BY q还取决于全家族排序与累积最小值，因此另将自身P变化和其他33个蛋白P变化按同样的两顺序分解。

## 4. 三个候选的归因结果

{table(c)}

仅加教育，三个候选q均约0.0330。仅加APOE，PDGFRB、NPTX1的q约0.1172，KLK6约0.0830；二者都加后分别约0.1065、0.1065、0.1047。这些模型各自按固定34家族作探索性BY，不把它们重选为新的主模型。

教育/APOE对两者都加时ΔQ的添加顺序平均贡献：

{table(shap)}

APOE单独加入后ΔQ在H点估计和V度量之间的分解：

{table(hv)}

PDGFRB和NPTX1的Q下降主要来自V改变；固定V看H点估计，Q的贡献反而上升。因此“APOE解释了这些蛋白的多少生物学信号”没有依据。KLK6的自身P从{kl.p_M0:.7f}变为{kl.p_new:.7f}，略微改善，q却从{kl.q_M0:.5f}升至{kl.q_new:.5f}：自身P的q贡献为{kl.q_change_own_p_contribution:.6f}，其他蛋白P变化贡献为{kl.q_change_other_proteins_p_contribution:.6f}。不能把KLK6的q增大叙述成自身病理响应减弱。

新增教育/APOE单独模型共70个真实设计，六情景各4000次、共168万独立模拟；固定原方差版本，全部通过本轮有限情景校准。它们仍不解决任意误差相关的识别问题。

## 5. 分子映射能说明什么

保留原始NULISA靶标/登录号、Soma特征和TMT聚合列。Soma对应NPTX1 X9256.78（SeqId 9256-78）、PDGFRB X3459.49（3459-49）、KLK6 X3450.4（3450-4）。这些是发布注释的一致性证据；当前没有NULISA抗体表位、Soma结合区域及原TMT肽段层证据。

官方UniProt规范序列及区域：NPTX1 432 aa，信号肽1–22、链23–432、pentraxin域226–428；PDGFRB 1106 aa，胞外33–532、跨膜533–553、胞内554–1106；KLK6 244 aa，信号肽1–16、活化前肽17–21、链22–244。这些是数据库结构注释，不是本数据已证实的实际测量对象。[NPTX1](https://www.uniprot.org/uniprotkb/Q15818/entry)、[PDGFRB](https://www.uniprot.org/uniprotkb/P09619/entry)、[KLK6](https://www.uniprot.org/uniprotkb/Q92876/entry)。

{table(molecular)}

除规范序列定位外，还下载UniProt 2026_03人类参考蛋白质组UP000005640的canonical＋isoform FASTA，共169,651条序列，对16个候选池相关肽段作精确及I/L等同检索。五个重点肽段在该限定参考集合中均未匹配其他基因。KLK6 LSELIQPLPLER、YTNWIQK同时匹配Q92876和Q92876-2；ESSQEQSSVVR仅在规范条目匹配。详细全部匹配清单保留。[参考蛋白质组及异构体定义](https://www.uniprot.org/help/proteome)。

结论只能到“指定参考集合中的基因来源特异性”。它不排除集合之外的序列、污染、转变干扰或非特异结合，也不能由内部胰酶肽段判定完整蛋白、裂解片段、成熟/酶活状态。PDGFRB的可溶胞外形式、NPTX1多聚状态及KLK6活化/剪切差异可以作为讨论中的待验证机制，不能写成已发生的原因。

## 6. 重点蛋白的直接正交证据

{table(direct)}

L04源字典称相对转变为标准化log2值，方法正文又提自然对数。本轮直接使用发布值，不再取log；两条有效转变等权平均为一个肽段，不把同肽段的两个碎片离子当成两个独立蛋白指标。对数底的固定正比例缩放不会改变相关、Wald或残差SD标准化系数。发布层面零/−4不套用TMT的缺失规则；不能从处理后的发布值重构原始L/H、LOD或QC CV。

重复记录按运行、子运行、进样顺序、原行号选择，完全不看关联强弱。完整三检测设计共同调整年龄、性别、Elecsys批次、N/S检测板、MRM板、运行/子运行及仪器。冻结门槛要求n≥80、残差df≥40和足够检测簇；所有同日重点肽段均因df不足而未进入正式相对约束检验。描述结果见`LEGACY_SAME_DATE_RELATIVE_DESCRIPTIVE.csv`，无正式P值。

对新增RID，先固定每人最早报告MRM日期，再判断同日pTau及协变量；排除全部571个Stage5 RID，得到95人。肽段回归pTau、年龄、性别及MRM/Elecsys技术因素；方差取板和运行CR1求和，正向单侧Student t的df=min(12,5)−1=4。方向来自原发现后已知的正向系数，属于事后单方向支持分析。五肽段共同BY，另作预先规定的去PINK敏感性；新增95人没有PINK，两个结果完全相同，不能算第二次复现。

{table(assoc)}

五肽段系数均正向。NPTX1两肽段原始单侧P<0.05；KLK6仅YTNWIQK原始P<0.05；全部BY q>0.05。该分析验证的问题是“肽段与pTau有无正向关联”，与“原三检测相对响应是否偏离共同方向”不同。不得以其正关联宣布原不一致模式获得独立验证。

十个新RID肽段设计各六情景4000次，共24万独立模拟，全部通过有限情景校准；方法在所检验零假设下明显保守，最大误拒绝率不超过1.25%，最大Wilson95上界约1.64%。保守性限制功效，不据此下调BY或改变推断方法。

## 7. Emory靶向MS的测量支持

V01为729临床样本行、719唯一RID/日期，10组重复；V02为833行、832唯一RID/日期，1组重复。两来源重叠不作为独立样本相加。REPLICATE1/2是实验室ID，不是两次强度；V02从第二个ID恢复真实板号，实际10板，与方法写9板的差异留档。跨发布按预定来源优先、同源按元数据顺序，未择优肽段/强度。

{table(support)}

共同三检测使用N/S/新增靶向MS以及对应共同设计，条件集合与原N/S/TMT不同，不能把P值差直接当成复现/失败。全部五肽段在固定BY家族下无相对约束偏离拒绝；这不证明测量一致。新RID仅84–87人，完整设计残差df<40，未作正式相对约束推断。

方法文献和本地方法的重标记肽段为肽段层内标，加入消化之后不能校准全部蛋白回收/消化损失，不是真实丰度金标准。原始导出未保留可恢复的AT+/AT−等QC池标签，因此不混池计算原始CV。作者方法的肽段/池CV另表提取并明确“来自方法，未由原始数据重算”。[Emory靶向CSF资源说明](https://doi.org/10.1038/s41597-023-02158-3)。

## 8. 四检测共同方向：可用的探索，不能替代主验证

此扩展在已看到Emory三检测结果后提出，对四个覆盖蛋白的全部五肽段执行。使用NULISA、Soma、原TMT和新增靶向MS，同一病例和四技术轴设计。允许N/S组内误差相关及两种MS组内误差相关，但跨亲和–MS误差相关仍设为零。在单一共同方向工作模型下，2×2跨组协方差块应秩1：

c_NT c_SM − c_NM c_ST = 0；R = c_NT c_SM / (c_NM c_ST)。

按四技术轴求和CR1、1df Wald；六类有限零假设情景各4000次、五设计共12万次，全数通过。该检验是本项目实现的条件性秩约束检验，四阶矩/tetrad思想不能被引用成已经发表验证过的相同有限簇算法。

{table(four)}

五肽段均未拒绝秩1，但不因此判等效。只有SOD1和NPTX2满足全部四个跨组协方差95%CI正向的信息门槛。NPTX2方向比0.870，95%CI 0.691–1.097；1.5倍工作容差的名义等效P=0.0120、BY q=0.1368。只能称相对精确的探索性共同方向支持，不能宣布BY等效、真实共同丰度已识别，更不能把NPTX2当作NPTX1验证。两档方向容差边界另作4万次校准，仍仅覆盖有限工作情景。

## 9. 共识的正面检验

定义d_S = H₁/(c_NT sqrt(v_S))，d_T = H₂/(c_NS sqrt(v_T))，单位为残差SD/每1SD pTau。联合影响函数包括b、c及v的不确定度；比值区间采用Fieller，避免弱分母时用窄的常规比值区间冒充精度。两条对比、双向边界共四个单侧检验，联合等效P取最大值；这是交并检验，必须全部边界成立。

{table(eq)}

门槛包含简单协方差结构点相容、三个交叉协方差95%CI全正、两个Fieller90区间有界以及本轮边界模拟通过。35合格设计（MAPT控制单列）共84万次边界模拟；正式探索家族仍为原34蛋白，不合格者P置1。11个合格者在所有工作容差下均无BY等效。这些事后工作容差不具有已确立的生物学阈值含义。[Fieller原理](https://doi.org/10.1111/j.2517-6161.1954.tb00159.x)、[BY多重检验](https://doi.org/10.1214/aos/1013699998)。

## 10. 三个重点如何进入论文

| 候选 | 可作为结果的内容 | 可讨论但尚未验证的解释 |
| --- | --- | --- |
| NPTX1 | 冻结主发现；简单单因子协方差点结构相容；直接两肽段序列映射；95名新增受试者两个肽段正向支持、BY未过 | 相对病理响应的检测差别有生物学或误差解释；不能声称具体多聚体/片段已识别 |
| PDGFRB | 冻结主发现；协方差存在Heywood点诊断；教育稳健、APOE条件敏感；直接靶向覆盖缺口明确 | 受体剪切/胞外形式是候选机制；无表位及原TMT肽段，不能定位 |
| KLK6 | 冻结主发现；原Soma–TMT残差相关约−0.246，简单正向单因子点结构不相容；直接三肽段；新增单肽名义支持，BY未过 | 基质/结合、不同形式、活化或技术误差均可能；不能由内肽判断酶活或裂解 |

主文结果优先呈现冻结主发现、四模型归因、三蛋白的注释/肽段证据矩阵；直接肽段结果作为有限正交支持附表/图。Emory其他四蛋白、四检测秩约束及操作性等效作为探索性补充，避免将不同问题的P值组合成新的发现。

论文核心表述应保持为“相对于指定共同方向模型的病理相关测量不一致”。若标题/摘要声称已分离或验证真实共同丰度、检测特异病理信号及其分子形式，当前证据不够。误差相关敏感性、APOE条件敏感性及没有外部复现必须影响结论措辞；可以写入局限，但不能只在局限中披露而在摘要保留更强结论。

## 11. 核查、交付与后续边界

独立复核{qa['passed']}/{qa['checks']}通过，包括140个归因模型独立OLS重拟合、原70模型复现、代数分解、原始3295个Emory缓存值、旧版直接肽段原值均值重建、十个新增单检测模型独立OLS/SUM2/t、完整参考蛋白质组独立I/L检索、Fieller/交并关系、四检测矩及BY基准核对。1010个模拟单元逐文件计数，共3,160,000次独立重复；原Stage1/4/5的2226项封存散列一致。小样本嵌套设计的数值秩复核记录见`NUMERICAL_QA_NOTE.md`。

本轮完成现有资料可执行的归因、映射、直接肽段支持和正式验证可评价性判断。未获得外部队列、表位、原TMT肽段或新湿实验，不把这些缺口标记为“已验证”。认知轨迹及后续临床转化属于下一分析问题，本轮未检验。

下一阶段优先整合主文证据、固定结果表述与图表，再根据论文对验证强度的要求寻求独立同日配对检测及PDGFRB直接覆盖。若要证明具体分子机制，需要N/S结合区域、原TMT肽段和针对形式/剪切位点的实验；仅增加同类型相关分析不足以完成识别。

完整文件索引见`README_研究结果入口.md`与Excel的INDEX页。普通结果不含个体表；`00_restricted`为本地复现数据，目录名不代表系统已配置访问控制。文献/方法依据见`REFERENCES_AND_METHOD_SOURCES.md`。
'''
 write('STAGE6_MAIN_REPORT_ZH.md',report)
 write('DECISION_MEMO_ZH.md','''# 本阶段决策备忘

本轮可执行分析已完成，原三项主发现与第五阶段规则保持不变。

1. 三候选仅加教育仍通过BY；加APOE后q均超过0.05。Q变化包含点估计和不确定度改变；KLK6自身P略改善，q升高主要来自家族其他P变化，不能一概称APOE解释掉信号。
2. 扩展旧档案后，找到NPTX1两肽段、KLK6三肽段；PDGFRB仍无直接靶向覆盖。五肽段在169651条人类参考序列中无跨基因匹配，但不证明表位或蛋白活化/裂解形式。
3. 直接同日重叠90–94人、残差df20–24，未通过冻结df40门槛，正式相对模式复现未执行。95名主样本之外的同ADNI受试者中，NPTX1两肽段正向，最小BY q=0.0563；无BY验证成功，不是外部队列。
4. 未显著不能判共识：34家族仅11个满足等效评价条件，在±0.2/0.35/0.5 SD工作容差下均无BY等效。
5. 其他四蛋白的四检测探索中，NPTX2方向比0.870（95%CI 0.691–1.097），但1.5倍工作容差BY q=0.137；不能代替NPTX1验证。

论文可以围绕条件模型下的可观察测量不一致继续推进。真实共同丰度、检测特异机制、外部复现和具体分子形式目前未被证明，结论强度需相应限定。下一步优先形成结果证据链与论文图表；更强机制结论需要额外检测/表位/肽段证据。
''')
 log('Chinese research report and decision memo built from checked results')

if __name__=='__main__':run()
