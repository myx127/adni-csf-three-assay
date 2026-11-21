from config import *
import math,importlib.metadata,platform
from openpyxl import Workbook
from openpyxl.styles import Font,PatternFill,Alignment
from openpyxl.utils import get_column_letter

def mdtable(df):
    def fmt(v):
        if v is None or (isinstance(v,float) and not np.isfinite(v)):return '—'
        if isinstance(v,(bool,np.bool_)):return '是' if v else '否'
        if isinstance(v,(float,np.floating)):return f'{v:.6g}'
        return str(v).replace('|',' / ').replace('\n',' ')
    lines=['| '+' | '.join(map(str,df.columns))+' |','| '+' | '.join(['---']*len(df.columns))+' |']
    lines+=['| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None)]
    return '\n'.join(lines)

def write(name,body):(ROOT/name).write_text(body.strip()+'\n',encoding='utf-8')

def run():
    qa=json.loads((RESULTS/'SCIENTIFIC_QA_SUMMARY.json').read_text(encoding='utf-8'))
    assert qa['failed']==0,'Science QA must pass before final report is built'
    core=pd.read_csv(RESULTS/'THREE_CORE_CANDIDATE_EVIDENCE_MATRIX.csv')
    model=pd.read_csv(RESULTS/'CORE_CANDIDATE_MODEL_EVIDENCE.csv')
    assoc=pd.read_csv(RESULTS/'INDEPENDENT_2019_NEW_RID_PEPTIDE_ASSOCIATIONS.csv')
    repeats=pd.read_csv(RESULTS/'CORE_PEPTIDE_TECHNICAL_REPEATABILITY.csv')
    qc=pd.read_csv(RESULTS/'MRM2019_PQC24_RAW_LH_VARIATION.csv')
    geom=pd.read_csv(RESULTS/'MRM2014_THREE_FOUR_ASSAY_DESIGN_AUDIT.csv')
    triang=pd.read_csv(RESULTS/'MRM2014_COLLECTION_CRF_DATE_TRIANGULATION.csv')
    overlap=pd.read_csv(RESULTS/'SUBJECT_OVERLAP_AND_INDEPENDENCE_AUDIT.csv')
    offsets=pd.read_csv(RESULTS/'CORE_PRIOR_ERROR_COVARIANCE_SENSITIVITY.csv')
    effects=pd.read_csv(RESULTS/'INDEPENDENT_PLATFORM_EFFECTS.csv')
    fieldtables={}
    table_core=core[['gene','M0_q_BY_shared34','APOE_q_BY_shared34','BOTH_q_BY_shared34','M0_q_BY_H34','APOE_q_BY_H34','BOTH_q_BY_H34']].rename(columns={'gene':'候选','M0_q_BY_shared34':'同向q：M0','APOE_q_BY_shared34':'同向q：APOE','BOTH_q_BY_shared34':'同向q：二者','M0_q_BY_H34':'偏离q：M0','APOE_q_BY_H34':'偏离q：APOE','BOTH_q_BY_H34':'偏离q：二者'})
    table_assoc=assoc[['gene','peptide','n','b_P_per_residual_SD','p_positive_direction','q_BY5']].rename(columns={'gene':'基因','peptide':'肽段','n':'人数','b_P_per_residual_SD':'标准化系数','p_positive_direction':'单侧P','q_BY5':'BY5 q'})
    table_rep=repeats[['release','gene','peptide','n_complete_pairs','pearson_r','pearson_pair_bootstrap95_low','pearson_pair_bootstrap95_high','Lin_CCC_descriptive']].rename(columns={'release':'实验','gene':'基因','peptide':'肽段','n_complete_pairs':'完整重复对','pearson_r':'Pearson r','pearson_pair_bootstrap95_low':'区间下限','pearson_pair_bootstrap95_high':'区间上限','Lin_CCC_descriptive':'描述性CCC'})
    table_geom=geom.drop_duplicates(['gene','assay_design'])[['gene','assay_design','n','rank_X','residual_df','depletion_column_G','max_hat_leverage']].rename(columns={'gene':'基因','assay_design':'检测设计','n':'人数','rank_X':'设计秩','residual_df':'残差df','depletion_column_G':'去高丰度柱数','max_hat_leverage':'最大帽子值'})
    table_pep=pd.read_csv(RESULTS/'CORE_PEPTIDE_SEQUENCE_REGION_SPECIFICITY.csv')[['gene','peptide','canonical_positions_1based','canonical_end_1based','EXACT_reference_sequence_matches','IL_EQUIVALENT_reference_sequence_matches']].rename(columns={'gene':'基因','peptide':'肽段','canonical_positions_1based':'起点','canonical_end_1based':'终点','EXACT_reference_sequence_matches':'精确匹配序列数','IL_EQUIVALENT_reference_sequence_matches':'I/L等同匹配数'})
    eq=pd.read_csv(S7/'results/CONSENSUS_INFORMATION_BOTTLENECKS.csv')
    csvout(RESULTS/'PRIOR_EQUIVALENCE_ALL_MARGINS_FAMILY34.csv',eq.loc[eq.in_family34])
    scientific_summary=f"{qa['checks']}/{qa['checks']}项科学检查通过；前8阶段{qa['prior_sealed_entries']}个封存条目、{qa['prior_unique_paths']}条唯一文件路径未改变"
    report=f"""
# 第9阶段：核心候选的证据整合与验证

执行日期：2026-10-28。核心候选固定为原主分析发现的NPTX1、PDGFRB、KLK6。本阶段结果知情地完成证据整合、独立计算复核、原始质控恢复、技术重复评估和验证可行性核查。前8阶段保持只读。

## 1. 对课题决策的影响

本阶段实质性补上了NPTX1/KLK6肽段来源、测量质控和技术重复证据，并独立复核全部105个三检测设计。三检测共同病理关联的14蛋白结果和原三项主模型偏离发现均复现。没有新增“严格相对模式验证成功”或“独立外部队列复现成功”。

新增的关键设计结论是：2014年MRM的N/S/MRM条件性设计为n144/141、残差df78/75；加入原TMT后，人数虽不变，设计秩升至125/124，残差df降至19/17。四检测验证无法沿用此前三检测几何资格；仅补齐日期仍不能使四检测达到冻结df40门槛。

论文可以继续围绕“共同病理关联与条件模型下相对响应偏离的共存”组织。NPTX1、KLK6最适合做核心案例；PDGFRB保留为条件性偏离候选。真实共同丰度、检测特异机制、AD特异性和因果通路尚不能由这些结果确认。

## 2. 三候选证据矩阵

各列分别沿用固定34蛋白家族BY。新同向关联与原偏离不是同一假设，不按两个家族各5%直接宣称全研究总体5%FDR。第8阶段已计算的合并68问题校正亦被独立复核，完整数值在CORE_CANDIDATE_MODEL_EVIDENCE.csv。

{mdtable(table_core)}

NPTX1与KLK6在三个模型中均有BY34同向关联支持；其相对响应偏离仅在原M0主模型达到BY。PDGFRB的同向关联q=0.054805835，超过0.05；展示时不因四舍五入将其升级为发现。三个核心候选在原RAW_P敏感性中，NPTX1/KLK6的偏离q约0.0311，PDGFRB约0.4183。

APOE相关q变化不等于APOE因果解释。既有第6阶段代数归因显示，KLK6自身P在APOE模型略改善，但全家族排序导致q增大；PDGFRB/NPTX1还涉及方差改变。保留此区分。

原工作容差±0.2、±0.35、±0.5残差SD的等效结果仍为0项BY成立。共同关联不等于响应等效，等效也不等于真实丰度被识别。

![全家族证据](figures/01_family34_integrated_evidence.png)

## 3. 独立重算验证了什么

新实现使用SVD构建设计空间，重新计算条件病理系数、残差协方差及9维联合影响函数，再按三个技术轴求和CR1得到方差。未调用前阶段统计估计函数。

- 全部35蛋白×M0/APOE/BOTH=105设计重新拟合；MAPT仍单列，不加入34家族。
- 独立复核双向同向关联IUT、原H约束、各自BY34和合并68校正。
- 315个蛋白/模型/检测效应的残差SD标准化及完整联合不确定度得到核验。
- 从2019原始发布重建两转变等权肽段值，独立复现五项95人分析及BY5。
- {scientific_summary}。

前阶段区间使用舍入临界值1.959964；本阶段绘图区间使用精确标准正态0.975分位数。最大差异约3.8×10⁻⁹，核验时按各自临界值解析对齐；点估计、方差、P与q没有因此改变。

上述是计算复核，不是新的生物学重复。本阶段未新增病理检验校准模拟；沿用已完成校准的统计规则。新增40,000次成对自助抽样仅用于10个技术重复相关区间，不能计作40,000个独立实验或新增病理验证。

![核心候选系数](figures/02_core_candidates_association_forest.png)

## 4. 从旧档案恢复了新的质控信息

进一步读取2019 Caprion压缩档内原始Supplementary data工作簿，而不只读取聚合CSV。附件提供774个样本注释，包括750个发布记录和24个过程质控PQC记录；RARC标注为有重复的研究样本。

750个发布Sample Name均唯一匹配附件，十个目标转变的发布归一化值与附件值最大差异约5.3×10⁻¹⁵，缺失位置一致。NPTX1 LENLEQYSR的两个转变各有5项缺失，未用插补代替。

原方法及附件图例明确：
- 0：通过峰质控。
- 1/2：低于检测限标志；原附件保留峰值以维持变异。
- 3/4：缺失或歧义，最终表设缺失。
- 预归一化对数为log2((Light+1)/Heavy)，不是log2(Light/Heavy+1)。

95人分析的950个转变测量全部标志0；24个PQC的240个目标转变测量也全部标志0。PQC原始Light/Heavy比值CV：KLK6约4.9%–8.0%，NPTX1约13.7%–15.6%。这说明所用目标测量有可追溯峰质控和过程稳定性支持，不构成所有样本的绝对准确度保证或其他平台精度估计。

PQC是混合池的过程重复，24条记录不代表24名独立供者。原始峰面积已恢复，但尚无全部色谱、峰形、信噪比和完整标准曲线供独立重新判峰。稳定偏差以及亲和平台误差也不能由MRM池CV排除。既有五肽段关联结果未因新质控改变。

![原始质控](figures/05_raw_transition_process_QC.png)

## 5. 核实的技术重复结果

2014年：16组RID/访视重复对，共32条记录，均有Sample_Duplicated=Yes。303行与原工作簿的320特征在既有6有效数字末位1单位容差内逐行核验通过；这是有公开精度限定的身份匹配，不是原始字节或数值完全相同。

2019年：20个重复RID/日期组另经同一SampleCode基础编号、.01/.02重复标签、RARC组核实，形成20对盲重复。10对跨仪器、13对跨运行/子运行、18对跨MS板。NPTX1 LENLEQYSR仅19个完整重复对，其他肽段20对。

在已发布归一化对数肽段尺度计算相关、Lin CCC和成对差；对每个肽段/实验单独做4000次成对自助抽样，并报告去一对相关范围。

{mdtable(table_rep)}

重复证据并非所有肽段同样精确。KLK6 ESSQEQSSVVR在2014年的相关约0.682，区间很宽；2019年KLK6三肽段相关约0.70–0.72。完整差值、CCC、RMS与去一对结果均保留。没有选择相关最高的肽段替换原先五肽段家族。

相关依赖重复集合中的样本间变异。2014和2019对象、实验和处理不同，不能仅按r高低推断哪个实验更准确。重复对只能评估重复性，无法排除同样本稳定偏差，也不能将原N/S/TMT的相关误差唯一分解为真实丰度及检测特异病理作用。

![技术重复](figures/03_peptide_technical_repeatability.png)

## 6. 95名主样本之外受试者：关联结果与批次影响

这95名受试者与Stage5全部571个RID零重叠，但来自同一ADNI母队列，属于受试者不重叠的单检测肽段支持。方向来自此前已知的正向关联，沿用原五肽段单侧检验、MRM板/运行SUM2 CR1及最小簇数减1的Student t自由度4。

{mdtable(table_assoc)}

五肽段系数均为正，BY5均未通过。最小q仍为NPTX1 FQLTFPLR的0.0563267；不能把“接近0.05”写成验证成功，也没有事后缩减为两基因或挑选单肽段重新校正。

对每个肽段逐一删除12个MRM板、5个运行，每次重新构建相同协变量模型，共85个描述性设计。全部可估计且系数保持正向；NPTX1标准化系数范围约0.436–0.601，KLK6约0.280–0.506。删除后的显著性未用于选择结果，也未作为新的发现。

95人中53人也在2014年MRM的287人中。两个MRM发布来自不同实验，但证据并不在受试者层面完全独立，不将其P值按独立实验直接合并。

![批次影响](figures/06_newRID_batch_influence.png)

## 7. 2014年日期链：进一步核查后仍未关闭

本阶段检查25张含样本标识的本地/附件表。精确匹配仅出现在2014自身Sample Info：LIMS_ID和External_Reference各303行；其他表未检出可连到实际采样日期的直接标识。2014 LIMS与2019 Sample Name基础编号无重叠，不能用看似相同的命名格式强连。

另核查CSFMETH、LOCLAB记录及官方数据字典：VISDATE可以来自assessment，也可以回退到Registry日期。它不天然等于该MRM标本实际采样日。287个MRM受试者/访视中，两张表各有274人与MASTER候选日匹配；这提高访视层面的日期一致性描述，仍缺“MRM样本→标本/分装→采样日”的直接证据。旧ADNI1记录中相应采样Yes/No字段为空，不能把没有Yes解释为没有采样。

新私有模板MRM2014_SPECIMEN_DATE_LINK_REQUEST_TEMPLATE_RESTRICTED.csv列出303条记录所需桥接字段，供用户后续通过获授权渠道核实。本阶段没有联系任何人或发送该表。

## 8. 三检测和四检测：分别审计真实交集

重复记录沿用固定处理日/进样顺序选择；不按病理效应选择。设计共同调整年龄、性别、Elecsys批次及相应检测技术因素，2014年另调整处理日和去高丰度柱固定效应。

{mdtable(table_geom)}

N/S/MRM：NPTX1两肽段均n144、秩66、df78；KLK6三肽段均n141、秩66、df75。N/S/MRM处理日簇数11/36/15，几何门槛表面达到n80/df40/G6，但源标本日期、2014实际设计误差率和2个去高丰度柱的潜在依赖仍未解决。最大帽子值接近1，不能把几何门槛通过当作推断校准通过。

N/S/TMT/MRM：NPTX1两肽段均n144、秩125、df19；KLK6三肽段均n141、秩124、df17。完整技术设计下均未通过原df40门槛。需要更多独立配对信息，或另行论证、事先规定并校准的新设计；目前没有执行2014病理系数或相对约束检验。

两个去高丰度柱的固定效应可调整均值差异，但不会自动消除任意柱内相关。仅2簇不足以给常规簇渐近推断提供普遍保证。原2019同日N/S/MRM设计仍仅df20–24，原严格相对模式检验未执行。

![设计几何](figures/04_MRM2014_three_four_assay_geometry.png)

## 9. 分子身份与解释层级

重新取得三项官方UniProt条目并保留版本/下载记录，规范序列与第6阶段一致。独立扫描冻结的人类参考蛋白质组canonical+isoform共169,651条序列，五肽段的精确及I/L等同匹配均未检出其他基因来源。

{mdtable(table_pep)}

NPTX1 FQLTFPLR位于pentraxin结构域，LENLEQYSR位于另一内部区域。KLK6三肽段位于成熟链内部；其中两条可同时匹配参考KLK6的两个异构序列。这些信息支持指定参考集合中的基因来源，不能证实实际测到完整蛋白、特定异构体、剪切片段或酶活状态。

PDGFRB规范序列标注胞外33–532、跨膜533–553、胞内554–1106。当前靶向面板没有直接P09619受体肽段；PDGF-BB配体不能替代PDGFRB受体。NULISA实际抗体识别区、Soma适配体结合区和原TMT逐肽段证据仍未取得。[UniProt NPTX1](https://www.uniprot.org/uniprotkb/Q15818/entry)、[PDGFRB](https://www.uniprot.org/uniprotkb/P09619/entry)、[KLK6](https://www.uniprot.org/uniprotkb/Q92876/entry)。

![分子肽段位置](figures/07_core_molecular_peptide_map.png)

## 10. 已明确的外部数据入口与访问状态

从原2025 TMT方法中定位到原始质谱Synapse编号syn59804727；未登录访问返回HTTP403，尚未取得原肽段/PSM矩阵。原始质谱文件本身也不自动等于可直接使用的逐肽段定量和样本映射。

外部Emory2022项目syn30549757及子文件公开元数据可读取，明确列出Master.3platform.2fluid.zip（syn30549894）、三平台RData（syn30550193）及分析流程PDF。逐文件下载均返回未认证403，未下载个体数据。该独立母队列采用Soma/Olink/TMT，三检测完整人数上限35，低于原冻结n80门槛，获得文件后仍需重新审计目标覆盖、病理指标和技术元数据。[原论文](https://doi.org/10.1186/s13195-022-01113-5)、[项目](https://www.synapse.org/Synapse:syn30549757)。

公开元数据登记、既有文献和UniProt注释不能计作外部个体级复现。完整请求状态、实体信息和最小缺失材料已留档；未进行账户操作或作者联系。

## 11. 识别假设与原敏感性如何进入结论

三项候选在此前125个工作协方差偏移点中均存在名义P跨0.05，数值范围分别见CORE_PRIOR_ERROR_COVARIANCE_SENSITIVITY.csv。该范围既不是数据校准的真实误差相关范围，也不是连续参数的严格稳健界；它说明约束偏离对测量误差假设敏感。

新增MRM重复和PQC提供的是该实验的可重复误差信息，未覆盖亲和平台，也不能识别所有样本稳定偏差。因此：
- 可以报告“在指定协变量、技术调整和误差结构工作假设下的可观察相对响应约束偏离”。
- 可以报告NPTX1/KLK6的三检测同向病理关联与M0偏离共存，并列展示APOE及原pTau尺度敏感性。
- 分子形式差异、基质/结合、复合物状态和技术偏差均保留为候选解释。
- 不把这些结果写成因果作用、AD特异信号、平台可互换或真实蛋白丰度/检测特异病理作用已唯一分离。

## 12. 下一步最小验证工作

1. 提供2014 LIMS/实验室编号到specimen/aliquot及实际采样日期的交叉表。三检测分支随后需独立完成2014设计校准；四检测分支还需先解决df19/17的信息不足。
2. 获取原TMT逐肽段/PSM、蛋白分组、修饰和GIS/plex/channel映射，优先核对三个核心候选。公开原始文件入口已明确，实际数据尚不可读。
3. 获取NULISA抗体及Soma识别区域/形式响应资料；对PDGFRB补直接受体检测。
4. 如具备实验条件，NPTX1可考察不同区域肽段与复合物/构象响应，KLK6增加活化或酶活特异测量，PDGFRB区分胞外/胞内区域。实验需能产生区分候选解释的观测，而不仅再测一次同基因总信号。
5. 独立队列应先核实同人配对、候选覆盖、独立病理和足够技术信息，再冻结验证问题。当前可以开始论文结果整合，不需要等待尚未取得的材料才整理已有发现。

## 13. 完成状态

本阶段分析执行与可复核交付完成，现有数据支持的计算复核、源质控、重复性、方向诊断、分子映射及验证资格审计均已实施。严格相对模式复现、独立外部队列验证和机制识别仍未完成，已逐项登记；本报告没有把它们写成验证成功。

入口：关键决策摘要.md、核心候选证据整合与验证_结果与审计.xlsx、字段说明与结果导航.md、可复现与运行说明.md。个人标识、逐样本原始峰值和日期请求模板仅放00_restricted，该名称本身不代表设置了访问控制。
"""
    write('STAGE9_MAIN_REPORT_ZH.md',report)
    write('关键决策摘要.md',f"""
# 第9阶段关键决策

1. 核心结论复现：14个三检测同向关联信号在M0/APOE/BOTH保持；NPTX1/KLK6的同向关联与M0相对响应偏离并存。PDGFRB同向q=0.0548，仍未通过BY。
2. 新恢复2019年24个PQC及20对有明确RARC标签的盲重复；2014年16对也经标志核实。95人分析的950个目标转变均原始QC0；PQC原始Light/Heavy CV约4.9%–15.6%。
3. 95人五肽段结果独立复现，85次删板/运行描述均正向，但BY5仍无验证成功；最小q=0.0563。95人与2014 MRM重叠53人，不能按独立人群合并。
4. 2014年日期桥接仍未关闭；数据字典显示VISDATE可回退为Registry日期。四检测设计新审计：NPTX1 n144/df19、KLK6 n141/df17，均未达冻结df40，补日期本身不足以启动四检测验证。
5. 原TMT公开入口syn59804727与外部Emory个体文件已定位，但未认证下载返回403。肽段/表位/独立外部复现缺口保留。
6. {scientific_summary}。交付含中文总报告、Excel、7组PNG/PDF、代码、全部表和SHA256。

论文可围绕共同病理关联与条件性相对反应差异继续推进；尚不能宣称真实共同丰度、具体分子形式或外部模式复现已被证明。
""")
    write('论文方法与结果建议_EN.md',"""
# Suggested manuscript wording

## Methods: evidence integration and measurement support

Following the discovery analyses, we integrated evidence for the three candidates selected from the original relative-response analysis (NPTX1, PDGFRB and KLK6). Same-direction pathology association, relative-response deviation and equivalence within working margins were treated as separate questions. The original 34-protein families and their BY adjustments were retained, and the previously specified combined 68-question sensitivity analysis was independently reproduced. This integration was informed by prior results and was not a prospectively registered discovery analysis.

We independently reconstructed the three-assay models using an SVD representation of the nuisance space. Covariance estimation used the sum of three one-way CR1 covariance matrices, rather than a CGM or CR2 estimator. Coefficient and residual-moment uncertainty was propagated jointly for relative-response constraints and residual-SD-standardized effects. The same participants and pathology scaling were used across the prespecified covariate models.

For orthogonal measurement support, we recovered raw light/heavy peak areas, transition quality flags and sample annotations from the archived 2019 Caprion supplementary workbook. Released sample names and normalized transition values were verified against the workbook. Peptides retained the original equal weighting of two finite transitions. The previously specified five-peptide association family in 95 participants disjoint from all 571 main-analysis participants was independently refitted without changing its one-sided direction, covariance estimator or multiplicity adjustment. Leave-one-plate and leave-one-run coefficients were descriptive diagnostics, with no significance-based selection.

Technical replicate pairs were verified separately within each release. The 2014 data contained 16 pairs explicitly flagged as duplicated. The 2019 pairs additionally required matching sample-code bases, replicate labels .01/.02, and RARC annotation. Pearson correlation, descriptive concordance correlation, paired differences and 4,000 paired-bootstrap correlation intervals were calculated for each peptide and release. These analyses assessed repeatability of released measurements, not external biological replication or absolute accuracy.

## Results: integrated evidence

Fourteen proteins retained BY-adjusted evidence of same-direction pTau association in each of the main, APOE-adjusted and education/APOE-adjusted models. NPTX1 and KLK6 also exhibited relative-response deviation in the original main model. The deviation evidence was covariate-sensitive, and did not retain BY significance after APOE adjustment. PDGFRB retained its original main-model deviation finding but did not meet the BY criterion for shared association (q=0.0548). No BY-supported equivalence finding was established at the previously evaluated working margins.

All 950 target-transition measurements used in the 95-participant analysis passed the archived peak-quality criterion. The ten target transitions also passed quality checks in 24 process-control measurements, with raw light-to-heavy ratio CVs ranging from 4.9% to 15.6%. Twenty blinded replicate pairs in the 2019 experiment were verified using sample annotations; one NPTX1 peptide had 19 complete pairs. These data provided additional measurement-quality support, while between-pair heterogeneity and potentially stable measurement bias remained.

All five subject-disjoint peptide–pTau coefficients were positive, and all 85 leave-one-batch coefficient diagnostics retained a positive direction. No peptide passed the original BY5 criterion; the smallest q-value remained 0.0563. These findings were therefore reported as limited orthogonal association support rather than successful replication of the original relative-response pattern.

## Results: verification limits

No specimen-specific date bridge for the 2014 MRM experiment was established. Clinical VISDATE fields could use registry-date fallback, and were insufficient to verify aliquot-level collection dates. Conditional three-assay designs had 144 or 141 participants and 78 or 75 residual degrees of freedom. Adding the original TMT measurements increased the full design rank, leaving only 19 or 17 residual degrees of freedom, below the frozen threshold of 40. No new 2014 candidate pathology association or relative-response test was fitted.

Public metadata identified external three-platform files and the original ADNI TMT raw-data repository, but unauthenticated file access was denied. No independent external individual-level replication or mechanistic identification was completed.

## Interpretation

Shared pathology association and relative-response heterogeneity can coexist. The findings describe observable conditional measurement behavior under specified modeling and error-structure assumptions. They do not uniquely identify true protein abundance, establish assay interchangeability, demonstrate causality or Alzheimer-specificity, or identify the molecular forms responsible for assay differences. Technical repeatability within MRM does not exclude stable assay bias or identify the complete cross-platform error covariance.
""")
    write('缺失材料最小清单与验证执行条件.md',"""
# 缺失材料与执行条件

## 2014 MRM

私有请求模板：00_restricted/MRM2014_SPECIMEN_DATE_LINK_REQUEST_TEMPLATE_RESTRICTED.csv。

需要逐行保留source_row、RID、VISCODE、LIMS_ID、External_Reference，并补充specimen_id、parent_specimen_id、aliquot_id、actual_collection_date/time、日期来源和实验室/生物库核实标志。若希望声称同一分装，需另提供与NULISA/Soma/TMT的父标本和分装链，不能只给RID/访视或进样日期。

三检测分支应在日期证据合格后冻结实际交集、重复规则和依赖轴，并完成实际设计校准。现有2个去高丰度柱不能只用固定效应调整后假定剩余依赖消失。四检测分支还需要解决df19/17信息不足。

## 原TMT

原方法指向syn59804727，匿名访问403。需要可分析的逐肽段/PSM定量、肽序列与修饰、唯一/共享肽、蛋白分组及归属规则、样本/日期、plex/channel和GIS映射。仅取得数千原始质谱文件不等于已经完成肽段重建或分子归因。

## 亲和检测

三个目标的NULISA抗体对识别区域、Soma适配体结合区、交叉反应、变性/复合物/剪切或修饰响应，以及必要的竞争和回收资料。基因或UniProt注释相同不等于识别同一完整分子。

## 直接验证

NPTX1保留两个不同区域肽段，KLK6保留三肽段并考虑活化/活性特异指标，PDGFRB补P09619受体直接肽段和区域信息。具体实验需在知道验证效应之前规定主要对比和判据，避免挑选最有利的肽段/区域。

## 外部Emory2022

公开项目syn30549757已核实，目标文件syn30549894和syn30550193未认证下载403。需要在获授权后取得个体矩阵、病例/采样信息、独立pTau和技术批次。已发表三平台完整人数上限35，低于当前冻结n80，可先审计描述性测量支持，不能自动沿用正式复现标签。

本阶段未申请数据权限、联系作者或发送任何模板。模板含个体标识，仅供用户按既有数据授权范围使用。目录名称00_restricted不构成技术访问控制。
""")
    write('参考文献与来源说明.md',"""
# 方法与原始资料

执行与检索日期：2026-10-28。项目内新数值来自原始发布或冻结结果的本阶段计算；外部资料仅用于方法、注释和数据可取得性，不计作个体级验证。

1. Nichols T, et al. Valid conjunction inference with the minimum statistic. NeuroImage. 2005;25:653–660. https://doi.org/10.1016/j.neuroimage.2004.12.005 。原文：https://www.fil.ion.ucl.ac.uk/spm/doc/papers/NicholsBrett.pdf 。支持所有组成命题均成立的合取逻辑；不为本项目有限簇边际P值提供普遍有效性保证。
2. Benjamini Y, Yekutieli D. The control of the false discovery rate in multiple testing under dependency. Ann Stat. 2001;29:1165–1188. https://doi.org/10.1214/aos/1013699998 。沿用固定家族BY。
3. Fieller EC. Some problems in interval estimation. JRSS B. 1954;16:175–185. https://doi.org/10.1111/j.2517-6161.1954.tb00159.x 。既有工作等效中比值区间的背景。
4. Schuirmann DJ. A comparison of the two one-sided tests procedure and the power approach for assessing equivalence of average bioavailability. J Pharmacokinet Biopharm. 1987. https://doi.org/10.1007/BF01068419 。工作容差等效与未拒绝零假设的区别。
5. MacKinnon JG, Nielsen MØ, Webb MD. Cluster-robust inference: A guide to empirical practice. J Econometrics. 2023;232:272–299. https://doi.org/10.1016/j.jeconom.2022.04.001 。有限簇、高杠杆和依赖结构；SUM3不是CGM/CR2。
6. Wilson EB. Probable inference, the law of succession, and statistical inference. JASA. 1927. https://doi.org/10.1080/01621459.1927.10502953 。既有有限模拟二项率区间。
7. ADNI/FNIH. 2014 Biomarkers Consortium MRM Data Primer. 本地原PDF路径及SHA见FINAL_LOCAL_SOURCE_MANIFEST.csv。项目说明：https://adni.bitbucket.io/reference/csfmrm.html ；对应原方法：https://adni.bitbucket.io/reference/docs/CSFMRM/04%20Apr%202014%20Bimarkers%20Consortium%20MRM%20Data%20Primer.pdf 。
8. ADNI/FNIH. Biomarkers Consortium CSF Proteomics Project, Data Primer. https://adni.bitbucket.io/reference/docs/CSFPROTEOMICS/BC_MRM_Data_Primer_-_Final.pdf 。实际计算使用本地csfproteomics_20190403.zip及其Supplementary data工作簿。峰质控代码与样本配对定义另见evidence/MRM2019_SUPPLEMENT_LEGEND.txt。
9. UniProt Q15818/NPTX1：https://www.uniprot.org/uniprotkb/Q15818/entry ；P09619/PDGFRB：https://www.uniprot.org/uniprotkb/P09619/entry ；Q92876/KLK6：https://www.uniprot.org/uniprotkb/Q92876/entry 。官方REST条目版本、时间和摘要见UNIPROT_REFRESH_SOURCE_MANIFEST.csv。结构注释不是实际结合区域的实验证明。
10. 冻结人类参考蛋白质组UP000005640 canonical+isoform，169,651条。原始gzip及SHA来自第6阶段；本阶段重新逐序列检索，未声称对所有未知序列普遍唯一。参考集合定义：https://www.uniprot.org/help/proteome 。
11. Dammer EB, et al. Multi-platform proteomic analysis of Alzheimer's disease cerebrospinal fluid and plasma reveals network biomarkers associated with proteostasis and the matrisome. Alzheimers Res Ther. 2022;14:174. https://doi.org/10.1186/s13195-022-01113-5 。项目：https://www.synapse.org/Synapse:syn30549757 。本阶段只获取公开元数据，个体文件下载403。
12. ADNI_EMORY_CSF_TMT_MS_Methods_20251007.pdf、ADNI_CruchagaLab_Methods_CSF_SomaScan7K_Proteomic_Data_20230817.pdf、ADNI_BSHRI_PLA_CSF_NULISA_CNS_METHODS_20260116.pdf。实际本地PDF路径、SHA和提取文本见ORIGINAL_ASSAY_METHOD_SOURCES.csv及evidence目录。TMT原始入口：https://www.synapse.org/Synapse:syn59804727 ，未认证访问403。
13. ADNI DATADIC_15May2026.csv：本阶段逐字段核对LOCLAB/CSFMETH的VISDATE和采样标志，见COLLECTION_CRF_FIELD_DEFINITIONS.csv。Registry回退语义不能替代标本采样日期。

部分网页需要JavaScript或出现访问检查；本地原始文件和官方REST快照分别留档。未使用未取得内容作为实际个体级分析依据。
""")
    write('字段说明与结果导航.md',"""
# 结果导航与字段解释

## 优先阅读

- THREE_CORE_CANDIDATE_EVIDENCE_MATRIX.csv：三个核心候选的汇总，含三模型共同关联/偏离q、直接肽段、质控、重复、验证状态。
- CORE_CANDIDATE_MODEL_EVIDENCE.csv：三候选×三模型；FAMILY34_INTEGRATED_EVIDENCE_ALL_MODELS.csv：全部34×三模型，MAPT不在家族。
- CORE_PEPTIDE_TECHNICAL_REPEATABILITY.csv：两个实验的五肽段技术重复。
- MRM2019_TARGET_TRANSITION_RAW_QC.csv与MRM2019_PQC24_RAW_LH_VARIATION.csv：原始转变QC和PQC比值CV。
- MRM2014_THREE_FOUR_ASSAY_DESIGN_AUDIT.csv：N/S/MRM与N/S/TMT/MRM分别审计，不混用df。
- INDEPENDENT_2019_NEW_RID_PEPTIDE_ASSOCIATIONS.csv、NEW_RID_LEAVE_ONE_BATCH_SUMMARY.csv：原五肽段结果复核和方向诊断。
- CANDIDATE_VALIDATION_STATUS_BY_LAYER.csv与FINAL_VALIDATION_GAP_AND_ACTION_REGISTRY.csv：哪些证据已支持、哪些验证未完成。

## 关键字段

q_BY_shared34：双向三检测同向关联的BY34结果；q_BY_H34：原相对响应零约束的BY34结果。q_BY_combined68_shared/H：将两类问题合为68个问题的既有敏感性，不能与独立家族q混称一个总体FDR。

shared_supported/H_deviation_supported按各自q≤0.05判定。BOTH_UNRESOLVED表示两问题均未获支持，不表示等效。SHARED_ONLY_H_UNRESOLVED表示共同关联获支持、偏离未定，不表示测量互换。

M0_equivalence_eligible与M0_q_equiv_working_margin05属于主模型旧等效分支；不是APOE/BOTH重新做的等效，数值在三模型行重复仅为背景。其他两个工作容差见PRIOR_EQUIVALENCE_ALL_MARGINS_FAMILY34.csv。

scaled_b单位是残差SD/每1个主分析pTau全局SD。scaled_SE_full_joint包括估计残差SD的不确定度。95人肽段区间沿用原分析条件于残差SD的定义，不能与全联合区间混称完全相同估计。

n_complete_pairs是有效重复对数，2019 LENLEQYSR为19，其余2019为20，2014为16。pearson_pair_bootstrap95是4000次成对抽样的百分位区间；Lin_CCC为样本矩描述。二者不证明绝对准确度。repeat_error_SD_equal_independent_error_working只在两次误差等方差、独立时有误差SD含义，不作为跨平台可靠性定量。

raw_LH_CV_percent_all_positive_no_flag_exclusion使用24个PQC的正且有限Light/Heavy比值，原始QC全0，因此未另按标志排除。它是池/实验的描述性CV，不是临床误差或全部平台CV。rawlog_formula中+1加在Light上、除Heavy前。

geometry_n80_df40_G6只描述几何门槛；strict_validation_eligible还需源标本日期、依赖假设及推断校准。所有2014行该字段均False。candidate_pathology_effect_fitted均False。

source_record_number/duplicate_pair_number是匿名顺序号，不是新增独立样本数。个体RID、样本编号、逐样本峰值与请求模板在00_restricted；Excel只汇总聚合结果。

SCIENTIFIC_QA_CHECKS.csv为本阶段检查，PRIOR_ALL_STAGE_SEAL_CHECKS_FINAL.csv为逐条旧封存核验；封存条目数和唯一文件路径数分别报告。WORKBOOK_SHEET_INDEX.csv列出全部CSV是否入Excel，过长单元格表保留CSV_ONLY，避免截断。

完整字段列表随各CSV表头保留；Excel INDEX提供来源文件、行列数和存储方式。
""")
    figindex=pd.read_csv(RESULTS/'FIGURE_INDEX.csv')
    figreview=[{'figure':r.figure,'PNG_and_PDF_present':(ROOT/r.PNG).exists() and (ROOT/r.PDF).exists(),'visual_review':'REVIEWED; forest legend and q labels separated; molecular legend added','time':now()} for r in figindex.itertuples()]
    csvout(RESULTS/'FIGURE_REVIEW_RECORD.csv',figreview)
    # Build a complete, indexed workbook; no truncation of long metadata fields.
    preferred=['THREE_CORE_CANDIDATE_EVIDENCE_MATRIX.csv','CORE_CANDIDATE_MODEL_EVIDENCE.csv','FAMILY34_INTEGRATED_EVIDENCE_ALL_MODELS.csv','CANDIDATE_VALIDATION_STATUS_BY_LAYER.csv','INDEPENDENT_2019_NEW_RID_PEPTIDE_ASSOCIATIONS.csv','CORE_PEPTIDE_TECHNICAL_REPEATABILITY.csv','MRM2019_PQC24_RAW_LH_VARIATION.csv','MRM2019_TARGET_TRANSITION_RAW_QC.csv','MRM2014_THREE_FOUR_ASSAY_DESIGN_AUDIT.csv','MRM2014_COLLECTION_CRF_DATE_TRIANGULATION.csv','NEW_RID_LEAVE_ONE_BATCH_SUMMARY.csv','FINAL_VALIDATION_GAP_AND_ACTION_REGISTRY.csv']
    files=[RESULTS/f for f in preferred]+[p for p in sorted(RESULTS.glob('*.csv')) if p.name not in preferred and p.name not in ['WORKBOOK_SHEET_INDEX.csv','DELIVERY_QA_CHECKS.csv']]
    wb=Workbook();idx=wb.active;idx.title='INDEX';idx.append(['sheet','source_CSV','rows','columns','storage','note'])
    index=[]
    for number,p in enumerate(files,1):
        d=pd.read_csv(p,dtype={c:str for c in pd.read_csv(p,nrows=0).columns if 'seed' in c.lower()})
        too_long=any(d[c].dropna().astype(str).str.len().gt(30000).any() for c in d)
        if len(d)>20000 or too_long:
            sheet='CSV_ONLY';note='Retained as CSV to avoid long-cell truncation or excessive sheet size'
        else:
            sheet=f'S{number:02d}_'+p.stem[:25];sheet=sheet[:31];note=''
            ws=wb.create_sheet(sheet);ws.append(list(d.columns));ws.freeze_panes='C2';ws.auto_filter.ref=f'A1:{get_column_letter(len(d.columns))}{len(d)+1}'
            for row in d.itertuples(index=False,name=None):
                output=[]
                for v in row:
                    if isinstance(v,np.generic):v=v.item()
                    if isinstance(v,float) and not math.isfinite(v):v=None
                    if isinstance(v,int) and abs(v)>2**53:v=str(v)
                    output.append(v)
                ws.append(output)
            for cell in ws[1]:cell.font=Font(color='FFFFFF',bold=True);cell.fill=PatternFill('solid',fgColor='23536B');cell.alignment=Alignment(wrap_text=True,vertical='center')
            ws.row_dimensions[1].height=42
            for j,c in enumerate(d.columns,1):
                ws.column_dimensions[get_column_letter(j)].width=14 if c in ['gene','accession','model','assay'] else 23
                for row in ws.iter_rows(min_row=2,min_col=j,max_col=j):
                    if isinstance(row[0].value,float):row[0].number_format='0.000000E+00' if row[0].value!=0 and abs(row[0].value)<.001 else '0.000000'
                    if 'seed' in c.lower():row[0].number_format='@'
        entry={'sheet':sheet,'source_CSV':p.name,'rows':len(d),'columns':len(d.columns),'storage':'CSV_ONLY' if sheet=='CSV_ONLY' else 'EXCEL_AND_CSV','note':note}
        index.append(entry);idx.append(list(entry.values()))
    for cell in idx[1]:cell.font=Font(color='FFFFFF',bold=True);cell.fill=PatternFill('solid',fgColor='23536B')
    idx.freeze_panes='A2'
    for col,width in [('A',34),('B',65),('C',12),('D',12),('E',23),('F',70)]:idx.column_dimensions[col].width=width
    bookpath=ROOT/'核心候选证据整合与验证_结果与审计.xlsx';wb.save(bookpath)
    csvout(RESULTS/'WORKBOOK_SHEET_INDEX.csv',index)
    versions={name:importlib.metadata.version(name) for name in ['numpy','pandas','scipy','matplotlib','openpyxl','pypdf','requests']}
    save_json(RESULTS/'ENVIRONMENT.json',{'python':sys.version,'platform':platform.platform(),'packages':versions})
    write('README.md',f"""
# 第9阶段交付入口

已完成现有数据可实施的核心候选证据整合、独立计算复核、MRM质控/重复性分析、分子映射和验证资格审计。严格相对模式验证、外部个体复现及具体机制识别尚未完成。

- [关键决策摘要](关键决策摘要.md)
- [中文总报告](STAGE9_MAIN_REPORT_ZH.md)
- [Excel结果与审计](核心候选证据整合与验证_结果与审计.xlsx)
- [论文英文方法/结果建议](论文方法与结果建议_EN.md)
- [字段与导航](字段说明与结果导航.md)
- [缺失材料与执行条件](缺失材料最小清单与验证执行条件.md)
- [参考文献与来源](参考文献与来源说明.md)
- [可复现说明](可复现与运行说明.md)

{scientific_summary}。全部结果CSV在results，7组PNG/PDF在figures。个体标识和原始目标峰值只在00_restricted，不应与聚合报告一起作为公共附件发送。

封存后默认运行code/run_all.ps1只读校验；重新分析应在新副本进行。
""")
    write('可复现与运行说明.md',"""
# 可复现与运行

## 默认只读核验

在PowerShell执行本目录下code/run_all.ps1，或python -X utf8 code/verify_delivery.py。默认只核对SHA256及完成状态，不重跑模拟、不更新下载、不覆盖封存结果。科学QA及交付QA结果分别保存在results对应CSV/JSON。

## 分析依赖与顺序

Python环境及库版本见results/ENVIRONMENT.json。无需R，也未安装新依赖。所有统计实现均在本阶段代码，未导入旧阶段会写文件的脚本。

在项目根目录下的新同级副本内依次执行：
bootstrap.py → independent_analysis.py → audit_metadata.py → audit_replicates_qc.py → audit_molecular.py → triangulate_dates.py → audit_public_sources.py → integrate_evidence.py → make_figures.py → validate_science.py → build_delivery.py → validate_delivery.py → seal_delivery.py。

重新计算应使用code/reproduce_in_copy.py创建新副本；它保留旧阶段只读来源，复制本阶段代码、协议与公开快照，运行前核验输入。公开来源默认使用冻结缓存，网络更新必须另存新分析分支。部分旧来源是用户本地绝对路径，移机后须先恢复相同授权数据位置并核对SHA。

先前未取得的外部个体文件不属于分析输入；没有通过替代渠道绕过403。PUBLIC_DATA_DOWNLOAD_AUDIT.json保留当时访问状态，默认重算不重新申请权限。

## 计算定义

三检测同向p+=max(P_k+)；p-=max(P_k-)；p_shared=min(1,2min(p+,p-))，34家族BY；合并68敏感性另列。原H约束继续包含全部系数、残差方差和协方差的联合不确定度，SUM3 CR1不称CGM或CR2。

95人肽段分析沿用五肽段家族、单侧正方向、SUM2 CR1和最小簇数减1的Student t。删板/运行只报告系数方向和变幅。重复相关区间每个肽段/实验4000次固定种子成对抽样，共40000次区间抽样，非新增病理零假设模拟。

2014身份核验在320特征保持缺失一致、6有效数字末位1单位容差内完成。2014病理关联未拟合，所有三/四检测几何表均为访视日期推断下的条件性审计。

## 文件组织

results：聚合统计、科学QA、源文件摘要、旧封存核验和工作簿索引。evidence：官方条目、公开元数据、方法提取和附件图例。00_restricted：个体标识、匹配记录、目标原始峰和模板。后者目录名不设置ACL。

SHA256_MANIFEST.csv覆盖交付文件，排除自身、临时文件和__pycache__。完成状态明确区分执行完成与严格验证未完成。科学图PNG/PDF已逐图检查；报告引用保持相对路径。
""")
    log(f'Delivery built:main report,decision memo,English manuscript text,field guide,{len(index)} indexed tables,{len(wb.sheetnames)} workbook sheets')

if __name__=='__main__':
    assert_unsealed()
    run()
