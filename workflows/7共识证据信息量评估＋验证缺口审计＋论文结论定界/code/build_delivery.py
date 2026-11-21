"""Summarize existing fixed results. No new clinical effects are fitted."""
from config import *
import argparse, importlib.metadata
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

PURPOSES={
 'CONSENSUS_INFORMATION_BOTTLENECKS':'逐蛋白、逐容差的信息瓶颈和冻结证据状态',
 'INFORMATION_BOTTLENECK_COUNTS':'34蛋白各容差瓶颈计数',
 'EVIDENCE_STATE_COUNTS':'冻结偏离、等效及未定状态计数',
 'FIELLER_PRECISION_AND_DENOMINATOR_INFORMATION':'Fieller集合、分母信息及精度',
 'CONDITIONAL_PRECISION_EXPANSION_SCENARIOS':'固定当前点和V/k的精度情景；不是功效',
 'STAGE6_EQUIVALENCE_REPRODUCTION':'105项原等效P复现',
 'POWER_SIMULATION_SUMMARY':'真等效DGP的门槛和等效成功率',
 'BOUNDARY_SIMULATION_SUMMARY':'新增边界单元两个阈值的统计量',
 'POWER_WORST_SCENARIO_INFORMATION':'每模板在三DGP中的最低成功率及最低Wilson下界',
 'POWER_TARGET_ATTAINMENT_COUNTS':'达到80%/90%要求的模板数',
 'POWER80_90_INDEPENDENT_PANEL_REQUIREMENTS':'有限网格首次达到80%/90%的独立面板倍数',
 'REAL_TEMPLATE_BOUNDARY_CHECKS':'各边界单元门槛前后错误等效及Wilson区间',
 'REAL_TEMPLATE_BOUNDARY_DECISIONS':'15模板×2端点×3容差的有限名义边界核查',
 'STAGE6_BOUNDARY_EXTREME_TAIL_PRECISION':'原边界模拟在小alpha阈值的重计数及MC精度',
 'POWER_ENGINE_INDEPENDENT_CHECKS':'扩展设计、投影、SUM3及IUT独立计算核查',
 'INDEPENDENT_SIMULATION_RECOUNTS':'720个新增模拟文件独立重计数',
 'VALIDATION_GAP_REGISTRY':'10个验证缺口、已补证据、所需资料及主张影响',
 'EXTERNAL_VALIDATION_OPPORTUNITIES':'外部资源可用性与适用性；没有新个体级复现',
 'EVIDENCE_INTERPRETATION_BOUNDARIES':'共同变化、等效、病理关联的用途边界',
 'MRM2014_CROSSWALK_AUDIT_COUNTS':'早期严格误差比较和候选日期；最终身份以FINAL表为准',
 'MRM2014_FINAL_IDENTITY_STATUS':'发布精度下的最终唯一身份及技术连接',
 'MRM2014_LINKAGE_STATUS':'2014技术和日期连接最终简表',
 'MRM2014_INFERRED_VISIT_DESIGN_SCREEN':'假定日期桥接成立的设计几何；无候选效应拟合',
 'MRM2014_RELEASE_SCALE_IDENTITY_DIAGNOSTICS':'发布尺度诊断；最终匹配没有使用仿射变换',
 'MRM2014_RELEASE_PRECISION_PROFILE':'六有效数字发布精度及版本差异',
 'MRM2014_SAMPLEINFO_FIELDS':'工作簿字段聚合覆盖',
 'EXTERNAL_EMORY2022_FILTERED_GENE_PLATFORM_COVERAGE':'公开过滤表基因覆盖；非个体完整度',
 'EXTERNAL_EMORY2022_CANDIDATE_ANNOTATION_HITS':'三候选在两种MS处理注释表的命中',
 'EXTERNAL_EMORY2022_ANNOTATION_SHEETS':'实际读取的两张公开注释表',
 'INPUT_SOURCE_CHECKS':'31个既有原始输入SHA核查',
 'ADDITIONAL_METADATA_INPUTS':'新增登记的3个只读本地元数据源',
 'METADATA_ARCHIVE_MEMBERS':'元数据压缩档成员',
 'REFERENCE_CODE_ORIGIN':'4份原方法代码来源及逐字节一致性',
 'PRIOR_SEALED_FILE_CHECKS':'前期封存文件核查',
 'PRIOR_SEALED_FILE_CHECKS_FINAL':'3469个前期封存条目最终独立重查',
 'INDEPENDENT_QA_CHECKS':'独立科学计算、来源和边界核查',
 'FIGURE_INDEX':'7组PNG/PDF路径及中文图注'}

FIELDS={
 'accession':'UniProt登录号；非受试者标识',
 'in_family34':'是否属于冻结34项家族；不缩为可评价子集',
 'margin':'沿用阶段6已知结果后的工作容差0.2/.35/.5',
 'delta_S / delta_T':'Soma/TMT相对NULISA病理响应偏离，单位为相应检测残差SD/病理SD',
 'equivalence_eligible_Stage6':'沿用原结构、协方差信息、Fieller和有限校准门槛',
 'nominal_equiv_eligible':'通过评价门槛且原IUT P≤.05',
 'BY_equivalence_Stage6':'冻结34家族BY等效q≤.05；本轮全False',
 'frozen_exact_zero_deviation_BY':'阶段5精确零偏离BY发现；不是宽容差不等效的判据',
 'point_structure_compatible':'正向单因子点结构且独特方差非负；不证明误差独立',
 'threshold':'NOMINAL或BY_RANK1_SUFFICIENT；后一项不是实际BY排序功效',
 'one_sided_alpha':'单侧阈值；Fieller置信水平1−2alpha',
 'bounded':'Fieller是否为有限有界区间；False不替换成窄区间',
 'k_independent_panels':'总的独立未来面板倍数；新个体和新技术簇',
 'gamma_common_P':'预设共同因子病理系数0/.55；非观测偏离事后真值',
 'future_n_scenario':'模板人数×面板倍数；不是无条件样本量建议',
 'rate':'包含原评价门槛的成功比例',
 'raw_success_count':'未加评价门槛的IUT成功次数',
 'gate_success_count':'同时通过门槛及IUT的次数',
 'worst_scenario_rate':'三预设DGP的最低成功比例',
 'worst_Wilson95_lower':'三DGP的Wilson95下界最小值；非所有情景同时95%保证',
 'both_endpoint_finite_boundary_checks_passed':'仅说明该模板/容差的SLOPES3、gamma.55、k1/8名义边界通过',
 'boundary_checked_only_SLOPES3_gamma_point55_at_k1_and8':'边界覆盖范围说明，不表示当前行条件已单独边界模拟',
 'fixed_cluster_variance_fraction_scenario_not_estimated':'假定固定簇方差份额0/.25/.5，非实际ICC',
 'required_full_independent_information_multiplier':'固定theta、V/k下所需信息倍数；inf表示此情景不能有限达到',
 'required_subject_multiplier_in_floor_scenario':'假定簇方差下限时的条件人数倍数；不是实际功效',
 'strict_same_day_replication_eligible':'2014均False；同访视不替代MRM标本日期证据',
 'fresh_candidate_association_fitted':'2014均False，本轮无新候选效应拟合',
 'geometry_passes_S6_n80_df40_G6':'仅数量/秩/簇数门槛，不证明推断有效',
 'n_Wilson_lower_meets_target':'各DGP下界均≥目标的模板数',
 'status':'按所属表解释；NOT_REACHED不表示任何更大规模均不可能'}

def table(path):return pd.read_csv(path,dtype={'seed':str})
def read(name):return table(RESULTS/(name+'.csv'))
def put(name,t):(ROOT/name).write_text(t.strip()+'\n',encoding='utf-8')
def md(df):
 def val(x):
  if pd.isna(x):return '—'
  if isinstance(x,(float,np.floating)):return ('∞' if x>0 else '−∞') if not np.isfinite(x) else str(int(x)) if x==int(x) else f'{x:.4g}'
  return str(x).replace('|','/').replace('\n',' ')
 return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(val(x) for x in row)+' |' for row in df.itertuples(index=False,name=None))

def historical_labels():
 p=RESULTS/'MRM2014_CROSSWALK_AUDIT_COUNTS.csv';d=pd.read_csv(p)
 d.audit_item=d.audit_item.replace({
 'L01_rows_with_at_least_one_identical_vector':'INITIAL_abs_tol_1e-5_L01_rows_with_any_match',
 'L01_rows_with_unique_identical_vector':'INITIAL_abs_tol_1e-5_L01_rows_with_unique_match',
 'L01_rows_with_multiple_identical_vectors':'INITIAL_abs_tol_1e-5_L01_rows_with_multiple_matches'})
 csvout(p,d)

def documents():
 historical_labels()
 qa=json.loads((RESULTS/'INDEPENDENT_QA_SUMMARY.json').read_text(encoding='utf-8'))
 info=read('CONSENSUS_INFORMATION_BOTTLENECKS')
 family=info.loc[info.in_family34&info.margin.eq(.5)]
 eligible=family.loc[family.equivalence_eligible_Stage6]
 nominal=eligible.loc[eligible.nominal_equiv_eligible,'gene'].tolist()
 counts=read('POWER_TARGET_ATTAINMENT_COUNTS')
 show=counts.loc[counts.gamma_common_P.eq(.55)&counts.margin.eq(.5)&counts.target.eq(.8),['k','threshold','median_worst_scenario_rate','n_Wilson_lower_meets_target']].sort_values(['k','threshold'])
 show.columns=['独立面板倍数','检验阈值','15模板最差情景率的中位数','达80%要求的模板数']
 req=read('POWER80_90_INDEPENDENT_PANEL_REQUIREMENTS')
 req=req.loc[req.gamma_common_P.eq(.55)&req.margin.eq(.5)&req.target_probability.eq(.8)&req.threshold.eq('BY_RANK1_SUFFICIENT')&req.first_tested_k_with_all3_scenario_Wilson_lower_above_target.notna(),['gene','first_tested_k_with_all3_scenario_Wilson_lower_above_target','future_n_scenario_at_that_k']].sort_values(['first_tested_k_with_all3_scenario_Wilson_lower_above_target','gene'])
 req.columns=['模板基因','网格内首次达80%的k','条件总人数']
 zero_upper=read('STAGE6_BOUNDARY_EXTREME_TAIL_PRECISION').zero_event_upper_if_zero_events.iloc[0]
 gaps=read('VALIDATION_GAP_REGISTRY')[['gap_id','issue','status','minimum_next_evidence_or_action']]
 report=f"""
# 共识证据信息量评估＋验证缺口审计报告

交付日期：2026-10-27。基于已知阶段1–6结果的审计与条件设计研究，不是前瞻注册；没有增加新的真实候选−病理效应检验。前期模型、发现家族、BY、容差和门槛均保持不变。

## 对研究决策的影响

1. 原34蛋白家族仍为0个BY等效，不能推断所有蛋白不一致。15项具有正向单因子点结构，11项通过完整等效评价门槛；最宽工作容差±0.5下6项名义等效、0项BY等效。正面共识证据受结构、精度和多重比较共同限制。
2. 真等效条件模拟显示当前规模经常证明不了等效。γ=.55、±.5时，15模板最差情景检出率的中位数为名义36.25%、保守BY足够条件0%；独立设计扩展8倍后为99.05%和80.35%。这些是明确假设下的设计情景，不预测实际发现。
3. 2014 MRM的303条发布记录已能在六有效数字发布精度下唯一连接原工作簿和技术信息。若采样日期桥接成立，NPTX1可有144人/df78、KLK6可有141人/df75。**MRM特异采样日期仍缺直接证据，不能列作复现完成。**
4. 外部Emory2022公开注释确有PDGFRB/KLK6的Soma/Olink/MS基因级三平台覆盖，提供新线索；三平台完整队列上限35，每基因完整人数和肽段身份未核实，尚不能按原冻结规则正式复现。

## 共识定义及证据层次

共同变化、病理响应等效和病理关联是三个层次。共同变化表示检测随人变化具有共同结构；响应等效表示检测灵敏度校正后的病理响应偏离足够小；病理关联则要求该共同方向确随独立病理指标变化。前两项不能自行证明真实完整蛋白浓度已识别。

令N/S/T为NULISA、SomaScan、TMT，b为同一调整模型中标准化独立pTau系数，c为残差跨检测协方差，v为残差方差：

H₁=bS·cNT−bN·cST；H₂=bT·cNS−bN·cST。

dS=H₁/(cNT·√vS)，dT=H₂/(cNS·√vT)。

等效评价dS/dT是否均位于既定容差，而非三个原始系数相同、仅相关为正或P均不显著。任意相关测量误差下，三次单检测仍不能唯一分离真实丰度与检测特异病理作用。本轮不生成个人F/D评分，不把分子形式或裂解机制称为已证实。

精确零偏离与宽容差等效是两个证据轴：很小但可检出的非零偏离，可以同时满足实用等效，不应强制互斥。本项目当前无同时成立实例：原主分析3项偏离发现（PDGFRB、NPTX1、KLK6），等效0，其余31项未定。3项仍需沿用阶段5的APOE调整敏感性说明，不能说所有设定均稳健。

## 真实数据的信息瓶颈

沿用35个冻结设计（34家族＋MAPT控制），同一年龄/性别及技术调整、SUM3 CR1和完整θ=(bN,bS,bT,vN,vS,vT,cNS,cNT,cST)联合影响函数；105项既有等效P全部复现。

34项中19项点结构不相容，另4项有点结构但缺足够的正协方差/分母信息，合计23项门槛不合格。这些状态都不等同生物学不一致。11项可评价蛋白：{', '.join(eligible.gene)}。

{md(read('INFORMATION_BOTTLENECK_COUNTS').rename(columns={'margin':'工作容差','main_bottleneck':'主要瓶颈','n_proteins':'蛋白数'}))}

±.5下名义等效为：{', '.join(nominal)}。±.35仅1项，±.2为0。名义证据不能改称BY共识发现；三档均为阶段6已知结果后的工作容差，没有已建立的生物学最小重要差异依据，不据结果扩大容差挑发现。

Fieller表提供90%置信集合、分母信号/SE和严格阈值集合。无界集合表示分母信息不足，不替换为delta窄区间。可覆盖区间的最小对称容差仅作精度描述，不是新选阈值。

## 条件精度与真正等效功效

精度外推保持当前θ不动，仅假设完整联合方差V/k下降。对容差m及正分母D，要求mD−H>0及mD+H>0，信息倍数由四个不等式和正协方差/Fieller门槛共同决定。点已超容差、结构不相容或分母不合格者不会得到虚假的有限成功样本量。这不是功效。

另假定固定簇方差份额f=0/.25/.5，V(r)/V=f+(1−f)/r。达到信息倍数k需r=(1−f)/(1/k−f)；分母不正时有限簇内扩样不能达到。f未从数据估计，仅说明旧簇内加人不能自动等同新增独立批次。

功效使用真等效DGP：15个残差协方差点可形成正向单因子模板，λ由跨协方差获得，ψ=v−λ²非负。其余19家族蛋白不截负ψ、不借其他蛋白补齐。15模板全部在主家族，MAPT没有适用模板。设b=λγ，γ=0/.55，真实H为0，不把观测偏离当真值。γ0也说明等效可以没有病理关联。

三DGP为IID高斯、三技术轴随机年龄斜率高斯、相应t5。斜率SD.15为固定压力情景，共同及检测误差系数独立，投影后二阶矩归一化到模板C。k1/2/4/8为总共k份独立未来面板，每份新个体和新技术簇；病理系数跨面板共同估计，nuisance面板内调整，df=k(n−rankZ)−1，CR1按新增簇重算。

15×4×3×2×2000=720000次功效重复。三档容差/两阈值共享同一重复，不重复计数。λ/ψ及协方差模板被条件性固定，没有积分其估计不确定性；DGP有限，不是通用样本量计算。

名义阈值.05。α*=.05/(34H₃₄)={ALPHA_BY_SUFFICIENT:.10g}，是无论其他P排序如何都足以通过原BY的保守条件。实际BY允许较后排序发现用更宽阈值，故**α*成功概率不是实际BY排序功效，也不是必需P值**。没有虚构跨蛋白相关，也没有缩小34家族。

下表为γ.55、±.5。先对每模板取三DGP最低成功比例，再对15模板取中位数。达标要求每个DGP的Wilson95下界均≥.80；这是单元MC不确定性描述，不是所有情景同时95%保证。

{md(show)}

保守BY足够条件下，有限网格中7个模板达到80%要求：

{md(req)}

例如NPTX1的k4/2216人，只表示若未来真为该等效模板、γ.55、±.5且新技术簇独立增加，网格首次达到要求的是4份设计；不能推断真实NPTX1现有偏离随扩样消失。未在8份内达到不证明任何更大规模均不可能。

## 有限边界校准的范围

另做15模板×k端点1/8×3容差×2对比×正负边界×2000=720000次独立重复，固定SLOPES3和γ.55。目标检测均值加±m√v，使其偏离恰处于边界，另一偏离为0。

90个模板−端点−容差组合均通过预设名义核查：四边界单元原始误等效Wilson95上界≤.06、失败率≤.005。最高原始误等效率.85%，最高Wilson上界1.357%，表现保守。只覆盖这些有限情景；k2/4、IID/t5及γ0没有新增独立边界证明。

原阶段6边界文件亦在α*重计数。每单元2000次，零事件Wilson95上界仍为{zero_upper:.6f}，高于α*={ALPHA_BY_SUFFICIENT:.6f}。零次不能当真实概率0，也不能精确认证这个极小尾部。本轮不升级普遍P值有效性或误差识别结论。

## 2014历史MRM：技术连接闭合，采样日期未闭合

L01发布CSV有303记录、287基线RID；L02整合工作簿有LIMS_ID、External_Reference、处理日、去高丰度柱和进样序号。CSV-Export经实验室键连接技术表，肽段矩阵列头对应LIMS_ID。两个2014发布为同一实验，不能计两次独立验证；2019为另一实验族。

早期全320肽段绝对差≤1e−5无匹配。发布精度检查显示L01为六有效数字，L02保留更多小数；精确六位舍入仍有每行0–6个肽段末位略异，最大原始差约.541末位单位，不能说浮点存储机制已证实。最终要求同RID/访视、缺失一致、全320值差不超过六有效数字末位1单位、且唯一匹配。**规则在精度检查后披露，不冒充事先注册。**最终匹配不使用仿射变换，不依据pTau关联挑记录。

303/303均唯一恢复技术记录，共15处理日、2去高丰度柱。这是来源精度下一致，不能说数值逐位完全相同。早期严格匹配0结果明确标为INITIAL；最终以FINAL_IDENTITY_STATUS为准。

287个RID/基线访视在另一标志物MASTER中有唯一采样日候选；149人与阶段5同访视且候选日相同。仍无MRM自身LIMS/标本号到采样日的直接桥接。本地ADNIMERGE及生物样本档案检查未补足。**MRM特异采样日期可直接核实数量仍为0。**

假定访视桥接成立，按处理日/进样序号确定技术重复优先级，可得：

| 候选 | 独立肽段 | 条件有效人数 | rankX | df | N/S/MRM日簇 |
| --- | --- | --- | --- | --- | --- |
| NPTX1 | FQLTFPLR、LENLEQYSR | 各144 | 66 | 78 | 11/36/15 |
| KLK6 | ESSQEQSSVVR、LSELIQPLPLER、YTNWIQK | 各141 | 66 | 75 | 11/36/15 |

以上仅设计几何，没有新候选效应拟合。虽超过n80/df40/各轴6簇门槛，仍有最大杠杆约1的饱和技术行；两去高丰度柱是否形成跨日依赖需核实，2014实际设计未有限簇校准。不能将几何通过写作推断/复现通过。

2019同日90–94人/df20–24仍仅描述用途。其95个主样本外新RID属于同一ADNI，单肽段−pTau关联不是三检测相对模式或真正外部队列复现。PDGF-BB配体不能验证PDGFRB受体。

## 外部来源审计

[Dammer等2022](https://doi.org/10.1186/s13195-022-01113-5)补充表4/5已下载并登记SHA，仅审计平台/基因注释。去高丰度及未去高丰度MS过滤表一致：PDGFRB、KLK6列于Soma/Olink/MS；NPTX1列于Soma/MS，Olink未列。未列不表示面板绝无靶标；基因注释不代替受体特异肽段、蛋白分组或个体完整度。

该研究36人，三平台完整队列上限35，不是各基因完整人数。论文数据入口为[Synapse](https://www.synapse.org/3platformEmory)。无NULISA、规模低于冻结n80；本轮未重析其个体配对矩阵，仅为外部测量支持机会。

[Kowal等2026预印本](https://doi.org/10.64898/2026.09.02.26361977)为67例混合神经退行疾病的Soma11K/NULISA双平台，声明可合理请求数据；无第三质谱，候选/pTau覆盖未核实，本轮未外发申请。[AMP-PDRD官方资源](https://amp-pdrd.org/data)有CSF Olink及MS/肽段元数据线索；PD队列不能默认验证AD特异假设，本轮未账号操作或取得配对个体数据。

## 缺口与下一步顺序

{md(gaps)}

先追踪2014的LIMS_ID/External_Reference→specimen_id/aliquot_id→collection_date/time和去高丰度/重复流程，可能激活141–144人现有配对资源。

日期核实后，冻结2014候选、肽段/重复规则、实际技术依赖和检验家族，校准该设计的有限簇推断及识别敏感性，再拟合候选效应。若日期取不到，应保留访视层面的探索用途；另定义误配敏感性，不能改称同日验证。

正面共识的新研究需同步增加独立个体及技术批次，在结果前确定有科学依据的容差。模板功效帮助限定设计范围，不能据此选“必然共识”蛋白。原TMT肽段、结合区域、分装和技术重复仍需资料或实验。

## 当前论文证据边界

可报告：冻结调整下三项病理相关相对响应偏离；等效证据受结构、分母精度及多重比较限制；历史MRM技术连接补齐并定位了可执行的资料需求。APOE敏感性和识别条件须与结果并列。

尚不支持：所有34项均不一致；已分离真实浓度或特定裂解/活化形式；三候选已严格同日/外部模式复现；PDGFRB已由PDGF-BB验证；0 BY等效证明没有共同信号。形式机制保留为待检验解释。

## 可复核性与交付

独立科学核查{qa['passed']}/{qa['checks']}通过，包括35个真实设计独立SVD重拟合及SUM3、105项等效、720文件144万次重计数、303行全320肽段及技术元数据。3469个前期封存条目未改。最终数量见results/INDEPENDENT_QA_SUMMARY.json及交付校验表。

见[汇总工作簿](共识信息量与验证缺口_汇总.xlsx)、[字段导航](字段说明与结果导航.md)、[参考依据](参考文献与来源说明.md)。figures提供7组PNG/矢量PDF。00_restricted包含个体连接、缓存及含公开个体附表的外部原始下载；目录名本身不设置系统权限。

封存后code/run_all.ps1默认只读查SHA。重建在旁置副本执行，见[运行说明](可复现与运行说明.md)。
"""
 put('共识证据信息量评估与验证缺口审计报告.md',report)
 put('关键决策摘要.md',f"""
# 下一步决策摘要

2026-10-27。两项审计完成；没有新增真实候选−病理效应拟合。

1. **0个BY等效不能解释成全部不一致。** 34项中15项点结构相容，11项完整可评价；±.5工作容差下6项名义等效、0项BY。原3项偏离发现仍有APOE调整敏感性。
2. **若真等效，现规模也常证明不了。** γ.55、±.5时15模板最差情景率中位数为名义36.25%、保守BY足够条件0%；达到各DGP的Wilson95下界≥80%的模板为3/15和0/15，8份独立设计为10/15和7/15。α*不是实际BY排序功效，旧批次内加人不能直接等同独立面板。
3. **2014 MRM值得先补资料。** 303记录已在发布精度下唯一连接，恢复15处理日/2去高丰度柱；若日期桥接成立，NPTX1可144人/df78，KLK6可141人/df75。MRM特异日期直接证据仍为0，未新拟合效应或校准此设计。
4. **PDGFRB有外部覆盖新线索，尚未验证。** Emory2022过滤注释列PDGFRB/KLK6于Soma/Olink/MS；人数上限35，无NULISA，每基因完整数/肽段未核实。NPTX1的Olink仅为过滤列表未列。

下一步先取2014实验室号→标本/分装号→采样日期及处理流程；成功后冻结实际设计、核实柱依赖并校准，再检验候选。正面共识需要更强独立信息和结果前确定的科学容差，不通过放宽现有容差补发现。

科学核查{qa['passed']}/{qa['checks']}通过；720新增单元、144万次模拟；3469个前期封存条目未改。见[完整报告](共识证据信息量评估与验证缺口审计报告.md)。
""")
 nav=['# 字段说明与结果导航','','2026-10-27。CSV为UTF-8 BOM。空值为不适用/未计算；inf为规定条件下无有限倍数，不是缺失。','','## 表级导航','']
 for p in sorted(RESULTS.glob('*.csv')):
  if p.name in ['DELIVERY_QA_CHECKS.csv','WORKBOOK_SHEET_INDEX.csv']:continue
  nav.append(f'- [{p.name}](results/{p.name})：{PURPOSES.get(p.stem,"辅助来源/诊断记录")}。')
 nav+=['','## 关键字段','','| 字段 | 含义 |','| --- | --- |']+[f'| {k} | {v} |' for k,v in FIELDS.items()]
 nav+=['','## 读取注意','','先读瓶颈和缺口登记，再读精度/功效。INITIAL_abs_tol_1e-5为早期严格比较；最终身份见FINAL_IDENTITY_STATUS和IDENTITY_ROUNDING_RULE，不用中间失败覆盖最终状态。','',
 '独立重复按file去重。720个新增文件各2000次，共144万；各文件可有多容差/阈值评价。bounded=False且上下界空表示未形成有界区间，不能当0。同访视和设计几何合格不是已核实同日复现。']
 put('字段说明与结果导航.md','\n'.join(nav))
 versions={p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','statsmodels','matplotlib','openpyxl','pypdf','pyarrow']}
 save_json(RESULTS/'SOFTWARE_VERSIONS.json',{'python':sys.version,'packages':versions,'date':'2026-10-27'})
 log('Chinese report, decision memo and navigation generated from fixed results')

def workbook():
 paths=[p for p in RESULTS.glob('*.csv') if p.name not in ['DELIVERY_QA_CHECKS.csv','WORKBOOK_SHEET_INDEX.csv']]
 priority=['CONSENSUS_INFORMATION_BOTTLENECKS','INFORMATION_BOTTLENECK_COUNTS','EVIDENCE_STATE_COUNTS','FIELLER_PRECISION_AND_DENOMINATOR_INFORMATION','CONDITIONAL_PRECISION_EXPANSION_SCENARIOS','POWER_TARGET_ATTAINMENT_COUNTS','POWER80_90_INDEPENDENT_PANEL_REQUIREMENTS','POWER_WORST_SCENARIO_INFORMATION','VALIDATION_GAP_REGISTRY','MRM2014_FINAL_IDENTITY_STATUS','MRM2014_INFERRED_VISIT_DESIGN_SCREEN','EXTERNAL_VALIDATION_OPPORTUNITIES']
 paths.sort(key=lambda p:(priority.index(p.stem) if p.stem in priority else len(priority),p.stem))
 records=[]
 for i,p in enumerate(paths,1):
  d=table(p);records.append({'sheet':f'S{i:02d}_{p.stem[:27]}','source':'results/'+p.name,'rows':len(d),'columns':len(d.columns),'purpose':PURPOSES.get(p.stem,'辅助来源/诊断记录')})
 index=pd.DataFrame(records);csvout(RESULTS/'WORKBOOK_SHEET_INDEX.csv',index)
 with pd.ExcelWriter(ROOT/'共识信息量与验证缺口_汇总.xlsx',engine='openpyxl') as writer:
  index.to_excel(writer,sheet_name='INDEX',index=False)
  for r in records:table(ROOT/r['source']).to_excel(writer,sheet_name=r['sheet'],index=False,inf_rep='inf')
  for ws in writer.book.worksheets:
   ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions;ws.sheet_view.zoomScale=85;ws.row_dimensions[1].height=44
   for cell in ws[1]:cell.fill=PatternFill('solid',fgColor='315D73');cell.font=Font(color='FFFFFF',bold=True);cell.alignment=Alignment(wrap_text=True,vertical='top')
   for row in ws.iter_rows(min_row=2):
    for cell in row:
     if isinstance(cell.value,str):cell.data_type='s';cell.alignment=Alignment(vertical='top',wrap_text=len(cell.value)<220)
     elif isinstance(cell.value,float):cell.number_format='0.000000E+00' if 0<abs(cell.value)<1e-4 else '0.0000000000'
   for j in range(1,ws.max_column+1):
    length=max((len(str(ws.cell(i,j).value or '')) for i in range(1,min(ws.max_row,55)+1)),default=14)
    ws.column_dimensions[get_column_letter(j)].width=min(52,max(14,length+2))
  ws=writer.book['INDEX'];ws.column_dimensions['B'].width=64;ws.column_dimensions['E'].width=85
  for i,r in enumerate(records,2):ws.cell(i,1).hyperlink=f"#'{r['sheet']}'!A1";ws.cell(i,1).font=Font(color='0563C1',underline='single')
 print({'workbook_data_sheets':len(records),'total_sheets':len(records)+1,'aggregate_rows':int(index.rows.sum())})

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--docs',action='store_true');parser.add_argument('--workbook',action='store_true');args=parser.parse_args()
 if args.docs:documents()
 if args.workbook:workbook()
