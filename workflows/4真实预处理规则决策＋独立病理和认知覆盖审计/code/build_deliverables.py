from __future__ import annotations
from audit_common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logging
logging.getLogger('fontTools.subset').setLevel(logging.ERROR)
from openpyxl import load_workbook
from openpyxl.styles import Font,PatternFill,Alignment
from openpyxl.formatting.rule import DataBarRule

def rules():
 items=[
 ('R01','索引','确定','571个RID唯一的精确同日三平台索引；不扩展到583个仅RID重叠者','同日不等于同一分装；Soma缺PTID，保留已核实的RID日期链','Stage1注册表'),
 ('R02','身份','确定','使用封存Stage1样本注册表和布尔可用性，不重新拼接近似蛋白样本','NULISA–TMT PTID一致；跨平台分装不能据此推断','Stage1身份审计'),
 ('R03','NULISA尺度','确定','已发布NPQ直接保留；不重复取log，不当作绝对浓度','保留QC_PASS及原PlateID','I01/I02/Stage1'),
 ('R04','NULISA重复','确定','同RID同日技术重复只合并已确认CSF重复，NPQ等权均值；保留重复来源','唯一已确认跨板CSF重复不在571人内；不能按RID合并纵向访视','Stage1'),
 ('R05','LOD标志','确定','NPQ数值与逐板LOD检出标志分开存储；TargetDetectability不作逐人检出','不以0、LOD或LOD/2替换低于LOD的NPQ','I01/I02/Stage1'),
 ('R06','低检出处理','候选需验证','数据层保留有效NPQ；数值完整路线作为首选评估路线，逐人LOD删除仅敏感性','低检出MME/YWHAZ/S100A12/BASP1不直接列为确证性高质量蛋白；保留数值不证明可定量','Stage1/Stage4交集'),
 ('R07','Soma尺度','确定','正的RFU作log2；缺失、非有限和非正值不取log、不填0','保留Soma原PlateId','I04/I06'),
 ('R08','Soma主映射','确定','43个单特征注释候选作为主候选池；逐特征报告','这是注释候选，并非已证实分子身份一致的最终Strict集','Stage1映射'),
 ('R09','多适配体','确定','26个多适配体保留逐特征支线；HBA1组合注释保留扩展支线','不默认median(log-z)为已验证统一测量；不依据病理效应选适配体','Stage1映射'),
 ('R10','TMT缺失','确定','源TMT -4置为缺失；只对正、有限样本强度取log2','-4规则仅限有文档依据的TMT字段，不推广到认知/Centiloid','I07/I08/Stage1'),
 ('R11','TMT GIS','确定','主预处理路线为log2(sample/GIS)，样本及对应蛋白GIS均须为正；raw log2为固定敏感性','不根据逐蛋白病理P值选择是否GIS；Stage1 68/70蛋白技术批次R²降低，仅描述技术改善','Stage1'),
 ('R12','TMT版本','外部来源待解释','使用SHA绑定的实际3907定量列版本；方法3447的差异保留来源问题','3904普通UniProt列＋3特殊列；不得编造过滤或版本原因','Stage1版本备忘'),
 ('R13','技术均值调整','候选需验证','保留NULISA/Soma/TMT三技术因子，联合协变量均值模型待仿真验证','三向稳健方差不能消除技术均值偏倚；不直接对每平台逐板残差化后套旧模型','Stage3/Stage4秩审计'),
 ('R14','尺度标准化','确定','模型输入保留上述尺度，作图可使用全样本统一固定仿射缩放并记录参数','不按A/T或诊断组分别z标准化；不将任意缩放当作识别假设','Stage2/Stage3'),
 ('R15','缺失','确定','每蛋白以同一组人同时满足三检测、主病理和协变量的完整病例计数','不插补病理/认知或跨平台缺失，不声称完整病例消除了选择偏倚','本轮覆盖'),
 ('R16','主病理锚点','确定数据路线','连续Elecsys pTau181，优先同一采样日；±90天为敏感性','567人均可同日配对；8/120技术范围，4人无精确值；模型中原值或log2函数形式仍需校准','E01/M02'),
 ('R17','病理批次','候选需验证','保留4个Elecsys批次；同日Batch9-only 528人作为重要敏感性设计','567人总路线需考虑病理检测批次的均值结构，不能假定第四来源没有误差','E01/M02'),
 ('R18','删失与比值','确定','技术限外或< >值保留范围标志；精确连续比值必须分子分母均精确有效','Aβ42限外仍可能判A类别；不能将外推值当精确测量','M02'),
 ('R19','A/T类别','确定用途','Aβ42<980、pTau>27仅作预注册操作性类别与覆盖描述','pTau27来源为临床进展研究；不是tau PET纤维tau阳性；pTau/Aβ42>.025是AD特征比值，不是独立T轴','阈值注册表'),
 ('R20','旧CSF','确定','AlzBio3与Elecsys分开；用已重标度字段，不混raw字段；MEDIAN与原始批次路线分列','MEDIAN代表的发布汇总机制未完全核实，不将其当独立新检测；不跨平台合并连续值','I10/M01'),
 ('R21','淀粉样PET','确定','按SCANDATE，主ROI QC1/2可用；QC2敏感性；跨示踪剂使用Centiloid','±90天335人，±365天338人；PROCESSDATE只作重处理版本排序','I11/M03'),
 ('R22','tau PET','确定用途限制','不同示踪剂分开；连续meta-temporal SUVR仅在各示踪剂内解释','±730天0人；多年后FTP只可作为另行设计的远期结局，不充当同期病理','I12/M04'),
 ('R23','认知版本','确定数据路线','现有2023版UWNPSYCHSUM的ADNI_MEM作为记忆终点候选；其他域作覆盖对照','文件名2026不代表分数更新至2026；ADSP-PHC不在已审计目录，取得后另行登记复审','I09/M05'),
 ('R24','认知时序','确定数据规则','基线±90天，随访>采样90天且限5年；±180天/3年/全部作为覆盖敏感性','至少基线＋1次随访568；基线＋2次且跨度≥2年420；后者不自动成为入组门槛','预先声明随访网格'),
 ('R25','年龄','确定','依据实际MM/YYYY出生信息及年份交叉核对；按月精度而非虚构精确生日','8人跨记录出生年/月冲突保守暂缺；563人年龄可用，均为月/年精度','I13实际格式审计'),
 ('R26','协变量','确定数据路线','基本年龄＋性别；认知扩展年龄＋性别＋教育(±365天)＋APOE ε4','性别、教育、APOE各571可用；两协变量集合联合均563。无实际协变量数值输出','I13/I14'),
 ('R27','诊断','确定用途','诊断类别只作覆盖与队列描述；±90天559、±365天571','DIAGNOSIS=3为dementia，不能直接改名AD；不默认纳入主调整集','I15字典'),
 ('R28','主蛋白门槛','候选需验证','按本轮逐蛋白n、技术设计秩、簇数和检出率进入下一轮真实设计仿真','80%/n≥100等仅描述；不据本轮人数直接冻结最终集合或多重检验家族','Stage1–4'),
 ('R29','2df与认知','确定用途限制','延续相对线性均值约束＋SUM3_CHI2为条件候选，需重验协变量完整模型','这是使用P估计方向的联合测量–病理检验；h1/h2不是个体F/D分数，不能直接放入认知斜率模型','Stage3'),
 ('R30','分析边界','确定','本轮仅字段/范围/类别/日期/覆盖/技术设计，无蛋白–病理或认知关联','独立病理指独立检测表和检测方法，不代表独立队列、独立分装或误差独立','本轮协议'),
 ('R31','病理构念重叠','确定标志','MAPT与pTau同属tau母蛋白测量轴，预先单列，并安排排除MAPT的敏感性','可作为近同构念对照候选；不能把其与pTau的关系视为跨通路独立验证。最终检验家族角色须在SAP预定','Stage1注释＋E01靶点')]
 df=pd.DataFrame(items,columns=['rule_id','topic','status','decision','qualification','evidence'])
 write_csv(RESULTS/'PREPROCESSING_DECISIONS.csv',df)
 save_json(RESULTS/'PREPROCESSING_RULES_v1.0.json',{'version':'1.0','date':'2026-10-27','scope':'Data rules and design candidates; not a final association SAP',
    'raw_protein_clinical_analysis_executed':False,'rules':df.to_dict('records')})
 actions=[
 ('N01','首要','完整均值模型','用真实缺失、三技术因子及病理批次结构模拟年龄/性别等协变量调整；验证相对线性均值约束及其影响函数','未完成；三簇方差不能代替均值去混杂'),
 ('N02','首要','真实设计功效/多重检验','逐蛋白n/19–36–65簇结构/秩亏风险进入仿真，预定科学可检测差异和多重检验家族','未完成；n≥100或80%覆盖不是充分资格'),
 ('N03','首要','测量识别边界','把残差误差相关、低LOD、非线性P函数和测量批次误差纳入敏感性；保留无法区分真实丰度/测量偏移的结论','Stage3结论仍适用，不能凭新覆盖解除'),
 ('N04','待核查','出生信息8人','用restricted标志中的RID回到原始记录/ADNI核对出生年/月冲突；无依据前暂缺，不人为选对分析有利记录','已定位数量，尚无权威真值'),
 ('N05','来源澄清','TMT 3447/3907','获取与该矩阵对应的方法/发布版本说明；现有文件哈希绑定并披露差异','外部来源未解释，不阻止本轮审计交付'),
 ('N06','认知方案','认知表示和随访模型','先定义可被识别的个体表示及交叉拟合，再验证混合模型与失访敏感性；避免直接套用h1/h2','尚无经验证的个体F/D分数'),
 ('N07','可选扩展','ADSP-PHC与远期PET','需要时下载登记PHC新版；tau PET作为多年后独立前瞻问题另设选择/时距模型','现有5年ADNI_MEM可继续设计，不必等待外部下载'),
 ('N08','SAP冻结前','正式口径与图表','依据校准/功效固定最终蛋白集、P函数形式、病理批次处理、敏感性和报告用语','不得只因有567人就宣布模型有效')]
 write_csv(RESULTS/'NEXT_ACTIONS.csv',pd.DataFrame(actions,columns=['id','priority','topic','action','state']))
 sources=[
 ('E01','Elecsys CSF','EXAMDATE','pTau181 / Aβ42 / Aβ40 / tTau','范围与缺失检查；保留BATCH','连续pTau主锚点；比值次级','M02'),
 ('I10','AlzBio3 CSF','DRAWDTE','重标度ABETA/TAU/PTAU','MEDIAN与原始检测批次分开','旧检测敏感性，不能混Elecsys','M01'),
 ('I11','淀粉样PET','SCANDATE','Centiloid / 发布A类别','qc_flag=1/2主ROI，=2敏感性','同期支持；不直接合并SUVR','M03'),
 ('I12','tau PET','SCANDATE','META_TEMPORAL_SUVR','qc_flag=1/2；示踪剂分开','同期覆盖不足；远期另立问题','M04'),
 ('I09','认知','EXAMDATE','ADNI_MEM/EF/EF2/LAN/VS','有限分值视有效；负值不自动缺失；唯一访视日','5年记忆随访候选','M05'),
 ('I13/I14','协变量','VISDATE/出生年月','年龄/性别/教育/APOE ε4','跨记录一致性、教育±365天、APUSABLE显式不通过排除','覆盖标志，不输出个体实际值','字典'),
 ('I15','诊断','EXAMDATE','CN / MCI / DEMENTIA','排除显式QC错误；同日冲突列明','描述，不默认主调整','字典')]
 write_csv(RESULTS/'CLINICAL_SOURCE_DECISIONS.csv',pd.DataFrame(sources,columns=['input_id','source','time_anchor','fields','validity','role','local_basis']))
 refs=[
 ('REF01','ADNI PET documentation','官方','https://adni.loni.usc.edu/quick-start-guide-asset/pet.html','SCANDATE/QC；Centiloid跨示踪剂；tau需分示踪剂'),
 ('REF02','ADNI UPENN CSF Biomarker Master reference','官方','https://adni.bitbucket.io/reference/upennbiomk_master.html','重标度与raw字段区别；DRAWDTE与RUNDATE'),
 ('REF03','Blennow et al. Scientific Reports 9:19024 (2019)','原始研究','https://www.nature.com/articles/s41598-019-54204-z','Elecsys pTau27用于临床进展风险分层；并非tau PET阈值'),
 ('REF04','Shaw et al. Ann Neurol (2009), CSF biomarker signature in ADNI','原始研究','https://pmc.ncbi.nlm.nih.gov/articles/PMC2696350/','AlzBio3原始平台阈值背景，不迁移到Elecsys'),
 ('REF05','ADNI Data Dictionary Search','官方','https://adni.loni.usc.edu/data-samples/data-dictionary-search/','字段编码与诊断含义；本地字典固定版本用于执行'),
 ('M01','UPENNBIOMK_MASTER_METHODS_20231109.txt','本地原始方法','results/M01_METHODS_EXTRACTED.txt','AlzBio3重标度和批次建议'),
 ('M02','UPENNBIOMK_ROCHE_ELECSYS_METHODS_20231109.txt','本地原始方法','results/M02_METHODS_EXTRACTED.txt','Aβ42 200–1700、tTau80–1300、pTau8–120；Aβ42<980、Aβ42/40<.061'),
 ('M03','ADNI_UCBerkeley_AmyloidPET_Methods_v3_2025-06-30.pdf','本地原始方法','results/M03_METHODS_EXTRACTED.txt','QC、Centiloid、发布A状态与18CL口径'),
 ('M04','ADNI_UCBerkeley_TauPET_Methods_2025-06-30.pdf','本地原始方法','results/M04_METHODS_EXTRACTED.txt','6mm、meta-temporal ROI、示踪剂、归一化已完成'),
 ('M05','ADNI_Methods_UWNPSYCHSUM_20231018.pdf (revised Nov2023)','本地原始方法','results/M05_METHODS_EXTRACTED.txt','最终2023更新；推荐ADSP-PHC；旧版分值有重算')]
 write_csv(RESULTS/'EVIDENCE_REGISTRY.csv',pd.DataFrame(refs,columns=['id','title','type','location','supports']))
 metrics=[('GIS_NUMERIC','三平台有效数值，Soma需候选内全部适配体有效，TMT需GIS有效；单特征43人为单适配体'),
 ('GIS_LOD','GIS_NUMERIC另要求NULISA所有重复都达对应逐板LOD，仅敏感性人数'),
 ('RAW_NUMERIC','TMT不作GIS分母约束；仍按正样本强度及Soma全部有效计数'),
 ('GIS_UNAMBIGUOUS','另要求Soma无组合accession注释；不是已证实分子身份'),
 ('n_same_day_ptau_age_sex','本蛋白三检测完整＋同日精确Elecsys pTau＋年龄性别有效'),
 ('n_same_day_ptau_fullcov','上列另含教育及APOE；本数据与基本协变量人数相同'),
 ('n_p_fullcov_cog_minimum','上列另含±90天记忆基线及5年内采样90天后至少1次记忆随访'),
 ('n_p_fullcov_cog_highinfo','上列认知条件加强为至少2次随访且基线至末访≥2年；不是正式纳入规则'),
 ('technical_design_rank','仅截距＋3个技术因子的加性虚拟变量矩阵秩；未含实际年龄、病理数值及其函数'),
 ('n_minus_technical_rank','n减上述技术秩；仅说明结构余量，不是最终模型自由度或有效样本量'),
 ('effective_clusters','(各簇人数和)^2/各簇人数平方和；衡量簇大小不均衡，不代替校准'),
 ('operational_AT','Elecsys Aβ42<980、pTau>27，同记录配对的操作性标签；不得表述为tau PET分期')]
 write_csv(RESULTS/'METRIC_DICTIONARY_ZH.csv',pd.DataFrame(metrics,columns=['metric','meaning']))
 meta=pd.read_csv(STAGE1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv',usecols=['accession','nulisa_target','annotation_strict_candidate'])
 meta['overlaps_primary_ptau_parent_protein']=meta.nulisa_target.eq('MAPT')
 meta['interpretation']=np.where(meta.overlaps_primary_ptau_parent_protein,'Tau parent-protein construct overlap; separate reporting and exclude-MAPT sensitivity','No identical parent-protein label flagged; not proof of error independence')
 write_csv(RESULTS/'CLINICAL_CONSTRUCT_OVERLAP_FLAGS.csv',meta)

def figures():
 plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':10,'pdf.fonttype':42})
 def save(fig,name):
  fig.savefig(FIGURES/(name+'.png'),dpi=220,bbox_inches='tight');fig.savefig(FIGURES/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
 df=pd.read_csv(RESULTS/'PATHOLOGY_WINDOW_COVERAGE.csv')
 fig,ax=plt.subplots(figsize=(9,4.8))
 for route,label,color in [('ELECSYS_PTAU','Elecsys pTau精确值','#264f87'),('ELECSYS_A4240_RATIO','Elecsys Aβ42/40精确比值','#977023'),('AMY_PET_CENTILOIDS','淀粉样PET Centiloid','#138a80'),('TAU_PET_ANY_TRACER','tau PET（任意示踪剂）','#b85050')]:
  part=df.loc[df.route.eq(route)];ax.plot(part.window_days,part.n_available,'o-',label=label,color=color)
 ax.set(xlabel='采样日前后允许窗口（天）',ylabel='可匹配人数（总索引571）',ylim=(-20,600),xticks=WINDOWS,title='同期覆盖：放宽窗口不能补足tau PET')
 ax.legend(loc='center right');ax.grid(alpha=.18);save(fig,'01_PATHOLOGY_WINDOWS')
 fig,ax=plt.subplots(figsize=(9,4.4))
 labels=['精确同日三平台','同日Elecsys pTau','再加年龄/性别','再加最低认知随访','若要求高信息认知随访']
 joins=pd.read_csv(RESULTS/'JOINT_COHORT_COVERAGE.csv').set_index('flag').n_participants
 vals=[571,567,int(joins.PTAU90_MINCOV),int(joins.PTAU90_MINCOV_COGMIN),int(joins.PTAU90_FULLCOV_COGHIGH)]
 bars=ax.barh(labels[::-1],vals[::-1],color=['#888fa2','#138a80','#286090','#3d78a2','#9aafc5'])
 for bar,n in zip(bars,vals[::-1]):ax.text(n+4,bar.get_y()+bar.get_height()/2,str(n),va='center')
 ax.set(xlim=(0,620),xlabel='人数；未叠加某一蛋白的缺失',title='临床与协变量交集：逐蛋白样本量仍需另算');save(fig,'02_COHORT_INTERSECTIONS')
 grid=pd.read_csv(RESULTS/'COGNITION_FOLLOWUP_GRID.csv').query('baseline_window_days==90 and horizon_years==5').pivot(index='minimum_total_visits',columns='minimum_span_years',values='n_participants')
 fig,ax=plt.subplots(figsize=(6.4,4.5));im=ax.imshow(grid.values,cmap='Blues',vmin=0,vmax=571)
 for i in range(3):
  for j in range(3):ax.text(j,i,str(grid.iloc[i,j]),ha='center',va='center',color='white' if grid.iloc[i,j]>400 else '#123')
 ax.set(xticks=range(3),xticklabels=grid.columns,yticks=range(3),yticklabels=grid.index,xlabel='要求基线至末访跨度（年）',ylabel='要求有效总访视次数',title='ADNI_MEM：±90天基线，5年随访覆盖网格')
 fig.colorbar(im,ax=ax,label='人数 / 571');save(fig,'03_COGNITIVE_GRID')
 df=pd.read_csv(RESULTS/'PER_PROTEIN_JOINT_COVERAGE.csv');df=df.loc[df.annotation_single_feature&df.route.eq('GIS_NUMERIC')].sort_values('n_same_day_ptau_age_sex')
 fig,ax=plt.subplots(figsize=(10,11.5));y=np.arange(len(df));ax.barh(y,df.n_three_assays,color='#d3dfec',label='三检测数值完整')
 ax.barh(y,df.n_same_day_ptau_age_sex,color='#326898',label='＋同日pTau、年龄/性别')
 ax.plot(df.n_p_fullcov_cog_highinfo,y,'o',ms=3.5,color='#b66d30',label='＋高信息认知随访（候选）')
 ax.set(yticks=y,yticklabels=df.gene,xlim=(0,595),xlabel='逐蛋白人数',title='43个单特征注释候选：真实联合覆盖（非最终蛋白集）');ax.legend(loc='lower right',fontsize=8);ax.grid(axis='x',alpha=.15)
 save(fig,'04_PER_PROTEIN_JOINT_COVERAGE')

def reports():
 qa=json.loads((RESULTS/'INDEPENDENT_QA_SUMMARY.json').read_text(encoding='utf-8'))
 p=pd.read_csv(RESULTS/'PER_PROTEIN_JOINT_COVERAGE.csv');main=p.loc[p.annotation_single_feature&p.route.eq('GIS_NUMERIC')]
 joins=pd.read_csv(RESULTS/'JOINT_COHORT_COVERAGE.csv').set_index('flag').n_participants
 thresholds=pd.read_csv(RESULTS/'THRESHOLD_REGISTRY.csv')
 def fmt(a):return str(int(a)) if float(a).is_integer() else str(a)
 def triple(col):return '–'.join(fmt(v) for v in [main[col].min(),main[col].median(),main[col].max()])
 text=f'''# Stage4 真实预处理规则决策与独立病理、认知覆盖审计

完成日期：2026-10-27。审计结论：已有数据足以继续设计“同日独立CSF pTau相关的跨检测相对测量模式”，并有淀粉样PET和记忆随访作为支持；目前不能把研究主张扩展为同期tau PET验证或已识别的真实丰度/测量偏差分解。

**最关键的决定。** 蛋白预处理采用NULISA已发布NPQ、Soma正RFU的log2、TMT正样本强度/GIS的log2；raw-TMT和逐人LOD删除保留敏感性支线。主候选池仍为43个单特征注释候选，26个多适配体与HBA1组合注释为扩展支线。主病理数据路线选择同日Elecsys pTau181；原值或log2的模型函数形式、病理批次调整、最终蛋白集和正式检验家族仍须真实设计仿真后冻结。

**完成范围与数据边界。** 原始15个输入加新增Elecsys及5份原始方法共21项来源散列均核对。先写v1.0协议，再有限读取临床字段；新增来源先登记再读取。出生信息格式和tau PET远时距扩展有单列v1.1补充，未改写原协议。未修改Stage1–3、原始数据、原SAP或旧论文。连续临床值仅临时转换为有效性、范围和预先文档化类别，未输出个体分数或分布；蛋白个体定量未与临床连接。此次关联模型数量为0。

“独立病理”是另外发布的临床检测表/检测方法，不是独立验证队列，也不证明同一分装、误差独立或可识别真实分子浓度。后续相对均值约束本身使用病理P估计残差共同方向，仍应称联合测量–病理分析。

| 571人中的来源/口径 | 同日 | ±90天 | ±365天 | 决定 |
|---|---:|---:|---:|---|
| Elecsys pTau181精确值 | 567 | 567 | 567 | 连续主病理锚点；扩大窗口无收益 |
| Elecsys Aβ42精确值 | 481 | 481 | 484 | 范围外不能当精确值 |
| Elecsys Aβ42操作性A类别 | 571 | 571 | 571 | 范围信息足够时仍可判类别 |
| Elecsys Aβ42/40精确比值 | 30 | 30 | 32 | 无法承担全队列连续A主分析 |
| Elecsys pTau/Aβ42精确比值 | 477 | 477 | 480 | AD特征比值，不是独立T轴 |
| AlzBio3 pTau（MEDIAN发布路线） | 571 | 571 | 571 | 与Elecsys分开，不能直接拼接 |
| 淀粉样PET Centiloid，QC1/2 | 63 | 335 | 338 | ±90天支持；±365天仅增3人 |
| tau PET，QC1/2任意示踪剂 | 0 | 0 | 0 | ±730天仍为0，不能作同期验证 |
| 诊断可解释类别 | 92 | 559 | 571 | dementia不自动等于AD |

以上均为输入可用性计数，非各蛋白最终分析人数。PET按采集日期SCANDATE匹配，处理日期只用于去除重处理重复。淀粉样PET跨示踪剂使用Centiloid；tau PET不同示踪剂不直接合并SUVR。发布的A状态与18CL是两个操作性定义，差异只描述为口径差异，见定义对照表。[ADNI官方PET说明](https://adni.loni.usc.edu/quick-start-guide-asset/pet.html)

**Elecsys具体口径。** pTau同日567人来自UPENNBIOMK9 528人、11批26人、13批8人、12批5人；Batch9-only应作为重要敏感性设计。原方法给出的Aβ42、tTau、pTau技术范围分别为200–1700、80–1300、8–120 pg/mL。全源Aβ42有511条超过技术上限的外推数字，本轮将其转为范围标志，不保留为精确连续值；这一511不是571人索引中的人数。逐索引覆盖必须见上表。比值必须分子分母都精确有效，不能用范围外Aβ42制造精确比值。相关依据见本地M02提取文本与源文件散列清单。

操作性Aβ42<980及pTau>27在同一Elecsys记录形成A/T标签，567人中A−T−189、A−T+47、A+T−134、A+T+197。27 pg/mL来自临床进展研究，不是统一的纤维tau PET阳性标准；本轮主指标是连续pTau，因此无需二分才能开展模型。[Blennow等原始研究](https://www.nature.com/articles/s41598-019-54204-z) 支持此阈值的用途限定。AlzBio3的192/23阈值单独登记，不能迁移至Elecsys；其重标度值和raw值是不同字段。[ADNI字段说明](https://adni.bitbucket.io/reference/upennbiomk_master.html)

**tau PET的真实限制。** 索引采样日期为2005-09-12至2013-01-18，而本地tau PET始于2015年。±3年仍0人，5年内27人、10年内109人、任意时间119人，均有FTP；最早有效扫描距索引1172天，任意时间最近扫描间隔中位数2176天。这些是多年后的潜在前瞻结局，受存活、回访和扫描选择影响，不能回填成采样时tau状态。MK6240无索引人匹配，PI2620仅1人且远隔5127天。

**认知与协变量。** ADNI_MEM基线±90天568人，均在采样90天后5年内至少有1次有效随访；若另要求至少2次随访且基线至末访≥2年，为420人。EF对应高信息量417人，EF2/LAN/VS均420人。次数/跨度组合网格全部交付，不能根据后续蛋白结果选择门槛。使用全部有基线与随访者进入适当纵向模型，可减少仅纳入长期回访者造成的选择；最终时间函数、随机效应与失访处理另验。

认知文件虽以2026命名，原方法说明ADNI_MEM等最终更新基于2023-10-18数据拉取，实际有效访视日期最晚2023-10-03。原方法建议今后使用ADSP-PHC；在已审计的Neuropsychological目录未找到该文件，未擅自切换终点。现有索引最晚2013年，5年随访窗口仍可用此固定版本设计；若取得PHC，须登记、复审版本和覆盖后再决定角色。认知及Centiloid的负数不能套用TMT的−4缺失规则。

性别、教育及APOE各571人可用；出生信息跨记录冲突8人保守暂缺，年龄可用563人，均为月/年精度。实际PTDOB为MM/YYYY、PTDOBYY为ISO日期字符串；只取后者年份作交叉核对，不虚构生日。教育取±365天内最近有效记录，APOE显式不合格排除、旧版空QC不自动判失败。年龄＋性别与加教育＋APOE两种集合联合可用数均563。同日pTau与任一协变量集交集560，再加最低认知随访557，再加高信息认知条件410。8人需回源核对，未经依据不能擅自挑选出生记录。

**逐蛋白人数和技术设计。** 43个单特征注释候选的最小–中位数–最大人数如下；所有门槛均为结构描述，不是正式筛选标准。

| 口径 | 最小–中位数–最大 |
|---|---|
| 三检测数值完整＋GIS | {triple('n_three_assays')} |
| ＋同日精确pTau | {triple('n_same_day_ptau')} |
| ＋年龄/性别（或完整协变量） | {triple('n_same_day_ptau_age_sex')} |
| ＋最低认知随访 | {triple('n_p_fullcov_cog_minimum')} |
| ＋高信息认知条件 | {triple('n_p_fullcov_cog_highinfo')} |
| 淀粉样PET±90天与本蛋白三检测交集 | {triple('n_amy_pet_90')} |

43候选中40个主路线人数≥100，余者S100A12 51、MME60、TEK91。仅截距＋三平台技术因子在571人设计中秩118，技术结构余量453；在560人pTau＋协变量集合中仍秩118，余量442。这里未放入实际年龄、性别数值、P或其函数，不是最终回归自由度。S100A12、MME、TEK逐蛋白技术设计余量仅1、8、30；不能不顾稀疏性直接套用大样本模型。

NULISA/Soma/TMT实际簇数为19/36/65，且逐蛋白缺失会改变簇数和不均衡程度。GIS不意味着批次偏倚消失。Stage3已表明原始技术混杂可在没有额外病理直接效应时产生接近100%的拒绝；稳健方差只处理不确定度，不能校正均值偏倚。因此本轮人数改善并未解除完整均值模型与真实设计再校准要求。

**识别和论文表述。** 当前2df相对线性均值约束检验仍为条件候选，不能宣称自动分离真实丰度和病理相关检测偏差；任意残差误差相关下的不可识别性仍存在。h1/h2是逐蛋白群体约束，不是个体F/D。认知延伸必须另行定义并验证个体表示、交叉拟合和纵向模型。MAPT与pTau同属tau母蛋白测量轴，已预先标记：单列结果并安排排除MAPT的敏感性，不把这种近同构念关系称作跨通路独立病理验证。论文主张应保留为可检验的跨检测相对模式及其与独立临床测量的关系，待识别假设与敏感性证据支持后再细化机制措辞。

**交付和复核。** 所有数据层决定见[31项规则](results/PREPROCESSING_DECISIONS.csv)及[机器可读规则](results/PREPROCESSING_RULES_v1.0.json)；数据来源、阈值、窗口、随访网格、逐蛋白人数、批次秩和后续动作均在[Excel总表](真实预处理与病理认知覆盖审计总表.xlsx)。[独立复核](results/INDEPENDENT_QA_SUMMARY.json)为{qa['passed']}/{qa['checks']}项通过、0项失败，重算E01/I09/I11/I12/I13的关键标志和全部逐蛋白交集，另外核对来源散列与禁入列。计数检验通过不等于关联模型已经校准。

个体RID、访视日期、类别和可用性仅放在00_restricted，不含原始连续临床值或基因型；该目录名称不代表另设操作系统访问权限。全部普通交付为汇总。详见[来源与依据](results/EVIDENCE_REGISTRY.csv)、[访问字段](results/FIELD_ACCESS_LOG.json)、[执行日志](logs/AUDIT_RUN.log)及最终SHA清单。

**下一步顺序。** 首先在真实缺失和批次结构下验证协变量完整均值模型、P的函数形式和Elecsys批次处理，再做逐蛋白功效/多重检验与最终蛋白集冻结。随后进入受限定的正式病理分析；认知表示与纵向延伸需要单独验证。8人出生信息与TMT3447/3907差异作为来源问题并行核查；无需为了增加人数放宽tau PET至多年或把Aβ42外推数字当精确值。
'''
 (ROOT/'STAGE4_AUDIT_REPORT_ZH.md').write_text(text,encoding='utf-8')
 decision='''# 真实预处理决策备忘

本轮已决定数据层转换与来源口径：NULISA NPQ直接保留，Soma正RFU取log2，TMT -4置缺失且正样本/GIS取log2；raw-TMT与LOD删人保留固定敏感性。技术标签全部保留。不能以病理关联或认知关联结果选择转换、适配体或分析窗口。

主病理数据路线选择同日Elecsys pTau181（567人），不扩大窗口换人数；Batch9-only（528人）作为重要设计敏感性。Aβ42的精确值与范围推导类别分开，Aβ42/40仅30人；AlzBio3单独分析，不能直接拼接。PET按采集日期，淀粉样跨示踪剂用Centiloid，tau示踪剂分开；同期tau PET不可用。

认知候选是当前固定2023版本ADNI_MEM，基线±90天，采样90天后至5年的随访；保留最低随访与高信息量两套人数。对认知及Centiloid不套用TMT的−4规则。年龄按出生月/年，8名出生信息冲突者暂缺；其余协变量覆盖完整。

仍未冻结的是最终蛋白资格门槛、P原值或log2函数、技术与病理批次的完整均值调整、相对均值约束影响函数、真实设计功效/多重检验和个体认知表示。它们影响统计有效性，不能由此次覆盖计数代替。三向SUM3_CHI2仅是Stage3的条件候选。逐平台预先板内残差化会改变共同方向，不直接套旧检验。

43个单特征注释候选不是最终Strict；26个多适配体不默认聚合，HBA1留扩展。低检出MME/YWHAZ/S100A12/BASP1需要测量质量与LOD敏感性验证；保留低LOD数值仅避免预先截断，不代表其是可靠定量。覆盖≥80%和n≥100只作描述，不能代替功效/秩/簇结构考核。

TMT3447/3907来源差异、8名出生信息和可选PHC新版获取均已列入下一步，不据猜测填补。未发出外部邮件。详见规则CSV的逐项证据与状态，及主报告中的人数、图表和识别限制。
'''
 (ROOT/'PREPROCESSING_DECISION_MEMO_ZH.md').write_text(decision,encoding='utf-8')
 readme='''# Stage4审计交付

先读[中文审计报告](STAGE4_AUDIT_REPORT_ZH.md)和[预处理决策备忘](PREPROCESSING_DECISION_MEMO_ZH.md)，逐表查阅[Excel总表](真实预处理与病理认知覆盖审计总表.xlsx)。

已完成真实规则、独立病理和认知覆盖、协变量冲突、逐蛋白布尔交集及技术设计秩审计。尚未进行蛋白–病理/认知关联，也未冻结完整模型与功效门槛。

results为汇总数据和QA；figures提供4组PNG/PDF；code含复现与独立验证；logs记录执行及修正；00_restricted仅含标识/日期/类别/可用性，无原始连续临床值。不要把restricted当公开附件。

复现：在本机PowerShell执行code/run_all.ps1。要求原输入、Stage1封存注册表/布尔掩码、Python及SOFTWARE_ENVIRONMENT.json中的依赖可用。运行会覆盖本Stage4生成物，不改写此前阶段和原始数据。根目录最终SHA清单不对自身取散列，排除Python缓存。
'''
 (ROOT/'README_审计总览.md').write_text(readme,encoding='utf-8')

def workbook():
 tables=[('规则决策','PREPROCESSING_DECISIONS'),('待处理事项','NEXT_ACTIONS'),('来源口径','CLINICAL_SOURCE_DECISIONS'),('阈值注册','THRESHOLD_REGISTRY'),
 ('病理窗口','PATHOLOGY_WINDOW_COVERAGE'),('临床协变量交集','JOINT_COHORT_COVERAGE'),('逐蛋白联合覆盖','PER_PROTEIN_JOINT_COVERAGE'),('43候选主路线',None),
 ('逐蛋白覆盖汇总','PER_PROTEIN_COVERAGE_SUMMARY'),('逐蛋白技术秩簇数','PER_PROTEIN_TECHNICAL_DESIGN'),('总体技术设计','COHORT_TECHNICAL_DESIGN'),
 ('认知随访网格','COGNITION_FOLLOWUP_GRID'),('认知域对照','COGNITIVE_DOMAIN_COVERAGE'),('认知日期结构','COGNITIVE_TIMING_STRUCTURE'),
 ('协变量覆盖','COVARIATE_AVAILABILITY'),('协变量冲突','COVARIATE_CONFLICT_COUNTS'),('协变量类别','COVARIATE_CATEGORY_COUNTS'),
 ('操作性类别计数','OPERATIONAL_CATEGORY_COUNTS'),('队列操作性AT','COHORT_OPERATIONAL_AT_COUNTS'),('逐蛋白操作性AT','PER_PROTEIN_OPERATIONAL_AT_COUNTS'),
 ('tau远期时距','TAU_PET_EXTENDED_TIMING'),('pTau检测批次','SELECTED_PTAU_ASSAY_BATCH_COUNTS'),('技术批次类别覆盖','TECHNICAL_BATCH_CATEGORY_OCCUPANCY'),
 ('淀粉样定义对照','AMY_RELEASE_VS_CL18_DEFINITIONS'),('源记录结构','SOURCE_ROW_STRUCTURE'),('源日期范围','SOURCE_TIMING_STRUCTURE'),
 ('重复冲突','DUPLICATE_AND_CONFLICT_AUDIT'),('数值状态非分布','CLINICAL_VALUE_STATUS_COUNTS'),('元数据类别','METADATA_CATEGORY_COUNTS'),
 ('出生字段格式','DEMOGRAPHIC_FORMAT_AUDIT'),('新增输入散列','ADDITIONAL_INPUT_MANIFEST'),('原输入散列检查',None),('指标字典','METRIC_DICTIONARY_ZH'),
 ('依据注册','EVIDENCE_REGISTRY'),('病理构念重叠','CLINICAL_CONSTRUCT_OVERLAP_FLAGS'),('独立QA','INDEPENDENT_QA_CHECKS')]
 path=ROOT/'真实预处理与病理认知覆盖审计总表.xlsx';counts=[]
 with pd.ExcelWriter(path,engine='openpyxl') as writer:
  for sheet,filename in tables:
   if sheet=='43候选主路线':
    frame=pd.read_csv(RESULTS/'PER_PROTEIN_JOINT_COVERAGE.csv');frame=frame.loc[frame.annotation_single_feature&frame.route.eq('GIS_NUMERIC')]
   elif sheet=='原输入散列检查':frame=pd.DataFrame(json.loads((RESULTS/'INPUT_HASH_CHECKS.json').read_text(encoding='utf-8')))
   else:frame=pd.read_csv(RESULTS/(filename+'.csv'))
   frame.to_excel(writer,sheet_name=sheet,index=False);counts.append({'sheet':sheet,'rows':len(frame),'columns':len(frame.columns)})
 book=load_workbook(path)
 for ws in book:
  ws.freeze_panes='C2';ws.auto_filter.ref=ws.dimensions
  for cell in ws[1]:cell.font=Font(bold=True,color='FFFFFF');cell.fill=PatternFill('solid',fgColor='2B557E');cell.alignment=Alignment(wrap_text=True,vertical='top')
  ws.row_dimensions[1].height=45
  for col in ws.columns:
   letter=col[0].column_letter;width=min(52,max(13,max(len(str(c.value or '')) for c in list(col)[:60])+2));ws.column_dimensions[letter].width=width
  for row in ws.iter_rows(min_row=2):
   for cell in row:cell.alignment=Alignment(vertical='top',wrap_text=ws.title in ['规则决策','待处理事项','来源口径','指标字典','依据注册'])
  if ws.title=='43候选主路线':
   for col in range(6,ws.max_column):ws.conditional_formatting.add(f'{ws.cell(2,col).column_letter}2:{ws.cell(2,col).column_letter}{ws.max_row}',DataBarRule(start_type='num',start_value=0,end_type='num',end_value=571,color='96BCD9'))
 book.save(path);write_csv(RESULTS/'WORKBOOK_SHEET_COUNTS.csv',counts)

if __name__=='__main__':
 rules();figures();reports();workbook();log('Chinese reports, rule registry, four PNG/PDF figure pairs, and workbook produced')
