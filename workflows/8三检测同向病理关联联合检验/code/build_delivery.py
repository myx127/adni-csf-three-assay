from config import *
import argparse,importlib.metadata
from openpyxl.styles import Font,PatternFill,Alignment
from openpyxl.utils import get_column_letter
PURPOSE={
 'PRIMARY_FAMILY34_SHARED_RESULTS':'主结果：原34蛋白、M0、双方向共同关联、原H及等效证据',
 'ALL35_THREE_MODEL_SHARED_RESULTS':'35蛋白×三模型全部参数、联合及不同家族校正',
 'PLATFORM_EFFECTS_ALL_MODELS':'三检测病理系数、SUM3标准误及含尺度不确定性的标准化区间',
 'MODEL_DISCOVERY_COUNTS':'各模型名义/BY/方向68/合并68发现计数',
 'EVIDENCE_STATE_COUNTS':'共同关联与原H偏离的四种允许共存状态',
 'APOE_SENSITIVITY_COMPARISONS':'同病例模型比较及自身P/其他33项P对q的代数贡献；非因果',
 'MAPT_CONTROL_ONLY':'同母蛋白构念对照；不进入34家族',
 'CALIBRATION_DECISIONS':'105设计逐项严格有限网格资格门槛',
 'CALIBRATION_MODEL_SUMMARY':'三模型校准概况',
 'COMPONENT_NULL_CALIBRATION':'各检测零系数单侧.025尾部检查',
 'NULL_CALIBRATION_SUMMARY':'25零假设组合×五DGP×三阈值全量评价；噪声复用',
 'POWER_SCENARIO_SUMMARY':'28预设同向效应组合的全量条件功效；非观察效应功效',
 'POWER_WORST_DGP_SUMMARY':'逐模板/效应的五DGP最低成功率和最低Wilson下界',
 'BALANCED_POWER_FAMILY_SUMMARY':'平衡效应的34模板功效概况',
 'SIMULATION_FILE_INVENTORY':'525个独立噪声文件及精确种子/SHA',
 'EFFECT_PATTERN_GRID':'25零假设及28功效系数组合，系数单位为模板残差SD/病理SD',
 'TEMPLATE_DESIGN_AUDIT':'105实际设计、人数、秩、簇及观测残差矩阵正定性',
 'ENGINE_INDEPENDENT_CHECKS':'1651项线性算子/原OLS方差/位置等变性/联合公式核查',
 'INDEPENDENT_QA_CHECKS':'独立科学计算核查',
 'INDEPENDENT_SIMULATION_RECOUNTS':'525文件全部53组合×三阈值的独立重计数',
 'PRIOR_PARAMETER_REPRODUCTION':'105设计与原系数/SE及H值复现',
 'REFERENCE_CODE_ORIGIN':'三份参考代码原封不动复制来源',
 'SOURCE_MANIFEST':'读取的原始缓存/聚合来源散列',
 'PRIOR_SEALED_FILE_CHECKS':'4323条前期封存核查',
 'PRIOR_SEALED_FILE_CHECKS_FINAL':'最终独立重查的前期封存条目',
 'FIGURE_INDEX':'7组科学图PNG/PDF及中文图注'}
CSV_ONLY={'NULL_CALIBRATION_SUMMARY.csv','POWER_SCENARIO_SUMMARY.csv'}
EXCLUDE={'WORKBOOK_SHEET_INDEX.csv','DELIVERY_QA_CHECKS.csv'}
def table(path):return pd.read_csv(path,dtype={'seed':str})
def read(name):return table(RESULTS/(name+'.csv'))
def put(name,body):(ROOT/name).write_text(body.strip()+'\n',encoding='utf-8')
def md(df):
 def val(x):
  if pd.isna(x):return '—'
  if isinstance(x,(float,np.floating)):return str(int(x)) if np.isfinite(x) and x==int(x) else f'{x:.5g}'
  return str(x).replace('|','/').replace('\n',' ')
 return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(val(x) for x in row)+' |' for row in df.itertuples(index=False,name=None))

def docs():
 allr=read('ALL35_THREE_MODEL_SHARED_RESULTS');main=read('PRIMARY_FAMILY34_SHARED_RESULTS');counts=read('MODEL_DISCOVERY_COUNTS')
 sets={m:set(allr.loc[allr.model.eq(m)&allr.shared_BY34_supported,'gene']) for m in MODELS};stable=set.intersection(*sets.values());assert len(stable)==14 and sets['M0']==stable
 genes=', '.join(sorted(stable));shared=main.loc[main.shared_BY34_supported]
 candidate=allr.loc[allr.accession.isin(CANDIDATES),['gene','model','p_shared_computed','q_BY_shared34','q_BY_combined68_shared','prior_H_q_BY_by_model','weakest_assay']]
 candidate=candidate.sort_values(['gene','model']);candidate.columns=['蛋白','模型','共同关联P','共同关联BY34q','合并68共同关联q','原H的BYq','最弱平台']
 power=read('BALANCED_POWER_FAMILY_SUMMARY');power=power.loc[power.model.eq('M0')&power.sign.eq(1)&power.threshold.isin(['NOMINAL','BY34_SUFFICIENT']),['magnitude','threshold','median_worst_rate','n_Wilson_lower_atleast80']]
 power.columns=['三个真实效应幅度','阈值','模板最差DGP率的中位数','各DGP下界均达80%的模板数']
 gates=read('CALIBRATION_MODEL_SUMMARY');qa=json.loads((RESULTS/'INDEPENDENT_QA_SUMMARY.json').read_text(encoding='utf-8'))
 oldineligible=shared.loc[~shared.prior_M0_equivalence_eligible,'gene'].tolist()
 scope=json.loads((RESULTS/'CALIBRATION_SCOPE.json').read_text(encoding='utf-8'))
 report=f"""
# 第八阶段：三检测同向病理关联联合检验

交付日期2026-10-28。已知阶段1–7结果后的新增方法分支；不是整体研究前瞻注册，不是独立队列复现。新方案与模拟引擎在本轮完整模拟及新增真实联合P计算前冻结。原34家族、原H发现及第七阶段等效结果保持原记录。

## 最重要的研究结果

1. M0基础模型34个蛋白均通过本轮评价门槛；18个名义共同关联，14个通过BY34。14个全部为正向三检测关联。APOE及教育/APOE模型仍是完全相同的14项，均通过各自BY34和方向性68项BY。
2. 14项为：{genes}。这是“三检测同向病理关联”的正面证据，不等于灵敏度校正后的响应比例等效，更不等于真实完整蛋白浓度已分离。
3. NPTX1和KLK6同时具有M0共同病理关联与原相对响应偏离证据。加入APOE后，两者共同关联仍保留，而原H偏离的BY证据不再达到.05；说明两种统计属性具有不同的调整敏感性。
4. PDGFRB的M0共同关联P=.005871、BY34 q=.054806，未达到主标准。APOE/BOTH的q约.07164/.08442。它仍是原M0相对偏离发现，但不能列入本轮BY共同关联清单。
5. M0合并34共同关联与原34H检验共68项校正，14项共同关联仍全部保留；APOE/BOTH下保留13项，NEFL不再达到这套更广家族标准。不同家族q必须标明，不能选择较有利者替代主标准。

{md(counts)}

上表n_original_H是各条件模型的既有H结果；只有M0的三项属于阶段5冻结原主发现。其他两模型的H结果是既有敏感性，不能追称独立复制。

## 新分支究竟检验什么

每个蛋白沿用同一组原三检测完整病例，独立病理为同日Elecsys log2pTau181，保留原临床560人的统一标准化。NULISA NPQ、Soma log2RFU及TMT log2样本/GIS不变。

三方程均为Y_k=α_k+b_kP+Zη_k+ε_k，其中Z含年龄、性别、NULISA/Soma/TMT技术因素及Elecsys批次。M0为基础模型，APOE与BOTH为固定扩展。每个蛋白三模型的个体索引相同。

H1为bN/bS/bT全部严格为正，或全部严格为负。H0是既非全部正亦非全部负，包括一个为0、一个真实反向、只有部分平台关联等。β是协变量条件下的线性投影关联，不是pTau的因果效应。

取SUM3 CR1系数方差子块，z_k=b_k/SE_k，使用本轮验证的正态参考。p_plus=max_k Φ(−z_k)，p_minus=max_k Φ(z_k)，p_shared=min[1,2min(p_plus,p_minus)]。

取最大值确保三检测全部提供证据；乘2处理两个候选方向，避免根据观察符号直接选单侧。若各分量单侧P在相应零假设下有效，则任意检测间依赖下，H0时两方向的错误概率之和不超过α；该结论以分量P的有效性为前提。[联合推断原理](https://www.fil.ion.ucl.ac.uk/spm/doc/papers/NicholsBrett.pdf)。

共同关联检验不需要载荷λ、独特方差ψ或正残差交叉协方差，因此不使用阶段6的等效结构/分母门槛。原等效门槛不合格但本轮支持共同关联的蛋白有：{', '.join(oldineligible)}。它们获得的是新问题的证据，没有被重新判定为响应等效。

## 固定多重检验

每模型34个蛋白级p_shared统一BY q≤.05；不先按方向、单项P、H发现或结构合格性筛选。MAPT单列。全部105个设计通过校准，没有因失败缩小家族。

固定附加结果包括：BH/Holm支持表；34蛋白×正/负共68个方向假设的BY；34共同关联＋34既有H的合并68项BY。主结果始终为新增分支的BY34。合并68用于检验更广主张集合的敏感性，不改变原冻结家族定义。

α34足够阈值={ALPHAS['BY34_SUFFICIENT']:.10g}，α68足够阈值={ALPHAS['BY68_SUFFICIENT']:.10g}。它们保证单项不论其他P排序都足以通过相应BY，不是必要P值，不是实际全BY排序功效。

## 核心候选

{md(candidate)}

NPTX1和KLK6在三个条件模型中均保留共同关联。其P/q降低只说明该统计量的证据变化，不能称为APOE增强了蛋白生物学效应。三个候选的M0最弱检测为Soma；KLK6扩展后最弱检测转为NULISA。PLATFORM_EFFECTS表保留每个平台的系数、原尺度SE、95区间及含残差SD估计不确定性的标准化效应。

APOE_SENSITIVITY_COMPARISONS进一步把q变化分成自身P和其他33项P的贡献，把计算P变化分成系数点估计与SE的两种替换顺序平均贡献。它们是代数分解，不是因果中介或“解释了多少百分比生物信号”。

## 与原测量模式结果整合

M0的34项分为：共同关联与H偏离兼有2项（NPTX1/KLK6）；共同关联、H偏离未定12项；H偏离、共同关联未定1项（PDGFRB）；两问题均未达到发现标准19项。

“H偏离未定”不等于等效；“共同关联未定”不等于无关联。第七阶段在三档工作容差下仍是0个BY响应等效，本轮没有改动容差、门槛或其结果。同一个蛋白可能兼有同向关联及响应比例差异，这正是本项目可讨论的多层测量证据。

## 校准及功效

105个实际设计（M0/APOE/BOTH各35，含单列MAPT）×五个DGP×每单元40000次独立噪声，共525单元、21000000次独立重复。25个零系数组合和28个功效组合复用每次噪声，通过固定X下OLS位置等变性精确平移均值；11.13亿组合评价不冒充11.13亿次独立模拟。

25个零组合包括全零、单零/双零其余强同向、一个反向/微弱反向的全部平台置换和正负镜像。28个功效组合包括平衡.10/.20/.35/.50/.75及单个最弱平台.10/.20/.35、其余.50。均为预设标准化真系数，不把观察效应当功效真值。

五DGP为IID高斯、随P异方差、NULISA轴随机年龄斜率、三轴随机年龄斜率及相应t5；斜率SD=.35。使用完整观测三维残差C的Cholesky，不要求单因子结构。完整投影后二阶矩归一到模板C。固定X、P、缺失子集；协方差模板是已知数据后的条件估计，没有积分其估计不确定性。

预定名义联合门槛为所有零组合Wilson95上界≤.06，零分量单侧.025上界≤.03，失败率≤.005；两个足够阈值处所有零组合Wilson95上界均≤对应alpha，才有正式资格。105/105均通过。

{md(gates)}

全网格最高名义误报率.355%，最高名义Wilson95上界约.4182%，反映在这些情景下偏保守。尾部Wilson上界/alpha最高约.9252。每单元零事件上界约.000096，较前期2000次的尾部精度改善。所有区间为单元区间，不是全网格同时95%保证，也不是任意误差分布和连续参数空间的P值有效性证明。

条件功效示例（M0、平衡正向效应；先取各模板五DGP最低成功率，再汇总34模板）：

{md(power)}

模板百分位不是置信区间，达到80%要求使用每个DGP的Wilson下界，不以观察效应倒算功效。完整不等效应、负向及其他条件模型结果保留在CSV；本轮没有虚构跨蛋白联合DGP来声称实际BY排序功效。

## 识别、验证和论文用途

可以报告：在既定调整、预处理和条件推断范围下，14项三检测同向pTau关联；其成员在APOE和教育/APOE敏感性保持；NPTX1/KLK6兼有基础模型相对响应偏离。

不能据此证明三平台测量相同分子形式、真实完整蛋白浓度或共同潜变量被唯一识别，也不能断言pTau造成蛋白变化。相同受试者、共享病理和检测资料并非三个独立队列；本轮不是独立复制。

原TMT肽段、结合区域、误差相关、分装链及2014MRM标本日期缺口仍按第七阶段保留。本轮没有新增2014候选效应或外部个体重分析。新增结果改善了观测层正面证据，不替代这些分子/实验验证。

## 可复核性及文件

1651项引擎独立核查及{qa['passed']}/{qa['checks']}项综合独立科学核查通过：105实际设计独立SVD/簇方差重拟合，三家族校正与statsmodels核对，525个文件全部53组合×三阈值的独立重计数，4323个前期封存条目未改。

见[主结果与审计工作簿](三检测同向病理关联_结果与审计.xlsx)、[字段导航](字段说明与结果导航.md)、[参考来源](参考文献与方法依据.md)、[运行说明](可复现与运行说明.md)。工作簿覆盖主要结果和聚合审计；两个全量模拟组合汇总仅保留CSV，INDEX明确标注。全部逐次合成噪声统计量在results/simulation_replicates；个体设计和矩仅00_restricted。
"""
 put('STAGE8_MAIN_REPORT_ZH.md',report)
 put('关键决策摘要.md',f"""
# 第八阶段关键结果

2026-10-28。三检测同向病理关联联合检验已完成。

- M0下34项中18项名义共同关联、14项BY共同关联，全部正向。APOE及教育/APOE调整仍为相同14项：{genes}。
- M0合并原H共68项校正，14项仍全部保留；两扩展模型保留13项，NEFL未过更广家族标准。
- NPTX1、KLK6兼有M0共同关联和原相对响应偏离；共同关联在APOE敏感性保留，原H的APOE独立证据仍不足。
- PDGFRB的M0共同关联q=.054806，未过主标准；不作舍入升级。
- 原0个BY响应等效保持不变。14项支持观测层共同病理关联，不代表响应比例等效、分子形式相同或独立复现。

525噪声单元、2100万次独立重复；25零/28功效组合共享重复，不重复计数。105设计通过预设有限门槛；{qa['passed']}项综合科学核查通过，4323条前期封存未改。

[完整报告](STAGE8_MAIN_REPORT_ZH.md)与[工作簿](三检测同向病理关联_结果与审计.xlsx)。
""")
 nav=['# 字段说明与结果导航','','CSV为UTF-8 BOM。原尺度b为每个全局病理SD的系数；scaled_b为b/相应检测残差SD。空值为不适用，不作0；MAPT无家族q。','','## 结果表','']
 for path in sorted(RESULTS.glob('*.csv')):
  if path.name in EXCLUDE:continue
  nav.append(f'- [{path.name}](results/{path.name})：{PURPOSE.get(path.stem,"来源/质量辅助记录")}。'+('全量CSV，不写入Excel。' if path.name in CSV_ONLY else ''))
 nav+=['','## 关键字段','',
 '- p_shared_computed：候选正态参考的双方向联合P；p_shared_for_BY：资格通过后用于正式家族的P。失败不会缩家族，正式P置1。',
 '- q_BY_shared34：新增共同关联34家族的主q；q_BY_direction_*68：正负方向共68假设的q；q_BY_combined68_*：新34共同关联＋旧34H合并q。不能互相替代或混称统一家族。',
 '- prior_H_q_BY_by_model：旧模型对应H q，M0是原主发现；其他是既有敏感性。prior_M0_eq_*仅属于原M0等效，不是APOE/BOTH新等效结果。',
 '- selected_direction为较小方向P的描述，弱或混合符号时仍可能有数值；只有相应发现字段才代表证据。directional68_supported另外检查方向主张。',
 '- weakest_assay为选定方向的最大分量P所在平台；min_directional_scaled_effect_descriptive是点估计描述，不是置信下界或个人分数。',
 '- calibrated/qualified通过有限固定DGP/效应网格，不代表任意误差下普遍有效。Wilson95均为单元区间，不是同时区间。',
 '- reps是独立噪声重复数。按file去重后525×40000=2100万；53组合/3阈值为复用评价，不能相加为独立重复。',
 '- 所有功效效应为预设真值，worst_DGP_rate取五DGP最低率；BY足够阈值概率不是实际全BY排序功效。',
 '- 所有seed在Excel中保存为字符串，避免超过15位整数精度损失。个体ID、P、X及影响函数只在00_restricted。']
 put('字段说明与结果导航.md','\n'.join(nav))
 save_json(RESULTS/'SOFTWARE_VERSIONS.json',{'python':sys.version,'packages':{x:importlib.metadata.version(x) for x in ['numpy','pandas','scipy','statsmodels','matplotlib','openpyxl','pypdf','pyarrow']},'date':'2026-10-28'})
 log('Chinese result report, decision memo and field navigation generated')

def workbook():
 paths=[p for p in RESULTS.glob('*.csv') if p.name not in EXCLUDE]
 priority=['PRIMARY_FAMILY34_SHARED_RESULTS','MODEL_DISCOVERY_COUNTS','ALL35_THREE_MODEL_SHARED_RESULTS','PLATFORM_EFFECTS_ALL_MODELS','APOE_SENSITIVITY_COMPARISONS','EVIDENCE_STATE_COUNTS','CALIBRATION_DECISIONS','CALIBRATION_MODEL_SUMMARY','BALANCED_POWER_FAMILY_SUMMARY']
 paths.sort(key=lambda p:(priority.index(p.stem) if p.stem in priority else len(priority),p.stem));records=[];counter=0
 for path in paths:
  frame=table(path);only=path.name in CSV_ONLY
  if not only:counter+=1
  records.append({'sheet':('CSV_ONLY' if only else f'S{counter:02d}_{path.stem[:27]}'),'source':'results/'+path.name,'rows':len(frame),'columns':len(frame.columns),'delivery':'CSV_ONLY' if only else 'EXCEL_AND_CSV','purpose':PURPOSE.get(path.stem,'来源/辅助审计')})
 index=pd.DataFrame(records);csvout(RESULTS/'WORKBOOK_SHEET_INDEX.csv',index)
 with pd.ExcelWriter(ROOT/'三检测同向病理关联_结果与审计.xlsx',engine='openpyxl') as writer:
  index.to_excel(writer,sheet_name='INDEX',index=False)
  for r in records:
   if r['delivery']!='CSV_ONLY':table(ROOT/r['source']).to_excel(writer,sheet_name=r['sheet'],index=False,inf_rep='inf')
  for ws in writer.book.worksheets:
   ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions;ws.sheet_view.zoomScale=85;ws.row_dimensions[1].height=44
   for cell in ws[1]:cell.fill=PatternFill('solid',fgColor='315D73');cell.font=Font(color='FFFFFF',bold=True);cell.alignment=Alignment(wrap_text=True,vertical='top')
   for row in ws.iter_rows(min_row=2):
    for cell in row:
     if isinstance(cell.value,str):cell.data_type='s';cell.alignment=Alignment(vertical='top',wrap_text=len(cell.value)<180)
     elif isinstance(cell.value,float):cell.number_format='0.000000E+00' if 0<abs(cell.value)<1e-4 else '0.0000000000'
   for j in range(1,ws.max_column+1):
    size=max((len(str(ws.cell(i,j).value or '')) for i in range(1,min(ws.max_row,40)+1)),default=14)
    ws.column_dimensions[get_column_letter(j)].width=min(52,max(14,size+2))
  ws=writer.book['INDEX'];ws.column_dimensions['B'].width=64;ws.column_dimensions['F'].width=85
  for i,r in enumerate(records,2):
   if r['delivery']!='CSV_ONLY':ws.cell(i,1).hyperlink=f"#'{r['sheet']}'!A1";ws.cell(i,1).font=Font(color='0563C1',underline='single')
 print({'workbook_data_sheets':counter,'workbook_total_sheets':counter+1,'full_simulation_tables_CSV_only':len(CSV_ONLY),'public_tables_indexed':len(records)})
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--docs',action='store_true');parser.add_argument('--workbook',action='store_true');args=parser.parse_args()
 if args.docs:docs()
 if args.workbook:workbook()
