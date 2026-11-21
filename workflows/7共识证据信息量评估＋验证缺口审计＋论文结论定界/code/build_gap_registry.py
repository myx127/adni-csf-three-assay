from config import *

def run():
 rows=[
  ('G01','2014 MRM发布与技术交叉表','CLOSED_WITH_PRECISION_QUALIFICATION','303行各320肽段均在6有效数字末位1单位内与唯一原工作簿记录一致；LIMS/15处理日/2去高丰度柱/进样序号已恢复','保留数值精度差异；两个发布不能计作两个实验','无需新临床效应拟合；独立复核全向量身份和元数据','提高下一步验证可执行性，不构成验证成功'),
  ('G02','2014 MRM精确采样日期','PARTIALLY_CLOSED_NOT_SPECIMEN_VERIFIED','287 RID/基线访视均有唯一MASTER采样日候选；149人与Stage5同访视且候选日期一致','缺MRM样本/LIMS到实际采样日期或带日期生物样本记录的直接桥接','LIMS_ID/External_Reference→specimen_id/aliquot_id→collection_date/time及访视；不要用运行日期替代','尚不能标为严格同日复现'),
  ('G03','2014潜在配对的设计与有限簇推断','CONDITIONAL_GEOMETRY_AVAILABLE_NOT_CALIBRATED','假定访视日期映射有效时，NPTX1两肽段n144/df78，KLK6三肽段n141/df75；N/S/处理日簇11/36/15','尚未验证日期假设和2014模型误拒绝；两去高丰度柱的依赖结构仍需审计','日期核实后冻结2014具体设计与重复规则，评估处理日/柱依赖并校准；然后才检验候选','潜在验证机会，不能把几何门槛通过当正式检验已合格'),
  ('G04','2019 MRM同日直接配对','OPEN_SAMPLE_INFORMATION','Stage6同日90–94人，df20–24，低于冻结df40；新RID95人的关联是另一问题','完整技术调整后信息不足','新增独立配对个体/技术簇，或新增经方法学验证的不同设计；现有冻结结果保持','仅描述性直接测量及有限单检测关联支持'),
  ('G05','PDGFRB直接检测','LOCAL_GAP_WITH_EXTERNAL_GENE_LEVEL_LEAD','本地靶向面板无P09619直接肽段；公开Emory2022注释表列PDGFRB于Soma/Olink/MS；PDGF-BB配体不代替受体','外部仅过滤后基因级覆盖；个体完整度、受体肽段/抗体识别区域未核实','取得该外部矩阵/受体肽段及同人配对、采样/批次；或新增验证过的P09619直接检测','发现外部覆盖线索，尚不能宣布PDGFRB正式验证成功'),
  ('G06','原TMT肽段及蛋白分组证据','OPEN_SOURCE_GRANULARITY','现有TMT压缩档仅聚合蛋白矩阵，精确肽段、序列覆盖与分组未提供','无法区分胞外/胞内肽段、异构体共享和特定裂解位置','肽段/PSM级矩阵、序列/修饰、蛋白组规则、唯一性、样本及批次映射','不能把PDGFRB剪切、NPTX1形式、KLK6活化当已证实机制'),
  ('G07','NULISA抗体/Soma结合区域','OPEN_REAGENT_INFORMATION','基因/登录号与适配体编号可核实；实际结合区域未提供','未知表位/结构或结合受基质/复合物影响情况','抗体对/适配体识别区域、交叉反应、形式/修饰响应、必要竞争/回收实验','共同注释不等于三检测观察相同完整分子'),
  ('G08','同一样本分装及误差相关','OPEN_EXPERIMENTAL_IDENTIFICATION','同RID/日期不证明同一分装；无足够同样本技术重复直接标定所有误差协方差','任意相关误差下真实共同丰度与检测特异病理作用不唯一','分装链、重复检测、前分析随机化、技术因素、必要标准/形式特异实验','保留条件性测量模型；本阶段信息量评估不解决识别'),
  ('G09','真正独立队列的三检测模式复现','OPEN_NOT_EXECUTED','登记外部三平台/双平台和PPMI资源线索；没有新的外部个体级复现结果','目标覆盖/同人同日、独立pTau、技术元数据、独立性和足够设计均需逐项核实','按EXTERNAL_VALIDATION_OPPORTUNITIES.csv逐项筛；获取合法可用数据后先审计/冻结/校准','文献、同ADNI新RID和另一检测均不能计为外部模式复现'),
  ('G10','病理相关共识的正面证据','OPEN_EVIDENCE_NOT_NEGATIVE_PROOF','34蛋白内11等效可评价；最宽容差6名义等效、0 BY；另有15个结构可作条件功效模板','共同方向/响应等效/病理关联须分层；容差为事后工作值','使用固定容差报告区间、功效情景和信息扩展假设；更强病理关联证据需另冻结问题','零BY等效不能写成没有共同信号或所有蛋白均不一致')]
 csvout(RESULTS/'VALIDATION_GAP_REGISTRY.csv',[dict(zip(['gap_id','issue','status','available_evidence','remaining_gap','minimum_next_evidence_or_action','claim_impact'],r)) for r in rows])
 external=[
  ('X01','Emory ADRC三平台研究2022','https://doi.org/10.1186/s13195-022-01113-5','36人，三平台完整者上限35；Soma/Olink/TMT，非NULISA','论文说明原数据/病例信息/代码存Synapse 3platformEmory；补充表仅审计注释','过滤后表中PDGFRB/KLK6见于三平台；NPTX1见于Soma/MS，Olink未列；未审计个体矩阵','同人三检测框架的外部测量支持线索，PDGFRB潜在覆盖','n35上限低于冻结n80；每基因完整人数、独立pTau、技术元数据、肽段及测量版本尚需核实','ANNOTATION_LEAD_NOT_READY_FOR_FROZEN_FORMAL_REPLICATION'),
  ('X02','NIH双平台CSF/血浆研究2026预印本','https://doi.org/10.64898/2026.09.02.26361977','67人；SomaScan11K＋NULISAseq，混合神经退行疾病','作者声明合理请求可获得数据，非本地现成数据','仅方法及可用性核查，未申请或下载个体矩阵','跨亲和检测测量支持线索','仅两平台；无第三正交检测；目标覆盖与独立pTau未核实；不能直接检验原三检测约束','METHOD_SUPPORT_ONLY_NOT_FORMAL_REPLICATION'),
  ('X03','AMP-PDRD/PPMI蛋白组资源','https://amp-pdrd.org/data','官方列CSF Olink及质谱/肽段/SDRF等元数据','平台提供注册用户访问工作区/数据；本轮无账号操作','官方方法和可用性核查，未取得配对个体数据','不同疾病队列中的检测结构可迁移性线索','PD队列不能默认验证AD特异pTau假设；同人交集、候选覆盖及日期/批次需核实','ACCESS_AND_APPLICABILITY_AUDIT_REQUIRED')]
 csvout(RESULTS/'EXTERNAL_VALIDATION_OPPORTUNITIES.csv',[dict(zip(['source_id','resource','primary_source','reported_design','availability','audit_completed','possible_use','limitations','readiness'],r)) for r in external])
 claims=[
  ('跨平台共同变化','可报告条件性共同方向及其结构/精度','不能等同真实浓度已分离；共享误差仍可能贡献'),
  ('病理响应等效','原固定34家族三容差均0 BY；不确定者保留未定','不能由未拒绝偏离推断等效；±0.2/0.35/0.5为工作容差'),
  ('精确零偏离与实用等效','二者是假设不同的两个证据轴，可并列报告','真实偏离零但幅度足够小，逻辑上可同时拒绝零且通过宽容差等效；不强行视作互斥'),
  ('病理相关共识信号','须另有病理关联证据，等效仅说明相对响应接近','gamma0模拟即使等效成立，也不表示存在病理关联'),
  ('功效与扩样','可报告15个既定协方差模板的真等效情景检出率、独立面板倍数和MC区间','不能将alpha*概率当实际BY功效；不能把V/k用于仅增加旧簇内人数'),
  ('2014验证机会','技术映射已闭合，149同访视日期候选，141–144潜在有效配对','MRM采样日未核实，推断未校准，尚未做新候选关联或正式复现')]
 csvout(RESULTS/'EVIDENCE_INTERPRETATION_BOUNDARIES.csv',[dict(zip(['topic','permitted_statement','limit'],r)) for r in claims])
 log('Ten validation gaps, three external resource leads and evidence interpretation boundaries registered; no new external validation claimed')

if __name__=='__main__':run()
