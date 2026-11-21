from config import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logging,importlib.metadata,sys,platform
logging.getLogger('fontTools.subset').setLevel(logging.ERROR)
from openpyxl import load_workbook
from openpyxl.styles import Font,PatternFill,Alignment
from openpyxl.formatting.rule import ColorScaleRule

def tables():
 main=pd.read_csv(RESULTS/'PRIMARY_ALL43_RESULTS.csv');variants=pd.read_csv(RESULTS/'ALL_FIXED_VARIANT_RESULTS.csv')
 family=main.loc[main.in_frozen_family].copy();hits=family.loc[family.BY_discovery].sort_values('p_value')
 cov=pd.read_csv(RESULTS/'COMMON_COVARIANCE_COMPATIBILITY_DIAGNOSTICS.csv')
 grid=pd.read_csv(RESULTS/'COVARIANCE_OFFSET_SENSITIVITY_SUMMARY.csv')
 evidence=hits[['accession','gene','n','p_value','q_BY','q_BH','p_Holm']].copy()
 for variant in ['RAW_P','RAW_TMT','LOD','BATCH9','FULLCOV','AMY']:
  q=variants.loc[variants.variant.eq(variant),['accession','n','p_value','q_BY']].rename(columns={'n':variant+'_n','p_value':variant+'_p','q_BY':variant+'_q_BY'})
  evidence=evidence.merge(q,on='accession',how='left')
 evidence=evidence.merge(cov,on=['accession','gene']).merge(grid,on='accession')
 csvout(RESULTS/'DISCOVERY_EVIDENCE_MATRIX.csv',evidence)
 summary=[]
 for label,frame in [('MAIN',family)]+[(v,g.loc[g.in_frozen_family]) for v,g in variants.groupby('variant')]:
  summary.append({'analysis':label,'m_frozen_family':len(frame),'n_BY05':int(frame.q_BY.le(.05).sum()),'n_BH05':int(frame.q_BH.le(.05).sum()),
   'min_p':frame.p_value.min(),'min_n':int(frame.n.min()),'median_n':float(frame.n.median()),'max_n':int(frame.n.max()),
   'main_discoveries_retained_BY':';'.join(hits.loc[hits.accession.isin(frame.loc[frame.q_BY.le(.05),'accession']),'gene'])})
 csvout(RESULTS/'ANALYSIS_COMPARISON_SUMMARY.csv',summary)
 screen=pd.read_csv(RESULTS/'PRE_EFFECT_PROTEIN_SCREEN.csv')
 powerstatus=screen[['accession','gene','single_annotation','role','n','design_status','exclusion_reasons']].copy()
 powerstatus['power_simulated']=powerstatus.design_status.eq('DESIGN_QUALIFIED')
 powerstatus['interpretation']=np.where(powerstatus.power_simulated,'Fixed synthetic-noise conditional power, six effects and two target assays','No primary-model power estimate: prespecified quality/design or mapping exclusion; not evidence of no biological effect')
 csvout(RESULTS/'POWER_EVALUABILITY_ALL70.csv',powerstatus)
 stress=pd.read_csv(RESULTS/'IDENTIFICATION_STRESS_RESULTS.csv')
 csvout(RESULTS/'STRESS_SCENARIO_SUMMARY.csv',stress.groupby('scenario').agg(n_proteins=('accession','size'),min_rate=('rejection_rate','min'),median_rate=('rejection_rate','median'),max_rate=('rejection_rate','max'),max_wilson_upper=('wilson95_high','max'),n_failed=('n_failed','sum')).reset_index())
 refs=[
 ('REF01','Cameron, Gelbach & Miller (2011). Robust Inference With Multiway Clustering. JBES29:238–249','https://cameron.econ.ucdavis.edu/research/JBESpaper2009version.pdf','原始论文作者预印本；经典CGM包含交集项。本文SUM3不是CGM，不能混称。'),
 ('REF02','MacKinnon, Nielsen & Webb. Cluster-robust inference: A guide to empirical practice. Journal of Econometrics (2023)','https://doi.org/10.1016/j.jeconom.2022.04.001','固定效应不自动消除簇内相关；杠杆及簇不均衡影响推断。'),
 ('REF03','Pustejovsky & Tipton (2018). Small-sample methods for cluster-robust variance estimation and hypothesis testing in fixed effects models','https://arxiv.org/abs/1601.01981','小簇/固定效应推断需专门处理；本轮未将自定义版本冒称CR2。'),
 ('REF04','Benjamini & Yekutieli (2001). The control of the false discovery rate in multiple testing under dependency. Ann Stat29:1165–1188','https://www.math.tau.ac.il/~ybenja/depApr27.pdf','BY依赖稳健调整以有效边际P为前提，不修复模型错误。'),
 ('REF05','ADNI PET documentation','https://adni.loni.usc.edu/quick-start-guide-asset/pet.html','采集日期匹配、QC和Centiloid跨示踪剂口径；本地文件版本见Stage4。')]
 csvout(RESULTS/'METHODS_EVIDENCE_REGISTRY.csv',pd.DataFrame(refs,columns=['id','title','source','supports_and_limits']))
 definitions=[
 ('主问题','检验病理相关的三检测线性均值系数与残差交叉协方差能否满足两条相对约束；不是单一真实丰度/偏差的识别。'),
 ('H1/H2','H1=bS*cNT−bN*cST；H2=bT*cNS−bN*cST。群体约束，不是个体评分。'),
 ('b_per_global_P_SD','同一完整模型中的病理P系数；P在临床560人中统一标准化。'),
 ('b_over_residual_SD','b除该检测模型残差SD；跨检测描述量，包含分母估计不确定度的95%区间，未作多重校正。'),
 ('q_BY','预先冻结家族的BY校正值；主34蛋白，淀粉样次级33蛋白，各为单独家族。'),
 ('q_BH/p_Holm','并列支持；不用于替代BY重定义主要发现。'),
 ('条件性发现','按冻结模型/检验规则成立；不表示对所有协变量、误差相关或测量机制都稳健。'),
 ('NULISA检测80%门槛','本阶段在实际蛋白效应读取前明确采用，需571人整体和该蛋白分析子集同时满足。'),
 ('rank_X/residual_df','完整共同设计的实际秩与n减秩，非虚拟变量名义列数。'),
 ('V0/V1/V2','V0三一向CR1协方差和；V1另乘(n−1)/df；V2=2V1。开发选V0，确认独立。'),
 ('条件功效','固定真实P/X/缺失/簇，σ=.6/.8/1.2、三向斜率.35等指定合成参数下的功效；不是实际蛋白可靠性估计或观察效应倒算功效。'),
 ('BY_FIRST_DISCOVERY','按.05/(m*H_m)评估单蛋白首发现阈值，不能直接视为整个BY程序功效。'),
 ('delta80_grid','所测试效应网格中首次达到80%的位置；不插值，不证明未测效应的MDE。'),
 ('协方差工作偏移网格','rho乘观察残差SD乘积的125个离散偏移；含SD导数；非已校准的检测误差相关范围。'),
 ('协方差相容性','正交叉协方差与非负独特方差是简单单因子/独立误差解释的点诊断；不作事后筛除或新发现标准。'),
 ('LOCO','每蛋白每技术因子删去H影响能量最大的一个簇，重新拟合；诊断性、未独立再校准，不能挑选有利删除结果。'),
 ('GIS与raw等价','同一人群且完整plex固定效应时，减去逐蛋白plex常数被固定效应完全吸收；不是独立复现。'),
 ('认知/tau PET','本轮未创建个体F/D、未分析认知斜率；无同期tau PET不能以多年后状态回填。')]
 csvout(RESULTS/'METRIC_DICTIONARY_ZH.csv',pd.DataFrame(definitions,columns=['term','definition']))
 save_json(RESULTS/'SOFTWARE_ENVIRONMENT.json',{'python':sys.version,'platform':platform.platform(),'packages':{p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','statsmodels','pyarrow','openpyxl','matplotlib','Pillow']},'BLAS_threads_per_worker':1,'simulation_workers':6})

def figures():
 plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':10,'pdf.fonttype':42})
 def save(fig,name):
  fig.savefig(FIGURES/(name+'.png'),dpi=220,bbox_inches='tight');fig.savefig(FIGURES/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
 main=pd.read_csv(RESULTS/'PRIMARY_FROZEN_FAMILY_RESULTS.csv').sort_values('p_value');hits=main.loc[main.BY_discovery]
 fig,ax=plt.subplots(figsize=(9,10));y=np.arange(len(main));color=np.where(main.BY_discovery,'#b55742','#315f89')
 ax.barh(y,-np.log10(main.q_BY.clip(lower=1e-10)),color=color);ax.invert_yaxis();ax.axvline(-np.log10(.05),ls='--',color='#555',lw=1)
 ax.set(yticks=y,yticklabels=main.gene,xlabel='−log10(BY q)，虚线为q=0.05',title='冻结主家族34蛋白：3个条件性相对约束偏离');ax.grid(axis='x',alpha=.15);save(fig,'01_PRIMARY_BY_RESULTS')
 conf=pd.read_csv(RESULTS/'CONFIRM_SIMULATION_SUMMARY.csv');conf=conf.loc[conf.version.eq(0)&conf.alpha.eq(.05)]
 worst=conf.groupby('design_id').agg(max_rate=('rate','max'),max_upper=('wilson95_high','max')).reindex(main.accession)
 fig,ax=plt.subplots(figsize=(9,4.8));xx=np.arange(len(worst));ax.plot(xx,worst.max_rate,'o',label='每蛋白六场景最大拒绝率',ms=4,color='#245c8c');ax.plot(xx,worst.max_upper,'_',label='对应最大Wilson95%上界',color='#c18039')
 ax.axhline(.05,ls='--',c='#666',label='名义5%');ax.axhline(.06,ls=':',c='#a33',label='预定校准审核上界6%');ax.set(xlabel='主家族蛋白（顺序同主结果图）',ylabel='合成兼容零场景拒绝率',ylim=(0,.065),title='独立确认：实际设计校准偏保守，不是普遍有效性的证明');ax.legend(fontsize=8,ncol=2);save(fig,'02_REAL_DESIGN_NULL_CALIBRATION')
 fig,axes=plt.subplots(1,len(hits),figsize=(11,4.6),sharex=True)
 if len(hits)==1:axes=[axes]
 for ax,(_,r) in zip(axes,hits.iterrows()):
  for j,(label,col) in enumerate(zip(['NULISA','Soma','TMT/GIS'],['N','S','T'])):
   b=r[col+'_b_over_residual_SD'];lo=r[col+'_scaled_b_CI95_low'];hi=r[col+'_scaled_b_CI95_high'];ax.errorbar(b,j,xerr=[[b-lo],[hi-b]],fmt='o',capsize=3,color=['#2d648e','#cb893d','#158779'][j])
  ax.axvline(0,ls=':',c='#888');ax.set(yticks=range(3),yticklabels=['NULISA','Soma','TMT/GIS'],title=f'{r.gene}\nn={int(r.n)}, BY q={r.q_BY:.3f}',xlabel='每1SD log2-pTau的系数 / 残差SD');ax.invert_yaxis()
 fig.suptitle('三检测调整后系数：区间为逐参数95%，不作单独发现判定',y=1.04);fig.tight_layout();save(fig,'03_THREE_DISCOVERY_COEFFICIENTS')
 variants=pd.read_csv(RESULTS/'ALL_FIXED_VARIANT_RESULTS.csv');cols=['MAIN','RAW_P','RAW_TMT','LOD','BATCH9','FULLCOV','AMY'];arr=[]
 for col in cols:
  part=main.set_index('accession').q_BY if col=='MAIN' else variants.loc[variants.variant.eq(col)].set_index('accession').q_BY
  arr.append(part.reindex(main.accession).to_numpy())
 values=np.array(arr).T;plot=np.ma.masked_invalid(-np.log10(np.maximum(values,1e-4)));fig,ax=plt.subplots(figsize=(9.5,10.8));im=ax.imshow(plot,cmap='YlOrRd',vmin=0,vmax=2,aspect='auto')
 ax.set(xticks=range(len(cols)),xticklabels=['主分析','原值pTau','raw-TMT','逐人LOD','Batch9','教育/APOE','淀粉样次级'],yticks=range(len(main)),yticklabels=main.gene,title='固定敏感性：颜色为−log10(BY q)，圆点为q≤0.05')
 plt.setp(ax.get_xticklabels(),rotation=35,ha='right')
 for i in range(len(main)):
  for j in range(len(cols)):
   if np.isfinite(values[i,j]) and values[i,j]<=.05:ax.plot(j,i,'o',ms=4,color='#20324a')
 fig.colorbar(im,ax=ax,label='−log10(BY q)；灰/空缺为未进入该家族或设计不足');save(fig,'04_FIXED_SENSITIVITY_MATRIX')
 power=pd.read_csv(RESULTS/'POWER_SUMMARY_BY_EFFECT.csv');fig,axes=plt.subplots(1,2,figsize=(10.5,4.8),sharey=True)
 for target,ax in zip([2,3],axes):
  for label,legend,color in [('NOMINAL05','名义5%','#28638e'),('BY_FIRST_DISCOVERY','BY首发现阈值','#bc7935')]:
   p=power.loc[power.target.eq(target)&power.threshold.eq(label)].sort_values('delta_SD');x=p.delta_SD.to_numpy();ax.plot(x,p.median_power,'o-',color=color,label=legend);ax.fill_between(x,p.min_power.to_numpy(),p.max_power.to_numpy(),alpha=.15,color=color)
  ax.axhline(.8,c='#555',ls='--',lw=1);ax.set(xlabel='预设直接偏离δ（零偏离残差SD单位）',title=f'偏离位于M{target}（'+('Soma' if target==2 else 'TMT')+'）',ylim=(0,1.03));ax.legend(fontsize=8)
 axes[0].set_ylabel('条件功效；线为34蛋白中位数，带为跨蛋白范围');fig.suptitle('真实设计＋固定合成噪声参数：不使用观察效应倒算功效');fig.tight_layout();save(fig,'05_CONDITIONAL_POWER')
 profile=pd.read_csv(RESULTS/'COVARIANCE_OFFSET_SENSITIVITY_SUMMARY.csv').set_index('accession').reindex(main.accession)
 fig,ax=plt.subplots(figsize=(9,10));y=np.arange(len(main));ax.hlines(y,profile.p_grid_min.clip(lower=1e-12),profile.p_grid_max,color='#a8b4c1',lw=2);ax.scatter(main.p_value,y,c=color if isinstance(color,str) else '#315f89',s=18,label='零偏移主P')
 ax.axvline(.05,c='#a44',ls='--',lw=1);ax.set(xscale='log',yticks=y,yticklabels=main.gene,xlabel='点式P；灰线为预设125点工作网格最小至最大值',title='协方差假设敏感性：灰线不是置信区间或连续稳健界');ax.invert_yaxis();save(fig,'06_COVARIANCE_WORKING_GRID')

def reports():
 main=pd.read_csv(RESULTS/'PRIMARY_FROZEN_FAMILY_RESULTS.csv');hits=main.loc[main.BY_discovery].sort_values('p_value')
 sap=json.loads((RESULTS/'MAIN_ANALYSIS_SAP_FROZEN.json').read_text(encoding='utf-8'));qa=json.loads((RESULTS/'INDEPENDENT_QA_SUMMARY.json').read_text(encoding='utf-8'))
 conf=pd.read_csv(RESULTS/'INDEPENDENT_CONFIRMATION_DECISIONS.csv');power=pd.read_csv(RESULTS/'POWER_SUMMARY_BY_EFFECT.csv')
 evidence=pd.read_csv(RESULTS/'DISCOVERY_EVIDENCE_MATRIX.csv');comparison=pd.read_csv(RESULTS/'ANALYSIS_COMPARISON_SUMMARY.csv').set_index('analysis')
 hitlines='\n'.join(f'| {r.gene} | {int(r.n)} | {r.Q:.3f} | {r.p_value:.6g} | {r.q_BY:.5f} |' for r in hits.itertuples())
 senslines='\n'.join(f'| {r.gene} | {r.RAW_P_q_BY:.4f} | {r.BATCH9_q_BY:.4f} | {r.FULLCOV_q_BY:.4f} | {r.AMY_q_BY:.4f} | {r.p_grid_max:.3f} |' for r in evidence.itertuples())
 powerlines='\n'.join(f'| M{r.target} | {r.delta_SD:g} | {100*r.median_power:.1f}% | {100*r.min_power:.1f}–{100*r.max_power:.1f}% |' for r in power.loc[power.threshold.eq('BY_FIRST_DISCOVERY')&power.delta_SD.isin([.35,.5,.75])].itertuples())
 report=f'''# Stage5：真实设计完整模型校准、逐蛋白功效与主分析

本轮已完成预定病理主分析。结论是：在冻结的年龄/性别、三平台技术批次和Elecsys批次调整模型下，34蛋白中PDGFRB、NPTX1、KLK6满足BY q≤0.05；这些是有明确假设条件的跨检测相对线性均值约束偏离，不能解释为已经识别的病理特异检测偏差或真实丰度分解。扩展教育/APOE后均不再通过BY，淀粉样PET次级家族无发现，协方差工作偏移网格也使三者的点式P跨越0.05。论文应完整报告此不稳定性。

完成日期2026-10-27。原始数据、Stage1–4、旧SAP和论文未改写；在本Stage5创建单独的新协议与冻结SAP。用户本轮明确授权定量主分析。首次蛋白定量连接发生在SAP冻结及独立确认之后，连接和个体值仅在00_restricted。没有创建个体F/D，没有分析认知斜率或诊断预测。

**先看结果及边界。** 主家族采用同日Elecsys pTau181的log2，并在临床完整560人中统一标准化。每蛋白使用三检测、病理和协变量共同完整的同一组人；主家族人数范围{int(main.n.min())}–{int(main.n.max())}，中位数{main.n.median():g}。NULISA保持NPQ，Soma用log2 RFU，TMT用log2(sample/GIS)。不按结果选适配体、变换或匹配窗口。

| 主发现蛋白 | n | 2df统计量Q | 原始P | BY q（m=34） |
|---|---:|---:|---:|---:|
{hitlines}

BH另有CCL2通过，BY未通过；主发现不据此扩大。Holm支持上述三项。MAPT预先作为与pTau同母蛋白构念对照单列，未纳入34蛋白家族，其主模型P约0.356；不能把这一检验当作检测所有tau形式差异的通用工具。

| 主发现 | 原值pTau q | Batch9-only q | 加教育/APOE q | 淀粉样次级q | 工作协方差网格最大点式P |
|---|---:|---:|---:|---:|---:|
{senslines}

逐人LOD敏感性中三项与主分析相同，说明这些候选的相关记录没有因该规则被移除，不能据此声称低LOD问题已被全面解决。raw-TMT的相同结果也不是独立稳健性证据：同一分析人群中，逐蛋白GIS是plex常数，完整plex固定效应将其完全吸收；35个合格设计的同人群残差等价检查全部通过。存在GIS缺失的人群不同仍须分别报告。

原值pTau下PDGFRB不再过BY，而NPTX1/KLK6保留；不同函数形式对应不同的线性投影问题，不能根据哪个显著更有利来改主方案。加入教育/APOE是预设的更强条件模型，三项都不再过BY。淀粉样PET对应另一病理轴、较小人数和不同时间匹配；没有发现不能称为证明不存在关联，但也不能称为独立支持。

**共同方向的解释限制。** KLK6的调整后Soma–TMT残差相关约−0.246，而另外两对为正；样本交叉协方差符号不支持一个简单实单因子＋互不相关误差的解释。PDGFRB在该简单模型的点分解中出现Soma独特方差略为负的Heywood现象，提示边界/模型不相容，不能据此给出稳定“共识丰度”。NPTX1的点协方差结构可容纳该简单解释，但这不证明误差独立。所有这些都是事后解释诊断，未改变家族、主P或发现标准；负独特方差点估计也不是单独的显著性检验。

各蛋白每技术轴删除H影响最大的一个簇后，主发现的P和校正值均已报告。删除各自影响最大的Soma簇时，三项的诊断性BY约0.076；说明临界发现对样本/影响配置敏感。这些删除结果未做独立再校准，仅作影响诊断，不替换主结果，也不选择最有利的删簇结果。

**如何冻结模型。** 在实际蛋白效应读取前，从43个单特征注释候选筛到35个设计合格蛋白：n≥200、完整设计残差df≥100、每技术轴≥10簇、NULISA在571人及分析子集中检出均≥80%。80%门槛在本阶段明确采用；此前Stage1/4的80%人数仅是描述。MME、UCHL1、SNAP25、YWHAZ、S100A12、BASP1、PDLIM5和TEK被排除，原因逐项保留。MAPT单列后主发现家族34。26个多适配体及HBA1组合注释未默认聚合，不进入本轮正式主家族。

三检测共用X=[截距、P、年龄、性别、NULISA板、Soma板、TMT plex、Elecsys批次]；避免每平台分别预处理后改变约束。年龄保留月/年精度，沿用8人出生冲突暂缺。完整X以实际矩阵秩而非名义列数处理，记录杠杆及P的偏杠杆。病理批次保留，Batch9-only是预设敏感性；教育/APOE扩展协变量不改变主家族。

令b为三个P系数，c为同一完整投影后的交叉协方差，c以n−rank(X)作尺度修正。主检验为H1=bS*cNT−bN*cST，H2=bT*cNS−bN*cST。影响函数同时包括b和c的估计不确定度；技术均值由完整X处理，三向方差不能代替均值调整。此法使用P估计共同方向，是联合测量–病理检验，不是完全病理盲的F/D构造。

V0是三个一向CR1协方差之和；V1另乘(n−1)/df；V2=2V1。开发按V0→V1→V2选择全部合格设计/场景通过的最小方差版本，得到V0。V0不含经典CGM交集扣除项，不能称为CGM；也未称为CR2。多向簇理论与实践背景见[原始多向簇论文](https://cameron.econ.ucdavis.edu/research/JBESpaper2009version.pdf)、[簇推断实践指南](https://doi.org/10.1016/j.jeconom.2022.04.001)及[小样本固定效应研究](https://arxiv.org/abs/1601.01981)。本轮实际校准结果本身才是所用自定义版本的范围证据。

**校准证据。** 开发210000次，独立确认840000次，35设计×6场景。场景含不等精度、异方差、固定效应后仍存在的19板随机斜率、三向随机斜率及t5重尾。确认V0最大名义5%拒绝率{100*conf.max_null_rate.max():.3f}%，最大Wilson95%上界{100*conf.max_null_wilson95_upper.max():.3f}%，无计算失败；全部通过预定上界6%的审核。145个额外不同敏感性/次级设计另有3480000次确认，210个支线记录中209个设计合格并通过，1个淀粉样设计因结构不足保留排除状态。完全相同设计按指纹复用，不能把复用计算重复计为新模拟。

主/支线还报告.005、.001尾部计数与区间。有限模拟不是任意真实误差结构下的P均匀性证明，零次拒绝也不表示尾概率为0。BY的依赖稳健性仍以有效边际P为前提，不能修复模型错设或测量相关误差。[Benjamini–Yekutieli原始论文](https://www.math.tau.ac.il/~ybenja/depApr27.pdf) 对应的是多重检验调整，而非本轮测量识别的证明。

**逐蛋白条件功效。** 840000次模拟，35设计×M2/M3×6效应×2000次。固定实际P、X、缺失子集和簇；合成噪声σ=.6/.8/1.2，共同T=.55P+U，三向随机斜率=.35，零偏离时预期回归残差SD归一为1。没有使用观察效应倒算功效，也未宣称这些噪声参数已代表每个实际蛋白的可靠性。下表为34个主家族蛋白在BY首发现阈值下的功效；这是单蛋白工作阈值，不是完整BY程序的家族功效。

| 偏离检测 | δ（模拟残差SD） | 中位功效 | 跨蛋白范围 |
|---|---:|---:|---:|
{powerlines}

因此名义5%功效高不等于经多重检验发现能力充分。完整逐蛋白表给出名义、Bonferroni、BY首发现阈值及Monte Carlo区间；δ80只在测试网格定位，不外推。未检出不构成跨检测一致的证据。

**识别压力验证。** 另有560000次预设压力模拟。弱共同信号下零拒绝，体现信息不足；本组固定误差相关和固定非线性场景也未表现为大量拒绝，说明两约束不能检出全部这些机制。遗漏技术均值的拒绝率约43.15%–80.325%，说明方差处理不能挽救均值混杂。不能将所有δ=0情景一概当作可观测H=0；误差相关、非线性和遗漏变量可改变被检验约束。

125点协方差偏移网格包含三个off-diagonal偏移及残差SD的估计导数。三项主发现均存在不显著的网格点，最大点式P见上表。该网格未由实际误差相关数据校准，不是连续稳健界或识别集合；也不能因选中某个网格点不显著就证明没有效应。它说明机制归因依赖未被额外信息识别的假设，Stage2–3观测等价限制仍然成立。

**独立复核与文件。** [独立QA](results/INDEPENDENT_QA_SUMMARY.json)为{qa['passed']}/{qa['checks']}项通过，35个主设计均用独立完整OLS/约束/簇方差公式重新拟合，BH/BY/Holm与statsmodels核对，保留模拟逐文件重数。总合成重复5930000次，其中开发、确认、支线、功效及压力分别计数，无将重放计为新重复。输入及此前关键产物散列、490项重建掩码、模型梯度和仿射不变性均核对；个体值仅在restricted。

请从[Excel总表](真实设计模型校准_逐蛋白功效与主分析总表.xlsx)查看全部70候选状态、43主候选完整结果、34家族、模拟、功效、敏感性、次级、解释诊断及QA。图提供PNG/PDF。[原始新协议](STAGE5_PROTOCOL_v1.0.md)、[冻结主SAP](results/MAIN_ANALYSIS_SAP_FROZEN.json)、[首次蛋白访问记录](results/FIRST_PROTEIN_VALUE_ACCESS.json)和最终SHA清单形成执行链。出生冲突8人、TMT3447/3907发布差异及分装身份仍是来源限制，没有据猜测补齐或发送外部邮件。

**论文当前可以写到哪里。** 可以报告预先限定的跨检测相对线性均值约束、真实设计下的校准和功效、3个条件性pTau相关候选及完整不稳健性证据。NPTX1在协方差点结构和原值P支线中相对更易解释，但仍受教育/APOE、影响簇和误差假设限制；不能把其提升为已验证机制。PDGFRB对P函数与扩展协变量敏感；KLK6首先表现为基本协方差结构不相容。现有结果不足以宣称普遍分离共识丰度与病理特异检测误差。进一步机制归因需要额外分子/重复检测或独立样本证据，不能靠继续挑模型消除这些限制。
'''
 (ROOT/'STAGE5_MAIN_REPORT_ZH.md').write_text(report,encoding='utf-8')
 methods='''# 完整模型与约束检验的实现

先用共同Z（截距、协变量及所有预定技术因子）对P残差化，rP=MZ P，a=rP/(rP′rP)。对每个检测使用相同X=[Z,P]，b=a′Y，E=MX Y。cjk=Σeij eik/df，df=n−rank(X)。这是自由度尺度修正，不保证复杂簇相关下无偏。

两约束为h=(bS*cNT−bN*cST，bT*cNS−bN*cST)。b的逐人贡献为a_i e_i；cjk的逐人贡献为(eij eik−平均残差乘积)/df。用解析Jacobian合成h贡献。每一技术轴按簇聚合并中心化，以G/(G−1)作CR1，三项相加为V0。Q=h′V0^−1h，以2df chi-square尾概率计算。无负特征值修补，病态矩阵保留失败。开发预定还比较V1=V0*(n−1)/df、V2=2V1；全局选择V0并独立确认。

独立加权经验扰动同时重拟合所有 nuisance 系数与c，验证影响函数；另用独立OLS、FWL、数值Jacobian、缩放不变性及完整簇求和验证。工作协方差偏移cjk−rho_jk*sqrt(cjj*ckk)的Jacobian包括对两个方差分量的导数，不将观察残差SD视为已知常数。

H=0是两个可观测相对线性均值约束，不等于任意尺度、载荷、非线性或个体评分的全面测量不变性。任意跨检测相关误差下，真实丰度、特异响应和误差相关不能凭三个单次检测唯一分离。简单单因子点诊断要求交叉协方差符号相容，且λ²分解不产生负独特方差；这些诊断不改变主发现家族。

功效固定真实P/X及缺失子集，以预定参数生成合成Y。不同蛋白的功效差异在这里来自设计几何、人数和簇，不代表已估计每个蛋白真实测量可靠性。BY首发现阈值只是单蛋白参照，不能与完整BY程序功效混淆。所有模型和家族在实际蛋白定量连接前固定；后加入的协方差相容性/影响诊断明确标为解释扩展，未重选主模型。
'''
 (ROOT/'METHOD_IMPLEMENTATION_ZH.md').write_text(methods,encoding='utf-8')
 english='''# Methods and results wording for review

We tested two observable relative linear-mean restrictions across three CSF assays, using the same nuisance design for each assay. The primary predictor was standardized log2 Elecsys CSF pTau181 from the same collection date. Models included age, sex, NULISA plate, SomaScan plate, TMT plex, and Elecsys assay batch. A prespecified quality screen preceded examination of protein–pathology effects. MAPT was retained separately as a construct-overlap control. The primary discovery family comprised 34 proteins.

The restrictions were h1=bS*cNT−bN*cST and h2=bT*cNS−bN*cST, where b denotes the adjusted pathology coefficient and c the residual cross-covariance. Their joint influence function included uncertainty in both quantities. A sum of three one-way CR1 covariance estimates was selected using synthetic development simulations and checked with independent simulations on the actual design. This estimator is neither the inclusion–exclusion CGM estimator nor a CR2 procedure. Benjamini–Yekutieli adjustment was applied within the frozen primary family; its interpretation remains conditional on valid marginal p-values.

PDGFRB, NPTX1, and KLK6 met the primary BY criterion (q approximately 0.029). All three failed the same BY criterion after additional adjustment for education and APOE ε4. No discoveries were obtained in the separately defined amyloid-PET family. The residual covariance pattern for KLK6 was incompatible at the point-estimate level with a simple single-factor model with mutually uncorrelated errors; PDGFRB exhibited a slightly negative estimated unique-variance component under that model. These findings and the covariance-offset sensitivity prevent attribution to uniquely identified pathology-specific measurement bias. We report conditional relative-pattern signals, not a validated separation of true protein abundance and assay-specific effects.

Power calculations used prespecified synthetic effect and noise parameters on the real covariate, missingness, and cluster design; they were not calculated from observed effects. No individual latent consensus/disagreement scores or cognitive slopes were estimated.
'''
 (ROOT/'METHODS_RESULTS_WORDING_EN.md').write_text(english,encoding='utf-8')
 readme='''# Stage5交付入口

先读[中文主报告](STAGE5_MAIN_REPORT_ZH.md)，查阅[Excel总表](真实设计模型校准_逐蛋白功效与主分析总表.xlsx)及[方法实现](METHOD_IMPLEMENTATION_ZH.md)。英文方法/结果措辞单独提供，避免将条件性结果写成已识别的测量机制。

本轮完成：真实共同设计、盲于蛋白效应的筛选和SAP冻结、开发/独立确认、逐蛋白条件功效、pTau主分析、五项固定敏感性、淀粉样次级分析、协方差/影响解释诊断及独立复核。认知表示/斜率和多年后tau PET未纳入此轮主问题。

results为汇总与模拟记录；00_restricted含个体连接、临床设计、变换蛋白及模型残差，仅作本地复现，不是公开附件。figures为六组PNG/PDF。根目录SHA清单排除自身和缓存。原始数据及Stage1–4保持封存。

运行code/run_all.ps1可按冻结设计复算主结果、功效汇总和QA；不会重新选择家族。模拟引擎的DEV/CONFIRM/POWER、run_variants.py及stress_analysis.py可完整生成相应重复，随机种子和64次分块均已固定。完整模拟若重跑，请在本Stage5目录副本中执行并保留原封存件，因为并行单元完成顺序可改变汇总CSV排序和散列，数值应按单元键比较；不要重新宣称这是第一次盲法分析。冻结记录不因复现重置。

依赖版本见results/SOFTWARE_ENVIRONMENT.json。8人出生冲突与TMT3447/3907解释仍是已披露来源限制，不代表本轮规定分析没有完成。
'''
 (ROOT/'README_研究结果入口.md').write_text(readme,encoding='utf-8')

def workbook():
 tablelist=[('主结果34','PRIMARY_FROZEN_FAMILY_RESULTS'),('43候选含排除','PRIMARY_ALL43_RESULTS'),('70候选盲态筛选','PRE_EFFECT_PROTEIN_SCREEN'),
 ('3发现证据矩阵','DISCOVERY_EVIDENCE_MATRIX'),('模型比较汇总','ANALYSIS_COMPARISON_SUMMARY'),('固定敏感性全部','ALL_FIXED_VARIANT_RESULTS'),('淀粉样次级','AMYLOID_SECONDARY_RESULTS'),
 ('开发方差选择','DEVELOPMENT_VARIANCE_SELECTION'),('独立确认决定','INDEPENDENT_CONFIRMATION_DECISIONS'),('开发模拟全表','DEV_SIMULATION_SUMMARY'),('主设计确认全表','CONFIRM_SIMULATION_SUMMARY'),
 ('支线设计','VARIANT_DESIGN_SCREEN'),('支线校准决定','VARIANT_CALIBRATION_DECISIONS'),('支线确认全表','VARIANT_CALIBRATION_SUMMARY'),
 ('70蛋白功效状态','POWER_EVALUABILITY_ALL70'),('逐蛋白条件功效','PER_PROTEIN_CONDITIONAL_POWER'),('功效效应汇总','POWER_SUMMARY_BY_EFFECT'),('80功效网格','PER_PROTEIN_DELTA80_GRID'),
 ('识别压力全部','IDENTIFICATION_STRESS_RESULTS'),('压力汇总','STRESS_SCENARIO_SUMMARY'),('工作协方差125点','COVARIANCE_OFFSET_GRID125'),('协方差敏感性汇总','COVARIANCE_OFFSET_SENSITIVITY_SUMMARY'),
 ('共同协方差解释','COMMON_COVARIANCE_COMPATIBILITY_DIAGNOSTICS'),('最大影响簇删除','TOP_INFLUENCE_CLUSTER_OMISSION_DIAGNOSTICS'),('GIS固定效应等价','GIS_FIXED_EFFECT_EQUIVALENCE_QA'),
 ('输入散列核对','SOURCE_HASH_CHECKS'),('重建掩码QA','PROTEIN_RECONSTRUCTION_QA'),('公式预检','PREFLIGHT_QA'),('梯度预检','DIAGNOSTIC_PREFLIGHT_QA'),
 ('独立主结果重算','INDEPENDENT_MAIN_NUMERIC_VALIDATION'),('模拟逐文件重数','INDEPENDENT_REPLICATE_RECOUNTS'),('独立QA','INDEPENDENT_QA_CHECKS'),('指标解释','METRIC_DICTIONARY_ZH'),('方法原始依据','METHODS_EVIDENCE_REGISTRY')]
 path=ROOT/'真实设计模型校准_逐蛋白功效与主分析总表.xlsx';counts=[]
 with pd.ExcelWriter(path,engine='openpyxl') as writer:
  for title,file in tablelist:
   df=pd.read_csv(RESULTS/(file+'.csv'));df.to_excel(writer,sheet_name=title,index=False);counts.append({'sheet':title,'rows':len(df),'columns':len(df.columns),'source':file+'.csv'})
 book=load_workbook(path)
 for ws in book:
  ws.freeze_panes='C2';ws.auto_filter.ref=ws.dimensions;ws.row_dimensions[1].height=55
  for cell in ws[1]:cell.font=Font(bold=True,color='FFFFFF');cell.fill=PatternFill('solid',fgColor='285679');cell.alignment=Alignment(wrap_text=True,vertical='top')
  for col in ws.columns:
   vals=list(col);ws.column_dimensions[vals[0].column_letter].width=min(45,max(13,max(len(str(c.value or '')) for c in vals[:45])+2))
  for row in ws.iter_rows(min_row=2):
   for cell in row:cell.alignment=Alignment(vertical='top',wrap_text=ws.title in ['指标解释','方法原始依据'])
  if ws.title in ['主结果34','3发现证据矩阵']:
   for c in range(1,ws.max_column+1):
    if ws.cell(1,c).value=='q_BY':
     letter=ws.cell(1,c).column_letter;ws.conditional_formatting.add(f'{letter}2:{letter}{ws.max_row}',ColorScaleRule(start_type='num',start_value=0,start_color='E8AC8F',end_type='num',end_value=1,end_color='FFFFFF'))
 book.save(path);csvout(RESULTS/'WORKBOOK_SHEET_COUNTS.csv',counts)

if __name__=='__main__':
 tables()
 if '--tables-only' not in sys.argv:figures();reports();workbook();log('Chinese report, English bounded wording, six figure pairs, and workbook generated')
