from __future__ import annotations
import json
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from openpyxl import Workbook
from openpyxl.styles import PatternFill,Font,Alignment
from openpyxl.utils import get_column_letter
from cal_config import *

LABELS={'IID_CHI2':'IID chi2','CR1_CHI2':'CR1 chi2','CR1_HT_F':'HT-type F','WILD_SCORE_T':'Wild score',
        'SUM3_CHI2':'PSD sum chi2 *','SUM3_MIN_F':'PSD sum min-F','CGM3_CHI2':'CGM3 chi2','PIGEONHOLE_CENTERED':'Product bootstrap'}
def read(name):return pd.read_csv(RESULTS/name)
def pct(x):return f'{x:.2%}'
def md(frame):
 def cell(x):return str(x).replace('|','/').replace('\n',' ')
 return '\n'.join(['| '+' | '.join(frame.columns)+' |','| '+' | '.join(['---']*len(frame.columns))+' |']+
     ['| '+' | '.join(map(cell,r))+' |' for r in frame.itertuples(index=False,name=None)])
def save_figure(fig,stem):
 for ext in ['png','pdf']:fig.savefig(FIGURES/(stem+'.'+ext),dpi=190,bbox_inches='tight')
 plt.close(fig)

def figures():
 confirmation=read('CONFIRMATION_METHOD_COMPARISON.csv')
 fig,axes=plt.subplots(1,2,figsize=(15,8.4))
 for ax,dimension in zip(axes,[8,2]):
  q=confirmation.loc[confirmation.dimension_df.eq(dimension)&confirmation.scenario_id.str.startswith('C')]
  rates=q.pivot(index='scenario_id',columns='method',values='rejection_rate').reindex(columns=METHODS)
  fail=q.pivot(index='scenario_id',columns='method',values='failure_rate').reindex(columns=METHODS)
  matrix=rates.to_numpy().copy();matrix[fail.to_numpy()>.01]=np.nan
  cmap=plt.colormaps['YlOrRd'].copy();cmap.set_bad('#d6d9dc')
  im=ax.imshow(matrix,vmin=0,vmax=.30,cmap=cmap,aspect='auto')
  ax.set_xticks(range(8),[LABELS[m] for m in METHODS],rotation=40,ha='left',fontsize=8);ax.xaxis.tick_top()
  ax.set_yticks(range(10),rates.index);ax.set_title(f'{dimension}-df restrictions | independent confirmation',pad=90,fontsize=11)
  for i in range(10):
   for j in range(8):
    text=f'F{fail.iloc[i,j]:.0%}' if fail.iloc[i,j]>.01 else f'{rates.iloc[i,j]:.1%}'
    ax.text(j,i,text,ha='center',va='center',fontsize=8,color='white' if np.isfinite(matrix[i,j]) and matrix[i,j]>.20 else 'black')
 fig.subplots_adjust(top=.73,bottom=.12,left=.06,right=.94,wspace=.16)
 cax=fig.add_axes([.955,.2,.012,.48]);fig.colorbar(im,cax=cax,label='Null rejection rate (color capped at 30%)')
 fig.text(.06,.025,'R=2,000 per cell; all planned replicates in denominator. Grey = >1% failure. F denotes failure rate.\n* Selected in development before confirmation. Low rejection alone does not establish adequate power or true molecular identity.',fontsize=9)
 save_figure(fig,'01_INDEPENDENT_NULL_CALIBRATION')
 fig,ax=plt.subplots(figsize=(8.5,5.8));x=np.arange(4)
 for dimension,offset,color in [(8,-.18,'#d78045'),(2,.18,'#286d97')]:
  q=confirmation.loc[confirmation.dimension_df.eq(dimension)&confirmation.scenario_id.str.startswith('A')&confirmation.method.eq('SUM3_CHI2')].sort_values('scenario_id')
  y=q.rejection_rate.to_numpy();err=np.maximum(0,np.vstack([y-q.wilson95_lower.to_numpy(),q.wilson95_upper.to_numpy()-y]))
  ax.bar(x+offset,y,.34,color=color,label=f'{dimension}-df / PSD sum chi2',yerr=err,capsize=3)
  for xx,yy in zip(x+offset,y):ax.text(xx,yy+.025,f'{yy:.1%}',ha='center',fontsize=9)
 ax.axhline(.8,color='.5',ls='--',lw=1,label='0.80 reference')
 ax.set(xticks=x,xticklabels=['A01: M2 offset .35','A02: M3 offset .50','A03: M3 loading .40','A04: offset .50 / crossed'],ylabel='Rejection proportion',ylim=(0,1.10),title='Independent power checks | R=2,000 per cell')
 ax.tick_params(axis='x',labelsize=8);ax.legend(frameon=False,fontsize=9);ax.grid(axis='y',alpha=.15)
 fig.text(.10,.025,'Artificial effects and noise scales; not actual ADNI protein power. The 2-df test does not cover all loading changes.',fontsize=8.5)
 fig.subplots_adjust(bottom=.18,top=.87);save_figure(fig,'02_CONTROL_POWER_TRADEOFF')
 profile=read('ASSUMPTION_GRID_ENVELOPE_SUMMARY.csv');rho=np.array([-.4,-.2,0,.2,.4])
 fig,axes=plt.subplots(1,2,figsize=(10.8,4.8))
 for dimension,color in [(8,'#c87840'),(2,'#286d97')]:
  q=profile.loc[profile.dimension_df.eq(dimension)&profile.method.eq('SUM3_CHI2')&profile.scenario_id.isin([f'R{i}' for i in range(5)])].sort_values('scenario_id')
  axes[0].plot(rho,q.zero_assumption_rejection,'o-',color=color,label=f'{dimension}-df, assume kappa=0')
  axes[0].plot(rho,q.grid_envelope_rejection,'s--',color=color,label=f'{dimension}-df, discrete max-P')
  axes[1].plot(rho,q.classification_flip_rate,'o-',color=color,label=f'{dimension}-df')
 axes[0].set(ylabel='Rejection proportion',title='No extra relative effect: fixed correlated errors')
 axes[1].set(ylabel='Fraction with any grid-dependent classification',title='Sensitivity is not mechanism identification')
 for ax in axes:ax.set_xlabel('Generating M1/M2 error correlation');ax.grid(alpha=.15);ax.legend(frameon=False,fontsize=8)
 fig.text(.07,.025,'Grid: kappa12 = [-0.4, -0.3, ..., 0.4] x 0.48 synthetic units. Discrete envelope only; not a continuous bound.\nA flip can start from a negative or a positive result; it is not the fraction of significant findings that fail replication.',fontsize=8.5)
 fig.subplots_adjust(bottom=.25,top=.85,wspace=.3);save_figure(fig,'03_ERROR_CORRELATION_SENSITIVITY')
 q=read('BOOTSTRAP_RESOLUTION_SUMMARY.csv');fig,axes=plt.subplots(1,2,figsize=(10.5,4.8))
 for ax,method in zip(axes,['WILD_SCORE_T','PIGEONHOLE_CENTERED']):
  a=q.loc[q.method.eq(method)].set_index('scenario_id').reindex(['C02','C03','C05','A01']);x=np.arange(4)
  ax.bar(x-.17,a.original_rejection,.32,label='B499');ax.bar(x+.17,a.high_B_rejection,.32,label='B1999')
  ax.set(xticks=x,xticklabels=a.index,ylim=(0,1.05),title=LABELS[method],ylabel='Rejection fraction in same 500 datasets')
  ax.legend(frameon=False);ax.grid(axis='y',alpha=.15)
  if max(a.original_rejection.max(),a.high_B_rejection.max())==0:
   ax.text(.5,.5,'0 / 500 rejections at both B levels\nin all four cells',transform=ax.transAxes,ha='center',color='.4')
 fig.text(.08,.025,'Same confirmation datasets replayed, not new independent replicates. More bootstrap draws do not repair structural assumptions.',fontsize=8.5)
 fig.subplots_adjust(bottom=.20,top=.85,wspace=.25);save_figure(fig,'04_BOOTSTRAP_RESOLUTION_REPLAY')

def method_registry():
 details=[
 ('IID_CHI2','个体独立chi-square','基线，忽略共享批次依赖'),('CR1_CHI2','单向CR1 chi-square','旧基线'),
 ('CR1_HT_F','Hotelling型有限簇F','候选近似，非完整CR2/HTZ，部分少簇情景未通过'),
 ('WILD_SCORE_T','中心化wild簇score-t','线性化零均值score，不是完整受限潜变量bootstrap'),
 ('SUM3_CHI2','三单向CR1之和＋chi-square','正半定保守候选；不作交集相减；仅在已验证条件下推荐'),
 ('SUM3_MIN_F','正半定和＋最小簇F','很保守，需单独看功效'),('CGM3_CHI2','三向容斥chi-square','非正定/奇异明确保留为失败，不作PSD修补'),
 ('PIGEONHOLE_CENTERED','乘积权重全矩bootstrap','完整重算矩或加权OLS，约束居中，交叉稀疏结构下可能很保守')]
 d=pd.DataFrame(details,columns=['method','中文名称','限制']);csv(RESULTS/'METHOD_DICTIONARY_ZH.csv',d);return d

def dictionary():
 rows=[('dimension_df','2为线性相对均值限制，8为组均值＋协方差限制；不是同一估计对象'),
 ('rejection_rate','拒绝数/全部计划重复，失败仍在分母'),('failure_rate','无法产生有效P的比例，不能将其当阴性或低功效'),
 ('wilson95_lower/upper','Monte Carlo比例的Wilson95区间，不是临床效应区间'),('n_replicates','每格人工重复数；重放不增加独立重复总数'),
 ('confirmation_null_gate_pass','冻结候选在全部声明范围零情景的Wilson上限≤7.5%，失败≤1%'),
 ('strong_power_gate_pass','预设A04强效应功效Wilson下限≥80%；不能外推真实效应或其他候选蛋白'),
 ('grid_envelope_rejection','离散kappa网格上max-P≤0.05的比例；不是连续区间严格上界'),
 ('classification_flip_rate','同一人工数据在网格内既有P≤0.05又有P>0.05，不限于原来阳性的样本'),
 ('known_true_fixed_kappa_rejection','仅模拟可知且真固定kappa恰在网格时报告；不是现实可用校正'),
 ('oracle_remove_quadratic_rejection','人工已知T移除q*T²后的CR1_HT_F比例；不适用于现实未知T'),
 ('bootstrap_B','开发/敏感性199，确认499，分辨率重放1999；bootstrap有plus-one修正'),
 ('seed','文本保存的64位种子；phase/replicate/数据与bootstrap子流相互独立'),
 ('C01–C10','兼容零情景；C03/04/05/07/08含簇级P分配且与误差独立'),
 ('A01–A04','人工相对偏移/载荷效应；2维均值检验不能覆盖所有载荷改变'),
 ('M1/M2/M3','人工平台，噪声组合不代表真实检测平台精度排序')]
 d=pd.DataFrame(rows,columns=['字段或概念','解释']);csv(RESULTS/'METRIC_DICTIONARY_ZH.csv',d);return d

def workbook():
 qa=json.loads((RESULTS/'COMPLETED_VALIDATION_QA.json').read_text(encoding='utf-8'))
 checks=pd.DataFrame(qa['checks']);checks['detail']=checks.detail.map(lambda x:json.dumps(x,ensure_ascii=False) if isinstance(x,(dict,list)) else str(x))
 notes=pd.DataFrame([
 ('结果定位','2df相对均值＋SUM3_CHI2通过本轮条件确认；不是正式SAP授权或真实丰度分解。'),
 ('保守性','确认零情景最高拒绝1.6%，显著低于5%；强A04功效92.35%，中等A01为75.85%。'),
 ('8df原方案','控制较好但强A04功效23.5%；原结果保留，低维检验另有补充方案和新种子。'),
 ('样本','116,000个独立人工重复；24,000敏感性重放及2,000分辨率重放不另算。'),
 ('行数','全部推断结果含重放1,136,000行，Parquet保存；Excel为汇总。'),
 ('识别限制','相关误差、非线性及选择性缺失仍需假设；校准不清除未处理技术混杂。'),
 ('F与wild','此处为非线性约束的候选改编，不是文献中完整CR2或完整受限模型wild bootstrap。'),
 ('真实数据边界','仅使用聚合人数/批次规模及冻结代码，未读个体蛋白、临床/病理或restricted。'),
 ('正式多重检验','B499路径最小P0.002，不能直接满足43蛋白BH首门槛；SUM3_CHI2为解析P。'),
 ],columns=['项目','说明'])
 sheets=[('阅读说明',notes),('冻结候选独立确认',read('INDEPENDENT_CONFIRMATION_GATES.csv')),
 ('两维八维确认对照',read('CONFIRMATION_METHOD_COMPARISON.csv')),('原8df全部实验',read('ALL_EXPERIMENT_SUMMARY.csv')),
 ('补充2df开发与确认',read('MEAN_FOLLOWUP_SUMMARY.csv')),('识别敏感性',read('SENSITIVITY_METHOD_SUMMARY.csv')),
 ('相关误差网格包络',read('ASSUMPTION_GRID_ENVELOPE_SUMMARY.csv')),('Bootstrap分辨率',read('BOOTSTRAP_RESOLUTION_SUMMARY.csv')),
 ('原开发选择依据',read('DEVELOPMENT_SELECTION_SCORES.csv')),('补充开发选择依据',read('MEAN_SELECTION_SCORES.csv')),
 ('主实验种子计划',read('EXPERIMENT_GRID.csv')),('补充种子计划',read('MEAN_FOLLOWUP_GRID.csv')),('人工情景',read('SCENARIO_REGISTRY.csv')),
 ('方法定义',method_registry()),('字段解释',dictionary()),('实现与复现QA',checks),('输入哈希',read('INPUT_SOURCE_HASH_CHECKS.csv'))]
 wb=Workbook();wb.remove(wb.active);schema=[]
 for name,frame in sheets:
  if 'seed' in frame.columns:frame=frame.copy();frame['seed']=frame.seed.map(str)
  ws=wb.create_sheet(name);ws.append(list(frame.columns))
  for row in frame.itertuples(index=False,name=None):ws.append([None if pd.isna(v) else v.item() if isinstance(v,np.generic) else v for v in row])
  ws.freeze_panes='C2';ws.auto_filter.ref=ws.dimensions;ws.sheet_view.zoomScale=85;ws.row_dimensions[1].height=42
  for c in ws[1]:c.font=Font(bold=True,color='FFFFFF');c.fill=PatternFill('solid',fgColor='24556E');c.alignment=Alignment(wrap_text=True,vertical='center')
  for j,col in enumerate(frame.columns,1):
   width=19 if pd.api.types.is_numeric_dtype(frame[col]) else min(70,max(24,len(str(col))))
   if name in ['阅读说明','字段解释']:width=32 if j==1 else 112
   ws.column_dimensions[get_column_letter(j)].width=width
   if pd.api.types.is_float_dtype(frame[col]):
    fmt='0.00%' if any(k in str(col) for k in ['rate','power','rejection','upper95','lower95','wilson','fraction']) else '0.0000'
    for cells in ws.iter_cols(min_col=j,max_col=j,min_row=2):
     for c in cells:c.number_format=fmt
  if name in ['阅读说明','字段解释']:
   for row in ws.iter_rows(min_row=2):
    ws.row_dimensions[row[0].row].height=42
    for c in row:c.alignment=Alignment(wrap_text=True,vertical='top')
  schema.append({'sheet':name,'rows':len(frame),'columns':len(frame.columns)})
 wb.properties.title='有限簇校准与识别假设敏感性验证';wb.save(ROOT/'有限簇校准与识别敏感性总表.xlsx');csv(RESULTS/'WORKBOOK_SCHEMA.csv',schema)

def reports():
 confirmation=read('CONFIRMATION_METHOD_COMPARISON.csv');gates=read('INDEPENDENT_CONFIRMATION_GATES.csv');profile=read('ASSUMPTION_GRID_ENVELOPE_SUMMARY.csv')
 rows=[]
 for sid in [f'C{i:02d}' for i in range(1,11)]+['A01','A02','A03','A04']:
  row={'情景':sid,'说明':next(s['label_zh'] for s in calibration_scenarios() if s['scenario_id']==sid)}
  for d,label in [(8,'8df保守联合检验'),(2,'2df保守均值检验')]:
   q=confirmation.loc[confirmation.dimension_df.eq(d)&confirmation.scenario_id.eq(sid)&confirmation.method.eq('SUM3_CHI2')].iloc[0]
   row[label]=f'{q.rejection_rate:.2%} [{q.wilson95_lower:.2%}, {q.wilson95_upper:.2%}]'
  rows.append(row)
 table=pd.DataFrame(rows);csv(RESULTS/'CONFIRMATION_KEY_RESULTS_ZH.csv',table)
 sensitivity=read('SENSITIVITY_METHOD_SUMMARY.csv');srows=[]
 for sid in ['R0','R2','R4','Q3','SAT','SEL2','RAW','GIS','OFF2','OFF3','RHOV']:
  q=sensitivity.loc[sensitivity.dimension_df.eq(2)&sensitivity.scenario_id.eq(sid)&sensitivity.method.eq('SUM3_CHI2')].iloc[0]
  p=profile.loc[profile.dimension_df.eq(2)&profile.scenario_id.eq(sid)&profile.method.eq('SUM3_CHI2')].iloc[0]
  srows.append({'情景':sid,'2df保守拒绝率':pct(q.rejection_rate),'离散误差网格包络拒绝率':pct(p.grid_envelope_rejection),'网格内存在分类翻转':pct(p.classification_flip_rate)})
 s_table=pd.DataFrame(srows);csv(RESULTS/'SENSITIVITY_KEY_RESULTS_ZH.csv',s_table)
 qa=json.loads((RESULTS/'COMPLETED_VALIDATION_QA.json').read_text(encoding='utf-8'))
 report=f'''# 有限簇推断校准与识别假设敏感性验证报告

日期：2026-10-27（Asia/Shanghai）。状态：**计算、独立确认及实现核验完成；获得一个范围收窄的条件候选，尚未冻结正式SAP。**

最可行的后续主分析候选是**2维线性相对均值限制＋三个单向CR1方差之和＋chi-square（SUM3_CHI2）**。它在全部预设确认零情景的最高拒绝率为1.60%（Wilson上限2.25%），无方法失败；强A04人工效应检出率92.35%（95%区间91.10%–93.44%）。这是偏保守的控制，不是精准校准到5%，也不是实际ADNI功效。

原8df保守联合检验的强A04检出率仅23.50%，控制通过但功效不足。有限簇F/wild score没有在全部不利情景通过；不能择取其中表现好的格用于推荐。

## 1. 设计与输入边界

主方案[CALIBRATION_PROTOCOL_v1.0.md](CALIBRATION_PROTOCOL_v1.0.md)包含74,000次独立人工重复。原8df开发暴露控制与功效取舍后，另写[低维补充方案](LOW_DIMENSIONAL_FOLLOWUP_PROTOCOL_v1.0.md)，使用新开发种子和新确认种子增加42,000次。合计**116,000次独立人工重复**。

24,000次2df敏感性数据为原敏感性种子重放；B1999实验重放确认集2,000次，均未加算独立次数。含重放共有1,136,000行推断结果，保存在Parquet。比较8条推断路径，全部失败留存。

仅使用聚合n=571、19/36/65批次规模及前阶段冻结代码；个体、P、T、误差和到批次的对应均重新生成。未读个体蛋白、临床/病理或restricted文件。误差参数不代表真实平台排名。原Stage1/2文件与原SAP保持不动。

开发→冻结选择→独立确认的顺序在代码中强制执行，确认种子与开发、前阶段及其他格不混用。本轮还让P在部分情景按簇分配且与误差独立，避免仅用大量簇内分组差异掩盖少簇风险；这不是实际病理分组结构。

## 2. 推断路径与估计对象

- 8df：组均值与非对角协方差联合稳定性，允许共同均值/方差和平台误差方差变化。
- 2df：P对三平台的线性斜率是否与残差协方差识别的共同方向平行，包含方向估计不确定性。它不检验全部非线性病理关系或全部载荷变化。
- CR1_HT_F是Hotelling型候选近似，**不是**已经为此非线性模型证明有效的完整CR2/HTZ。
- WILD_SCORE_T是中心化簇影响函数的studentized符号bootstrap，**不是**完整受限潜变量模型的全流程bootstrap。
- Pigeonhole路径独立抽各向簇权重后相乘，完整重算组矩或加权OLS，再对约束居中；不是只扰动P值，也不是完整潜变量参数零模型bootstrap。
- SUM3_CHI2将三个单向CR1相加，不作交集相减，因此正半定且可能很保守。只按本轮适用范围判断，不声称任意多向依赖下均有保证。
- 旧CGM3容斥路径保留为失败基线，没有特征值修补、静默删失败或伪逆。

开发规则：范围内每格零情景Wilson上限≤8%、失败≤1%，再按A01/A02平均检出率排序。两类估计对象均在查看各自确认数据前选定SUM3_CHI2。确认规则上限≤7.5%、失败≤1%；强A04功效区间下限≥80%。这些规则不能事后跳过19簇或不平衡情景。

## 3. 独立确认结果

每格R=2,000，方括号为Monte Carlo Wilson95区间。比例分母包括失败，并非仅成功样本。

{md(table)}

2df候选通过ONE_WAY与CROSSED两个声明范围的控制及强效应门槛，但中等A01仅75.85%，载荷效应A03为0；**不能把它说成覆盖所有不一致模式的万能检验**。8df保守方法不能因零情景很少拒绝就称充分有效。

其他方法的确认结果、失败和完整开发依据均在[全部确认对照](results/CONFIRMATION_METHOD_COMPARISON.csv)及[冻结确认门槛](results/INDEPENDENT_CONFIRMATION_GATES.csv)。没有在确认失败后改选另一方法声称盲验证成功。

## 4. 识别假设敏感性

24个新种子情景包含固定/变化误差相关、固定二次及饱和响应、MCAR/LOD/地板替代/选择性缺失、批次混杂、GIS后共享噪声、弱信号和不同平台偏移。下表是已确认的2df保守方法；完整2df/8df各路径另外保留。

{md(s_table)}

最明确的警示是RAW：没有额外病理直接项，只有未处理技术混杂，2df保守检验仍拒绝99.90%；离散误差包络后仍为87.20%。**控制随机簇依赖不等于去除了系统技术混杂。**GIS压力情景的低拒绝也不能被外推为所有蛋白都应使用GIS。

相关误差网格为kappa12=[−0.4,−0.3,...,0.4]×0.48人工单位。max-P只对该离散网格求包络；不是连续区间的严格上确界，也不是实际检测误差已校准的范围。固定真kappa恰在网格、且单点推断有效时才有相应条件含义；变化的相关误差和其他模型错误不能自动被覆盖。

“网格内存在分类翻转”可以从原本阴性开始，不等于“原本阳性结果失效的比例”。例如R0有27.70%的数据在某些网格点改变分类，但假设0时拒绝仅1%。图表与原表明确保留这个分母区别。

固定二次响应在对称P分布下可能不改变2df所针对的线性斜率，因此低拒绝不是证明任意非线性均可忽略。人工已知T移除q*T²仅是信息基准；不能用有误差的真实读数直接逆变换冒充同样校正。弱信号、低完整人数和低功效仍不能被P>0.05认证为稳定关系。

前阶段的观测等价证明仍成立：不同潜在信号、误差相关和平台直接项可产生同一观测分布。校准成功不会新增识别信息，不能据此声称分离真丰度或proteoform机制。

## 5. Bootstrap分辨率检查

开发/敏感性B199、确认B499，plus-one修正，统一P≤0.05；重放四个确认情景各500个原有数据集，将B提高到1999。详见[分辨率汇总](results/BOOTSTRAP_RESOLUTION_SUMMARY.csv)。提高B不自动修复检验对象或非线性/相关误差假设。

B499路径最小P=0.002，高于43蛋白BH首门槛0.05/43；这些bootstrap路径若用于正式多重检验，还需更高B或事前规定的序贯MC。当前条件候选SUM3_CHI2为解析P，不受这个数值分辨率限制，但仍需后续多蛋白验证。

## 6. 对论文的建议

**建议主推断先收窄为“独立病理指标相关的跨平台线性相对均值偏离”**，将2df保守方法列为候选。8df关系重塑可作为方法开发/描述或预设支持分析；本轮不支持把低功效版本直接当主要发现工具。

正式实施前仍须完成：真实技术处理规则、独立病理及随访覆盖、按实际每蛋白缺失结构的功效、多重检验和整体流程bootstrap/敏感性。低检出、信号弱、误差假设敏感或技术混杂明显的蛋白保留为不可有力评价/解释不确定。

这项主张应保持为测量关系的条件性证据。若要归因于分子形式、裂解或表位，需要匹配肽段、复测或其他独立识别信息；不能因为本轮有限簇控制通过就升级机制结论。

## 7. 实现验证与交付

实现/冻结/计划/种子核验{qa['passed']}/{qa['total']}通过，包括wild快速公式对显式studentization、全权重矩/加权OLS对逐组计算、plus-one、Wilson边界、全部逐格人数及源哈希。核验通过不代表所有方法通过校准。

补充代码一次验证接口修改被冻结检查拦截，确认文件当时为0；重放14个开发文件后全部SHA相同，重新冻结才运行确认。独立Wilson复核器的0/全拒绝端点求根错误已修正，生产P值和区间未更改。详见[修正记录](logs/FAILURES_AND_CORRECTIONS.md)。

- [Excel总表](有限簇校准与识别敏感性总表.xlsx)、[下一步方法决定](NEXT_METHOD_DECISION_MEMO.md)
- [校准图](figures/01_INDEPENDENT_NULL_CALIBRATION.png)、[控制与功效](figures/02_CONTROL_POWER_TRADEOFF.png)、[误差相关敏感性](figures/03_ERROR_CORRELATION_SENSITIVITY.png)、[分辨率重放](figures/04_BOOTSTRAP_RESOLUTION_REPLAY.png)，均附PDF。
- [访问与冻结记录](results/ACCESS_AND_PROTOCOL_RECORD.json)、[完成验证](results/COMPLETED_VALIDATION_QA.json)、[交付核验](results/DELIVERY_QA.json)、[输入输出SHA](CALIBRATION_INPUT_OUTPUT_SHA256.csv)。
- `results/replicates`保存逐重复Parquet；`code/run_all.ps1`重建流程。日志/Excel/PDF时间和哈希会随重跑变化，统计数值应复现。

方法背景：[Cameron等2008](https://doi.org/10.1162/rest.90.3.414)、[Pustejovsky与Tipton](https://arxiv.org/abs/1601.01981)、[Owen2007](https://doi.org/10.1214/07-AOAS122)、[Owen与Eckles](https://arxiv.org/abs/1106.2125)。本轮具体非线性限制及其结果是本课题改编的数值验证，不是称这些论文已证明其在ADNI的适用性。
'''
 (ROOT/'FINITE_CLUSTER_AND_ASSUMPTION_REPORT_ZH.md').write_text(report,encoding='utf-8')
 memo='''# 下一步方法决定

本轮计算与确认已完成；正式SAP和蛋白集尚未冻结。

| 对象 | 判断 | 用途与限制 |
| --- | --- | --- |
| 2df相对均值＋SUM3_CHI2 | 条件候选，偏保守 | 通过本轮范围内控制与强效应门槛；只覆盖线性相对均值，不能证明真丰度或全部载荷变化 |
| 8df联合限制＋SUM3_CHI2 | 控制通过但功效不足 | A04仅23.5%；不建议本版本作为主要发现工具 |
| 本轮F/wild score改编 | 全范围确认不足 | 不按有利情景挑选使用，也不冒充完整CR2或受限模型bootstrap |
| Product bootstrap/最小簇F | 过于保守 | 保留完整结果，不因少拒绝就认定充分有效 |
| 原CGM3 | 失败保留 | 不用修补矩阵隐藏不稳定性 |
| 相关误差网格 | 敏感性描述 | 离散范围、非现实校准值，不是连续识别集合或机制证明 |

下一阶段先确定真实技术处理与独立病理覆盖，再按实际每蛋白结构验证2df候选的功效和多重检验。正式主问题可收窄为病理相关的跨平台线性相对均值偏离；若继续主张全面分离共同丰度与病理特异测量机制，仍需额外识别信息。
'''
 (ROOT/'NEXT_METHOD_DECISION_MEMO.md').write_text(memo,encoding='utf-8')
 (ROOT/'README_阅读顺序.md').write_text('''# 有限簇推断校准＋识别假设敏感性验证

完成116,000次独立人工重复、独立确认与敏感性重放。2df相对均值＋保守正半定方差是条件候选；8df原方案功效不足。没有读取个体蛋白或临床/病理文件。

1. [中文结果报告](FINITE_CLUSTER_AND_ASSUMPTION_REPORT_ZH.md)
2. [Excel总表](有限簇校准与识别敏感性总表.xlsx)
3. [方法决定](NEXT_METHOD_DECISION_MEMO.md)
4. [主方案](CALIBRATION_PROTOCOL_v1.0.md)与[低维补充方案](LOW_DIMENSIONAL_FOLLOWUP_PROTOCOL_v1.0.md)
5. [交付QA](results/DELIVERY_QA.json)与[哈希清单](CALIBRATION_INPUT_OUTPUT_SHA256.csv)

正式SAP尚未冻结。校准通过不等于已去除技术混杂或已识别真实丰度/分子机制。
''',encoding='utf-8')

def main():
 figures();workbook();reports();log('Reports, workbook and four PNG/PDF pairs built')
if __name__=='__main__':main()
