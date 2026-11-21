from config import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
plt.rcParams.update({'font.family':'Microsoft YaHei','font.size':10,'axes.unicode_minus':False,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
COL={'N':'#2979b9','S':'#dd7c21','T':'#329b70'}
LABEL={'N':'NULISA','S':'SomaScan','T':'TMT'}
GG=list(GENES.values());inventory=[]

def save(fig,name,description):
 fig.savefig(FIGURES/(name+'.png'),bbox_inches='tight');fig.savefig(FIGURES/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
 for ext in ['png','pdf']:inventory.append({'figure_id':name,'file':name+'.'+ext,'description':description,'sha256':sha(FIGURES/(name+'.'+ext))})

def run():
 e=pd.read_csv(RESULTS/'PLATFORM_EFFECTS_ALL27.csv');r=pd.read_csv(RESULTS/'ALL42_DIRECTIONAL_EXPLANATION_RESULTS.csv');h=pd.read_csv(RESULTS/'PAIRED_CONTRAST_ESTIMATES_ALL42.csv')
 # Pointwise intervals; the paired test, not a visual overlap rule, addresses change.
 fig,axs=plt.subplots(1,3,figsize=(14,4.8),sharey=True)
 for ax,g in zip(axs,GG):
  for j,a in enumerate(COL):
   d=e.loc[e.gene.eq(g)&e.assay.eq(a)].set_index('model').loc[['NO_P','LOG_P','RAW_P','RCS_P']];x=np.arange(4)+(j-1)*.11
   ax.errorbar(x,d.gamma,yerr=1.959964*d.SE_full_joint,fmt='o-',capsize=3,color=COL[a],label=LABEL[a],lw=1.4,markersize=4)
  ax.axhline(0,color='.6',lw=.8);ax.set_xticks(range(4),['无pTau','log pTau','原值pTau','样条pTau']);ax.set_title(g);ax.grid(axis='y',alpha=.15)
 axs[0].set_ylabel('记忆年变化 / 原冻结蛋白残差SD');axs[-1].legend(fontsize=9);fig.suptitle('01 同一样本、固定尺度：加入pTau如何改变各检测的认知关联',y=1.03)
 fig.text(.5,-.025,'点及逐项95%区间；不是同时置信区间。候选来自同一结局的前轮筛选。',ha='center',fontsize=9);fig.tight_layout();save(fig,'F01_PLATFORM_COEFFICIENTS','Same-case fixed-scale platform coefficients; pointwise intervals.')
 fig,axs=plt.subplots(1,3,figsize=(14,4.6),sharey=True)
 for ax,g in zip(axs,GG):
  for j,cname in enumerate(['S_MINUS_N','T_MINUS_N']):
   d=h.loc[h.gene.eq(g)&h.model.isin(PRIMARY_COMPARISONS)&h.contrast.eq(cname)].set_index('model').loc[PRIMARY_COMPARISONS];y=np.arange(3)+(j-.5)*.17
   ax.errorbar(d.estimate,y,xerr=1.959964*d.SE_full_joint,fmt='o',capsize=3,color=['#bb5b26','#397d72'][j],label=['Soma − NULISA','TMT − NULISA'][j])
  ax.axvline(0,color='.5',lw=.8);ax.set_yticks(range(3),['log变化','原值变化','样条变化']);ax.invert_yaxis();ax.set_title(g);ax.set_xlabel('加入pTau前后，跨检测系数差异的变化')
  rr=r.loc[r.gene.eq(g)&r.model.isin(PRIMARY_COMPARISONS)].set_index('model')
  text='；'.join(m.split('_')[0]+': q9='+format(rr.loc[m,'q_BY_family'],'.3g') for m in PRIMARY_COMPARISONS);ax.text(.02,-.24,text,transform=ax.transAxes,fontsize=8)
 handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,1.04),ncol=2,fontsize=9);fig.suptitle('02 直接配对检验：检验差异本身是否发生改变',y=1.12);fig.text(.5,-.08,'逐项95%区间；q9只整理本轮范围，未处理前轮同结局候选筛选。',ha='center',fontsize=9);fig.tight_layout();save(fig,'F02_PAIRED_CHANGE','Paired conditioning changes with full cross-model covariance; q values describe selected scope.')
 fig,ax=plt.subplots(figsize=(13,4.1));vals=np.array([[r.loc[r.gene.eq(g)&r.model.eq(m),'q_BY_family'].iloc[0] for m in MODELS] for g in GG]);im=ax.imshow(-np.log10(np.clip(vals,1e-8,1)),cmap='YlOrRd',vmin=0,vmax=3,aspect='auto')
 for i in range(3):
  for j in range(9):
   rr=r.loc[r.gene.eq(GG[i])&r.model.eq(MODELS[j])].iloc[0];ax.text(j,i,format(vals[i,j],'.3g')+(' *' if rr.supported_family else '')+(' GATE' if not rr.calibration_qualified else ''),ha='center',va='center',fontsize=9,color='white' if -np.log10(max(vals[i,j],1e-8))>1.9 else 'black')
 ax.set_xticks(range(9),MODELS,rotation=35,ha='right');ax.set_yticks(range(3),GG);ax.set_title('03 九种条件下的跨检测差异：BY33范围内q值');fig.colorbar(im,ax=ax,label='−log10(q33)');fig.text(.1,-.035,'* q33≤0.05且校准合格；发现后的定向敏感性分析，未校正前轮候选筛选。',fontsize=9);fig.tight_layout();save(fig,'F03_MODEL_DIFFERENCE_MATRIX','27 model-difference tests, q within secondary33; gate and selection boundary explicit.')
 d=pd.read_csv(RESULTS/'OVB_ASSAY_DECOMPOSITION.csv');m=pd.read_csv(RESULTS/'FWL_WEIGHTED_MOMENT_EXPLANATION.csv');fig,axs=plt.subplots(2,3,figsize=(13,7.8))
 for k,g in enumerate(GG):
  dd=d.loc[d.gene.eq(g)&d.new_model.eq('LOG_P')].set_index('assay').loc[list(COL)];axs[0,k].bar(list(LABEL.values()),dd.after_minus_before,color=list(COL.values()));axs[0,k].axhline(0,color='.5',lw=.8);axs[0,k].set_title(g);axs[0,k].set_ylabel('条件化带来的系数增量')
  mm=m.loc[m.gene.eq(g)].set_index('assay').loc[list(COL)];axs[1,k].bar(list(LABEL.values()),mm.weighted_partial_corr_X_pTau,color=list(COL.values()));axs[1,k].set_ylim(0,1);axs[1,k].set_ylabel('蛋白–pTau加权偏相关')
 fig.suptitle('04 精确代数分解：pTau关联项与蛋白–pTau耦合共同决定系数变化',y=1.015);fig.text(.5,-.018,'β无P − β有P = δ × π；多项样条时对各基函数求和。这是回归恒等式，不是因果中介比例。',ha='center',fontsize=9);fig.tight_layout();save(fig,'F04_EXACT_SUPPRESSION_DECOMPOSITION','Exact OVB identity and FWL partial correlations; not mediation.')
 curves=pd.read_csv(RESULTS/'ADJUSTED_PTAU_CURVES.csv');fig,axs=plt.subplots(2,3,figsize=(14,7.4))
 for k,g in enumerate(GG):
  for model,style in [('LOG_P','--'),('RCS_P','-')]:
   d=curves.loc[curves.gene.eq(g)&curves.model.eq(model)]
   for a in COL:axs[0,k].plot(d.pTau_pg_per_ml,d[a+'_protein_component_in_fixedSD'],style,color=COL[a],label=LABEL[a]+(' 样条' if style=='-' else ' log'),lw=1.5)
   axs[1,k].plot(d.pTau_pg_per_ml,d.memory_change_component_relative_median,style,color='#725198',label=model)
  for row in [0,1]:
   axs[row,k].set_xscale('log');axs[row,k].set_xticks([10,20,40,60],['10','20','40','60']);axs[row,k].xaxis.set_minor_locator(NullLocator());axs[row,k].axhline(0,lw=.6,color='.6');axs[row,k].set_xlabel('pTau（pg/mL；对数轴）')
  axs[0,k].set_title(g);axs[0,k].set_ylabel('蛋白分量 / 固定SD');axs[1,k].set_ylabel('记忆年变化分量');axs[1,k].legend(fontsize=8)
 axs[0,0].legend(fontsize=7,ncol=2);fig.suptitle('05 时间信息加权的pTau曲线：原协变量调整后，相对pTau中位数',y=1.015);fig.text(.5,-.018,'描述性曲线；不是条件效应、因果轨迹或个体预测。',ha='center',fontsize=9);fig.tight_layout();save(fig,'F05_ADJUSTED_PTAU_CURVES','Weighted adjusted protein and outcome profiles; linear vs restricted cubic spline.')
 info=pd.read_csv(RESULTS/'EXPOSURE_INFORMATION_AND_LEVERAGE.csv');fig,axs=plt.subplots(2,3,figsize=(14,7.5));mm=['NO_P','LOG_P','BG2_NO_P','BG2_LOG_P','BG3_NO_P','BG3_LOG_P']
 for k,g in enumerate(GG):
  for a in COL:
   dd=e.loc[e.gene.eq(g)&e.assay.eq(a)].set_index('model').loc[mm];ii=info.loc[info.gene.eq(g)&info.assay.eq(a)].set_index('model').loc[mm]
   axs[0,k].plot(range(6),dd.gamma,'o-',color=COL[a],label=LABEL[a],ms=4);axs[1,k].plot(range(6),ii.fraction_of_original_NO_P_information,'o-',color=COL[a],ms=4)
  axs[0,k].set_title(g);axs[0,k].axhline(0,color='.6',lw=.7);axs[0,k].set_ylabel('认知关联系数（固定SD）');axs[1,k].set_ylabel('蛋白残差信息 / 原无P模型')
  for row in [0,1]:axs[row,k].set_xticks(range(6),mm,rotation=50,ha='right');axs[row,k].grid(axis='y',alpha=.15)
 axs[0,0].legend(fontsize=8);fig.suptitle('06 跨蛋白背景调整：效应变化与信息损失需要同时解释',y=1.02);fig.text(.5,-.02,'BG2：Soma+TMT背景；BG3再加NULISA。代理背景可能包含生物信号，不能当作实测脑脊液总蛋白。',ha='center',fontsize=9);fig.tight_layout();save(fig,'F06_BACKGROUND_AND_INFORMATION','Background proxy sensitivity and residual exposure information.')
 toy=pd.read_csv(RESULTS/'SINGLE_TRUE_SIGNAL_SYNTHETIC_REPRODUCTION.csv');grid=pd.read_csv(RESULTS/'SINGLE_TRUE_SIGNAL_ERROR_VARIANCE_GRID.csv');fig,axs=plt.subplots(1,2,figsize=(13,4.8))
 for a in COL:
  dd=toy.loc[toy.assay.eq(a)].sort_values('pTau_adjusted');axs[0].plot([0,1],dd.analytical_standardized_coefficient,'o-',color=COL[a],label=LABEL[a]);axs[0].scatter([0,1],dd.synthetic_OLS_coefficient,color=COL[a],marker='x',s=50)
 axs[0].set_xticks([0,1],['未调整P','调整P']);axs[0].set_ylabel('合成模型标准化系数');axs[0].set_title('一个真实蛋白信号，不同误差方差');axs[0].legend(fontsize=9)
 scatter=axs[1].scatter(grid.NO_P_S_MINUS_N,grid.P_S_MINUS_N,c=np.log10(grid.theta_S/grid.theta_N),cmap='coolwarm',s=24,alpha=.8);axs[1].axhline(0,color='.7');axs[1].axvline(0,color='.7');axs[1].set_xlabel('无P的Soma − NULISA差异');axs[1].set_ylabel('有P的Soma − NULISA差异');axs[1].set_title('125组预设误差方差的解析结果');fig.colorbar(scatter,ax=axs[1],label='log10(θS/θN)')
 fig.suptitle('07 反例：条件化后出现检测差异，并不要求存在不同分子形式',y=1.025);fig.text(.5,-.025,'F=0.8P+U；Y=−0.06P+0.10U+η；三检测载荷均为1。纯合成示例，未拟合ADNI误差方差。',ha='center',fontsize=9);fig.text(.5,-.08,'左图圆点：解析值；叉号：20万合成人OLS。右图部分点重叠，因为θT不影响Soma−NULISA对比。',ha='center',fontsize=9);fig.tight_layout();save(fig,'F07_SINGLE_SIGNAL_COUNTEREXAMPLE','Analytical and synthetic counterexample; not fitted to ADNI.')
 cal=pd.read_csv(RESULTS/'NULL_CALIBRATION_WORST_SUMMARY.csv');lo=pd.read_csv(RESULTS/'LEAVE_ONE_TECHNICAL_CLUSTER_POINTS.csv');fig,axs=plt.subplots(1,2,figsize=(13.5,5.2));xx=np.arange(42)
 for name,marker in [('NOMINAL','o'),('BY9_SUFFICIENT','s'),('BY33_SUFFICIENT','^'),('BY42_SUFFICIENT','x')]:
  dd=cal.loc[cal.alpha_name.eq(name)].set_index('test_id').loc[r.test_id];axs[0].plot(xx,dd.max_Wilson95_upper/dd.alpha,marker,ms=3,label=name)
 axs[0].axhline(1,color='red',lw=.8);axs[0].axhline(1.2,color='.6',lw=.7,ls='--');axs[0].set_xlabel('42个固定检验（索引）');axs[0].set_ylabel('最坏零假设Wilson上限 / α');axs[0].set_title('有限情景校准');axs[0].legend(fontsize=7)
 for k,g in enumerate(GG):
  dd=lo.loc[lo.gene.eq(g)];axs[1].scatter(dd.CHANGE_S_MINUS_N,k+np.linspace(-.18,.18,len(dd)),s=10,alpha=.6,label=g)
  point=r.loc[r.gene.eq(g)&r.model.eq('LOG_SHIFT'),'S_MINUS_N'].iloc[0];axs[1].plot(point,k,'kD',ms=5)
 axs[1].set_yticks(range(3),GG);axs[1].axvline(0,color='.6',lw=.8);axs[1].set_xlabel('Soma − NULISA：加入log pTau后的差异增量');axs[1].set_title('逐一剔除技术簇的点估计（无新增P值）');fig.suptitle('08 校准范围与技术簇影响',y=1.02);fig.text(.5,-.02,'每检验×5情景×60,000噪声；7个零均值结构。黑菱形为全样本；点估计诊断不提供额外显著性。',ha='center',fontsize=9);fig.tight_layout();save(fig,'F08_CALIBRATION_AND_CLUSTER_INFLUENCE','Finite-scenario calibration and descriptive leave-one-cluster point refits.')
 csvout(RESULTS/'FIGURE_INVENTORY.csv',inventory);log(f'{len(inventory)//2} figure groups generated as PNG and vector PDF.')
if __name__=='__main__':assert_unsealed();run()
