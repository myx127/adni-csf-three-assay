from config import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
import textwrap
plt.rcParams.update({'font.family':'Microsoft YaHei','axes.unicode_minus':False,'pdf.fonttype':42,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
COL={'N':'#0072B2','S':'#D55E00','T':'#009E73'}
def qfmt(q):return '<0.001' if q<.001 else f'{q:.3f}'
def run():
 formal=pd.read_csv(RESULTS/'FORMAL204_CLINICAL_RESULTS.csv');effects=pd.read_csv(RESULTS/'PLATFORM_CLINICAL_EFFECTS_ALL144.csv');con=pd.read_csv(RESULTS/'PAIRWISE_CLINICAL_CONTRASTS_ALL144.csv');sens=pd.read_csv(RESULTS/'PRESPECIFIED_SENSITIVITY_RESULTS.csv');dec=pd.read_csv(RESULTS/'CALIBRATION_QUALIFICATION_DECISIONS.csv')
 idx=[]
 def save(fig,name,title,caption):
  caption += ' N=NULISA，S=SomaScan，T=TMT-MS。校准未合格者仅作描述，正式资格见结果表。'
  fig.suptitle(title,fontsize=15,y=.99)
  fig.text(.02,.014,caption,fontsize=8,va='bottom',wrap=True)
  fig.tight_layout(rect=[0,.06,1,.96])
  for ext in ['png','pdf']:fig.savefig(FIGURES/(name+'.'+ext),dpi=180,bbox_inches='tight',facecolor='white')
  plt.close(fig);idx.append({'figure_id':len(idx)+1,'basename':name,'title':title,'caption':caption,'png':'figures/'+name+'.png','pdf':'figures/'+name+'.pdf'})
 main=formal.loc[formal.branch.eq('MEM5Y')].copy()
 fig,ax=plt.subplots(figsize=(14,9));labels=[]
 for i,g in enumerate(FOCUS14):
  e=effects.loc[effects.gene.eq(g)&effects.branch.eq('MEM5Y')]
  r=main.loc[main.gene.eq(g)&main.test.eq('SAME')].iloc[0];d=main.loc[main.gene.eq(g)&main.test.eq('DIFF')].iloc[0]
  labels.append(f'{g}  n={r.n_people}  q共同={qfmt(r.q_BY_within_family)} / q差异={qfmt(d.q_BY_within_family)}')
  for a,off in [('N',-.21),('S',0),('T',.21)]:
   z=e.loc[e.assay.eq(a)].iloc[0];ax.errorbar(z.standardized_effect,i+off,xerr=1.95996398454*z.standardized_SE_full_joint,fmt='o',ms=4,color=COL[a],capsize=2,label=a if i==0 else None)
 ax.set_yticks(range(14),labels);ax.invert_yaxis();ax.axvline(0,color='#555',lw=1);ax.grid(axis='x',alpha=.2);ax.set_xlabel('ADNI-MEM分数单位/年，每1个固定蛋白残差SD');ax.legend(title='检测 N / S / T',ncol=3)
 save(fig,'01_focus14_memory_progression_forest','14个共同病理信号：三检测记忆进展关联','区间为点95%置信区间；q来自完整主68项BY家族。正值表示条件下更正向的认知年变化，不代表因果保护。')
 fig,ax=plt.subplots(figsize=(12,9));labs={'S_MINUS_N':('S−N',-.21,'#0072B2'),'T_MINUS_N':('T−N',0,'#D55E00'),'T_MINUS_S':('T−S',.21,'#009E73')}
 for i,g in enumerate(FOCUS14):
  for name,(label,off,color) in labs.items():
   z=con.loc[con.gene.eq(g)&con.branch.eq('MEM5Y')&con.contrast.eq(name)].iloc[0];ax.errorbar(z.difference_standardized_effect,i+off,xerr=1.95996398454*z.SE_full_joint,fmt='o',ms=4,color=color,capsize=2,label=label if i==0 else None)
 ax.set_yticks(range(14),FOCUS14);ax.invert_yaxis();ax.axvline(0,color='#444');ax.grid(axis='x',alpha=.2);ax.set_xlabel('标准化临床关联系数之差：ADNI-MEM单位/年');ax.legend(ncol=3)
 save(fig,'02_focus14_paired_assay_contrasts','14个共同病理信号：配对检测差异','区间为点95%区间，已保留检测方程和蛋白定标不确定度；跨0不构成临床等效证据。正式结论使用2df DIFF及主68项BY。')
 cols=[('MEM5Y','SAME'),('MEM5Y','DIFF'),('MEM_BASE','SAME'),('MEM_BASE','DIFF'),('EF5Y','SAME'),('EF5Y','DIFF')]
 mat=np.array([[formal.loc[formal.gene.eq(g)&formal.branch.eq(b)&formal.test.eq(t),'q_BY_within_family'].iloc[0] for b,t in cols] for g in FOCUS14])
 fig,ax=plt.subplots(figsize=(11,9));im=ax.imshow(np.minimum(-np.log10(np.maximum(mat,1e-12)),4),aspect='auto',cmap='Blues',vmin=0,vmax=4)
 for i in range(14):
  for j in range(6):
   ax.text(j,i,qfmt(mat[i,j]),ha='center',va='center',color='white' if mat[i,j]<.01 else '#222',fontsize=9)
   if mat[i,j]<=.05:ax.add_patch(Rectangle((j-.48,i-.48),.96,.96,fill=False,edgecolor='#D55E00',lw=1.8))
 ax.set_yticks(range(14),FOCUS14);ax.set_xticks(range(6),['记忆进展\n共同','记忆进展\n差异','基线记忆\n共同','基线记忆\n差异','执行进展\n共同','执行进展\n差异']);fig.colorbar(im,ax=ax,label='−log10(q)')
 save(fig,'03_focus14_clinical_evidence_matrix','14个共同病理信号的临床证据矩阵','格内为完整家族BY q：记忆进展主68项，基线记忆与执行进展合并次136项。边框表示q≤0.05；子集不重新校正。')
 integ=pd.read_csv(RESULTS/'INTEGRATED_PATHOLOGY_CLINICAL_FAMILY34.csv').set_index('gene');genes=list(FAMILY34.values());fields=['M0_pathology_same_q_BY34','MEM5Y_SAME_q_BY_within_family','MEM5Y_DIFF_q_BY_within_family','MEM_BASE_SAME_q_BY_within_family','EF5Y_SAME_q_BY_within_family']
 mm=np.array([[integ.loc[g,c] for c in fields] for g in genes]);fig,ax=plt.subplots(figsize=(9,15));im=ax.imshow(np.minimum(-np.log10(np.maximum(mm,1e-12)),4),cmap='Blues',vmin=0,vmax=4,aspect='auto')
 ax.set_yticks(range(34),[g+(' *' if g in FOCUS14 else '') for g in genes]);ax.set_xticks(range(5),['既有pTau共同','记忆进展共同','记忆进展差异','基线记忆共同','执行进展共同'],rotation=25,ha='right')
 for i in range(34):
  for j in range(5):ax.text(j,i,qfmt(mm[i,j]),ha='center',va='center',fontsize=8,color='white' if mm[i,j]<.01 else '#222')
 fig.colorbar(im,ax=ax,shrink=.5,label='−log10(q)')
 save(fig,'04_full34_pathology_clinical_context','原34蛋白完整范围中的病理与临床证据','*为固定14个共同病理关联蛋白。各列使用已声明的不同检验家族；不能将颜色交集解释为因果路径或全项目FDR。')
 fig,axs=plt.subplots(1,4,figsize=(18,10),sharey=True,sharex=True)
 for ax,b,title in zip(axs,['MEM5Y']+SENSITIVITIES,['主模型','不调整pTau','高信息随访','个体等权']):
  for i,g in enumerate(FOCUS14):
   for a,off in [('N',-.2),('S',0),('T',.2)]:
    z=effects.loc[effects.gene.eq(g)&effects.branch.eq(b)&effects.assay.eq(a)].iloc[0];ax.errorbar(z.standardized_effect,i+off,xerr=1.95996398454*z.standardized_SE_full_joint,fmt='o',ms=3,color=COL[a],lw=.8)
  ax.axvline(0,color='#555',lw=1);ax.set_title(title);ax.grid(axis='x',alpha=.2);ax.set_xlabel('记忆分数单位/年')
 axs[0].set_yticks(range(14),FOCUS14);axs[0].invert_yaxis()
 save(fig,'05_focus14_prespecified_sensitivities','14个共同病理信号的固定敏感性','所有分支使用相同蛋白定标单位；无pTau保持同病例。敏感性BY28仅作诊断，不替代主68项结论。区间为点95%。')
 follow=pd.read_csv(RESULTS/'FOLLOWUP_INFORMATION_FAMILY34.csv');ff=follow.loc[follow.branch.eq('MEM5Y')].set_index('gene').loc[FOCUS14]
 fig,axs=plt.subplots(1,2,figsize=(13,9),sharey=True);yy=np.arange(14)
 axs[0].errorbar(ff.median_span_years,yy,xerr=np.array([ff.median_span_years-ff.span_q25,ff.span_q75-ff.median_span_years]),fmt='o',color='#0072B2',capsize=3)
 axs[0].set_yticks(yy,[f'{g} (n={int(ff.loc[g,"n_people"])})' for g in FOCUS14]);axs[0].invert_yaxis();axs[0].set_xlabel('实际随访跨度中位数及四分位区间（年）');axs[0].axvline(5,color='#999',ls='--')
 axs[1].barh(yy,100*ff.last_visit_ge4p5years_fraction,color='#56B4E9');axs[1].set_xlabel('末次访视达到4.5年的参与者比例（%）')
 save(fig,'06_focus14_followup_information','五年观察窗口内的实际随访信息','观察窗口最长五年，不表示所有人均随访五年。时间权重贡献和其他20蛋白的完整描述见CSV。')
 fig,axs=plt.subplots(1,2,figsize=(13,6));x=np.arange(len(dec));colors=['#0072B2' if g in FOCUS14 else '#888888' for g in dec.gene]
 axs[0].scatter(x,dec.max_nominal_Wilson95_upper,c=colors,s=13);axs[0].axhline(.06,color='#D55E00',ls='--');axs[0].set_ylabel('相关零假设最坏Wilson95上界');axs[0].set_xlabel('固定设计序号');axs[0].set_title('名义0.05检验，资格界0.06')
 axs[1].scatter(x,dec.max_tail_Wilson_over_alpha,c=colors,s=13);axs[1].axhline(1,color='#D55E00',ls='--');axs[1].set_ylabel('小尾部Wilson上界/对应alpha');axs[1].set_xlabel('固定设计序号');axs[1].set_title('BY充分阈值的小尾部校准')
 save(fig,'07_actual_design_finite_calibration','实际设计的有限情景校准','蓝色为14重点蛋白。资格还要求边际单侧检验和数值失败率合格；有限情景合格不代表所有误差或失访机制有效。')
 pow=pd.read_csv(RESULTS/'PRIMARY_WORKING_EFFECT_POWER.csv');cols=[(d,a) for d in [.02,.05,.1] for a in ['NOMINAL','BY68_SUFFICIENT']];mm=np.array([[pow.loc[pow.gene.eq(g)&np.isclose(pow.working_effect_abs,d)&pow.alpha_name.eq(a),'worst_DGP_and_sign_rate'].iloc[0] for d,a in cols] for g in FOCUS14])
 fig,ax=plt.subplots(figsize=(11,9));im=ax.imshow(mm,vmin=0,vmax=1,cmap='YlGnBu',aspect='auto');ax.set_yticks(range(14),FOCUS14);ax.set_xticks(range(6),[f'δ={d:g}\n'+('名义' if a=='NOMINAL' else 'BY充分阈值') for d,a in cols])
 for i in range(14):
  for j in range(6):ax.text(j,i,f'{100*mm[i,j]:.1f}%',ha='center',va='center',color='white' if mm[i,j]>.65 else '#222',fontsize=9)
 fig.colorbar(im,ax=ax,label='最坏工作情景及方向下的通过率')
 save(fig,'08_focus14_conditional_working_power','共同记忆关联的条件工作效应功效','δ单位为ADNI-MEM/年/蛋白残差SD，不是临床MCID；BY充分阈值通过率是实际BY排序功效的充分条件下界。')
 csvout(RESULTS/'FIGURE_INDEX.csv',idx)
 log('Eight figure groups generated in PNG and PDF.')
if __name__=='__main__':assert_unsealed();run()

