from config import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
from scipy.stats import norm
for p in [Path('C:/Windows/Fonts/msyh.ttc'),Path('C:/Windows/Fonts/simhei.ttf')]:
 if p.exists():font_manager.fontManager.addfont(str(p));plt.rcParams['font.family']=font_manager.FontProperties(fname=str(p)).get_name();break
plt.rcParams.update({'axes.unicode_minus':False,'font.size':10,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
GENES=['NPTX1','KLK6','PDGFRB'];ASSAYS=['N','S','T'];COLORS={'N':'#185e9d','S':'#c37b19','T':'#278966'}
FIGINDEX=[]
def save(fig,name,description):
 fig.savefig(FIGURES/(name+'.png'),dpi=220,bbox_inches='tight',facecolor='white');fig.savefig(FIGURES/(name+'.pdf'),bbox_inches='tight',facecolor='white');plt.close(fig)
 FIGINDEX.append({'figure':name,'description':description,'PNG':'figures/'+name+'.png','PDF':'figures/'+name+'.pdf'})

def run():
 effects=pd.read_csv(RESULTS/'PLATFORM_CLINICAL_EFFECTS_ALL43.csv');formal=pd.read_csv(RESULTS/'FORMAL32_CLINICAL_RESULTS.csv')
 fig,axes=plt.subplots(1,3,figsize=(13,5.8),sharey=True)
 for ax,gene in zip(axes,GENES):
  part=effects.loc[effects.gene.eq(gene)&effects.branch.eq('MEM5Y')].set_index('assay').loc[ASSAYS]
  for i,assay in enumerate(ASSAYS):
   r=part.loc[assay];ax.errorbar(r.standardized_effect,i,xerr=[[r.standardized_effect-r.standardized_CI95_low],[r.standardized_CI95_high-r.standardized_effect]],fmt='o',capsize=3,color=COLORS[assay])
  ax.axvline(0,color='grey',lw=.8);ax.set_xlim(-.025,.09);ax.set_title(f'{gene} | n={int(part.n_people.iloc[0])}',fontweight='bold');ax.grid(axis='x',alpha=.2)
  ff=formal.loc[formal.gene.eq(gene)&formal.branch.eq('MEM5Y')].set_index('test')
  ax.text(.02,-.30,f'同向关联：P={ff.loc["SAME","p_computed"]:.4g}，BY6 q={ff.loc["SAME","q_BY_within_family"]:.3f}\n检测间差异：BY6 q={ff.loc["DIFF","q_BY_within_family"]:.3f}',transform=ax.transAxes,fontsize=10)
  ax.set_xlabel('ADNI_MEM年变化系数 / 蛋白残差SD')
 axes[0].set_yticks(range(3),['NULISA','SomaScan','TMT']);axes[0].invert_yaxis()
 fig.suptitle('主分析：三检测均指向更正向的记忆年变化，尚无BY6发现\n区间为点估计95%，未作多重校正；调整同日pTau及共同协变量',y=1.01)
 fig.subplots_adjust(wspace=.25,bottom=.27,top=.80)
 save(fig,'01_primary_memory_progression_forest','主要6检验与9个观测平台效应；点区间不等于家族发现')
 contrasts=pd.read_csv(RESULTS/'PAIRWISE_CLINICAL_CONTRASTS_ALL43.csv')
 fig,axes=plt.subplots(1,3,figsize=(13,5.8),sharey=True)
 for ax,gene in zip(axes,GENES):
  part=contrasts.loc[contrasts.gene.eq(gene)&contrasts.branch.eq('MEM5Y')].set_index('contrast').loc[['S_MINUS_N','T_MINUS_N','T_MINUS_S']]
  for i,r in enumerate(part.itertuples()):
   ax.plot([r.CI95_within_design_3pair_Bonf_low,r.CI95_within_design_3pair_Bonf_high],[i,i],color='#bfc5ca',lw=4)
   ax.errorbar(r.difference_standardized_effect,i,xerr=[[r.difference_standardized_effect-r.CI95_pointwise_low],[r.CI95_pointwise_high-r.difference_standardized_effect]],fmt='o',color='#185e9d',capsize=3)
  ax.axvline(0,color='grey',lw=.8);ax.set_xlim(-.085,.065);ax.set_title(gene,fontweight='bold');ax.set_xlabel('检测间年变化关联系数之差');ax.grid(axis='x',alpha=.2)
 axes[0].set_yticks(range(3),['Soma − NULISA','TMT − NULISA','TMT − Soma']);axes[0].invert_yaxis()
 fig.suptitle('未获得跨检测认知关联差异证据，也未检验临床等效\n蓝色为点估计95%区间；浅灰为同设计三对比Bonferroni区间',y=1.02);fig.tight_layout(rect=(0,0,1,.90))
 save(fig,'02_primary_paired_assay_contrasts','三检测直接配对差异及区间；不以未拒绝判等效')
 fig,axes=plt.subplots(1,3,figsize=(14,6),sharey=True)
 branches=[('MEM_BASE','基线记忆'),('EF5Y','执行功能年变化'),('AMY_CL','淀粉样PET Centiloid')]
 for ax,(branch,label) in zip(axes,branches):
  for offset,assay in zip([-.18,0,.18],ASSAYS):
   part=effects.loc[effects.branch.eq(branch)&effects.assay.eq(assay)].set_index('gene').loc[GENES]
   ax.errorbar(part.standardized_effect,np.arange(3)+offset,xerr=[part.standardized_effect-part.standardized_CI95_low,part.standardized_CI95_high-part.standardized_effect],fmt='o',color=COLORS[assay],capsize=2,label=assay)
  ax.axvline(0,color='grey',lw=.8);ax.set_title(label);ax.set_xlabel('结局单位 / 蛋白残差SD');ax.grid(axis='x',alpha=.2)
 axes[0].set_yticks(range(3),GENES);axes[0].invert_yaxis()
 handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,['NULISA','SomaScan','TMT'],loc='lower center',ncol=3,frameon=False)
 fig.suptitle('次级临床与病理关联：无BY26发现\n全部条件于pTau；不同结局横轴单位不同，不能直接比较效应大小',y=1.02);fig.tight_layout(rect=(0,.08,1,.90))
 save(fig,'03_secondary_clinical_pathology_effects','基线记忆、EF进展及独立淀粉样CL；没有新增校正后发现')
 spec=[('MEM5Y','主5年窗口'),('MEM_NO_PTAU','不加pTau'),('MEM_RAW_P','原值pTau'),('MEM3Y','3年窗口'),('MEM_HIGHINFO','高信息随访'),('MEM_EQUAL_PERSON','个体等权'),('MEM_BATCH9','仅Batch9'),('MEM_BASE_ADJUSTED','条件于基线记忆'),('MEM_QUADRATIC','共同二次时间'),('MEM_PET_MAIN','PET子集＋CL')]
 fig,axes=plt.subplots(1,3,figsize=(15,9),sharey=True)
 for ax,gene in zip(axes,GENES):
  for offset,assay in zip([-.18,0,.18],ASSAYS):
   part=effects.loc[effects.gene.eq(gene)&effects.assay.eq(assay)].set_index('branch').loc[[s for s,_ in spec]]
   ax.errorbar(part.standardized_effect,np.arange(len(spec))+offset,xerr=[part.standardized_effect-part.standardized_CI95_low,part.standardized_CI95_high-part.standardized_effect],fmt='o',ms=3,color=COLORS[assay],alpha=.85,capsize=2,label=assay)
  ax.axvline(0,color='grey',lw=.8);ax.set_title(gene,fontweight='bold');ax.set_xlabel('记忆年变化系数 / 固定蛋白残差SD');ax.grid(axis='x',alpha=.2)
 axes[0].set_yticks(range(len(spec)),[s for _,s in spec]);axes[0].invert_yaxis()
 handles,_=axes[0].get_legend_handles_labels();fig.legend(handles,['NULISA','SomaScan','TMT'],loc='lower center',ncol=3,frameon=False)
 fig.suptitle('预设敏感性：保持蛋白单位，病例与条件集合变化分别标记\n全部敏感性均无BY6诊断发现；PET子集与加入CL的影响另作描述拆分',y=1.01);fig.tight_layout(rect=(0,.06,1,.92))
 save(fig,'04_memory_prespecified_sensitivity','9类预设敏感性和主模型；不重选模型')
 cond=pd.read_csv(RESULTS/'PTAU_CONDITIONING_ALGEBRAIC_DECOMPOSITION.csv');pet=pd.read_csv(RESULTS/'PET_SUBSET_VS_CL_POINT_DECOMPOSITION.csv')
 fig,axes=plt.subplots(2,3,figsize=(14,8.5))
 for col,gene in enumerate(GENES):
  ax=axes[0,col]
  for assay in ASSAYS:
   r=cond.loc[cond.gene.eq(gene)&cond.assay.eq(assay)].iloc[0];ax.plot([0,1],[r.standardized_b_without_ptau,r.standardized_b_with_ptau],'-o',color=COLORS[assay],label=assay)
  ax.set_xticks([0,1],['不加pTau','加入pTau']);ax.set_title(gene,fontweight='bold');ax.axhline(0,color='grey',lw=.7);ax.set_ylabel('记忆年变化关联系数');ax.grid(axis='y',alpha=.2)
  ax=axes[1,col]
  for assay in ASSAYS:
   r=pet.loc[pet.gene.eq(gene)&pet.assay.eq(assay)].iloc[0];ax.plot(range(3),[r.full_main_standardized_b,r.PET_subset_same_covariates_b,r.PET_subset_plus_CL_b],'-o',color=COLORS[assay])
  ax.set_xticks(range(3),['全主样本','PET子集\n相同协变量','PET子集\n再加CL']);ax.axhline(0,color='grey',lw=.7);ax.set_ylabel('记忆年变化关联系数');ax.grid(axis='y',alpha=.2)
 handles,_=axes[0,0].get_legend_handles_labels();fig.legend(handles,['NULISA','SomaScan','TMT'],loc='lower center',ncol=3,frameon=False)
 fig.suptitle('结果知情的描述性条件模型拆分：没有新增假设检验\n上：同人群加/不加pTau；下：PET样本选择与CL调整分步展示；均非因果中介',y=1.02);fig.tight_layout(rect=(0,.07,1,.92))
 save(fig,'05_descriptive_conditioning_decomposition','pTau代数关系和PET子集/CL点估计拆分，非因果或新发现')
 import pickle
 ts=pickle.load((PRIVATE/'ANALYSIS_TEMPLATES_RESTRICTED.pkl').open('rb'));data=pickle.load((PRIVATE/'ASSEMBLED_CLINICAL_RESTRICTED.pkl').open('rb'))
 reach=pd.read_csv(RESULTS/'FOLLOWUP_REACH_COUNTS.csv')
 fig,axes=plt.subplots(1,2,figsize=(12,5.5));colors=['#185e9d','#c37b19','#278966']
 for gene,acc,color in zip(GENES,CANDIDATES,colors):
  part=reach.loc[reach.gene.eq(gene)];axes[0].plot(part.year,part.n_with_last_visit_ge_year_minus_half/part.denominator_main,'o-',label=gene,color=color)
  t=ts['MEM5Y__'+acc];w=np.sort([data['trajectories'][(i,'ADNI_MEM',5)]['w'] for i in t['ids']])[::-1];axes[1].plot(np.r_[0,np.arange(1,len(w)+1)/len(w)],np.r_[0,np.cumsum(w)/sum(w)],color=color,label=gene)
 axes[0].set_xticks(range(1,6));axes[0].set_ylim(0,1.03);axes[0].set_xlabel('观察年份（末次访视达到该年−0.5年）');axes[0].set_ylabel('主样本比例');axes[0].set_title('中位跨度约3年；4.5年以上约12%')
 axes[1].plot([0,1],[0,1],'--',color='grey',lw=.8);axes[1].axvline(.1,color='grey',lw=.7);axes[1].set_xlabel('按时间信息从高到低的个体比例');axes[1].set_ylabel('累计时间信息比例');axes[1].set_title('前10%个体约贡献27%时间信息')
 for ax in axes:ax.legend(frameon=False);ax.grid(alpha=.2)
 fig.suptitle('五年是观察窗口；重复访视不能当成独立受试者数',y=1.02);fig.tight_layout()
 save(fig,'06_followup_information_and_selection','随访长度和时间信息集中度；不假定所有人均满5年')
 cal=pd.read_csv(RESULTS/'CALIBRATION_QUALIFICATION_DECISIONS.csv');cal=cal.loc[cal.formal].copy()
 order=[b+'__'+a for b in FORMAL_BRANCHES for a in CANDIDATES if b!='MEM_APOE_INTERACTION' or a=='P09619'];cal=cal.set_index('design_id').loc[order]
 labels=[r.gene+' '+{'MEM5Y':'记忆进展','MEM_BASE':'基线记忆','EF5Y':'EF进展','AMY_CL':'淀粉样CL','MEM_AMY_INTERACTION':'CL交互','MEM_APOE_INTERACTION':'APOE交互'}[r.branch] for r in cal.itertuples()]
 y=np.arange(len(cal));colors=['#4e9a78' if x else '#cf6b61' for x in cal.qualified]
 fig,axes=plt.subplots(1,2,figsize=(13,9),sharey=True)
 axes[0].barh(y,100*cal.max_nominal_false_rate,color=colors)
 axes[0].plot(100*cal.max_nominal_Wilson95_upper,y,'k|',ms=10);axes[0].axvline(5,color='grey',ls='--');axes[0].axvline(6,color='red',ls=':');axes[0].set_xlabel('最坏名义误拒绝率及Wilson上界（%）');axes[0].set_title('名义上界要求≤6%')
 axes[1].barh(y,cal.max_tail_Wilson_over_alpha,color=colors);axes[1].axvline(1,color='red',ls=':');axes[1].set_xlabel('极小尾部Wilson上界 / 对应α');axes[1].set_title('BY充分阈值尾部比值要求≤1')
 axes[0].set_yticks(y,labels);axes[0].invert_yaxis()
 fig.suptitle('正式16设计：13合格；3个淀粉样交互设计未通过\n未通过的6项假设保留在次26家族中，正式P置1',y=1.02);fig.tight_layout()
 save(fig,'07_formal_design_calibration','名义与小尾部独立资格；自由度过关不等于推断合格')
 power=pd.read_csv(RESULTS/'PRIMARY_WORKING_EFFECT_POWER.csv');fig,axes=plt.subplots(1,3,figsize=(13,5.6),sharey=True)
 for ax,gene in zip(axes,GENES):
  for an,col,label in [('NOMINAL','#185e9d','名义α=0.05'),('BY6_SUFFICIENT','#c37b19','BY6第一步充分门槛'),('BY32_SUFFICIENT','#278966','BY32第一步充分门槛')]:
   part=power.loc[power.gene.eq(gene)&power.alpha_name.eq(an)&power.pattern.str.endswith('_1')].copy();part['mag']=part.pattern.str.extract(r'BALANCED_([0-9.]+)_')[0].astype(float);part=part.sort_values('mag');ax.plot(part.mag,part.minimum_conditional_power,'o-',color=col,label=label)
  ax.set_title(gene);ax.set_xticks([.02,.05,.1]);ax.set_ylim(-.02,1.03);ax.set_xlabel('三检测相同的工作效应（MEM单位/年/蛋白SD）');ax.grid(alpha=.2)
 axes[0].set_ylabel('最坏DGP下条件通过率');handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,frameon=False)
 fig.suptitle('条件功效：中小共同效应的信息有限\nBY曲线是第一步充分门槛通过率，不是实际排序后的BY功效；工作效应不是临床阈值',y=1.03);fig.tight_layout(rect=(0,.10,1,.91))
 save(fig,'08_conditional_working_effect_power','固定设计、工作效应下的功效；非观测效应事后功效')
 csvout(RESULTS/'FIGURE_INDEX.csv',FIGINDEX)
 log(f'{len(FIGINDEX)} scientific figure groups generated in PNG/PDF')
if __name__=='__main__':
 assert_unsealed();run()
