from config import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm
import logging
logging.getLogger('fontTools').setLevel(logging.ERROR)
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'ps.fonttype':42,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
records=[]
COLORS={'N':'#297fa0','S':'#bd7b42','T':'#438c66'}
def load(name):return pd.read_csv(RESULTS/(name+'.csv'))
def save(fig,name,caption):
 for ext in ['png','pdf']:fig.savefig(FIGURES/(name+'.'+ext),dpi=300,bbox_inches='tight',facecolor='white')
 plt.close(fig);records.append({'figure':name,'PNG':'figures/'+name+'.png','PDF':'figures/'+name+'.pdf','caption':caption})
def primary():return load('PRIMARY_FAMILY34_SHARED_RESULTS').sort_values('q_BY_shared34')

def forest():
 order=primary().gene.tolist();f=load('PLATFORM_EFFECTS_ALL_MODELS');f=f.loc[f.model.eq('M0')&f.in_family34];fig,ax=plt.subplots(figsize=(10.5,13))
 for j,a in enumerate(['N','S','T']):
  part=f.loc[f.assay.eq(a)].set_index('gene').loc[order];y=np.arange(len(order))+(j-1)*.20
  ax.errorbar(part.scaled_b,y,xerr=np.array([part.scaled_b-part.scaled_CI95_low,part.scaled_CI95_high-part.scaled_b]),fmt='o',markersize=3.8,capsize=1.6,color=COLORS[a],label={'N':'NULISA','S':'SomaScan','T':'TMT-MS'}[a],linewidth=1)
 ax.axvline(0,color='.4',lw=.9);ax.set(yticks=range(len(order)),yticklabels=order,xlabel='P-tau coefficient / assay residual SD (95% interval)',title='Three observed pathology associations | same cases, same covariates');ax.invert_yaxis();ax.legend(frameon=False,loc='lower right')
 fig.text(.5,.01,'M0; 34-protein family. Scaling uncertainty is included. Intervals are pointwise, not simultaneous.\nA same-direction association does not establish proportional response or latent abundance.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.06,1,1])
 save(fig,'01_platform_association_forest','34蛋白三检测标准化pTau系数及含尺度估计不确定性的区间；共同关联不等于响应等效。')

def evidence():
 f=primary();fig,ax=plt.subplots(figsize=(10.5,6.4));labels={'SHARED_ASSOCIATION_AND_RELATIVE_DEVIATION':'Shared + relative deviation','SHARED_ASSOCIATION_ONLY_DEVIATION_UNRESOLVED':'Shared; deviation unresolved','RELATIVE_DEVIATION_SHARED_ASSOCIATION_UNRESOLVED':'Deviation; shared unresolved','UNRESOLVED':'Both unresolved'}
 colors=['#9b5b82','#3f8d89','#c08a46','#adb5bc']
 for (state,label),color in zip(labels.items(),colors):
  part=f.loc[f.evidence_state.eq(state)];ax.scatter(-np.log10(part.q_BY_shared34.clip(lower=1e-300)),-np.log10(part.prior_H_q_BY_by_model.clip(lower=1e-300)),s=55,label=f'{label} (n={len(part)})',color=color,edgecolor='white',linewidth=.5)
 for gene in ['NPTX1','KLK6','PDGFRB','MDH1','NRGN','PGK1']:
  r=f.loc[f.gene.eq(gene)].iloc[0];x=-np.log10(r.q_BY_shared34);y=-np.log10(r.prior_H_q_BY_by_model);ax.annotate(gene,(x,y),xytext=(6,8 if gene!='PDGFRB' else -16),textcoords='offset points',fontsize=9)
 for fn in [ax.axvline,ax.axhline]:fn(-np.log10(.05),color='#ad5757',ls='--',lw=1)
 ax.set(xlabel='Shared association: -log10(BY q, family 34)',ylabel='Original relative deviation: -log10(BY q)',title='Shared association and relative-response deviation can coexist');ax.margins(.12);ax.legend(frameon=False,loc='upper right',fontsize=9)
 fig.text(.5,.015,'M0; each axis retains its own 34-test family. Combined-68 results are reported separately.\nUnresolved deviation is not equivalence. Coincident points may overlap; all proteins are in the CSV.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.08,1,1])
 save(fig,'02_shared_and_deviation_evidence','共同关联与原相对偏离双证据轴；12共享关联、2两者兼有、1仅偏离、19未定，分别沿用各自家族。')

def sensitivity():
 allr=load('ALL35_THREE_MODEL_SHARED_RESULTS');order=primary().gene.tolist();matrix=allr.loc[allr.in_family34].pivot(index='gene',columns='model',values='q_BY_shared34').reindex(order)[MODELS]
 fig,ax=plt.subplots(figsize=(7.8,12.5));im=ax.imshow((-np.log10(matrix.clip(lower=1e-300))).clip(upper=6),cmap='YlGnBu',vmin=0,vmax=6,aspect='auto')
 for i in range(len(matrix)):
  for j in range(3):
   q=matrix.iloc[i,j];ax.text(j,i,f'{q:.2g}',ha='center',va='center',fontsize=8,color='white' if q<1e-3 else '.15')
   if q<=.05:ax.add_patch(plt.Rectangle((j-.49,i-.48),.98,.96,fill=False,edgecolor='#a86242',lw=1))
 ax.set(xticks=range(3),xticklabels=['M0','+ APOE','+ Education/APOE'],yticks=range(len(matrix)),yticklabels=matrix.index,title='Shared association: fixed covariate sensitivity');fig.colorbar(im,ax=ax,label='-log10(q), color capped at 6',fraction=.05,pad=.03)
 fig.text(.5,.015,'Numbers are actual BY q values in each fixed 34-protein family; borders mark q <=0.05.\nSame individuals and same P-tau scale in all models. A q change is not causal mediation.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.055,1,1])
 save(fig,'03_shared_association_APOE_sensitivity','同一病例三模型q热图，三模型相同14项通过；色阶截断仅为显示，文字保留实际q。')

def calibration():
 null=load('NULL_CALIBRATION_SUMMARY');f=null.loc[null.threshold.eq('NOMINAL')].copy();f['null_group']=f.pattern.str.extract(r'^(ZERO_ALL|ONE_ZERO|TWO_ZERO|MIXED_STRONG|MIXED_NEAR_ZERO)')
 g=f.groupby('null_group').agg(rate=('rate','max'),upper=('Wilson95_high','max'));fig,axes=plt.subplots(1,2,figsize=(13,5))
 axes[0].bar(range(len(g)),100*g.rate,color='#4b8c9e');axes[0].scatter(range(len(g)),100*g.upper,color='#b07744',label='Max cellwise Wilson upper')
 axes[0].set(xticks=range(len(g)),xticklabels=['Mixed near zero','Mixed strong','One zero','Two zero','All zero'],ylabel='False shared-association rejection (%)',title='Nominal alpha=5%; displayed range below 0.6%',ylim=(0,.6));axes[0].tick_params(axis='x',rotation=18);axes[0].legend(frameon=False,fontsize=8)
 for j,threshold in enumerate(['BY34_SUFFICIENT','BY68_SUFFICIENT']):
  part=null.loc[null.threshold.eq(threshold)].copy();part['ratio']=part.Wilson95_high/part.alpha;s=part.groupby('model').ratio.max().reindex(MODELS)
  axes[1].bar(np.arange(3)+(j-.5)*.32,s,width=.32,label=threshold.replace('_SUFFICIENT',''),color=['#4b8c9e','#b58d63'][j])
 axes[1].axhline(1,color='#ad5757',ls='--');axes[1].set(xticks=range(3),xticklabels=MODELS,ylabel='Maximum Wilson upper / tested alpha',title='Strict small-tail finite-grid checks',ylim=(0,1.15));axes[1].legend(frameon=False,fontsize=9)
 fig.text(.5,.015,'105 designs x 5 DGPs x 40000 independent noise draws; all 25 null configurations reuse each draw.\nIntervals are cellwise. Passing finite grids does not prove validity for every error distribution or parameter value.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.12,1,1])
 save(fig,'04_calibration_and_small_tail_checks','部分零假设及小尾部校准；525噪声单元、每单元4万次；所有有限门槛通过，不作全参数保证。')

def power():
 f=load('POWER_WORST_DGP_SUMMARY');f=f.loc[f.in_family34&f.pattern.str.startswith('BALANCED')&f.sign.eq(1)]
 fig,axes=plt.subplots(1,2,figsize=(12,5))
 for ax,threshold,title in zip(axes,['NOMINAL','BY34_SUFFICIENT'],['Nominal alpha=0.05','BY34 sufficient threshold']):
  for model,color in zip(MODELS,['#317f9c','#a57547','#568b67']):
   g=f.loc[f.model.eq(model)&f.threshold.eq(threshold)].groupby('magnitude').worst_DGP_rate
   med=g.median();ax.plot(med.index,med,'o-',label=model,color=color)
   if model=='M0':ax.fill_between(med.index,g.quantile(.1),g.quantile(.9),color=color,alpha=.15)
  ax.axhline(.8,color='.5',ls='--',lw=.8);ax.set(xlabel='True coefficient in each assay (residual SD / P SD)',ylabel='Minimum success probability across 5 fixed DGPs',ylim=(0,1.03),title=title);ax.legend(frameon=False)
 fig.text(.5,.015,'Balanced positive-association scenarios, conditional on 34 covariance templates. Curves: template median.\nShading: M0 10th-90th template percentiles, not confidence bands. This is not actual BY-ranking power.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.12,1,1])
 save(fig,'05_conditional_shared_association_power','真实同向关联条件功效，五DGP最差率及模板变异；不使用观察效应作真值，不是全BY排序功效。')

def weakest():
 f=primary();f['strength']=[max(min(r.N_z,r.S_z,r.T_z),-max(r.N_z,r.S_z,r.T_z)) for r in f.itertuples()];f=f.sort_values('strength',ascending=False)
 fig,ax=plt.subplots(figsize=(9.5,12));ax.barh(range(len(f)),f.strength,color=f.weakest_assay.map(COLORS))
 ax.axvline(norm.isf(.025),color='.4',ls='--',label='Nominal bidirectional threshold')
 ax.axvline(norm.isf(ALPHAS['BY34_SUFFICIENT']/2),color='#ad5757',ls=':',label='Conservative BY34 sufficient threshold')
 ax.axvline(0,color='.5',lw=.7);ax.set(yticks=range(len(f)),yticklabels=f.gene,xlabel='max(min(z), -max(z)); negative = mixed observed signs',title='The weakest assay determines shared-association evidence');ax.invert_yaxis();ax.legend(frameon=False,fontsize=9,loc='lower right')
 fig.text(.5,.015,'Bar color: weakest assay in the selected direction (NULISA blue, Soma brown, TMT green).\nThe dotted line is sufficient for BY regardless of other P values; it is not the actual rank-specific cutoff.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.06,1,1])
 save(fig,'06_weakest_platform_information','最弱平台的方向性z信息；混合观察方向不能通过共同关联；BY足够线不等于实际排序阈值。')

def candidates():
 platform=load('PLATFORM_EFFECTS_ALL_MODELS');allr=load('ALL35_THREE_MODEL_SHARED_RESULTS');fig,axes=plt.subplots(3,2,figsize=(13,12))
 for row,gene in enumerate(['NPTX1','KLK6','PDGFRB']):
  ax=axes[row,0];labels=[]
  for i,(model,assay) in enumerate([(m,a) for a in ['N','S','T'] for m in MODELS]):
   r=platform.loc[platform.gene.eq(gene)&platform.model.eq(model)&platform.assay.eq(assay)].iloc[0];ax.errorbar(r.scaled_b,i,xerr=np.array([[r.scaled_b-r.scaled_CI95_low],[r.scaled_CI95_high-r.scaled_b]]),fmt={'M0':'o','APOE':'^','BOTH':'s'}[model],color=COLORS[assay],markersize=5,capsize=2);labels.append(assay+' / '+model)
  ax.axvline(0,color='.5',lw=.8);ax.set(yticks=range(9),yticklabels=labels,title=gene+': observed pathology coefficients',xlabel='Coefficient / residual SD (95% interval)');ax.invert_yaxis()
  part=allr.loc[allr.gene.eq(gene)].set_index('model').loc[MODELS];ax=axes[row,1]
  for j,col,label,color in [(0,'q_BY_shared34','Shared association','#398b89'),(1,'prior_H_q_BY_by_model','Relative deviation','#9e6f8e')]:
   ax.bar(np.arange(3)+(j-.5)*.30,-np.log10(part[col].clip(lower=1e-300)),width=.30,label=label,color=color)
  ax.axhline(-np.log10(.05),color='#ad5757',ls='--');ax.set(xticks=range(3),xticklabels=MODELS,ylabel='-log10(BY q)',title=gene+': two distinct evidence axes')
  if row==0:ax.legend(frameon=False,fontsize=9)
 fig.text(.5,.01,'Each question retains its own fixed family; combined-68 results are in the tables.\nShared association does not prove response equivalence, true abundance, molecular form, or independent replication.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.05,1,1])
 save(fig,'07_core_candidate_evidence','NPTX1/KLK6共同关联在APOE敏感性保留，原偏离证据具有模型依赖；PDGFRB未过共同关联BY。')

if __name__=='__main__':
 forest();evidence();sensitivity();calibration();power();weakest();candidates()
 csvout(RESULTS/'FIGURE_INDEX.csv',records);log('Seven scientific figures exported as PNG/vector-PDF pairs')
