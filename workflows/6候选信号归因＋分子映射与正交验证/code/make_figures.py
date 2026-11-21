from config import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import logging
logging.getLogger('fontTools').setLevel(logging.ERROR)
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'ps.fonttype':42,'font.size':10,
 'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','savefig.facecolor':'white'})
COLORS={'NPTX1':'#247BA0','PDGFRB':'#9963A5','KLK6':'#DD8452'}
GENES=['PDGFRB','NPTX1','KLK6'];inventory=[]

def save(fig,name,caption):
 for ext in ['png','pdf']:fig.savefig(FIGURES/(name+'.'+ext),dpi=300,bbox_inches='tight')
 plt.close(fig);inventory.append({'figure':name,'PNG':'figures/'+name+'.png','PDF':'figures/'+name+'.pdf','caption':caption})

def covariates():
 f=pd.read_csv(RESULTS/'COVARIATE_FOUR_MODEL_RESULTS.csv');fig,axes=plt.subplots(1,3,figsize=(12.4,4.6),sharey=True);names=['M0','EDU','APOE','BOTH']
 for ax,gene in zip(axes,GENES):
  x=f.loc[f.gene.eq(gene)].set_index('model').reindex(names);ax.bar(np.arange(4),x.Q,color=COLORS[gene],alpha=.85,width=.68)
  for i,r in enumerate(x.itertuples()):ax.text(i,r.Q+.35,f'q={r.q_BY_exploratory:.3f}',ha='center',fontsize=9)
  ax.set(xticks=np.arange(4),xticklabels=['Base','+Education','+APOE','+Both'],title=f'{gene} (n={int(x.n.iloc[0])})',ylim=(0,18));ax.tick_params(axis='x',rotation=20)
 axes[0].set_ylabel('Joint relative-constraint statistic Q')
 fig.suptitle('Covariate attribution on the same cases and fixed 34-protein family',fontsize=14,y=1.02)
 fig.text(.5,-.02,'Labels: BY q values. Added models are postdiscovery explanations; primary findings remain frozen.',ha='center',fontsize=9)
 fig.tight_layout();save(fig,'01_covariate_attribution','相同病例四模型Q及固定34家族BY q；新增模型为事后解释，不替换主分析。')

def decomposition():
 h=pd.read_csv(RESULTS/'H_V_AND_B_C_DECOMPOSITION.csv');q=pd.read_csv(RESULTS/'BY_Q_OWN_VS_FAMILY_DECOMPOSITION.csv');fig,axes=plt.subplots(1,2,figsize=(12,4.9))
 parts=[(axes[0],h.loc[h.comparison.eq('M0_to_APOE')],'delta_Q_H_point_contribution','delta_Q_variance_contribution','Change in Q','H point','Variance metric'),
  (axes[1],q.loc[q.model.eq('APOE')],'q_change_own_p_contribution','q_change_other_proteins_p_contribution','Change in BY q','Own P value','Other family P values')]
 for ax,data,c1,c2,yl,l1,l2 in parts:
  data=data.set_index('gene').reindex(GENES);x=np.arange(3);ax.bar(x-.18,data[c1],.36,label=l1,color='#3B82A0');ax.bar(x+.18,data[c2],.36,label=l2,color='#D89545')
  total=data[c1]+data[c2]
  ax.plot(x,total,'ko',label='Total change',markersize=5);ax.axhline(0,color='.5',linewidth=.8);ax.set(xticks=x,xticklabels=GENES,ylabel=yl);ax.legend(loc='best',frameon=False,fontsize=9)
 axes[0].set_title('APOE-only: H and V change together');axes[1].set_title('APOE-only: q depends on the entire family')
 fig.text(.5,-.01,'Exact two-order Shapley decompositions. Algebraic model attribution; not causal mediation or a biological percentage.',ha='center',fontsize=9)
 fig.tight_layout();save(fig,'02_Q_and_BY_decomposition','APOE单独加入后的精确代数分解；KLK6的q升高由其他家族P值改变主导，不能称自身信号变弱。')

def sequences():
 regions=pd.read_csv(RESULTS/'UNIPROT_REGIONS_AND_PROCESSING.csv');peps=pd.read_csv(RESULTS/'LEGACY_PEPTIDE_CANONICAL_MAPPING.csv');lengths={'PDGFRB':1106,'NPTX1':432,'KLK6':244}
 fig,axes=plt.subplots(3,1,figsize=(13.8,9),gridspec_kw={'hspace':.75})
 for ax,gene in zip(axes,GENES):
  L=lengths[gene];ax.add_patch(Rectangle((1,.65),L-1,.18,color='#d4dae0'))
  for r in regions.loc[regions.gene.eq(gene)].itertuples():
   if r.feature_type=='Signal':ax.add_patch(Rectangle((r.start,.65),r.end-r.start+1,.18,color='#93bfe5'))
   if r.feature_type=='Propeptide':ax.add_patch(Rectangle((r.start,.65),r.end-r.start+1,.18,color='#efa24d'))
   if r.feature_type=='Transmembrane':ax.add_patch(Rectangle((r.start,.62),r.end-r.start+1,.24,color='#bc4c59'));ax.text((r.start+r.end)/2,.95,'TM',ha='center',fontsize=9)
  for i,r in enumerate(peps.loc[peps.gene.eq(gene)].itertuples()):
   start=int(r.canonical_start_1based);end=int(r.canonical_end_1based);yy=.3 if i%2==0 else -.02
   ax.plot([start,end],[yy,yy],linewidth=5,color='#2d926c',solid_capstyle='butt');ax.plot([(start+end)/2]*2,[yy+.02,.64],color='#2d926c',linewidth=.7)
   ax.text((start+end)/2,yy-.13,f'{r.peptide}  {start}-{end}',ha='center',va='top',fontsize=10)
  if gene=='PDGFRB':
   ax.text(280,.48,'Extracellular 33-532',ha='center',fontsize=10);ax.text(815,.48,'Cytoplasmic 554-1106',ha='center',fontsize=10)
   ax.text(L/2,.15,'No direct peptide coverage in audited targeted panels',ha='center',color='#a43f4c',fontsize=11)
  ax.set(xlim=(-.02*L,1.03*L),ylim=(-.46,1.15),yticks=[],title=f'{gene} | canonical length {L} aa',xlabel='Canonical amino-acid position')
  for side in ['left','right','top']:ax.spines[side].set_visible(False)
 fig.suptitle('What is mapped: sequence regions and released legacy MRM peptides',fontsize=15,y=.98)
 fig.text(.5,.006,'Blue: signal peptide; orange: activation peptide; red: transmembrane; green: measured tryptic peptide.\nNo NULISA/Soma binding epitope or original TMT peptide-level evidence is available. Internal peptides do not establish intact/active species.',ha='center',fontsize=9)
 save(fig,'03_molecular_sequence_mapping','官方规范序列与旧版直接肽段位置。表位、实际蛋白形式和活化状态仍未知；PDGFRB未被直接靶向覆盖。')

def coverage():
 matrix=np.array([[0,0,0],[0,2,3],[0,2,3]])
 fig,axes=plt.subplots(1,2,figsize=(12.5,4.8),gridspec_kw={'width_ratios':[1,1.3]})
 ax=axes[0];ax.imshow(matrix,cmap='Blues',vmin=0,vmax=3,aspect='auto')
 for i in range(3):
  for j in range(3):ax.text(j,i,str(matrix[i,j]),ha='center',va='center',color='white' if matrix[i,j]>=2 else '.2',fontsize=17)
 ax.set(xticks=range(3),xticklabels=GENES,yticks=range(3),yticklabels=['Emory 2022/2023','2014 MRM: audit only','2019 Caprion MRM'],title='Direct peptide coverage (distinct sequences)')
 f=pd.read_csv(RESULTS/'LEGACY_DIRECT_DESIGN_SCREEN.csv');f=f.loc[f.cohort.eq('S5_SAME_DATE')];ax=axes[1];labels=[r.gene+' '+r.peptide for r in f.itertuples()];y=np.arange(len(f))
 ax.barh(y-.17,f.n,.32,label='Exact-date matched n',color='#438caa');ax.barh(y+.17,f.residual_df,.32,label='Residual df',color='#db9257');ax.axvline(40,linestyle='--',color='#ae3346',label='Frozen df >=40 gate')
 ax.set(yticks=y,yticklabels=labels,xlim=(0,115),title='Direct relative-pattern design is insufficient');ax.invert_yaxis();ax.legend(frameon=False,fontsize=8,loc='lower right')
 fig.text(.5,-.01,'2014 aliases L01/L02 count as one experiment; exact-date/batch crosswalk is unresolved. 2019 is a distinct experiment, same ADNI cohort.',ha='center',fontsize=9)
 fig.tight_layout();save(fig,'04_direct_coverage_and_design','原三项发现的直接靶向覆盖和同日设计门槛；覆盖不等于推断合格。')

def legacy_associations():
 f=pd.read_csv(RESULTS/'LEGACY_NEW_RID_PEPTIDE_PTAU_ASSOCIATION.csv');f=f.loc[f.cohort.eq('NEW_RID_CLINICAL')].reset_index(drop=True);fig,ax=plt.subplots(figsize=(11.8,5.2));y=np.arange(5)
 for i,r in enumerate(f.itertuples()):
  ax.errorbar(r.b_P_per_residual_SD,i,xerr=[[r.b_P_per_residual_SD-r.CI_per_residual_SD_conditional_on_SD_low],[r.CI_per_residual_SD_conditional_on_SD_high-r.b_P_per_residual_SD]],fmt='o',color=COLORS[r.gene],capsize=3,markersize=7)
  ax.text(1.15,i,f'P+={r.p_positive_direction:.4f}   BY q={r.q_BY5_single_assay:.3f}',va='center',fontsize=10)
 ax.axvline(0,color='.55',linestyle='--');ax.set(yticks=y,yticklabels=[r.gene+' | '+r.peptide for r in f.itertuples()],xlim=(-.15,1.85),xlabel='Peptide-pTau coefficient per residual SD / 1 SD pTau',title='95 subjects disjoint from all 571 Stage5 subjects | same ADNI cohort');ax.invert_yaxis()
 fig.text(.5,-.015,'SUM(plate CR1, run CR1); one-sided positive-direction t test (df=4); fixed 5-peptide BY.\n95% intervals scale the raw coefficient CI by observed residual SD, conditional on that SD. This tests a different question from relative-pattern replication.',ha='center',fontsize=9)
 fig.tight_layout();save(fig,'05_new_RID_direct_peptide_associations','95名主样本之外的ADNI受试者、五肽段正向pTau关联。图中标准化区间以观测残差SD为条件；全部BY q>0.05。')

def four_correlations():
 f=pd.read_csv(RESULTS/'FOUR_ASSAY_RESIDUAL_CORRELATIONS.csv');pairs=list(f[['gene','peptide']].drop_duplicates().itertuples(index=False,name=None));names=['NULISA','Soma','TMT','targetMS'];fig,axes=plt.subplots(1,5,figsize=(16.2,4.2))
 for ax,(gene,pep) in zip(axes,pairs):
  part=f.loc[f.gene.eq(gene)&f.peptide.eq(pep)];R=np.eye(4)
  for r in part.itertuples():
   i=names.index(r.assay_i);j=names.index(r.assay_j);R[i,j]=R[j,i]=r.residual_correlation_after_same_four_assay_design
  im=ax.imshow(R,vmin=-1,vmax=1,cmap='RdBu_r')
  for i in range(4):
   for j in range(4):ax.text(j,i,f'{R[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if abs(R[i,j])>.6 else '.15')
  ax.set(xticks=range(4),xticklabels=['N','S','T','M'],yticks=range(4),yticklabels=['N','S','T','M'],title=gene+'\n'+pep)
 fig.colorbar(im,ax=axes.ravel().tolist(),fraction=.02,pad=.02);fig.suptitle('Same-design residual correlations | four other genes, not the three discoveries',fontsize=14,y=1.04)
 fig.text(.5,.02,'N: NULISA | S: Soma | T: original TMT | M: Emory targeted MS. Same subjects and common four-assay adjustment; descriptive residual structure.',ha='center',fontsize=9)
 save(fig,'06_four_assay_residual_correlations','四种其他蛋白五肽段的四检测残差相关；同病例共同设计，不作为原三项发现验证。')

def ratios():
 f=pd.read_csv(RESULTS/'FOUR_ASSAY_CROSSBLOCK_RESULTS.csv');fig,ax=plt.subplots(figsize=(11.3,5.3))
 ax.axvspan(1/1.5,1.5,color='#b6d9cc',alpha=.4,label='Exploratory working tolerance x1.5');ax.axvline(1,color='.4',linestyle='--')
 for i,r in enumerate(f.itertuples()):
  if np.isfinite(r.direction_ratio):
   strong=r.cross_block_cov_CI95_all_positive;ax.errorbar(r.direction_ratio,i,xerr=[[r.direction_ratio-r.direction_ratio_CI95_low],[r.direction_ratio_CI95_high-r.direction_ratio]],fmt='o' if strong else 'o',mfc='#2d8468' if strong else 'white',color='#2d8468' if strong else '#999999',capsize=3)
  else:ax.text(.12,i,'Ratio undefined: negative cross-block covariance',va='center',fontsize=9)
  ax.text(5.7,i,f'q(eq x1.5)={f.iloc[i]["q_BY5_equiv_factor_1.5"]:.3f}',va='center',fontsize=9)
 ax.set_xscale('log');ax.set(xlim=(.1,14),yticks=range(len(f)),yticklabels=[r.gene+' | '+r.peptide for r in f.itertuples()],xlabel='R = (c_NT x c_SM) / (c_NM x c_ST), with 95% CI',title='Exploratory cross-block direction | no BY-supported equivalence');ax.invert_yaxis();ax.legend(frameon=False,loc='lower left',fontsize=9)
 fig.text(.5,-.015,'Filled: all four cross-block covariance 95% CIs positive. Hollow: insufficient cross-block information.\nZero cross-block error covariance remains an assumption. Non-rejection of a tetrad is not evidence of equivalence.',ha='center',fontsize=9)
 fig.tight_layout();save(fig,'07_four_assay_direction_ratios','跨亲和-MS协方差块方向比值及95%CI；NPTX2区间较窄，但1.5倍工作容差的BY q=0.137。')

def response_equivalence():
 eq=pd.read_csv(RESULTS/'OPERATIONAL_EQUIVALENCE_ALL35.csv');eligible=eq.loc[eq.in_family34&eq.margin_working_SD.eq(.5)&eq.eligible_for_operational_equivalence,['accession','gene']]
 ci=pd.read_csv(RESULTS/'RELATIVE_RESPONSE_FIELLER_INTERVALS.csv');ci=ci.loc[ci.margin_working_SD.eq(.5)];fig,axes=plt.subplots(1,2,figsize=(12,7),sharey=True)
 for ax,contrast,title in zip(axes,['S_relative_to_N','T_relative_to_N'],['Soma relative to NULISA','TMT relative to NULISA']):
  ax.axvspan(-.5,.5,color='#c2daca',alpha=.5);ax.axvline(0,color='.4',linewidth=.8)
  for i,r in enumerate(eligible.itertuples()):
   p=ci.loc[ci.accession.eq(r.accession)&ci.contrast.eq(contrast)].iloc[0];ax.errorbar(p.delta_point,i,xerr=[[p.delta_point-p.fieller90_lower],[p.fieller90_upper-p.delta_point]],fmt='o',color='#386f96',markersize=4,capsize=2)
  ax.set(title=title,xlabel='Relative-response contrast (residual SD / 1 SD pTau)',yticks=range(len(eligible)),yticklabels=eligible.gene)
 axes[0].set_xscale('symlog',linthresh=.5);axes[0].set_xticks([-8,-2,-.5,0,.5]);axes[0].set_xticklabels(['-8','-2','-0.5','0','0.5'])
 axes[0].set_xlabel('Relative-response contrast (symlog beyond +/-0.5)')
 axes[0].invert_yaxis();fig.suptitle('11 eligible proteins | Fieller 90% intervals | 0 of 34 BY-equivalent',fontsize=14,y=1.01)
 fig.text(.5,-.01,'Shading: widest exploratory working margin (+/-0.5 SD). Other proteins fail information/structure gates.\nMargins 0.2 / 0.35 / 0.5 are postdiscovery working criteria, not established biological equivalence standards.',ha='center',fontsize=9)
 fig.tight_layout();save(fig,'08_operational_response_equivalence','11个符合结构和信息要求的蛋白、两个相对响应的Fieller90%区间；固定34家族在三档工作容差下均无BY等效。')

if __name__=='__main__':
 covariates();decomposition();sequences();coverage();legacy_associations();four_correlations();ratios();response_equivalence()
 csvout(RESULTS/'FIGURE_INDEX.csv',inventory);log(f'{len(inventory)} scientific figures exported as standalone PNG and vector PDF')
