from config import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logging
logging.getLogger('fontTools').setLevel(logging.ERROR)
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'ps.fonttype':42,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
records=[]
def save(fig,name,caption):
 for ext in ['png','pdf']:fig.savefig(FIGURES/(name+'.'+ext),dpi=300,bbox_inches='tight',facecolor='white')
 plt.close(fig);records.append({'figure':name,'PNG':'figures/'+name+'.png','PDF':'figures/'+name+'.pdf','caption':caption})

def bottlenecks():
 info=pd.read_csv(RESULTS/'CONSENSUS_INFORMATION_BOTTLENECKS.csv');info=info.loc[info.in_family34];fig,axes=plt.subplots(1,2,figsize=(13.2,4.8),gridspec_kw={'width_ratios':[1,1.5]})
 axes[0].barh(range(5),[34,15,11,6,0],color=['#536e82','#729eab','#8db79c','#bd9bc3','#cc7664'])
 axes[0].set(yticks=range(5),yticklabels=['Frozen family','Point structure compatible','All equivalence gates','Nominal equiv. (margin 0.5)','BY equiv. (margin 0.5)'],xlim=(0,38),title='Different reasons for a zero BY count');axes[0].invert_yaxis()
 for i,n in enumerate([34,15,11,6,0]):axes[0].text(n+.5,i,str(n),va='center')
 cats=['GATE_INELIGIBLE','POINT_OUTSIDE_WORKING_MARGIN','INTERVAL_PRECISION_INSUFFICIENT','NOMINAL_EQUIV_NOT_BY'];labels=['Gate ineligible','Point outside margin','Interval imprecision','Nominal only'];colors=['#aab5be','#d89561','#6b99ba','#ad8ab9'];left=np.zeros(3)
 for cat,label,color in zip(cats,labels,colors):
  counts=np.array([int((info.loc[info.margin.eq(m),'main_bottleneck']==cat).sum()) for m in MARGINS]);axes[1].barh(range(3),counts,left=left,color=color,label=label)
  for i,c in enumerate(counts):
   if c:axes[1].text(left[i]+c/2,i,str(c),ha='center',va='center',fontsize=10)
  left+=counts
 axes[1].set(yticks=range(3),yticklabels=['+/-0.2 SD','+/-0.35 SD','+/-0.5 SD'],xlabel='Proteins in the same fixed family (n=34)',title='Working-margin information bottlenecks');axes[1].invert_yaxis();axes[1].legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.18),ncol=2,fontsize=9)
 fig.tight_layout();save(fig,'01_consensus_evidence_bottlenecks','34蛋白家族的评价门槛及三档容差瓶颈；6个名义等效不等于BY等效，门槛不合格不等于真实不一致。')

def intervals():
 info=pd.read_csv(RESULTS/'CONSENSUS_INFORMATION_BOTTLENECKS.csv');genes=info.loc[info.in_family34&info.margin.eq(.5)&info.equivalence_eligible_Stage6,['accession','gene']]
 ci=pd.read_csv(RESULTS/'FIELLER_PRECISION_AND_DENOMINATOR_INFORMATION.csv');fig,axes=plt.subplots(1,2,figsize=(13,7),sharey=True)
 for ax,contrast,title in zip(axes,['S_relative_to_N','T_relative_to_N'],['Soma relative to NULISA','TMT relative to NULISA']):
  ax.axvspan(-.5,.5,color='#d6e7dc',alpha=.7);ax.axvline(0,color='.5',linewidth=.8)
  for i,g in enumerate(genes.itertuples()):
   part=ci.loc[ci.accession.eq(g.accession)&ci.contrast.eq(contrast)]
   for off,label,color in [(-.12,'NOMINAL','#337a9b'),(.12,'BY_RANK1_SUFFICIENT','#bd7b42')]:
    r=part.loc[part.threshold.eq(label)].iloc[0]
    if r.bounded:ax.plot([r.lower,r.upper],[i+off]*2,color=color,linewidth=1.6);ax.plot(r.delta_point,i+off,'o',color=color,markersize=4)
    else:ax.plot(r.delta_point,i+off,'x',color=color,markersize=5);ax.text(.98,i+off,'unbounded',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=7,color=color)
  ax.set_xscale('symlog',linthresh=.5);ax.set(title=title,yticks=range(len(genes)),yticklabels=genes.gene,xlabel='Relative-response contrast (symlog beyond +/-0.5)');ax.margins(x=.08)
 axes[0].invert_yaxis();fig.suptitle('Fieller confidence sets at two thresholds | actual observed data',y=1.02,fontsize=14)
 fig.text(.5,-.015,'Blue: 90% Fieller interval. Orange: confidence set for conservative BY-sufficient alpha*=0.000357.\nUnbounded sets are labeled, not replaced by narrow intervals. Shading: widest fixed working margin; no new tolerance selected.',ha='center',fontsize=9)
 fig.tight_layout();save(fig,'02_Fieller_precision_and_multiplicity','真实观测的Fieller90%与BY足够阈值置信集合，明确标示无界区间；工作容差保持不变。')

def power_heatmaps():
 f=pd.read_csv(RESULTS/'POWER_WORST_SCENARIO_INFORMATION.csv');f=f.loc[f.gamma_common_P.eq(.55)&f.in_family34];order=f[['accession','gene']].drop_duplicates().sort_values('gene');fig,axes=plt.subplots(2,3,figsize=(13.8,12),sharey=True)
 for row,threshold,title in [(0,'NOMINAL','Nominal alpha=0.05'),(1,'BY_RANK1_SUFFICIENT','Conservative BY-sufficient alpha*')]:
  for col,m in enumerate(MARGINS):
   part=f.loc[f.threshold.eq(threshold)&f.margin.eq(m)];mat=part.pivot(index='gene',columns='k_independent_panels',values='worst_scenario_rate').reindex(order.gene)[K_GRID];ax=axes[row,col];im=ax.imshow(mat,vmin=0,vmax=1,cmap='YlGnBu',aspect='auto')
   for i in range(len(mat)):
    for j in range(4):ax.text(j,i,f'{mat.iloc[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if mat.iloc[i,j]>.6 else '.15')
   ax.set(title=title+f'\nWorking margin +/-{m} SD',xticks=range(4),xticklabels=['1x','2x','4x','8x'],yticks=range(len(mat)),yticklabels=mat.index,xlabel='New independent complete design panels')
 fig.colorbar(im,ax=axes.ravel().tolist(),fraction=.02,pad=.025,label='Minimum gated success probability across 3 fixed DGPs')
 fig.suptitle('True relative-response equivalence: 15 real-covariance templates, common gamma=0.55',y=.98,fontsize=14)
 fig.text(.5,.015,'2000 repeats per cell; zero relative deviation imposed. Covariance templates are conditional estimates.\nAll technical cluster counts grow with k. Lower row is a sufficient rejection condition, not actual BY-ranking power.',ha='center',fontsize=10)
 fig.subplots_adjust(top=.90,bottom=.10,wspace=.30,hspace=.27,right=.87)
 save(fig,'03_true_equivalence_power_heatmaps','15个真实协方差模板、真实等效假设下的最差情景检出率；独立面板扩展与固定簇内扩样不同，下排不是实际BY排序功效。')

def targets():
 f=pd.read_csv(RESULTS/'POWER_TARGET_ATTAINMENT_COUNTS.csv');f=f.loc[f.gamma_common_P.eq(.55)&f.margin.eq(.5)&f.target.eq(.8)];fig,axes=plt.subplots(1,2,figsize=(12,4.8))
 for ax,threshold,title,color in zip(axes,['NOMINAL','BY_RANK1_SUFFICIENT'],['Nominal criterion','Conservative BY-sufficient criterion'],['#4385a3','#a8774d']):
  p=f.loc[f.threshold.eq(threshold)].set_index('k').reindex(K_GRID);ax.bar(range(4),p.n_Wilson_lower_meets_target,color=color)
  for j,r in enumerate(p.itertuples()):ax.text(j,r.n_Wilson_lower_meets_target+.3,str(r.n_Wilson_lower_meets_target),ha='center')
  ax.set(xticks=range(4),xticklabels=['1x','2x','4x','8x'],ylim=(0,16),title=title,ylabel='Templates meeting 80% criterion (of 15)',xlabel='Independent complete design panels')
 fig.suptitle('Widest fixed working margin (+/-0.5 SD): information needed under true equivalence',fontsize=14,y=1.01)
 fig.text(.5,-.015,'Criterion: Wilson 95% lower bound exceeds 0.80 in every one of 3 DGPs. Fixed gamma=0.55.\nThese are conditional scenario counts, not expected new discoveries or unconditional sample-size recommendations.',ha='center',fontsize=9);fig.tight_layout()
 fig.subplots_adjust(bottom=.20)
 save(fig,'04_independent_design_target_attainment','在最宽工作容差下达到80%检出要求的模板数；要求三个DGP的Wilson95下界均达标。')

def expansion():
 f=pd.read_csv(RESULTS/'CONDITIONAL_PRECISION_EXPANSION_SCENARIOS.csv');f=f.loc[f.in_family34&f.margin.eq(.5)&f.threshold.eq('BY_RANK1_SUFFICIENT')];fig,axes=plt.subplots(1,2,figsize=(12.5,4.9))
 labels=['Structure / point / floor fail','Already sufficient','Needs <=2x','Needs <=4x','Needs <=8x','Needs >8x'];colors=['#aeb8c0','#6ca27c','#75aeb7','#648eb4','#aa8db5','#cf926b']
 for floor,ax,title in [(0.,axes[0],'All variance can fall with independent information'),(.25,axes[1],'Illustrative 25% fixed-cluster variance floor')]:
  part=f.loc[f.fixed_cluster_variance_fraction_scenario_not_estimated.eq(floor)];counts=np.zeros(6,int)
  for r in part.itertuples():
   k=r.required_subject_multiplier_in_floor_scenario
   if not r.point_structure_compatible or not r.point_inside_both_margins or not np.isfinite(k):counts[0]+=1
   else:counts[1 if k<=1 else 2 if k<=2 else 3 if k<=4 else 4 if k<=8 else 5]+=1
  ax.barh(range(6),counts,color=colors);ax.set(yticks=range(6),yticklabels=labels,title=title,xlim=(0,35));ax.invert_yaxis()
  for i,n in enumerate(counts):ax.text(n+.3,i,str(n),va='center')
 fig.suptitle('Current point held fixed: conditional precision scenarios, NOT power',fontsize=14,y=1.02)
 fig.text(.5,-.01,'Fixed margin +/-0.5 SD. The 25% floor is an assumption, not estimated from the data.\nIncreasing people within existing clusters cannot be assumed to provide the same information gain as new independent panels.',ha='center',fontsize=9);fig.tight_layout()
 fig.subplots_adjust(bottom=.20)
 save(fig,'05_precision_expansion_assumptions','保持当前点估计的条件精度情景，比较全部方差随独立信息下降与假定25%簇方差下限；不是事后功效或样本量承诺。')

def mrm():
 fig,axes=plt.subplots(1,2,figsize=(13,5.2));labels=['Released measurement rows','Resolved technical records','Unique baseline subjects','Same-visit date candidates in S5','NPTX1 conditional valid cases','KLK6 conditional valid cases','Verified MRM specimen dates'];values=[303,303,287,149,144,141,0]
 axes[0].barh(range(7),values,color=['#527e98']*3+['#a18e65']*3+['#b55a5a']);axes[0].set(yticks=range(7),yticklabels=labels,title='2014 audit: identity closed, collection dates open',xlim=(0,340));axes[0].invert_yaxis()
 for i,n in enumerate(values):axes[0].text(n+3,i,str(n),va='center')
 dfs=[24,23,20,78,75];lab=['2019 NPTX1 peptide1','2019 NPTX1 peptide2','2019 KLK6 peptides','2014 NPTX1 (inferred date)','2014 KLK6 (inferred date)'];axes[1].barh(range(5),dfs,color=['#879eae']*3+['#bca068']*2);axes[1].axvline(40,color='#a13d4b',linestyle='--',label='Frozen residual df >=40')
 axes[1].set(yticks=range(5),yticklabels=lab,title='Potential design gain requires date verification',xlabel='Residual degrees of freedom',xlim=(0,100));axes[1].invert_yaxis();axes[1].legend(frameon=False,fontsize=9,loc='lower right')
 fig.text(.5,-.01,'2014: published-value precision reconciled over all 320 peptides; no new candidate associations fitted.\nInferred visit dates are not verified source-specimen dates. Conditional geometry does not establish inferential validity.',ha='center',fontsize=9);fig.tight_layout()
 save(fig,'06_2014_validation_gap_and_potential','2014技术映射闭合、日期缺口未闭合；潜在有效配对141–144和df75–78均以日期映射成立为条件，尚未正式推断。')

def external():
 f=pd.read_csv(RESULTS/'EXTERNAL_EMORY2022_FILTERED_GENE_PLATFORM_COVERAGE.csv');f=f.loc[f.sheet.eq('4.Overlap.CSF.MSdepleted')].set_index('gene').reindex(['PDGFRB','NPTX1','KLK6']);cols=['SomaScan_filtered_gene_present','Olink_filtered_gene_present','MS_filtered_gene_present'];M=f[cols].to_numpy(int);fig,axes=plt.subplots(1,2,figsize=(12.2,4.5),gridspec_kw={'width_ratios':[1,1.5]});ax=axes[0];ax.imshow(M,cmap='Blues',vmin=0,vmax=1,aspect='auto')
 for i in range(3):
  for j in range(3):ax.text(j,i,'Listed' if M[i,j] else 'Not listed',ha='center',va='center',color='white' if M[i,j] else '.2')
 ax.set(xticks=range(3),xticklabels=['SomaScan','Olink','TMT-MS'],yticks=range(3),yticklabels=M.shape[0] and list(f.index),title='External Emory 2022: filtered gene coverage')
 axes[1].axis('off');axes[1].text(0,.92,'What is now documented',fontsize=14,fontweight='bold');axes[1].text(0,.78,'PDGFRB and KLK6: 3-platform annotation\nNPTX1: Soma/MS; Olink not listed\nDepleted and undepleted MS lists agree',va='top',fontsize=11,linespacing=1.7)
 axes[1].text(0,.42,'What remains unresolved',fontsize=14,fontweight='bold');axes[1].text(0,.30,'Three-platform cohort upper bound: 35\nPer-gene complete n and peptide identity unverified\nNo NULISA; no external individual-level reanalysis',va='top',fontsize=11,linespacing=1.7)
 fig.text(.5,.015,'Source: Dammer et al., Alzheimer\'s Research & Therapy (2022), Supplementary Tables 4/5, DOI 10.1186/s13195-022-01113-5.\nPublished annotation is an opportunity for measurement support; it is not a completed replication result.',ha='center',fontsize=9);fig.subplots_adjust(wspace=.3,bottom=.2,top=.9)
 save(fig,'07_external_coverage_opportunity','公开Emory2022过滤后基因级覆盖；PDGFRB/KLK6三平台线索与NPTX1两平台线索，不等于个体完整度或正式外部复现。')

if __name__=='__main__':
 bottlenecks();intervals();power_heatmaps();targets();expansion();mrm();external();csvout(RESULTS/'FIGURE_INDEX.csv',records);log('Seven scientific figures generated as PNG/vectorPDF pairs')
