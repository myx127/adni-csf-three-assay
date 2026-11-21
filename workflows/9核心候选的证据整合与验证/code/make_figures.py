from config import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch
from matplotlib import font_manager
fonts=[Path('C:/Windows/Fonts/msyh.ttc'),Path('C:/Windows/Fonts/simhei.ttf')]
for p in fonts:
    if p.exists():font_manager.fontManager.addfont(str(p));plt.rcParams['font.family']=font_manager.FontProperties(fname=str(p)).get_name();break
plt.rcParams.update({'axes.unicode_minus':False,'pdf.fonttype':42,'ps.fonttype':42,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
COLORS={'M0':'#145DA0','APOE':'#D17B0F','BOTH':'#2A8C68'}
INDEX=[]
def save(fig,name,description):
    fig.savefig(FIGURES/(name+'.png'),dpi=220,bbox_inches='tight',facecolor='white')
    fig.savefig(FIGURES/(name+'.pdf'),bbox_inches='tight',facecolor='white')
    INDEX.append({'figure':name,'description':description,'PNG':str((FIGURES/(name+'.png')).relative_to(ROOT)),'PDF':str((FIGURES/(name+'.pdf')).relative_to(ROOT))})
    plt.close(fig)

def run():
    family=pd.read_csv(RESULTS/'FAMILY34_INTEGRATED_EVIDENCE_ALL_MODELS.csv')
    names=['NPTX1','KLK6','PDGFRB']
    order=names+sorted(set(family.loc[family.shared_supported,'gene'])-set(names))+sorted(set(family.gene)-set(names)-set(family.loc[family.shared_supported,'gene']))
    states=['BOTH_UNRESOLVED','SHARED_ONLY_H_UNRESOLVED','H_ONLY_SHARED_UNRESOLVED','SHARED_AND_H_DEVIATION']
    colors=['#e4e7ea','#7dc5d4','#f2bf76','#8870b5']
    A=np.array([[states.index(family.loc[family.gene.eq(g)&family.model.eq(m),'evidence_state'].iloc[0]) for m in MODELS] for g in order])
    fig,ax=plt.subplots(figsize=(8,11.5));ax.imshow(A,cmap=ListedColormap(colors),vmin=-.5,vmax=3.5,aspect='auto')
    ax.set_yticks(range(34),order);ax.set_xticks(range(3),['主模型 M0','＋APOE','＋教育/APOE'])
    symbols=['未定','共同关联','反应偏离','两者并存']
    for i in range(34):
        for j in range(3):ax.text(j,i,symbols[A[i,j]],ha='center',va='center',fontsize=9,color='white' if A[i,j]==3 else '#202020')
    for tick in ax.get_yticklabels():
        if tick.get_text() in names:tick.set_fontweight('bold')
    ax.set_title('34蛋白：共同病理关联与相对反应偏离分别评价\n各模型均沿用固定BY34；“未定”不等于等效',pad=14)
    ax.set_yticks(np.arange(-.5,34,1),minor=True);ax.grid(which='minor',axis='y',color='white',lw=1.5);ax.tick_params(which='minor',left=False)
    fig.subplots_adjust(bottom=.1)
    ax.legend(handles=[Patch(facecolor=c,label=l) for c,l in zip(colors,symbols)],loc='upper center',bbox_to_anchor=(.5,-.05),ncol=4,frameon=False)
    save(fig,'01_family34_integrated_evidence','34蛋白三模型的两类证据；不以未拒绝判等效')
    effects=pd.read_csv(RESULTS/'INDEPENDENT_PLATFORM_EFFECTS.csv')
    fig,axes=plt.subplots(1,3,figsize=(13,6.2),sharey=True)
    for ax,gene in zip(axes,names):
        for model,offset in zip(MODELS,[-.2,0,.2]):
            p=effects.loc[effects.gene.eq(gene)&effects.model.eq(model)].set_index('assay').loc[['N','S','T']]
            y=np.arange(3)+offset;ax.errorbar(p.scaled_b,y,xerr=[p.scaled_b-p.scaled_CI95_low,p.scaled_CI95_high-p.scaled_b],fmt='o',ms=4,capsize=2,color=COLORS[model],label=model)
        ax.axvline(0,color='grey',lw=.8);ax.set_title(gene,fontweight='bold');ax.grid(axis='x',alpha=.2)
        p=family.loc[family.gene.eq(gene)&family.model.eq('M0')].iloc[0]
        ax.text(.02,-.38,f'M0 同向关联 q={p.q_BY_shared34:.4g}\nM0 反应偏离 q={p.q_BY_H34:.4g}',transform=ax.transAxes,fontsize=10)
        ax.set_xlabel('每1 SD pTau的调整系数 / 残差SD')
    axes[0].set_yticks(range(3),['NULISA','SomaScan','TMT']);axes[0].invert_yaxis()
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.89),ncol=3,frameon=False)
    fig.suptitle('核心候选：同向病理关联可以与相对响应差异并存\n95%区间纳入系数与残差SD的联合不确定度',y=1.02)
    fig.subplots_adjust(wspace=.25,bottom=.3,top=.77)
    save(fig,'02_core_candidates_association_forest','三候选三个模型三检测的标准化病理系数与联合不确定度')
    rep=pd.read_csv(RESULTS/'CORE_PEPTIDE_TECHNICAL_REPEATABILITY.csv')
    peps=['FQLTFPLR','LENLEQYSR','ESSQEQSSVVR','LSELIQPLPLER','YTNWIQK']
    fig,ax=plt.subplots(figsize=(10,6))
    for rel,off,col in [('MRM2014',-.12,'#145DA0'),('MRM2019',.12,'#D17B0F')]:
        d=rep.loc[rep.release.eq(rel)].set_index('peptide').loc[peps]
        y=np.arange(5)+off;ax.errorbar(d.pearson_r,y,xerr=[d.pearson_r-d.pearson_pair_bootstrap95_low,d.pearson_pair_bootstrap95_high-d.pearson_r],fmt='o',capsize=3,color=col,label=rel)
        for yy,n in zip(y,d.n_complete_pairs):ax.text(1.025,yy,f'n={n}',va='center',fontsize=9,color=col)
    ax.set_xlim(-.15,1.13);ax.set_yticks(range(5),['NPTX1 '+p if i<2 else 'KLK6 '+p for i,p in enumerate(peps)]);ax.invert_yaxis();ax.axvline(0,lw=.7,color='grey');ax.axvline(1,lw=.7,color='grey')
    ax.set_xlabel('重复对Pearson相关及成对自助法95%区间（描述性）');ax.set_title('核实的技术重复：2014年16对，2019年20对\n不同重复集合的相关受样本间变异影响，不能直接排名检测准确度')
    ax.legend(loc='lower left',frameon=False);ax.grid(axis='x',alpha=.2)
    save(fig,'03_peptide_technical_repeatability','两发布家族各肽段的技术重复相关及4000次成对自助法区间')
    geom=pd.read_csv(RESULTS/'MRM2014_THREE_FOUR_ASSAY_DESIGN_AUDIT.csv').drop_duplicates(['gene','assay_design'])
    fig,ax=plt.subplots(figsize=(10,5.3))
    geom=geom.set_index(['gene','assay_design']).loc[[('NPTX1','N_S_MRM'),('NPTX1','N_S_TMT_MRM'),('KLK6','N_S_MRM'),('KLK6','N_S_TMT_MRM')]].reset_index()
    y=np.arange(4);ax.barh(y,geom.rank_X,color='#abb7c2',label='被模型设计占用的秩')
    ax.barh(y,geom.residual_df,left=geom.rank_X,color=['#78b9ba' if x>=40 else '#e6a096' for x in geom.residual_df],label='残差自由度')
    for yy,row in zip(y,geom.itertuples()):
        ax.text(row.rank_X/2,yy,f'秩 {row.rank_X}',ha='center',va='center')
        ax.text(row.rank_X+row.residual_df/2,yy,f'df={row.residual_df}',ha='center',va='center')
        ax.text(row.n+2,yy,f'n={row.n}',va='center')
    ax.set_xlim(0,161);ax.set_yticks(y,[g+'\n'+('N/S/MRM' if a=='N_S_MRM' else 'N/S/TMT/MRM') for g,a in zip(geom.gene,geom.assay_design)]);ax.invert_yaxis()
    ax.set_xlabel('观测人数 = 模型秩 + 残差自由度');ax.set_title('2014年MRM：增加TMT后，四检测设计df降至19/17\n冻结df门槛40；全部设计的标本日期和2014推断校准仍未完成')
    ax.legend(loc='lower left',bbox_to_anchor=(0,-.27),ncol=2,frameon=False)
    save(fig,'04_MRM2014_three_four_assay_geometry','三/四检测交集和设计秩独立重算；不把几何信息判为验证通过')
    qc=pd.read_csv(RESULTS/'MRM2019_PQC24_RAW_LH_VARIATION.csv')
    qc=qc.sort_values(['gene','peptide','transition'],ascending=[False,True,True]).reset_index(drop=True)
    fig,ax=plt.subplots(figsize=(10,6.4));y=np.arange(len(qc));vals=qc.raw_LH_CV_percent_all_positive_no_flag_exclusion
    ax.barh(y,vals,color=['#145DA0' if g=='NPTX1' else '#6c9e80' for g in qc.gene])
    ax.set_yticks(y,[g+' '+p+' '+t.split('.')[2] for g,p,t in zip(qc.gene,qc.peptide,qc.transition)]);ax.invert_yaxis()
    for yy,v in zip(y,vals):ax.text(v+.15,yy,f'{v:.1f}%',va='center')
    ax.set_xlim(0,18);ax.set_xlabel('24个过程质控样本：原始Light/Heavy比值CV（%）')
    ax.set_title('恢复原始质控：10个转变的24个PQC峰标志均为0\n95人分析的950个转变测量也均为0（通过原始峰质控）')
    ax.grid(axis='x',alpha=.2)
    save(fig,'05_raw_transition_process_QC','恢复PQC标签后的原始轻重比值CV；不作为临床精度认证')
    inf=pd.read_csv(RESULTS/'NEW_RID_LEAVE_ONE_BATCH_DESCRIPTIVE.csv')
    asc=pd.read_csv(RESULTS/'INDEPENDENT_2019_NEW_RID_PEPTIDE_ASSOCIATIONS.csv').set_index('peptide')
    fig,ax=plt.subplots(figsize=(10,6))
    for i,p in enumerate(peps):
        for axis,offset,color in [('MRM_PLATE',-.13,'#145DA0'),('MRM_RUN',.13,'#D17B0F')]:
            part=inf.loc[inf.peptide.eq(p)&inf.deleted_axis.eq(axis)&inf.status.eq('DESCRIPTIVE_ONLY')]
            ax.scatter(part.scaled_b,np.full(len(part),i+offset),s=25,alpha=.6,color=color,label=('删一检测板' if axis=='MRM_PLATE' else '删一运行') if i==0 else None)
        ax.scatter(asc.loc[p,'b_P_per_residual_SD'],i,c='black',marker='D',s=42,label='完整95人' if i==0 else None,zorder=4)
    ax.axvline(0,color='grey',lw=.8);ax.set_yticks(range(5),[('NPTX1 ' if i<2 else 'KLK6 ')+p for i,p in enumerate(peps)]);ax.invert_yaxis();ax.set_xlim(-.03,.7)
    ax.set_xlabel('调整后pTau系数 / 肽段残差SD（描述性）');ax.set_title('85次逐板/运行删除：五肽段关联方向保持正向\n未选择删除后P值；原BY5结论不变，不能计为新复现')
    ax.legend(loc='lower right',frameon=False);ax.grid(axis='x',alpha=.2)
    save(fig,'06_newRID_batch_influence','固定协变量重建85个删批次设计，仅描述系数方向和变幅')
    regions=pd.read_csv(RESULTS/'CORE_UNIPROT_REGIONS.csv');peptide=pd.read_csv(RESULTS/'CORE_PEPTIDE_SEQUENCE_REGION_SPECIFICITY.csv')
    versions=pd.read_csv(RESULTS/'UNIPROT_SEQUENCE_VERSION_AUDIT.csv').set_index('gene')
    fig,axes=plt.subplots(3,1,figsize=(12,7.5));types={'Signal':('#b7c4cf',.3),'Propeptide':('#e5b467',.3),'Domain':('#8bb9a6',.18),'Transmembrane':('#d77979',.3),'Topological domain':('#bed6e0',.3)}
    for ax,gene in zip(axes,['NPTX1','PDGFRB','KLK6']):
        length=versions.loc[gene,'length_aa'];ax.hlines(0,1,length,color='grey',lw=3)
        for row in regions.loc[regions.gene.eq(gene)&regions.feature_type.isin(types)].itertuples():
            if pd.notna(row.start) and pd.notna(row.end):
                color,height=types[row.feature_type];ax.barh(0,row.end-row.start+1,left=row.start,height=height,color=color,alpha=.9)
        ps=peptide.loc[peptide.gene.eq(gene)]
        for i,row in enumerate(ps.itertuples()):
            start=int(str(row.canonical_positions_1based).split('|')[0]);end=start+len(row.peptide)-1
            ax.plot([start,end],[.35,.35],color='#7a3880',lw=5);ax.annotate(row.peptide,(start,.35),xytext=(0,12+14*(i%2)),textcoords='offset points',fontsize=9,ha='center',rotation=0)
        ax.set_xlim(0,length+5);ax.set_ylim(-.45,1.05);ax.set_yticks([]);ax.set_title(f'{gene} | {int(length)} aa'+(' | 无直接MRM肽段覆盖' if gene=='PDGFRB' else ''),loc='left',fontweight='bold');ax.set_xlabel('规范序列氨基酸位置')
    fig.suptitle('肽段来源已核实；实际表位、完整蛋白形式与活性仍未知',y=1.01)
    fig.tight_layout(h_pad=1.5,rect=(0,.05,1,1))
    fig.legend(handles=[Patch(facecolor='#b7c4cf',label='信号肽'),Patch(facecolor='#e5b467',label='前肽'),Patch(facecolor='#bed6e0',label='拓扑区域'),Patch(facecolor='#8bb9a6',label='结构域'),Patch(facecolor='#d77979',label='跨膜区'),Patch(facecolor='#7a3880',label='MRM肽段')],loc='lower center',ncol=6,frameon=False)
    save(fig,'07_core_molecular_peptide_map','三蛋白规范序列区域及五条靶向肽段位置；不宣称已识别具体分子形式')
    csvout(RESULTS/'FIGURE_INDEX.csv',INDEX)
    log(f'{len(INDEX)} scientific figure groups generated as PNG/PDF')

if __name__=='__main__':
    assert_unsealed()
    run()
