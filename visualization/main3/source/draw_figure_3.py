"""Distinct cognitive-result views: evidence plane, cross-outcome map, assay profiles, log-q paths."""
from pathlib import Path
import json,logging,warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.text import Text
from PIL import Image
BASE=Path(__file__).resolve().parents[1]
DATA=json.loads((BASE/'audit/PLOT_INPUT.json').read_text(encoding='utf-8'))
W,H=180,230
INK='#293C4B';GRID='#E7EDF0';ROSE='#A45D80'
FILLS={'NULISA':'#D4977D','SomaScan':'#76ACA0','TMT-MS':'#92A5C3'}
EDGES={'NULISA':'#A3674F','SomaScan':'#3D7F73','TMT-MS':'#5B759C'}
SHAPES={'NULISA':'o','SomaScan':'s','TMT-MS':'D'}
PROTEIN_COLORS={'NPTX2':'#A2577D','TAFA5':'#6D6A9D','TIMP3':'#A28B56'}
PROTEIN_SHAPES={'NPTX2':'o','TAFA5':'^','TIMP3':'s'}
logging.getLogger('fontTools').setLevel(logging.ERROR)
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','font.family':'Arial',
 'figure.facecolor':'white','savefig.facecolor':'white','axes.unicode_minus':False,'axes.linewidth':.65,
 'xtick.major.size':2.7,'ytick.major.size':2.7,'xtick.major.width':.65,'ytick.major.width':.65})
T={
 'EN':dict(A='Shared memory and assay differences',B='NPTX2: secondary cognition',C='Memory associations adjusted for p-tau181',
      D='Memory differences in sensitivity analyses',shared='Shared p-tau',other='Other proteins',
      ax='Shared-memory evidence\n−log$_{10}$ q',ay='Between-assay difference evidence\n−log$_{10}$ q',
      none='No proteins met\nshared-direction\nq ≤ 0.05',overlap='24 proteins\ncoincide',
      bx='Baseline ADNI-MEM association\nper protein residual SD',by='Annual ADNI-EF change\nper protein residual SD',
      baseline='Baseline',ef='EF change',cy='ADNI-MEM units/year\nper protein residual SD',diff='Difference q',
      dy='Difference q',settings=['Primary\n68 tests','All outcomes\n204 tests','Longer follow-up\n28 tests','Equal weighting\n28 tests']),
 'ZH':dict(A='共同记忆关联与跨检测差异',B='NPTX2：次要认知结局',C='经p-tau181调整的记忆关联',D='敏感性分析中的记忆关联差异',
      shared='共同p-tau关联',other='其他蛋白',ax='共同记忆关联证据\n−log$_{10}$ q',ay='跨检测差异证据\n−log$_{10}$ q',
      none='没有蛋白达到\n共同方向\nq≤0.05',overlap='24个蛋白\n位置重合',
      bx='基线ADNI-MEM关联\n每蛋白残差SD',by='年度ADNI-EF变化\n每蛋白残差SD',baseline='基线记忆',ef='执行功能变化',
      cy='ADNI-MEM分/年\n每蛋白残差SD',diff='差异q',dy='差异q',
      settings=['主要分析\n68项检验','全部认知结局\n204项检验','较长随访\n28项检验','个体等权\n28项检验'])}

class Figure:
    def __init__(self,lang):
        self.lang=lang;self.t=T[lang];self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=180);self.axes={};self.primitives={}
    def font(self,size=8,bold=False):
        name=('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN' else ('msyhbd.ttc' if bold else 'msyh.ttc')
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/name),family='Arial' if self.lang=='EN' else 'Microsoft YaHei',size=size,weight='bold' if bold else 'normal')
    def text(self,x,y,s,size=8,bold=False,ha='left',color=INK,rotation=0):
        return self.fig.text(x/W,1-y/H,s,fontproperties=self.font(size,bold),color=color,ha=ha,va='center',linespacing=1.1,rotation=rotation,zorder=10)
    def panel(self,letter,title,x,y):
        self.text(x,y,letter,11.5,True);self.text(x+7,y,title,9.0,True)
    def axis(self,key,x,y,w,h):
        ax=self.fig.add_axes([x/W,1-(y+h)/H,w/W,h/H]);self.axes[key]=ax
        ax.spines[['top','right']].set_visible(False);ax.spines[['bottom','left']].set_color('#708692');ax.tick_params(colors=INK,pad=3,direction='out')
        return ax
    def ticks(self,ax,size=7.4):
        for s in ax.get_xticklabels()+ax.get_yticklabels():s.set_fontproperties(self.font(size));s.set_color(INK)
        ax.xaxis.get_offset_text().set_fontproperties(self.font(size));ax.yaxis.get_offset_text().set_fontproperties(self.font(size))
    def evidence_plane(self):
        self.panel('A',self.t['A'],3,7)
        for x,label,marker,color,filled in [(18,self.t['shared'],'D',ROSE,True),(56,self.t['other'],'o','#95AAB7',False)]:
            self.fig.lines.append(Line2D([x/W],[1-16.1/H],transform=self.fig.transFigure,marker=marker,markersize=3.8,
                markerfacecolor=color if filled else 'white',markeredgecolor=color,markeredgewidth=.7,linestyle='none'))
            self.text(x+3,16.1,label,7.4)
        ax=self.axis('A',16,25,66,49);cut=DATA['evidence_threshold']
        ax.set(xlim=(-.06,1.9),ylim=(-.10,3.5),xticks=[0,.5,1,1.5],yticks=[0,1,2,3])
        ax.axvspan(cut,1.9,color='#EDF5F1',alpha=.75,zorder=0);ax.axhspan(cut,3.5,color='#F8EEE7',alpha=.55,zorder=0)
        ax.axvline(cut,color='#9A8B95',linestyle=(0,(3,2)),linewidth=.75);ax.axhline(cut,color='#9A8B95',linestyle=(0,(3,2)),linewidth=.75)
        rows=[r for r in DATA['memory_evidence'] if r['Qualified']]
        for cls in [False,True]:
            group=[r for r in rows if r['Shared p-tau181']==cls]
            ax.scatter([r['Shared evidence'] for r in group],[r['Difference evidence'] for r in group],
                s=19 if cls else 22,marker='D' if cls else 'o',facecolor=ROSE if cls else 'white',
                edgecolors=ROSE if cls else '#91A6B3',linewidths=.7,zorder=5 if cls else 4)
        labels={'NPTX2':(.74,2.72),'TAFA5':(.13,3.18),'TIMP3':(.14,1.69)}
        for gene,pos in labels.items():ax.text(*pos,gene,fontproperties=self.font(7.8,True),color=ROSE,ha='left',va='center')
        ax.text(1.58,2.24,self.t['none'],fontproperties=self.font(7.2),color='#4F6E63',ha='center',va='center',linespacing=1.15)
        ax.annotate(self.t['overlap'],xy=(.015,.015),xytext=(.80,.39),textcoords='data',fontproperties=self.font(7.1),
            color=INK,ha='center',va='center',linespacing=1.12,arrowprops={'arrowstyle':'-','color':'#9CACB6','lw':.65},zorder=6)
        self.ticks(ax);self.text(49,85.4,self.t['ax'],7.6,False,'center');self.text(4.7,49.5,self.t['ay'],7.5,False,'center',rotation=90)
        self.primitives['A']={'points':rows,'withheld':['CST3'],'formal_family':68,'origin_overlap_count':24,
            'x_limits':[-.06,1.9],'y_limits':[-.1,3.5],'threshold':cut,'no_jitter':True,'no_shared_discoveries':True}
    def secondary(self):
        self.panel('B',self.t['B'],95,7)
        for x,assay in zip([99,127,157],DATA['assays']):
            self.fig.lines.append(Line2D([x/W],[1-16.1/H],transform=self.fig.transFigure,marker=SHAPES[assay],markersize=3.8,
                markerfacecolor=FILLS[assay],markeredgecolor=EDGES[assay],markeredgewidth=.6,linestyle='none'))
            self.text(x+3,16.1,assay,7.4,True,color=EDGES[assay])
        ax=self.axis('B',110,25,66,49)
        ax.set(xlim=(0,.62),ylim=(-.01,.18),xticks=[0,.2,.4,.6],yticks=[0,.05,.10,.15])
        ax.axhline(0,color='#9AAAB5',linestyle=(0,(3,2)),linewidth=.7);ax.axvline(0,color='#9AAAB5',linestyle=(0,(3,2)),linewidth=.7)
        for r in DATA['secondary_map']:
            assay=r['Assay'];x=r['Baseline effect'];y=r['EF change effect']
            ax.errorbar(x,y,xerr=[[x-r['Baseline CI lower']],[r['Baseline CI upper']-x]],
                yerr=[[y-r['EF change CI lower']],[r['EF change CI upper']-y]],fmt=SHAPES[assay],markersize=4.8,
                mfc=FILLS[assay],mec=EDGES[assay],mew=.65,color=EDGES[assay],elinewidth=.9,capsize=2.0,capthick=.7,zorder=4)
        qbase=next(r['q'] for r in DATA['secondary_tests'] if r['Outcome']=='Baseline memory' and r['Test']=='Assay differences')
        qef=next(r['q'] for r in DATA['secondary_tests'] if r['Outcome']=='Executive function change' and r['Test']=='Assay differences')
        ax.text(.035,.172,self.t['diff']+'\n'+self.t['baseline']+f': {qbase:.3f}\n'+self.t['ef']+f': {qef:.4f}',
             fontproperties=self.font(7.3),ha='left',va='top',color=INK,linespacing=1.2)
        self.ticks(ax);self.text(143,85.4,self.t['bx'],7.6,False,'center');self.text(98,49.5,self.t['by'],7.4,False,'center',rotation=90)
        self.primitives['B']={'points':DATA['secondary_map'],'tests':DATA['secondary_tests'],'family':136,
          'horizontal_and_vertical_intervals_are_marginal':True,'joint_ellipse_not_estimated':True,'x_limits':[0,.62],'y_limits':[-.01,.18]}
    def profiles(self):
        self.panel('C',self.t['C'],3,99)
        self.text(4.5,142,self.t['cy'],7.5,False,'center',rotation=90)
        for x,gene in zip([16,76,136],DATA['genes']):
            r=[r for r in DATA['memory_effects'] if r['Protein']==gene]
            n=r[0]['Participants'];q=DATA['primary_diff_tests'][gene]['Assay difference q']
            self.text(x+20.5,108,gene+f' · n={n}',8.2,True,'center')
            self.text(x+20.5,114.7,self.t['diff']+f'={q:.3g}',7.4,False,'center')
            ax=self.axis('C_'+gene,x,121,41,42)
            ax.set(xlim=(-.55,2.55),ylim=(-.05,.18),xticks=[0,1,2],yticks=[-.05,0,.05,.10,.15])
            ax.set_xticklabels(DATA['assays']);ax.axhline(0,color='#8C9DA8',linestyle=(0,(3,2)),linewidth=.75,zorder=1)
            for tick in [.05,.10,.15]:ax.axhline(tick,color=GRID,linewidth=.5,zorder=0)
            byassay={a:next(v for v in r if v['Assay']==a) for a in DATA['assays']}
            betas=[byassay[a]['Effect'] for a in DATA['assays']]
            ax.plot([0,1,2],betas,color='#AAB6BF',linewidth=.8,zorder=2)
            for i,assay in enumerate(DATA['assays']):
                v=byassay[assay];b=v['Effect']
                ax.errorbar(i,b,yerr=[[b-v['CI lower']],[v['CI upper']-b]],fmt=SHAPES[assay],markersize=4.4,
                    mfc=FILLS[assay],mec=EDGES[assay],mew=.65,color=EDGES[assay],elinewidth=.9,capsize=2.4,capthick=.75,zorder=4)
            self.ticks(ax,7.3)
            if gene!='NPTX2':ax.set_yticklabels([])
            self.primitives['C_'+gene]={'points':r,'joint_difference_q':q,'family':68,'x_is_categorical_not_time':True,'same_y_limits':[-.05,.18]}
    def robustness(self):
        self.panel('D',self.t['D'],3,180)
        ax=self.axis('D',16,191,161,27)
        ax.set(xlim=(-.15,3.75),ylim=(.00075,1),yscale='log',xticks=[0,1,2,3],yticks=[.001,.01,.1,1])
        ax.set_xticklabels(self.t['settings']);ax.set_yticklabels(['0.001','0.01','0.1','1'])
        ax.minorticks_off();ax.axhline(.05,color='#8D9AA6',linewidth=.8,linestyle=(0,(3,2)),zorder=1)
        for tick in [.001,.01,.1]:ax.axhline(tick,color=GRID,linewidth=.5,zorder=0)
        paths=[]
        for gene in DATA['genes']:
            r=[next(v for v in DATA['robustness'] if v['Protein']==gene and v['Analysis']==a) for a in DATA['robustness_analyses']]
            values=[v['q'] for v in r]
            ax.plot(range(4),values,color=PROTEIN_COLORS[gene],linewidth=1.1,marker=PROTEIN_SHAPES[gene],markersize=3.8,
                    markerfacecolor=PROTEIN_COLORS[gene],markeredgecolor='white',markeredgewidth=.4,zorder=4)
            ax.text(3.16,values[-1],gene,fontproperties=self.font(7.7,True),color=PROTEIN_COLORS[gene],ha='left',va='center')
            paths.append({'protein':gene,'values':r})
        self.ticks(ax,7.3);self.text(5.0,204.5,self.t['dy'],7.5,False,'center',rotation=90)
        self.primitives['D']={'paths':paths,'cutoff_q':.05,'log_limits':[.00075,1],
           'family_sizes':[68,204,28,28],'q_changes_not_formal_effect_change_tests':True,'protein_colours':PROTEIN_COLORS}
    def layout(self):
        self.fig.canvas.draw();renderer=self.fig.canvas.get_renderer();records=[];seen=set();overlap=[]
        for a in self.fig.findobj(Text):
            if not a.get_visible() or not a.get_text().strip() or id(a) in seen:continue
            seen.add(id(a));bb=a.get_window_extent(renderer);v=np.asarray([bb.x0,bb.y0,bb.x1,bb.y1])/self.fig.dpi*25.4
            row={'text':a.get_text(),'font_pt':float(a.get_fontsize()),'x0_mm':float(v[0]),'y0_mm':float(H-v[3]),'x1_mm':float(v[2]),'y1_mm':float(H-v[1])}
            assert row['x0_mm']>=.1 and row['x1_mm']<=W-.1 and row['y0_mm']>=.1 and row['y1_mm']<=H-.1,('Text clipping',self.lang,row)
            records.append(row)
        for i,a in enumerate(records):
            for b in records[i+1:]:
                xx=min(a['x1_mm'],b['x1_mm'])-max(a['x0_mm'],b['x0_mm']);yy=min(a['y1_mm'],b['y1_mm'])-max(a['y0_mm'],b['y0_mm'])
                if xx>.22 and yy>.22:overlap.append([a['text'],b['text'],round(xx,2),round(yy,2)])
        assert not overlap,('Text overlap',self.lang,overlap)
        return records
    def save(self):
        self.evidence_plane();self.secondary();self.profiles();self.robustness();labels=self.layout()
        title='Figure 3. Cognitive associations of CSF proteins across three assays'
        self.fig.savefig(BASE/f'Figure_3_{self.lang}.pdf',metadata={'Title':title,'Creator':'Paired ADNI CSF proteomic study'})
        self.fig.savefig(BASE/f'editable/Figure_3_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_3_{self.lang}.png',dpi=300);self.fig.savefig(BASE/f'Figure_3_{self.lang}_preview.png',dpi=180)
        tif=BASE/f'high_resolution/Figure_3_{self.lang}.tiff';self.fig.savefig(tif,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        with Image.open(tif) as im:
            rgb=im.convert('RGB');rgb.save(tif,compression='tiff_lzw',dpi=(600,600));rgb.close()
        (BASE/f'audit/PLOT_PRIMITIVES_{self.lang}.json').write_text(json.dumps(self.primitives,ensure_ascii=False,indent=2),encoding='utf-8')
        (BASE/f'audit/LAYOUT_{self.lang}.json').write_text(json.dumps(labels,ensure_ascii=False,indent=2),encoding='utf-8')
        return {'language':self.lang,'size_mm':[W,H],'panels':4,'minimum_font_pt':min(v['font_pt'] for v in labels),
          'text_overlap':False,'text_clipped':False,'evidence_points':33,'memory_profiles':3,'secondary_assay_points':3}

def main():
    warnings.filterwarnings('error',message='Glyph.*missing');outputs=[]
    for lang in ['EN','ZH']:
        f=Figure(lang);outputs.append(f.save());plt.close(f.fig)
    (BASE/'audit/LAYOUT_VALIDATION.json').write_text(json.dumps(outputs,ensure_ascii=False,indent=2),encoding='utf-8')
    (BASE/'audit/STYLE_SPECIFICATION.json').write_text(json.dumps({'size_mm':[W,H],'assay_fills':FILLS,'assay_edges':EDGES,
         'assay_shapes':SHAPES,'protein_colours_panel_D':PROTEIN_COLORS,'formats':['two-test evidence scatter','cross-outcome confidence crosshairs','categorical assay profiles','log-q sensitivity paths'],
         'no_heatmap_or_bubble_matrix':True,'outside_grey_explanatory_footnotes':False},ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(outputs))
if __name__=='__main__':main()
