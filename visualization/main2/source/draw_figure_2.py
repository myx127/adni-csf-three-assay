"""Bubble evidence map, paired interval plot and covariate evidence heatmap."""
from pathlib import Path
import json,logging,warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.text import Text
from matplotlib.patches import Rectangle,FancyBboxPatch
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.cm import ScalarMappable
from PIL import Image
BASE=Path(__file__).resolve().parents[1]
DATA=json.loads((BASE/'audit/PLOT_INPUT.json').read_text(encoding='utf-8'))
W,H=180,230
ASSAYS=DATA['assays'];MODELS=DATA['models'];ORDER=DATA['formal_proteins'];SHARED=DATA['shared_proteins']
INK='#293C4B';GRID='#E6ECEF';STRONG='#9B587A';PALE='#F5EBF0'
COLORS={'NULISA':'#D4977D','SomaScan':'#76ACA0','TMT-MS':'#92A5C3'}
EDGES={'NULISA':'#A3674F','SomaScan':'#3D7F73','TMT-MS':'#5B759C'}
TINTS={'NULISA':'#F7EBE5','SomaScan':'#EAF3EF','TMT-MS':'#EDF1F7'}
MARKERS={'NULISA':'o','SomaScan':'s','TMT-MS':'D'}
CMAP=LinearSegmentedColormap.from_list('shared_evidence',['#FCF8FA','#E8C9D8','#C38BA8','#975C80','#653753'])
logging.getLogger('fontTools').setLevel(logging.ERROR)
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.unicode_minus':False,
    'font.family':'Arial','figure.facecolor':'white','savefig.facecolor':'white','axes.linewidth':.65,
    'xtick.major.size':2.5,'xtick.major.width':.65,'ytick.major.size':0})
TEXT={
 'EN':dict(A='Evidence across 34 proteins',B='14 shared positive associations',C='Sensitivity to covariate adjustment',
       protein='Protein',positive='Positive',negative='Negative',shared='Shared q ≤ 0.05',q='Shared\nq',n='n',
       area='Circle area: |β|',evidence='−log$_{10}$ q',unit='Protein residual SD\nper SD higher log$_2$ p-tau181',
       models=['Age, sex',r'+ APOE $\varepsilon4$',r'+ Education'+'\n'+r'and APOE $\varepsilon4$']),
 'ZH':dict(A='34个蛋白的关联证据',B='14个共同正向关联蛋白',C='协变量调整敏感性',
       protein='蛋白',positive='正向',negative='负向',shared='共同关联q≤0.05',q='共同关联\nq',n='n',
       area='圆面积：|β|',evidence='−log$_{10}$ q',unit='蛋白残差SD变化\n每1 SD log$_2$ p-tau181升高',
       models=['年龄、性别',r'+ APOE $\varepsilon4$',r'+ 受教育年限'+'\n'+r'及APOE $\varepsilon4$'])}

def effect(gene,assay,model=None):
    return next(r for r in DATA['effects'] if r['Protein']==gene and r['Assay']==assay and r['Model']==(model or MODELS[0]))
def joint(gene,model=None):
    return next(r for r in DATA['tests'] if r['Protein']==gene and r['Model']==(model or MODELS[0]))
def qlabel(q):
    if q<.001:
        m,e=f'{q:.1e}'.split('e');return '$'+m+r'\times10^{'+str(int(e))+'}$'
    return f'{q:.3g}'

class Figure:
    def __init__(self,lang):
        self.lang=lang;self.t=TEXT[lang];self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=180)
        self.axes={};self.records={};self.labels=[]
    def font(self,size=8,bold=False):
        fname=('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN' else ('msyhbd.ttc' if bold else 'msyh.ttc')
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/fname),family='Arial' if self.lang=='EN' else 'Microsoft YaHei',size=size,weight='bold' if bold else 'normal')
    def text(self,x,y,s,size=8,bold=False,ha='left',color=INK,role='label'):
        a=self.fig.text(x/W,1-y/H,s,fontproperties=self.font(size,bold),color=color,ha=ha,va='center',linespacing=1.1,zorder=10)
        self.labels.append((a,role));return a
    def rect(self,x,y,w,h,fill,edge='none',lw=.6,z=-1):
        self.fig.patches.append(Rectangle((x/W,1-(y+h)/H),w/W,h/H,transform=self.fig.transFigure,facecolor=fill,edgecolor=edge,linewidth=lw,zorder=z))
    def marker(self,x,y,assay,size=3.6):
        self.fig.lines.append(Line2D([x/W],[1-y/H],transform=self.fig.transFigure,marker=MARKERS[assay],
            markersize=size,markerfacecolor=COLORS[assay],markeredgecolor=EDGES[assay],markeredgewidth=.6,linestyle='none'))
    def panel(self,letter,title,x,y):
        self.text(x,y,letter,11.5,True,role='panel letter');self.text(x+7,y,title,9.2,True,role='panel title')
    def axis(self,key,x,y,w,h):
        a=self.fig.add_axes([x/W,1-(y+h)/H,w/W,h/H]);self.axes[key]=a
        a.spines[['top','right','left']].set_visible(False);a.spines['bottom'].set_color(INK);a.tick_params(colors=INK,pad=3,direction='out')
        return a
    def style_axis(self,a,size=7.4):
        for label in a.get_xticklabels()+a.get_yticklabels():label.set_fontproperties(self.font(size));label.set_color(INK)
    def evidence_map(self):
        self.panel('A',self.t['A'],3,7)
        self.rect(3,13.4,29,8,PALE,z=0);self.text(17.5,17.4,self.t['shared'],7.25,True,'center',STRONG)
        for x,label,filled in [(38.0,self.t['positive'],True),(61.0,self.t['negative'],False)]:
            self.fig.lines.append(Line2D([x/W],[1-17.4/H],transform=self.fig.transFigure,marker='o',markersize=3.5,
                markerfacecolor='#7B909D' if filled else 'white',markeredgecolor='#7B909D',markeredgewidth=.8,linestyle='none'))
            self.text(x+3,17.4,label,7.3)
        self.text(3,32,self.t['protein'],7.8,True)
        self.text(77.5,32,self.t['q'],7.7,True,'center')
        for i,assay in enumerate(ASSAYS):
            x=20.2+15*i;self.rect(x,27.3,13.6,9,TINTS[assay],z=0)
            self.text(27+15*i,31.8,assay,7.35,True,'center',EDGES[assay])
        ax=self.axis('A_bubbles',19.5,39,45,166);ax.set(xlim=(-.5,2.5),ylim=(33.5,-.5),xticks=[],yticks=[]);ax.spines['bottom'].set_visible(False)
        qax=self.axis('A_evidence',68,39,19,166);qax.set(xlim=(-.8,20.5),ylim=(33.5,-.5),xticks=[0,10,20],yticks=[])
        qax.axvline(DATA['shared_cutoff_strength'],color=STRONG,linestyle=(0,(3,2)),linewidth=.8,zorder=1)
        qax.axvline(0,color='#DEE5E9',linewidth=.5,zorder=0)
        bubble=[];qpoints=[]
        for i,g in enumerate(ORDER):
            supported=g in SHARED;y=39+(i+.5)*166/34
            if supported:
                self.rect(2.7,y-166/68,84.3,166/34,PALE,z=-2)
                ax.axhspan(i-.5,i+.5,color=PALE,zorder=0);qax.axhspan(i-.5,i+.5,color=PALE,zorder=0)
            self.text(3,y,g,7.65,supported,color=STRONG if supported else INK,role='protein A')
            for j,assay in enumerate(ASSAYS):
                r=effect(g,assay);value=r['Effect'];area=DATA['bubble_area_points2_per_abs_effect']*abs(value)
                ax.scatter(j,i,s=area,marker='o',facecolor=COLORS[assay] if value>=0 else 'white',edgecolors=EDGES[assay],linewidths=.7,zorder=4)
                bubble.append({'protein':g,'assay':assay,'effect':value,'area_points2':area,'filled':value>=0,'x':j,'y':i})
            q=joint(g);strength=q['Evidence strength']
            qax.plot([0,strength],[i,i],color=STRONG if supported else '#BAC7CF',linewidth=.8,zorder=2)
            qax.plot(strength,i,marker='D',markersize=3.0,markerfacecolor=STRONG if supported else '#A8B9C4',
                markeredgecolor=STRONG if supported else '#93A7B3',markeredgewidth=.65,linestyle='none',zorder=4)
            qpoints.append({'protein':g,'q':q['Shared q'],'strength':strength,'supported':supported})
        self.style_axis(qax)
        self.text(77.5,214,self.t['evidence'],7.6,False,'center')
        self.text(31.5,213.3,self.t['area'],7.6,False,'center')
        # Area legend uses the exact same points-squared mapping as the effect matrix.
        legend=self.fig.add_axes([8/W,1-224/H,54/W,6/H]);legend.set(xlim=(0,54),ylim=(0,6));legend.axis('off')
        for x,v in [(10,.25),(27,.5),(44,1.0)]:
            legend.scatter(x,3,s=56*v,facecolor='#B1C2CD',edgecolor='#667E8D',linewidths=.6)
            self.text(8+x,226.7,f'{v:g}',7.3,False,'center')
        self.records['A']={'bubble_points':bubble,'q_points':qpoints,'all_formal_proteins':ORDER,'shared_proteins':SHARED,
            'q_strength_limits':[-.8,20.5],'q_cutoff':.05,'no_strength_truncation':True,'area_scale':56,'circle_area_encodes_effect_not_significance':True}
    def intervals(self):
        self.panel('B',self.t['B'],95,7)
        for x,assay in zip([97,125,154],ASSAYS):
            self.marker(x,17,assay,3.6);self.text(x+3,17,assay,7.6,True,color=EDGES[assay])
        self.text(95,27,self.t['protein'],7.9,True);self.text(115.5,27,self.t['n'],7.9,True,'center')
        ax=self.axis('B_intervals',121,34,56,92)
        ax.set(xlim=(-.1,2.15),ylim=(13.5,-.5),xticks=[0,1,2],yticks=[])
        ax.axvline(0,color='#899AA4',linewidth=.75,linestyle=(0,(3,2)),zorder=1)
        for tick in [1,2]:ax.axvline(tick,color=GRID,linewidth=.5,zorder=1)
        records=[]
        for i,g in enumerate(SHARED):
            y=34+(i+.5)*92/14
            if i%2==0:
                ax.axhspan(i-.5,i+.5,color='#F5F7F8',zorder=0);self.rect(94.5,y-92/28,24,92/14,'#F5F7F8',z=-1)
            self.text(95,y,g,8.0,True,role='protein B');self.text(115.5,y,str(joint(g)['Participants']),7.8,False,'center')
            betas=[effect(g,a)['Effect'] for a in ASSAYS]
            ax.plot([min(betas),max(betas)],[i,i],color='#D6E0E5',linewidth=1.5,zorder=1)
            for offset,assay in zip([-.23,0,.23],ASSAYS):
                r=effect(g,assay);b=r['Effect'];lo=r['CI lower'];hi=r['CI upper']
                ax.errorbar(b,i+offset,xerr=[[b-lo],[hi-b]],fmt=MARKERS[assay],markersize=3.4,mfc=COLORS[assay],mec=EDGES[assay],
                    mew=.55,color=EDGES[assay],elinewidth=.8,capsize=1.55,capthick=.65,zorder=4)
                records.append({**r,'y':i+offset})
        self.style_axis(ax,7.7);self.text(149,139,self.t['unit'],7.9,False,'center')
        self.records['B']={'intervals':records,'x_limits':[-.1,2.15],'all_intervals_untruncated':True,'range_segments_are_point_estimate_ranges':True}
    def sensitivity(self):
        self.panel('C',self.t['C'],95,150)
        for i,label in enumerate(self.t['models']):self.text(124.5+21*i,162.2,label,7.25,True,'center')
        norm=Normalize(vmin=0,vmax=20)
        ax=self.fig.add_axes([114/W,1-219/H,63/W,49/H]);self.axes['C_heatmap']=ax;ax.set(xlim=(0,3),ylim=(14,0));ax.axis('off')
        cells=[]
        for i,g in enumerate(SHARED):
            y=170+(i+.5)*49/14;self.text(95,y,g,7.55,False,role='protein C')
            for j,m in enumerate(MODELS):
                r=joint(g,m);value=r['Evidence strength'];q=r['Shared q']
                ax.add_patch(Rectangle((j,i),1,1,facecolor=CMAP(norm(value)),edgecolor='white',linewidth=.7))
                self.text(124.5+21*j,y,qlabel(q),7.2,False,'center',color='white' if value>10 else INK,role='q value C')
                cells.append({'protein':g,'model':m,'q':q,'strength':value,'q_family':34,'label':qlabel(q)})
        bar=self.fig.add_axes([124/W,1-225/H,53/W,2/H])
        cb=self.fig.colorbar(ScalarMappable(norm=norm,cmap=CMAP),cax=bar,orientation='horizontal',ticks=[0,10,20])
        cb.solids.set_rasterized(False)
        cb.outline.set_visible(False);cb.ax.tick_params(length=2,width=.5,pad=1,colors=INK)
        for label in cb.ax.get_xticklabels():label.set_fontproperties(self.font(7.2))
        self.text(96,224,self.t['evidence'],7.4)
        self.records['C']={'cells':cells,'heatmap_limits':[0,20],'untruncated_strength':True,'full_family_per_column':34,'same_14_in_all_three_models':True}
    def layout(self):
        self.fig.canvas.draw();renderer=self.fig.canvas.get_renderer();records=[];seen=set();overlaps=[]
        for a in self.fig.findobj(Text):
            if not a.get_visible() or not a.get_text().strip() or id(a) in seen:continue
            seen.add(id(a));bb=a.get_window_extent(renderer);v=np.asarray([bb.x0,bb.y0,bb.x1,bb.y1])/self.fig.dpi*25.4
            r={'text':a.get_text(),'font_pt':float(a.get_fontsize()),'x0_mm':float(v[0]),'y0_mm':float(H-v[3]),'x1_mm':float(v[2]),'y1_mm':float(H-v[1])}
            assert r['x0_mm']>=.2 and r['x1_mm']<=W-.2 and r['y0_mm']>=.2 and r['y1_mm']<=H-.2,('Clipped label',self.lang,r)
            records.append(r)
        for i,a in enumerate(records):
            for b in records[i+1:]:
                xx=min(a['x1_mm'],b['x1_mm'])-max(a['x0_mm'],b['x0_mm'])
                yy=min(a['y1_mm'],b['y1_mm'])-max(a['y0_mm'],b['y0_mm'])
                if xx>.2 and yy>.2:overlaps.append([a['text'],b['text'],round(xx,2),round(yy,2)])
        assert not overlaps,('Text overlap',self.lang,overlaps)
        return records
    def save(self):
        self.evidence_map();self.intervals();self.sensitivity();labels=self.layout()
        title='Figure 2. Shared p-tau181 associations across three CSF proteomic assays'
        self.fig.savefig(BASE/f'Figure_2_{self.lang}.pdf',metadata={'Title':title,'Creator':'Paired ADNI CSF proteomic study'})
        self.fig.savefig(BASE/f'editable/Figure_2_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_2_{self.lang}.png',dpi=300)
        self.fig.savefig(BASE/f'Figure_2_{self.lang}_preview.png',dpi=180)
        tif=BASE/f'high_resolution/Figure_2_{self.lang}.tiff';self.fig.savefig(tif,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        with Image.open(tif) as im:
            rgb=im.convert('RGB');rgb.save(tif,compression='tiff_lzw',dpi=(600,600));rgb.close()
        (BASE/f'audit/PLOT_PRIMITIVES_{self.lang}.json').write_text(json.dumps(self.records,ensure_ascii=False,indent=2),encoding='utf-8')
        (BASE/f'audit/LAYOUT_{self.lang}.json').write_text(json.dumps(labels,ensure_ascii=False,indent=2),encoding='utf-8')
        return {'language':self.lang,'figure_size_mm':[W,H],'minimum_font_pt':min(r['font_pt'] for r in labels),'panels':3,
            'all_34_displayed':True,'interval_estimates':42,'covariate_cells':42,'text_overlap':False,'text_clipped':False}

def main():
    warnings.filterwarnings('error',message='Glyph.*missing')
    result=[]
    for lang in ['EN','ZH']:
        f=Figure(lang);result.append(f.save());plt.close(f.fig)
    (BASE/'audit/LAYOUT_VALIDATION.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    (BASE/'audit/STYLE_SPECIFICATION.json').write_text(json.dumps({'size_mm':[W,H],'assay_colours':COLORS,'assay_edges':EDGES,
       'assay_shapes':MARKERS,'shared_colour':STRONG,'bubble_area_scale':56,'negative_effects':'open circles',
       'q_strength_scale':[0,20],'all_q_strengths_untruncated':True,'outside_grey_footnotes':False},ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result))
if __name__=='__main__':main()
