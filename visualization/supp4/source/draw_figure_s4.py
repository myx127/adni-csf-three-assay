"""Draw complete, aligned assay-specific effects and shared pathology evidence."""
from pathlib import Path
import json,sys,warnings
import logging
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.text import Text
from PIL import Image

logging.getLogger('fontTools').setLevel(logging.ERROR)
BASE=Path(__file__).resolve().parent.parent
AUDIT=BASE/'audit'
DATA=json.loads((AUDIT/'PLOT_INPUT.json').read_text(encoding='utf-8'))
W,H=180,232
PALETTE=['#9BCFC1','#A4C4D1','#B9A0CF','#E4A9CD','#F1C5D4','#C2D6B4']
COLORS={'NULISA':'#9BCFC1','SomaScan':'#A4C4D1','TMT-MS':'#B9A0CF'}
EDGES={'NULISA':'#467E6F','SomaScan':'#507F96','TMT-MS':'#7A5994'}
MARKERS={'NULISA':'o','SomaScan':'s','TMT-MS':'D'}
INK='#263C59';GRID='#E9EDF0';STRONG='#A95788'
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
    'axes.linewidth':.7,'xtick.major.width':.7,'ytick.major.width':.7,'xtick.major.size':3,'ytick.major.size':0,
    'axes.unicode_minus':False,'savefig.facecolor':'white','figure.facecolor':'white'})

TEXT={
 'EN':dict(figure='Supplementary Figure S4',D='Shared\nassociation',protein='Protein',n='n',beta='Standardized β',q='−log₁₀ q',
     unit='Protein residual SD per SD of log₂ p-tau181',supported='Shared q ≤ 0.05 (14 proteins)',cutoff='q = 0.05'),
 'ZH':dict(figure='补充图 S4',D='共同\n关联',protein='蛋白',n='n',beta='标准化β',q='−log₁₀ q',
     unit='每1 SD log₂ p-tau181对应的蛋白残差SD变化',supported='共同关联q≤0.05（14个蛋白）',cutoff='q = 0.05')}

class Figure:
    def __init__(self,lang):
        self.lang=lang;self.t=TEXT[lang];self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=160)
        self.axes={};self.primitives={}
    def font(self,size=9,bold=False):
        filename=('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN' else ('msyhbd.ttc' if bold else 'msyh.ttc')
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/filename),size=size,weight='bold' if bold else 'normal')
    def text(self,x,y,value,size=9,bold=False,ha='left'):
        return self.fig.text(x/W,1-y/H,value,fontproperties=self.font(size,bold),color=INK,ha=ha,va='center',linespacing=1.05)
    def header(self,letter,title,x,width,color):
        self.fig.patches.extend([
            Rectangle((x/W,1-23/H),width/W,11/H,transform=self.fig.transFigure,facecolor=color,alpha=.30,edgecolor='none',zorder=-1),
            Rectangle((x/W,1-23/H),5.8/W,11/H,transform=self.fig.transFigure,facecolor=color,edgecolor='none')])
        self.text(x+2.9,17.5,letter,11.5,True,'center');self.text(x+7,17.5,title,9.5,True)
    def axis(self,key,x,y,w,h):
        ax=self.fig.add_axes([x/W,1-(y+h)/H,w/W,h/H]);self.axes[key]=ax
        ax.spines[['top','right','left']].set_visible(False)
        ax.spines['bottom'].set_color(INK);ax.tick_params(colors=INK,direction='out',pad=3)
        ax.set_axisbelow(True)
        return ax
    def style(self,ax):
        for label in ax.get_xticklabels()+ax.get_yticklabels():label.set_fontproperties(self.font(9));label.set_color(INK)
        for label in [ax.xaxis.label,ax.yaxis.label]:label.set_fontproperties(self.font(9));label.set_color(INK)
    def draw(self):
        self.text(2,5.7,self.t['figure'],13,True)
        columns=[('A','NULISA',32,32),('B','SomaScan',66.5,32),('C','TMT-MS',101,32)]
        for letter,assay,x,width in columns:self.header(letter,assay,x-1,width+1,COLORS[assay])
        self.header('D',self.t['D'],137,41,'#C2D6B4')
        self.fig.patches.append(Rectangle((2/W,1-29.5/H),4.8/W,2.7/H,transform=self.fig.transFigure,
            facecolor='#F1C5D4',edgecolor='none',alpha=.62))
        self.text(8,28.2,self.t['supported'],8.8)
        self.fig.lines.append(Line2D([138/W,145/W],[1-28.2/H,1-28.2/H],transform=self.fig.transFigure,
            color=STRONG,linewidth=.85,linestyle=(0,(3,2))))
        self.text(147,28.2,self.t['cutoff'],8.8)
        self.text(2,35,self.t['protein'],9,True);self.text(29.3,35,self.t['n'],9,True,'right')
        yy=np.arange(34);supported=np.asarray(DATA['shared_supported'],bool)
        for letter,assay,x,width in columns:
            ax=self.axis(letter,x,40,width,174)
            e=DATA['assays'][assay];beta=np.asarray(e['beta']);lo=np.asarray(e['ci_low']);hi=np.asarray(e['ci_high'])
            for i,yes in enumerate(supported):
                if yes:ax.axhspan(i-.48,i+.48,color='#F1C5D4',alpha=.21,zorder=0)
            ax.axvline(0,color='#82909A',linewidth=.75,linestyle=(0,(3,2)),zorder=1)
            ax.errorbar(beta,yy,xerr=np.vstack([beta-lo,hi-beta]),fmt=MARKERS[assay],
                markersize=3.9,mfc=COLORS[assay],mec=EDGES[assay],mew=.55,color=EDGES[assay],elinewidth=.8,capsize=1.9,capthick=.65,zorder=3)
            ax.set(xlim=(-.65,2.15),ylim=(33.65,-.65),xticks=[-.5,0,1,2],yticks=[],xlabel=self.t['beta'])
            ax.xaxis.grid(True,color=GRID,linewidth=.5)
            self.primitives[letter]=dict(assay=assay,proteins=DATA['proteins'],beta=beta.tolist(),ci_low=lo.tolist(),ci_high=hi.tolist(),
                x_limits=[-.65,2.15],y=yy.tolist(),highlighted=supported.tolist())
        # Protein names and sample sizes share the exact y positions of all four panels.
        anchor=self.axes['A']
        for i,(gene,n) in enumerate(zip(DATA['proteins'],DATA['participants'])):
            ypx=anchor.transData.transform((0,i))[1];y_mm=H-ypx/self.fig.dpi*25.4
            if supported[i]:
                self.fig.patches.append(Rectangle((1/W,1-(y_mm+2.4)/H),29/W,4.8/H,
                    transform=self.fig.transFigure,facecolor='#F1C5D4',alpha=.21,edgecolor='none',zorder=-1))
            self.text(2,y_mm,gene,8.8,bool(supported[i]));self.text(29.3,y_mm,str(n),8.8,False,'right')
        q=np.asarray(DATA['shared_q']);strength=-np.log10(q)
        ax=self.axis('D',138,40,39,174)
        for i,yes in enumerate(supported):
            if yes:ax.axhspan(i-.48,i+.48,color='#F1C5D4',alpha=.21,zorder=0)
        ax.axvline(-np.log10(.05),color=STRONG,linestyle=(0,(3,2)),linewidth=.85,zorder=1)
        for i,value in enumerate(strength):
            ax.plot([0,value],[i,i],color='#E4A9CD' if supported[i] else '#A4C4D1',linewidth=.8,zorder=2)
        ax.scatter(strength[supported],yy[supported],s=20,marker='o',facecolor='#E4A9CD',edgecolors=STRONG,linewidths=.55,zorder=3)
        ax.scatter(strength[~supported],yy[~supported],s=17,marker='o',facecolor='white',edgecolors='#507F96',linewidths=.55,zorder=3)
        ax.set(xlim=(-1,21),ylim=(33.65,-.65),xticks=[0,5,10,15,20],yticks=[],xlabel=self.t['q'])
        ax.xaxis.grid(True,color=GRID,linewidth=.5)
        self.primitives['D']=dict(proteins=DATA['proteins'],q=q.tolist(),minus_log10_q=strength.tolist(),
            cutoff_q=.05,cutoff_minus_log10=float(-np.log10(.05)),highlighted=supported.tolist(),display_truncated=False)
        for axis in self.axes.values():self.style(axis)
        self.axes['D'].xaxis.label.set_fontproperties(FontProperties(fname='C:/Windows/Fonts/segoeui.ttf',size=9))
        self.fig.canvas.draw();return self
    def layout(self):
        self.fig.canvas.draw();renderer=self.fig.canvas.get_renderer();records=[];seen=set()
        for artist in self.fig.findobj(Text):
            if not artist.get_visible() or not artist.get_text().strip() or id(artist) in seen:continue
            seen.add(id(artist));box=artist.get_window_extent(renderer);mm=np.array([box.x0,box.y0,box.x1,box.y1])/self.fig.dpi*25.4
            assert mm[0]>=-.15 and mm[2]<=W+.15 and mm[1]>=-.15 and mm[3]<=H+.15,('text outside canvas',artist.get_text(),mm.tolist())
            records.append(dict(text=artist.get_text(),point_size=float(artist.get_fontsize()),x0_mm=float(mm[0]),y0_mm=float(H-mm[3]),x1_mm=float(mm[2]),y1_mm=float(H-mm[1])))
        overlaps=[]
        for i,a in enumerate(records):
            for b in records[i+1:]:
                xx=min(a['x1_mm'],b['x1_mm'])-max(a['x0_mm'],b['x0_mm']);yy=min(a['y1_mm'],b['y1_mm'])-max(a['y0_mm'],b['y0_mm'])
                if xx>.35 and yy>.35:overlaps.append((a['text'],b['text'],round(xx,2),round(yy,2)))
        assert not overlaps,('text overlap',self.lang,overlaps)
        return records
    def save(self):
        records=self.layout()
        for folder in ['editable','high_resolution']:(BASE/folder).mkdir(exist_ok=True)
        title='Supplementary Figure S4: Shared associations of CSF proteins with p-tau181 across three assays'
        self.fig.savefig(BASE/f'Figure_S4_{self.lang}.pdf',metadata={'Title':title,'Creator':'Study supplementary figures'})
        self.fig.savefig(BASE/f'editable/Figure_S4_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_S4_{self.lang}.png',dpi=300)
        tif=BASE/f'high_resolution/Figure_S4_{self.lang}.tiff';self.fig.savefig(tif,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        with Image.open(tif) as im:
            rgb=im.convert('RGB');rgb.save(tif,compression='tiff_lzw',dpi=(600,600));rgb.close()
        (AUDIT/f'LAYOUT_{self.lang}.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
        (AUDIT/f'PLOT_PRIMITIVES_{self.lang}.json').write_text(json.dumps(self.primitives,ensure_ascii=False,indent=2),encoding='utf-8')
        return dict(language=self.lang,width_mm=W,height_mm=H,minimum_font_pt=min(r['point_size'] for r in records),checked_text_elements=len(records),text_clipped=False,text_overlap=False)

def main():
    warnings.filterwarnings('error',message='Glyph.*missing')
    results=[]
    for lang in ['EN','ZH']:
        figure=Figure(lang).draw();results.append(figure.save());plt.close(figure.fig)
    (AUDIT/'LAYOUT_VALIDATION.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
    (AUDIT/'STYLE_SPECIFICATION.json').write_text(json.dumps(dict(palette=PALETTE,assay_colours=COLORS,assay_shapes=MARKERS,
        user_reference='Water green, pale blue, lavender, mauve, pale rose and sage: the user-specified Figure S4 palette.',
        explanatory_footnotes_outside_plot=False,figure_size_mm=[W,H],common_assay_effect_axis=[-.65,2.15],q_axis_truncated=False),ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(results,ensure_ascii=False),flush=True)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');main()
