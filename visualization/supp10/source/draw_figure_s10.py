"""Draw bilingual S14 as native vector heatmaps and conditional-power comparisons."""
from pathlib import Path
import json, logging, sys, warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from PIL import Image

logging.getLogger('fontTools').setLevel(logging.ERROR)
BASE = Path(__file__).resolve().parent.parent
AUDIT = BASE / 'audit'
DATA = json.loads((AUDIT/'PLOT_INPUT.json').read_text(encoding='utf-8'))
W, H = 180, 258
INK, BLUE, ROSE = '#263A46', '#527F94', '#A36872'
PALETTE = ['#FAF5E8','#D0E3DB','#92B3CD','#35594F','#B0C6A8','#D99DA3','#E8E2DB']
CMAP = LinearSegmentedColormap.from_list('conditional_power',
    [(0,'#FAF5E8'),(.35,'#D0E3DB'),(.70,'#92B3CD'),(1,'#35594F')],N=1024)
TEXT = {
 'EN':dict(figure='Supplementary Figure S10',protein='Protein',participants='n',
    A='Nominal threshold',B='BY first-step threshold',
    effect='Assumed |γ| in all three assays',
    focus='Shaded labels: shared p-tau181',power='Power (%)',
    C='One weaker association\nNominal threshold',D='One weaker association\nBY first-step threshold',
    dot='Protein',median='Median',
    configurations=['All 0.05','NULISA\nweaker','SomaScan\nweaker','TMT-MS\nweaker']),
 'ZH':dict(figure='补充图 S10',protein='蛋白',participants='人数',
    A='名义检验阈值',B='BY首步阈值',effect='三种检测共同假定的效应幅度 |γ|',
    focus='底色标记：共同p-tau181关联',power='功效（%）',
    C='单一检测关联较弱\n名义检验阈值',D='单一检测关联较弱\nBY首步阈值',
    dot='蛋白',median='中位数',
    configurations=['均为0.05','NULISA\n较弱','SomaScan\n较弱','TMT-MS\n较弱'])
}
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
    'axes.linewidth':.7,'xtick.major.width':.7,'ytick.major.width':.7,
    'xtick.major.size':2.5,'ytick.major.size':0,'axes.unicode_minus':False,
    'savefig.facecolor':'white','figure.facecolor':'white'})

def percentage_label(probability):
    if probability == 0: return '0'
    if probability*100 < .01: return '<0.01'
    return f'{probability*100:.2f}'

def contrast_ink(colour):
    def lum(rgb):
        v=np.array(rgb[:3]);v=np.where(v<=.04045,v/12.92,((v+.055)/1.055)**2.4)
        return float(v@[.2126,.7152,.0722])
    light=lum(colour);dark=lum(matplotlib.colors.to_rgba(INK))
    return 'white' if 1.05/(light+.05) > (light+.05)/(dark+.05) else INK

class Figure:
    def __init__(self,language):
        self.lang=language;self.t=TEXT[language]
        self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=160)
        self.primitives={};self.axes={}
    def font(self,size=8.8,bold=False,math=False):
        name='segoeui.ttf' if math else (('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN'
             else ('msyhbd.ttc' if bold else 'msyh.ttc'))
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/name),size=size,weight='bold' if bold else 'normal')
    def text(self,x,y,label,size=8.8,bold=False,ha='left',math=False,colour=INK):
        return self.fig.text(x/W,1-y/H,label,fontproperties=self.font(size,bold,math),color=colour,
                             ha=ha,va='center',linespacing=1.02)
    def box(self,x,y,width,height,colour,edge='none',lw=0):
        self.fig.patches.append(Rectangle((x/W,1-(y+height/2)/H),width/W,height/H,
            transform=self.fig.transFigure,facecolor=colour,edgecolor=edge,linewidth=lw,zorder=-2))
    def line(self,xs,ys,colour,lw=.7):
        self.fig.lines.append(Line2D(np.array(xs)/W,1-np.array(ys)/H,
            transform=self.fig.transFigure,color=colour,linewidth=lw,zorder=1))
    def header(self,key,x,y,width,fill,badge,two_line=False):
        self.box(x,y,width,10,fill);self.box(x,y,5.5,10,badge)
        self.text(x+2.75,y,key,11.2,True,'center')
        if two_line:self.text(x+6.8,y,self.t[key],9.2,True)
        else:
            self.text(x+6.8,y-2.2,self.t[key],9.2,True)
            self.text(x+6.8,y+2.3,'P ≤ 0.05' if key=='A' else 'P ≤ 1.53 × 10⁻⁴',8.8,math=True)
    def axis(self,key,x,top,width,height):
        ax=self.fig.add_axes([x/W,1-(top+height)/H,width/W,height/H],facecolor='none')
        ax.spines[['top','right','left']].set_visible(False);ax.spines['bottom'].set_color(INK)
        ax.tick_params(colors=INK,direction='out',pad=2);ax.set_axisbelow(True)
        self.axes[key]=ax;return ax
    def style(self,ax):
        for label in ax.get_xticklabels()+ax.get_yticklabels()+[ax.xaxis.label,ax.yaxis.label]:
            label.set_fontproperties(self.font());label.set_color(INK)
        ticks=ax.get_xticklabels()
        if ticks: ticks[0].set_horizontalalignment('left');ticks[-1].set_horizontalalignment('right')
    def draw(self):
        self.text(2,6,self.t['figure'],13,True)
        self.header('A',40,21,62,'#FAF5E8','#92B3CD')
        self.header('B',116,21,62,'#D0E3DB','#B0C6A8')
        self.text(109,34,self.t['effect'],9,ha='center')
        self.text(2,41,self.t['protein'],9,True);self.text(34,41,self.t['participants'],9,True,'right')
        rows={g:47+i*3.68 for i,g in enumerate(DATA['genes'])}
        for gene,y in rows.items():
            if gene in DATA['focus14']:self.box(1.5,y,34,3.68,'#D0E3DB')
            if gene=='CST3':self.box(1.5,y,34,3.68,'#E8E2DB')
            n=next(r['n'] for r in DATA['balanced'] if r['gene']==gene)
            self.text(2,y,gene+('†' if gene=='CST3' else ''))
            self.text(34,y,str(n),ha='right')
        for key,x,threshold in [('A',40,'NOMINAL'),('B',116,'BY68_SUFFICIENT')]:
            cells=[]
            for j,effect in enumerate([.02,.05,.1]):
                xc=x+(j+.5)*62/3
                self.text(xc,41,f'{effect:.2f}',9,True,'center')
                for gene,y in rows.items():
                    r=next(r for r in DATA['balanced'] if r['gene']==gene and r['effect']==effect and r['threshold']==threshold)
                    face=CMAP(r['rate']) if r['qualified'] else '#E8E2DB'
                    self.box(x+j*62/3,y,62/3,3.68,face,edge='white',lw=.3)
                    label=percentage_label(r['rate']) if r['qualified'] else '—'
                    self.text(xc,y,label,8.8,ha='center',colour=contrast_ink(matplotlib.colors.to_rgba(face)))
                    cells.append(dict(gene=gene,effect=effect,threshold=threshold,
                        probability=r['rate'],displayed_as_power=r['qualified'],label=label))
            self.primitives[key]=dict(cells=cells,scale=[0,1],sample_counts_from_original_designs=True,
                shared_pathology_row_labels=DATA['focus14'],masked_gene='CST3')
        # Shared continuous vector colour key for both upper panels.
        self.text(2,178.2,self.t['focus'],8.8)
        for i in range(256):
            col=CMAP((i+.5)/256)
            self.box(100+i*78/256,177.5,78/256+(.025 if i<255 else 0),3,col,edge=col,lw=.05)
        for value in [0,25,50,75,100]:
            x=100+78*value/100
            self.text(x,183.5,str(value),8.8,ha='left' if value==0 else ('right' if value==100 else 'center'))
        self.text(139,189.5,self.t['power'],8.8,ha='center')
        self.fig.lines.append(Line2D([3/W],[1-190/H],transform=self.fig.transFigure,
            linestyle='none',marker='o',ms=3.2,mfc='#92B3CD',mec=BLUE,mew=.45))
        self.text(7,190,self.t['dot'])
        self.line([39,39],[188.2,191.8],ROSE,1.3);self.text(43,190,self.t['median'])
        self.header('C',2,203,84,'#FAF5E8','#92B3CD',two_line=True)
        self.header('D',92,203,86,'#D0E3DB','#B0C6A8',two_line=True)
        qualified=[g for g in DATA['genes'] if g!='CST3']
        spread=np.linspace(-2.05,2.05,len(qualified))
        for key,x,labelx,threshold,limit,ticks in [
             ('C',28,2,'NOMINAL',75,[0,25,50,75]),
             ('D',122,94,'BY68_SUFFICIENT',1.2,[0,.4,.8,1.2])]:
            ax=self.axis(key,x,212,56,33)
            points=[];medians=[]
            for j,config in enumerate(DATA['configurations']):
                y=216+j*8
                self.text(labelx,y,self.t['configurations'][j],8.8)
                values=[]
                for gene,offset in zip(qualified,spread):
                    r=next(r for r in DATA['comparison'] if r['gene']==gene and
                           r['configuration']==config['key'] and r['threshold']==threshold)
                    value=r['rate']*100;values.append(value)
                    ax.plot(value,y+offset,linestyle='none',marker='o',ms=2.7,
                            mfc='#92B3CD',mec=BLUE,mew=.35,alpha=.75,clip_on=False)
                    points.append(dict(gene=gene,configuration=config['key'],threshold=threshold,
                                       probability=r['rate'],x_percent=value,vertical_offset_mm=float(offset)))
                median=float(np.median(values));ax.vlines(median,y-2.35,y+2.35,color=ROSE,lw=1.2,zorder=4)
                medians.append(dict(configuration=config['key'],median_percent=median))
            ax.set(xlim=(0,limit),ylim=(245,212),xticks=ticks,yticks=[],xlabel=self.t['power'])
            ax.set_xticklabels([str(t) if t!=0 else '0' for t in ticks])
            ax.xaxis.grid(True,color='#E4EBE7',linewidth=.5);self.style(ax)
            self.primitives[key]=dict(points=points,medians=medians,x_limits_percent=[0,limit],
                categorical_vertical_spread_only=True,quantitative_x_values_unchanged=True,
                excluded_unqualified_gene='CST3',protein_count=33)
        self.fig.canvas.draw();return self
    def layout(self):
        self.fig.canvas.draw();renderer=self.fig.canvas.get_renderer();records=[];seen=set()
        for artist in self.fig.findobj(Text):
            if not artist.get_visible() or not artist.get_text().strip() or id(artist) in seen:continue
            seen.add(id(artist));box=artist.get_window_extent(renderer)
            mm=np.array([box.x0,box.y0,box.x1,box.y1])/self.fig.dpi*25.4
            assert mm[0]>=-.15 and mm[2]<=W+.15 and mm[1]>=-.15 and mm[3]<=H+.15,('clipped text',self.lang,artist.get_text(),mm.tolist())
            records.append(dict(text=artist.get_text(),point_size=float(artist.get_fontsize()),
                x0_mm=float(mm[0]),y0_mm=float(H-mm[3]),x1_mm=float(mm[2]),y1_mm=float(H-mm[1])))
        overlaps=[]
        for i,a in enumerate(records):
            for b in records[i+1:]:
                xx=min(a['x1_mm'],b['x1_mm'])-max(a['x0_mm'],b['x0_mm'])
                yy=min(a['y1_mm'],b['y1_mm'])-max(a['y0_mm'],b['y0_mm'])
                if xx>.35 and yy>.35:overlaps.append((a['text'],b['text'],round(xx,2),round(yy,2)))
        assert not overlaps,('text overlap',self.lang,overlaps)
        return records
    def save(self):
        records=self.layout()
        for folder in ['editable','high_resolution']:(BASE/folder).mkdir(exist_ok=True)
        title='Supplementary Figure S10: Conditional power for detecting shared associations with annual memory change'
        self.fig.savefig(BASE/f'Figure_S10_{self.lang}.pdf',metadata={'Title':title,'Creator':'Study supplementary figures'})
        self.fig.savefig(BASE/f'editable/Figure_S10_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_S10_{self.lang}.png',dpi=300)
        tif=BASE/f'high_resolution/Figure_S10_{self.lang}.tiff'
        self.fig.savefig(tif,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        with Image.open(tif) as im:
            rgb=im.convert('RGB');rgb.save(tif,compression='tiff_lzw',dpi=(600,600));rgb.close()
        (AUDIT/f'LAYOUT_{self.lang}.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
        (AUDIT/f'PLOT_PRIMITIVES_{self.lang}.json').write_text(json.dumps(self.primitives,ensure_ascii=False,indent=2),encoding='utf-8')
        return dict(language=self.lang,width_mm=W,height_mm=H,minimum_font_pt=min(r['point_size'] for r in records),
                    checked_text_elements=len(records),text_clipped=False,text_overlap=False)

def main():
    warnings.filterwarnings('error',message='Glyph.*missing')
    results=[]
    for lang in ['EN','ZH']:
        figure=Figure(lang).draw();results.append(figure.save());plt.close(figure.fig)
    (AUDIT/'LAYOUT_VALIDATION.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
    (AUDIT/'STYLE_SPECIFICATION.json').write_text(json.dumps(dict(palette=PALETTE,
        choice='Cream-to-blue-to-deep-green conditional power, sage pathology labels and rose medians.',
        figure_size_mm=[W,H],panels=4,proteins=34,qualified_power_cells=198,masked_cells=6,
        comparison_points=264,identical_heatmap_colour_scale=True,separate_comparison_axes=True,
        numeric_power_units='Percent',quantitative_coordinates_not_jittered=True,
        loose_footer_notes=False),indent=2),encoding='utf-8')
    print(json.dumps(results),flush=True)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');main()
