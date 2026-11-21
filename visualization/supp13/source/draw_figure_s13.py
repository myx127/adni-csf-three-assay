"""Draw residual correlations, fixed-scale coefficient shifts and full information retention."""
from pathlib import Path
import json,logging,sys,warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from PIL import Image

logging.getLogger('fontTools').setLevel(logging.ERROR)
BASE=Path(__file__).resolve().parent.parent
AUDIT=BASE/'audit'
DATA=json.loads((AUDIT/'PLOT_INPUT.json').read_text(encoding='utf-8'))
W,H=180,258
INK='#263A46';GRID='#E9EEEB'
COLORS={'NULISA':'#FAD6C9','SomaScan':'#B0C6A8','TMT-MS':'#92B3CD'}
EDGES={'NULISA':'#B5806D','SomaScan':'#6F8C64','TMT-MS':'#547F9E'}
MARKERS={'NULISA':'o','SomaScan':'s','TMT-MS':'D'}
PALETTE=['#FAD6C9','#B0C6A8','#92B3CD','#FAF5E8','#D0E3DB','#E8E2DB','#527F94','#C289A0','#35594F']
CORR_CMAP=LinearSegmentedColormap.from_list('BlueCreamRose',['#527F94','#FAF5E8','#C289A0'])
INFO_CMAP=LinearSegmentedColormap.from_list('CreamSageForest',['#FAF5E8','#D0E3DB','#35594F'])
UPPER_ROWS=np.array([52,59,66,82,89,96,112,119,126.])
LOWER_ROWS=np.array([193,197.7,202.4,207.1,214.1,218.8,223.5,230.5,235.2])
TEXT={
 'EN':dict(figure='Supplementary Figure S13',A='Covariate-adjusted\ncorrelations',
   B='Memory associations before\nand after p-tau181 adjustment',
   C='Protein information retained across adjustment models',
   assay='Assay',corrheads=['Protein–\np-tau181','Protein–\nmemory','p-tau181–\nmemory'],
   before='Without\np-tau181',after='With log₂\np-tau181',r='Partial r',
   units='ADNI-MEM/year\nper protein residual SD',
   modelhead='p-tau181 and\nbackground adjustment',
   models=['None','Linear log₂','Linear untransformed','Spline log₂',
           'Two scores only','Two scores + linear log₂','Two scores + spline log₂',
           'Three scores only','Three scores + linear log₂'],
   retained='Information retained (%)'),
 'ZH':dict(figure='补充图 S13',A='共同协变量调整后的\n相关关系',
   B='p-tau181调整前后的\n记忆关联',
   C='不同调整模型保留的蛋白信息',
   assay='检测',corrheads=['蛋白–\np-tau181','蛋白–\n记忆','p-tau181–\n记忆'],
   before='不调整\np-tau181',after='调整对数\np-tau181',r='偏相关r',
   units='ADNI-MEM分/年\n每1蛋白残差SD',
   modelhead='p-tau181及\n背景调整',
   models=['均不加入','对数线性','原值线性','对数样条',
           '仅两背景','两背景及对数线性','两背景及对数样条',
           '仅三背景','三背景及对数线性'],
   retained='信息保留比例（%）')
}
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.linewidth':.7,
    'xtick.major.width':.7,'ytick.major.width':.7,'xtick.major.size':3,'ytick.major.size':0,
    'axes.unicode_minus':False,'savefig.facecolor':'white','figure.facecolor':'white'})

def luminance(rgb):
    a=np.asarray(rgb[:3]);linear=np.where(a<=.04045,a/12.92,((a+.055)/1.055)**2.4)
    return float(linear@np.array([.2126,.7152,.0722]))

def readable_text(rgb):
    lum=luminance(rgb);ink_lum=luminance(matplotlib.colors.to_rgba(INK))
    contrast_ink=(lum+.05)/(ink_lum+.05);contrast_white=1.05/(lum+.05)
    preferred=('white',contrast_white) if contrast_white>contrast_ink else (INK,contrast_ink)
    # Mid-tone cells need black rather than the softer figure text colour.
    return ('black',(lum+.05)/.05) if preferred[1]<4.5 else preferred

class Figure:
    def __init__(self,lang):
        self.lang=lang;self.t=TEXT[lang]
        self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=160)
        self.axes={};self.primitives={}
    def font(self,size=9,bold=False,math=False):
        name='segoeui.ttf' if math else (('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN' else ('msyhbd.ttc' if bold else 'msyh.ttc'))
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/name),size=size,weight='bold' if bold else 'normal')
    def text(self,x,y,label,size=9,bold=False,ha='left',math=False,rotation=0,color=INK):
        return self.fig.text(x/W,1-y/H,label,fontproperties=self.font(size,bold,math),color=color,ha=ha,va='center',
                             linespacing=1.06,rotation=rotation,rotation_mode='anchor')
    def box(self,x,y,width,height,fill,alpha=1,edge='none',lw=0):
        self.fig.patches.append(Rectangle((x/W,1-(y+height/2)/H),width/W,height/H,transform=self.fig.transFigure,
            facecolor=fill,alpha=alpha,edgecolor=edge,linewidth=lw,zorder=-1))
    def header(self,panel,x,y,width,badge,fill,tall=False):
        h=10 if tall else 6.4
        self.box(x,y,width,h,fill);self.box(x,y,5.5,h,badge)
        self.text(x+2.75,y,panel,11.2,True,'center')
        self.text(x+6.8,y,self.t[panel],9.2 if tall else 10,True)
    def style(self,ax):
        for label in ax.get_xticklabels()+ax.get_yticklabels()+[ax.xaxis.label,ax.yaxis.label]:
            label.set_fontproperties(self.font(9));label.set_color(INK)
    def draw(self):
        self.text(2,5.7,self.t['figure'],13,True)
        self.header('A',2,19,90,'#FAD6C9','#FAF5E8',True)
        self.header('B',98,19,80,'#92B3CD','#E8E2DB',True)
        self.text(3.5,32,self.t['assay'],8.8,True)
        for j,label in enumerate(self.t['corrheads']):self.text(33.5+21*j,32,label,8.8,True,'center')
        for x,filled,key in [(102,False,'before'),(143,True,'after')]:
            self.fig.lines.append(Line2D([x/W],[1-32/H],transform=self.fig.transFigure,linestyle='none',marker='o',
                markersize=4,markerfacecolor=INK if filled else 'white',markeredgecolor=INK,markeredgewidth=.7))
            self.text(x+3.5,32,self.t[key],8.8,math=self.lang=='EN' and '₂' in self.t[key])
        for gene,n,y in zip(DATA['proteins'],[553,557,552],[43,73,103]):
            self.box(2,y,176,5.4,'#D0E3DB',.5)
            self.text(3.5,y,f'{gene}  (n={n})',9,True)
        r=np.array([row['correlations'] for row in DATA['rows']])
        corr_labels=[];contrasts=[]
        for i,(row,y) in enumerate(zip(DATA['rows'],UPPER_ROWS)):
            self.text(3.5,y,row['assay'],8.8)
            printed=[]
            for j,v in enumerate(row['correlations']):
                rgba=CORR_CMAP((v+1)/2);color,ratio=readable_text(rgba);contrasts.append(ratio)
                self.box(23+21*j,y,21,6,rgba,edge='white',lw=.7)
                label=f'{v:.3f}'.replace('-','−');printed.append(label)
                self.text(33.5+21*j,y,label,8.8,ha='center',color=color)
            corr_labels.append(printed)
        self.primitives['A']=dict(order=[dict(protein=r['protein'],assay=r['assay']) for r in DATA['rows']],
            correlation_order=['Protein-p-tau181','Protein-memory slope','p-tau181-memory slope'],
            correlations=r.tolist(),displayed_text=corr_labels,colour_range=[-1,1],colours_truncated=False,
            common_covariates_only=True,third_variable_not_added_to_conditioning=True)
        ax=self.fig.add_axes([102/W,1-131/H,74/W,83/H],facecolor='none');self.axes['B']=ax
        ax.spines[['top','right','left']].set_visible(False);ax.spines['bottom'].set_color(INK)
        ax.tick_params(colors=INK,direction='out',pad=3);ax.set_axisbelow(True)
        for low,high in [(48,69.5),(78,99.5),(108,129.5)]:
            ax.vlines(0,low,high,color='#A6B3C0',linestyle=(0,(3,2)),linewidth=.75,zorder=1)
            for tick in [.05,.1,.15]:ax.vlines(tick,low,high,color=GRID,linewidth=.5,zorder=0)
        points=[]
        for row,y in zip(DATA['rows'],UPPER_ROWS):
            assay=row['assay'];before=row['before'];after=row['after']
            ax.plot([before,after],[y,y],color=EDGES[assay],linewidth=.9,zorder=2)
            ax.plot(before,y,linestyle='none',marker=MARKERS[assay],markersize=4.4,mfc='white',mec=EDGES[assay],mew=.85,zorder=4)
            ax.plot(after,y,linestyle='none',marker=MARKERS[assay],markersize=4.4,mfc=COLORS[assay],mec=EDGES[assay],mew=.65,zorder=5)
            points.append(dict(protein=row['protein'],assay=assay,before=before,after=after,y=float(y),fixed_original_scale=True))
        ax.set(xlim=(-.025,.15),ylim=(131,48),xticks=[0,.05,.1,.15],yticks=[],xlabel=self.t['units'])
        ax.set_xticklabels(['0','0.05','0.10','0.15']);self.style(ax)
        self.primitives['B']=dict(points=points,x_limits=[-.025,.15],point_estimates_only=True,
            before_model='No p-tau181',after_model='Linear log2 p-tau181',
            no_individual_time_trajectory=True,no_new_tests=True)
        cax=self.fig.add_axes([23/W,1-138/H,63/W,2.8/H]);self.axes['correlation_colourbar']=cax
        cb=self.fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(-1,1),cmap=CORR_CMAP),cax=cax,orientation='horizontal',ticks=[-1,0,1])
        cb.solids.set_rasterized(False);cb.solids.set_edgecolor('face');cb.outline.set_linewidth(.55)
        cb.outline.set_edgecolor(INK);cb.set_label(self.t['r'],labelpad=3);cax.tick_params(colors=INK,pad=2);self.style(cax)
        self.header('C',2,160,176,'#B0C6A8','#D0E3DB')
        self.text(2,182,self.t['modelhead'],8.8,True)
        for k,gene in enumerate(DATA['proteins']):
            self.box(51+42*k,171,42,5.6,'#FAF5E8')
            self.text(72+42*k,171,gene,9.2,True,'center')
        for k,col in enumerate(DATA['information']['columns']):
            self.text(58+14*k,183,col['assay'],8.8,ha='center',rotation=55)
        fractions=np.array(DATA['information']['fractions']);info_labels=[]
        for i,(model,y) in enumerate(zip(DATA['models'],LOWER_ROWS)):
            label=self.t['models'][i]
            self.text(2,y,label,8.8,math=self.lang=='EN' and '₂' in label)
            printed=[]
            for j,value in enumerate(fractions[i]):
                rgba=INFO_CMAP(value);color,ratio=readable_text(rgba);contrasts.append(ratio)
                self.box(51+14*j,y,14,4.7,rgba,edge='white',lw=.6)
                label=f'{100*value:.1f}';printed.append(label)
                self.text(58+14*j,y,label,8.8,ha='center',color=color)
            info_labels.append(printed)
        self.primitives['C']=dict(models=DATA['models'],columns=DATA['information']['columns'],
            fractions=fractions.tolist(),displayed_percentages=info_labels,colour_range=[0,100],colours_truncated=False,
            common_no_ptau_no_background_reference=True,no_new_tests=True)
        cax=self.fig.add_axes([115/W,1-246/H,62/W,2.8/H]);self.axes['information_colourbar']=cax
        cb=self.fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,100),cmap=INFO_CMAP),cax=cax,orientation='horizontal',ticks=[0,25,50,75,100])
        cb.solids.set_rasterized(False);cb.solids.set_edgecolor('face');cb.outline.set_linewidth(.55)
        cb.outline.set_edgecolor(INK);cb.set_label(self.t['retained'],labelpad=3);cax.tick_params(colors=INK,pad=2);self.style(cax)
        self.minimum_cell_contrast=min(contrasts)
        assert self.minimum_cell_contrast>=4.5,('low cell contrast',self.minimum_cell_contrast)
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
        title='Supplementary Figure S13: Regression relationships underlying changes in memory associations after p-tau181 adjustment'
        self.fig.savefig(BASE/f'Figure_S13_{self.lang}.pdf',metadata={'Title':title,'Creator':'Study supplementary figures'})
        self.fig.savefig(BASE/f'editable/Figure_S13_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_S13_{self.lang}.png',dpi=300)
        tif=BASE/f'high_resolution/Figure_S13_{self.lang}.tiff'
        self.fig.savefig(tif,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        with Image.open(tif) as im:
            rgb=im.convert('RGB');rgb.save(tif,compression='tiff_lzw',dpi=(600,600));rgb.close()
        (AUDIT/f'LAYOUT_{self.lang}.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
        (AUDIT/f'PLOT_PRIMITIVES_{self.lang}.json').write_text(json.dumps(self.primitives,ensure_ascii=False,indent=2),encoding='utf-8')
        return dict(language=self.lang,width_mm=W,height_mm=H,minimum_font_pt=min(r['point_size'] for r in records),
            checked_text_elements=len(records),minimum_cell_text_contrast=self.minimum_cell_contrast,text_clipped=False,text_overlap=False)

def main():
    warnings.filterwarnings('error',message='Glyph.*missing')
    results=[]
    for lang in ['EN','ZH']:
        figure=Figure(lang).draw();results.append(figure.save());plt.close(figure.fig)
    (AUDIT/'LAYOUT_VALIDATION.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
    (AUDIT/'STYLE_SPECIFICATION.json').write_text(json.dumps(dict(
        palette=PALETTE,assay_colours=COLORS,assay_shapes=MARKERS,
        correlation_gradient=['#527F94','#FAF5E8','#C289A0'],information_gradient=['#FAF5E8','#D0E3DB','#35594F'],
        choice='Continue assay colors; use blue/cream/rose for signed correlations and cream/sage/forest for retained information.',
        figure_size_mm=[W,H],panels=3,correlation_range=[-1,1],memory_coefficient_range=[-.025,.15],
        information_range_percent=[0,100],colours_truncated=False,colourbars_vector=True,
        loose_footer_notes=False,point_estimates_not_individual_trajectories=True),ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(results,ensure_ascii=False),flush=True)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');main()
