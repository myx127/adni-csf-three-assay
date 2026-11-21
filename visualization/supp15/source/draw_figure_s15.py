"""Draw the common-signal illustration, preserved synthetic check and full pairwise grid."""
from pathlib import Path
import json, logging, sys, warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.text import Text
from PIL import Image

logging.getLogger('fontTools').setLevel(logging.ERROR)
BASE=Path(__file__).resolve().parent.parent
AUDIT=BASE/'audit'
DATA=json.loads((AUDIT/'PLOT_INPUT.json').read_text(encoding='utf-8'))
W,H=180,258
INK='#263A46'
ASSAYS=['NULISA','SomaScan','TMT-MS']
FACES=['#FAD6C9','#B0C6A8','#92B3CD']
EDGES=['#AA785F','#6F8C64','#527F94']
PALETTE=['#FAD6C9','#B0C6A8','#92B3CD','#D0E3DB','#FAF5E8','#527F94','#C289A0','#E2B9C7']
CMAP=LinearSegmentedColormap.from_list('association_contrast',
    ['#527F94','#FAF5E8','#C289A0'],N=1024)
NORM=Normalize(-.075,.075)
TEXT={
 'EN':dict(figure='Supplementary Figure S15',
    A='Illustrative signal model',B='Population and synthetic\nassociations',
    C='Without p-tau181\nadjustment',D='With p-tau181\nadjustment',
    P='p-tau181\nP',L='Common signal\nL',Y='Outcome\nY',U='Other\ncomponent U',
    equal='Equal loadings',ylabel='Outcome association (γ)',
    xticks=['Without\nadjustment','With\nadjustment'],
    population='Population',synthetic='Synthetic sample',
    reference='Reference assay error variance',comparison='Comparison assay error variance',
    contrast='Association contrast'),
 'ZH':dict(figure='补充图 S15',
    A='共同信号的示例模型',B='总体解析关联\n与合成样本估计',
    C='不调整\np-tau181',D='调整\np-tau181',
    P='p-tau181\nP',L='共同信号\nL',Y='结局\nY',U='其他共同\n成分 U',
    equal='检测载荷相同',ylabel='结局关联系数（γ）',xticks=['不调整','调整'],
    population='总体解析值',synthetic='合成样本',
    reference='参照检测的误差方差',comparison='比较检测的误差方差',
    contrast='检测间关联对比')
}
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.linewidth':.7,
    'xtick.major.width':.7,'ytick.major.width':.7,'xtick.major.size':2.5,'ytick.major.size':2.5,
    'axes.unicode_minus':False,'savefig.facecolor':'white','figure.facecolor':'white'})

def contrast_ink(colour):
    def lum(rgb):
        v=np.array(rgb[:3]);v=np.where(v<=.04045,v/12.92,((v+.055)/1.055)**2.4)
        return float(v@[.2126,.7152,.0722])
    light=lum(colour);dark=lum(matplotlib.colors.to_rgba(INK))
    return 'white' if 1.05/(light+.05)>(light+.05)/(dark+.05) else INK

def label_contrast(value):
    if value==0:return '0'
    return f'{value:.3f}'.replace('-','−')

class Figure:
    def __init__(self,language):
        self.lang=language;self.t=TEXT[language]
        self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=160)
        self.primitives={}
    def font(self,size=8.8,bold=False,math=False):
        name='segoeui.ttf' if math else (('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN'
            else ('msyhbd.ttc' if bold else 'msyh.ttc'))
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/name),size=size,weight='bold' if bold else 'normal')
    def text(self,x,y,label,size=8.8,bold=False,ha='left',math=False,colour=INK,rotation=0):
        return self.fig.text(x/W,1-y/H,label,fontproperties=self.font(size,bold,math),color=colour,
            ha=ha,va='center',linespacing=1.04,rotation=rotation)
    def box(self,x,y,width,height,face,edge='none',lw=0):
        self.fig.patches.append(Rectangle((x/W,1-(y+height/2)/H),width/W,height/H,
            transform=self.fig.transFigure,facecolor=face,edgecolor=edge,linewidth=lw,zorder=-2))
    def line(self,xs,ys,colour=INK,lw=.7):
        self.fig.lines.append(Line2D(np.array(xs)/W,1-np.array(ys)/H,transform=self.fig.transFigure,
            color=colour,linewidth=lw,zorder=0))
    def arrow(self,x0,y0,x1,y1):
        self.fig.patches.append(FancyArrowPatch((x0/W,1-y0/H),(x1/W,1-y1/H),
            transform=self.fig.transFigure,arrowstyle='-|>',mutation_scale=7.5,
            color=INK,linewidth=.75,shrinkA=0,shrinkB=0,zorder=0))
    def header(self,key,x,y,width,face,badge):
        self.box(x,y,width,10,face);self.box(x,y,5.5,10,badge)
        self.text(x+2.75,y,key,11.2,True,'center');self.text(x+6.8,y,self.t[key],9.2,True)
    def legend(self,x,y,label,marker,face,edge,math=False):
        self.fig.lines.append(Line2D([x/W],[1-y/H],transform=self.fig.transFigure,linestyle='none',
            marker=marker,ms=4,mfc=face,mec=edge,mew=.8))
        self.text(x+3.5,y,label,8.8,math=math)
    def draw(self):
        self.text(2,6,self.t['figure'],13,True)
        self.header('A',2,23,86,'#FAF5E8','#B0C6A8')
        self.header('B',96,23,82,'#D0E3DB','#92B3CD')
        # Planar directed construction: arrows specify hypothetical data generation only.
        nodes={'P':(20,44,'#D0E3DB'),'L':(66,44,'#FAF5E8'),
               'Y':(20,80,'#E2B9C7'),'U':(66,80,'#D0E3DB')}
        for name,(x,y,fill) in nodes.items():
            self.box(x-13,y,26,14,fill,edge='#A6B3C0',lw=.55)
            self.text(x,y,self.t[name],8.8,ha='center')
        self.arrow(33,44,53,44);self.text(43,38.5,'0.8',8.8,ha='center')
        self.arrow(20,51,20,73);self.text(24,62,'−0.06',8.8)
        self.arrow(66,73,66,51);self.text(70,62,'1',8.8)
        self.arrow(53,80,33,80);self.text(43,85,'0.10',8.8,ha='center')
        self.line([79,88,88,18],[44,44,96,96])
        self.text(47,92,self.t['equal'],8.8,ha='center')
        for i,(name,x,face,edge) in enumerate(zip(ASSAYS,[18,45,72],FACES,EDGES)):
            self.arrow(x,96,x,103)
            self.box(x-13,110,26,14,face,edge=edge,lw=.6)
            self.text(x,107.8,name,8.8,ha='center')
            self.text(x,112.8,'L + e'+['₁','₂','₃'][i],8.8,ha='center',math=True)
        self.primitives['A']=dict(model=DATA['model'],assays=ASSAYS,
            edges=[['P','L',.8],['P','Y',-.06],['U','L',1],['U','Y',.1]],
            measurement_equations=['L + e1','L + e2','L + e3'],
            illustrative_data_generation_not_empirical_causal_graph=True)

        ax=self.fig.add_axes([110/W,1-116/H,67/W,79/H],facecolor='none')
        ax.spines[['top','right']].set_visible(False)
        ax.spines[['left','bottom']].set_color(INK);ax.tick_params(colors=INK,direction='out',pad=2)
        ax.set_axisbelow(True);ax.yaxis.grid(True,color='#E4EBE7',linewidth=.5)
        bpoints=[]
        for name,face,edge in zip(ASSAYS,FACES,EDGES):
            rr=[next(r for r in DATA['example'] if r['assay']==name and r['adjusted']==adjusted) for adjusted in [False,True]]
            population=[r['population'] for r in rr];synthetic=[r['simulation'] for r in rr]
            ax.plot([0,1],population,color=edge,lw=1.05,marker='o',ms=5.3,mfc='white',mec=edge,mew=.95)
            ax.plot([0,1],synthetic,linestyle='none',marker='x',ms=3.9,color=edge,mew=.95)
            for j,r in enumerate(rr):
                bpoints.append(dict(assay=name,adjusted=r['adjusted'],x=j,population=r['population'],
                    simulation=r['simulation'],assumed_error_variance=r['error_variance']))
        ax.set(xlim=(-.18,1.18),ylim=(0,.14),xticks=[0,1],yticks=[0,.04,.08,.12],ylabel=self.t['ylabel'])
        ax.set_xticklabels(self.t['xticks']);ax.set_yticklabels(['0','0.04','0.08','0.12'])
        for label in ax.get_xticklabels()+ax.get_yticklabels()+[ax.xaxis.label,ax.yaxis.label]:
            label.set_fontproperties(self.font());label.set_color(INK)
        ax.yaxis.labelpad=6
        for x,name,psi,face,edge in zip([98,126,154],ASSAYS,[3,.1,4],FACES,EDGES):
            self.legend(x,133,name+'\nψ = '+str(psi),'o',face,edge)
        self.legend(98,145,self.t['population'],'o','white','#527F94')
        self.legend(138,145,self.t['synthetic'],'x','white','#527F94')
        self.primitives['B']=dict(points=bpoints,population_markers='open circles',synthetic_markers='crosses',
            synthetic_participants=200000,assumed_error_variances=[3,.1,4],
            x_categories=['Without p-tau181 adjustment','With p-tau181 adjustment'],y_limits=[0,.14],
            fixed_preadjustment_standardization=True,no_inference_intervals_added=True)

        for key,x,which,face,badge in [('C',30,'unadjusted','#FAF5E8','#B0C6A8'),
                                     ('D',123,'adjusted','#D0E3DB','#92B3CD')]:
            self.header(key,2 if key=='C' else 95,158,86 if key=='C' else 83,face,badge)
            plotted=[]
            for i,ref in enumerate(DATA['grid_levels']):
                yc=181.5+i*11
                self.text(x-4,yc,str(ref),8.8,ha='right')
                for j,comp in enumerate(DATA['grid_levels']):
                    xc=x+(j+.5)*11
                    r=next(r for r in DATA['pairwise_grid'] if r['reference_error_variance']==ref and r['comparison_error_variance']==comp)
                    value=r[which];colour=CMAP(NORM(value))
                    self.box(x+j*11,yc,11,11,colour,edge='white',lw=.5)
                    label=label_contrast(value)
                    self.text(xc,yc,label,8.8,ha='center',colour=contrast_ink(colour))
                    plotted.append(dict(reference_error_variance=ref,comparison_error_variance=comp,value=value,
                        label=label,source_multiplicity=5))
            for j,comp in enumerate(DATA['grid_levels']):
                self.text(x+(j+.5)*11,235.5,str(comp),8.8,ha='center')
            self.text(x+27.5,242,self.t['comparison'],8.8,ha='center')
            self.text(x-25,203.5,self.t['reference'],8.8,ha='center',rotation=90)
            self.primitives[key]=dict(cells=plotted,contrast='comparison minus reference',
                colour_limits=[-.075,.075],reference='NULISA',applies_equally_to=['SomaScan minus NULISA','TMT-MS minus NULISA'],
                original_three_assay_grid_rows=125,unique_pairwise_settings=25,unused_third_error_variance_collapsed=True)
        self.text(2,251.3,self.t['contrast'],8.8)
        for i in range(256):
            col=CMAP((i+.5)/256)
            self.box(63+i*115/256,251.3,115/256+(.025 if i<255 else 0),2.6,col,edge=col,lw=.05)
        for value in [-.075,0,.075]:
            x=63+(value+.075)/.15*115
            self.text(x,256,str(value).replace('-','−'),8.8,ha='left' if value<0 else ('right' if value>0 else 'center'))
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
        title='Supplementary Figure S15: Illustrative simulations of measurement error and p-tau181 adjustment'
        self.fig.savefig(BASE/f'Figure_S15_{self.lang}.pdf',metadata={'Title':title,'Creator':'Study supplementary figures'})
        self.fig.savefig(BASE/f'editable/Figure_S15_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_S15_{self.lang}.png',dpi=300)
        tif=BASE/f'high_resolution/Figure_S15_{self.lang}.tiff'
        self.fig.savefig(tif,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        with Image.open(tif) as im:
            rgb=im.convert('RGB');rgb.save(tif,compression='tiff_lzw',dpi=(600,600));rgb.close()
        (AUDIT/f'LAYOUT_{self.lang}.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
        (AUDIT/f'PLOT_PRIMITIVES_{self.lang}.json').write_text(json.dumps(self.primitives,ensure_ascii=False,indent=2),encoding='utf-8')
        return dict(language=self.lang,width_mm=W,height_mm=H,minimum_font_pt=min(r['point_size'] for r in records),
            checked_text_elements=len(records),text_clipped=False,text_overlap=False)

def main():
    warnings.filterwarnings('error',message='Glyph.*missing');results=[]
    for lang in ['EN','ZH']:
        fig=Figure(lang).draw();results.append(fig.save());plt.close(fig.fig)
    (AUDIT/'LAYOUT_VALIDATION.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    (AUDIT/'STYLE_SPECIFICATION.json').write_text(json.dumps(dict(palette=PALETTE,
        choice='Peach, sage and blue assay illustration with a shared blue-cream-rose contrast scale.',
        figure_size_mm=[W,H],panels=4,heatmap_cells=50,population_coefficients=6,synthetic_coefficients=6,
        contrast_scale=[-.075,.075],same_colour_scale_for_both_matrices=True,
        all_125_original_grid_rows_preserved_in_source_data=True,
        illustrative_parameters_not_empirical_assay_ranking=True,loose_footer_notes=False),indent=2),encoding='utf-8')
    print(json.dumps(results),flush=True)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');main()
