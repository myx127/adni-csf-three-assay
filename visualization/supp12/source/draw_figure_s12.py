"""Draw direct paired changes, separating the original 9- and 33-test correction families."""
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
XLIM=[-.06,.1]
PALETTE=['#FAD6C9','#B0C6A8','#FAF5E8','#92B3CD','#D99DA3','#D0E3DB','#A6B3C0','#E8E2DB']
CONTRAST_COLORS={'SomaScan − NULISA':'#B0C6A8','TMT-MS − NULISA':'#92B3CD'}
EDGES={'SomaScan − NULISA':'#6F8C64','TMT-MS − NULISA':'#547F9E'}
MARKERS={'SomaScan − NULISA':'s','TMT-MS − NULISA':'D'}
INK='#263A46';ROSE='#A36872';GRID='#E9EEEB'
CMAP=LinearSegmentedColormap.from_list('ApricotRose',['#FAF5E8','#FAD6C9','#D99DA3'])
TEXT={
 'EN':dict(figure='Supplementary Figure S12',legendhead='Contrast\nchange',
   A='p-tau181 functional form',B='p-tau181 with background adjustment',
   families=['Nine comparisons (all shown)','33 comparisons (six shown)'],
   protein='Protein',n='n',term='Added p-tau181\nterm',background='Background\nscores',
   effect='Change in assay contrast',qhead='Joint change q',
   primary_rows=['Linear log₂','Linear,\nuntransformed','Spline log₂'],
   background_rows=['SomaScan +\nTMT-MS','All three\nassays'],
   units='Change in contrast\n(ADNI-MEM/year per protein residual SD)',q='−log₁₀ q'),
 'ZH':dict(figure='补充图 S12',legendhead='对比\n变化',
   A='p-tau181调整形式',B='背景调整后加入对数线性p-tau181',
   families=['9项家族（全部展示）','33项家族（展示6项）'],
   protein='蛋白',n='n',term='加入的p-tau181\n形式',background='背景\n分数',
   effect='检测间对比的变化',qhead='联合变化q',
   primary_rows=['对数线性','原值线性','对数样条'],
   background_rows=['SomaScan +\nTMT-MS','三种检测'],
   units='对比变化\n（ADNI-MEM分/年，每1蛋白残差SD）',q='−log₁₀ q')
}
PANEL_SPEC=[
 dict(panel='A',title_y=27,head_y=39,rows=[51,58.5,66,78,85.5,93,105,112.5,120],axis_top=46,axis_bottom=125,
      gene_centres=[58.5,85.5,112.5],group_bounds=[(46.5,70.5),(73.5,97.5),(100.5,124.5)],indices=[0,1,2]),
 dict(panel='B',title_y=149,head_y=160,rows=[173,181,194,202,215,223],axis_top=168,axis_bottom=228,
      gene_centres=[177,198,219],group_bounds=[(168.5,185.5),(189.5,206.5),(210.5,227.5)],indices=[3,4])
]
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.linewidth':.7,
    'xtick.major.width':.7,'ytick.major.width':.7,'xtick.major.size':3,'ytick.major.size':0,
    'axes.unicode_minus':False,'savefig.facecolor':'white','figure.facecolor':'white'})

def q_label(q):
    if q>=.001:return f'{q:.3g}'
    mantissa,exponent=f'{q:.2e}'.split('e')
    superscript=str(int(exponent)).translate(str.maketrans('0123456789-','⁰¹²³⁴⁵⁶⁷⁸⁹⁻'))
    return f'{mantissa} × 10{superscript}'

class Figure:
    def __init__(self,lang):
        self.lang=lang;self.t=TEXT[lang]
        self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=160)
        self.axes={};self.primitives={}
    def font(self,size=9,bold=False,math=False):
        name='segoeui.ttf' if math else (('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN' else ('msyhbd.ttc' if bold else 'msyh.ttc'))
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/name),size=size,weight='bold' if bold else 'normal')
    def text(self,x,y,label,size=9,bold=False,ha='left',math=False):
        return self.fig.text(x/W,1-y/H,label,fontproperties=self.font(size,bold,math),color=INK,ha=ha,va='center',linespacing=1.06)
    def box(self,x,y,width,height,fill,alpha=1,edge='none',lw=0):
        self.fig.patches.append(Rectangle((x/W,1-(y+height/2)/H),width/W,height/H,transform=self.fig.transFigure,
            facecolor=fill,alpha=alpha,edgecolor=edge,linewidth=lw,zorder=-1))
    def style(self,ax):
        for label in ax.get_xticklabels()+ax.get_yticklabels()+[ax.xaxis.label,ax.yaxis.label]:
            label.set_fontproperties(self.font(9));label.set_color(INK)
    def draw(self):
        self.text(2,5.7,self.t['figure'],13,True)
        self.text(2,15.5,self.t['legendhead'],8.8,True)
        for contrast,x in [('SomaScan − NULISA',35),('TMT-MS − NULISA',110)]:
            self.fig.lines.append(Line2D([x/W],[1-15.5/H],transform=self.fig.transFigure,linestyle='none',
                marker=MARKERS[contrast],markersize=4,markerfacecolor=CONTRAST_COLORS[contrast],
                markeredgecolor=EDGES[contrast],markeredgewidth=.6))
            self.text(x+3,15.5,contrast,9)
        for k,spec in enumerate(PANEL_SPEC):
            panel=spec['panel'];ty=spec['title_y'];hy=spec['head_y']
            self.box(2,ty,176,6.4,'#FAF5E8' if panel=='A' else '#E8E2DB')
            self.box(2,ty,5.5,6.4,'#FAD6C9' if panel=='A' else '#D99DA3')
            self.text(4.75,ty,panel,11.2,True,'center')
            self.text(9,ty,self.t[panel],9.5,True)
            self.text(177,ty,self.t['families'][k],8.8,False,'right')
            self.text(2,hy,self.t['protein'],8.8,True)
            self.text(36,hy,self.t['n'],8.8,True,'right')
            self.text(43,hy,self.t['term'] if panel=='A' else self.t['background'],8.8,True)
            self.text(111,hy,self.t['effect'],8.8,True,'center')
            self.text(162.5,hy,self.t['qhead'],8.8,True,'center')
            for j,(low,high) in enumerate(spec['group_bounds']):
                if j%2==0:self.box(2,(low+high)/2,176,high-low,'#D0E3DB',.22)
            order=[(gene,index) for gene in DATA['proteins'] for index in spec['indices']]
            rows=np.array(spec['rows'],float)
            top,bottom=spec['axis_top'],spec['axis_bottom']
            ax=self.fig.add_axes([80/W,1-bottom/H,62/W,(bottom-top)/H],facecolor='none')
            self.axes[panel]=ax
            ax.spines[['top','right','left']].set_visible(False);ax.spines['bottom'].set_color(INK)
            ax.tick_params(colors=INK,direction='out',pad=3);ax.set_axisbelow(True)
            for low,high in spec['group_bounds']:
                ax.vlines(0,low,high,color='#A6B3C0',linestyle=(0,(3,2)),linewidth=.8,zorder=1)
                for tick in [-.05,.05,.1]:ax.vlines(tick,low,high,color=GRID,linewidth=.5,zorder=0)
            groups=[]
            for j,contrast in enumerate(CONTRAST_COLORS):
                b=np.array([DATA['candidates'][gene]['contrasts'][contrast]['estimate'][index] for gene,index in order])
                lo=np.array([DATA['candidates'][gene]['contrasts'][contrast]['ci_low'][index] for gene,index in order])
                hi=np.array([DATA['candidates'][gene]['contrasts'][contrast]['ci_high'][index] for gene,index in order])
                yp=rows+(j-.5)*3
                ax.errorbar(b,yp,xerr=np.vstack([b-lo,hi-b]),fmt=MARKERS[contrast],markersize=4,
                    mfc=CONTRAST_COLORS[contrast],mec=EDGES[contrast],mew=.6,color=EDGES[contrast],elinewidth=.8,
                    capsize=1.6,capthick=.65,zorder=3)
                groups.append(dict(contrast=contrast,estimate=b.tolist(),ci_low=lo.tolist(),ci_high=hi.tolist(),y=yp.tolist()))
            ax.set(xlim=XLIM,ylim=(bottom,top),xticks=[-.05,0,.05,.1],yticks=[],xlabel=self.t['units'])
            ax.set_xticklabels(['−0.05','0','0.05','0.10']);self.style(ax)
            for gene,y in zip(DATA['proteins'],spec['gene_centres']):
                self.text(2,y,gene,9.2,True)
                self.text(36,y,str(DATA['candidates'][gene]['participants']),8.8,ha='right')
            q=[];labels=[]
            for (gene,index),y in zip(order,rows):
                label=self.t['primary_rows'][index] if panel=='A' else self.t['background_rows'][index-3]
                self.text(43,y,label,8.8,math=self.lang=='EN' and '₂' in label)
                value=DATA['candidates'][gene]['q'][index];q.append(value)
                self.box(147,y,31,5.4,CMAP(-np.log10(value)/8))
                if value<=.05:self.box(147,y,31,5.4,'none',edge=ROSE,lw=.85)
                label=q_label(value);labels.append(label)
                self.text(162.5,y,label,8.8,ha='center',math=True)
            self.primitives[panel]=dict(order=[dict(protein=gene,comparison=DATA['comparisons'][index]) for gene,index in order],
                participants=[DATA['candidates'][gene]['participants'] for gene,index in order],
                groups=groups,q=q,minus_log10_q=(-np.log10(q)).tolist(),displayed_q_text=labels,
                supported=[value<=.05 for value in q],original_family_size=9 if panel=='A' else 33,
                displayed_tests=len(order),other_tests_in_original_family=0 if panel=='A' else 27,
                x_limits=XLIM,q_colour_range=[0,8],q_colours_truncated=False,
                direction='Adjusted contrast minus reference contrast',cross_model_covariance_retained=True)
        cax=self.fig.add_axes([123/W,1-246.2/H,55/W,2.8/H]);self.axes['colourbar']=cax
        cb=self.fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,8),cmap=CMAP),cax=cax,orientation='horizontal',ticks=[0,2,4,6,8])
        cb.solids.set_rasterized(False);cb.solids.set_edgecolor('face')
        cb.outline.set_linewidth(.55);cb.outline.set_edgecolor(INK)
        cb.set_label(self.t['q'],labelpad=3);cax.tick_params(colors=INK,pad=2);self.style(cax)
        cax.xaxis.label.set_fontproperties(self.font(9,math=True))
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
        title='Supplementary Figure S12: Paired changes in between-assay memory association contrasts after p-tau181 adjustment'
        self.fig.savefig(BASE/f'Figure_S12_{self.lang}.pdf',metadata={'Title':title,'Creator':'Study supplementary figures'})
        self.fig.savefig(BASE/f'editable/Figure_S12_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_S12_{self.lang}.png',dpi=300)
        tif=BASE/f'high_resolution/Figure_S12_{self.lang}.tiff'
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
    (AUDIT/'STYLE_SPECIFICATION.json').write_text(json.dumps(dict(
        palette=PALETTE,contrast_colours=CONTRAST_COLORS,contrast_shapes=MARKERS,
        choice='Use sage squares for change in SomaScan minus NULISA and gray-blue diamonds for change in TMT-MS minus NULISA; rose for joint q.',
        figure_size_mm=[W,H],panels=2,common_change_axis=XLIM,original_family_sizes=[9,33],displayed_joint_tests=[9,6],
        q_colour_range=[0,8],q_colours_truncated=False,colourbar_vector=True,loose_footer_notes=False),ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(results,ensure_ascii=False),flush=True)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');main()
