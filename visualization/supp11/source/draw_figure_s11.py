"""Draw nine clinically interpretable adjustment models for three selected candidates."""
from pathlib import Path
import json, logging, sys, warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from PIL import Image

logging.getLogger('fontTools').setLevel(logging.ERROR)
BASE=Path(__file__).resolve().parent.parent
AUDIT=BASE/'audit'
DATA=json.loads((AUDIT/'PLOT_INPUT.json').read_text(encoding='utf-8'))
W,H=180,200
PALETTE=['#FAD6C9','#B0C6A8','#FAF5E8','#92B3CD','#D99DA3','#D0E3DB','#A6B3C0','#E8E2DB']
COLORS={'NULISA':'#FAD6C9','SomaScan':'#B0C6A8','TMT-MS':'#92B3CD'}
EDGES={'NULISA':'#B5806D','SomaScan':'#6F8C64','TMT-MS':'#547F9E'}
MARKERS={'NULISA':'o','SomaScan':'s','TMT-MS':'D'}
INK='#263A46';ROSE='#A36872';GRID='#E9EEEB'
CMAP=LinearSegmentedColormap.from_list('ApricotRose',['#FAF5E8','#FAD6C9','#D99DA3'])
ROWS=np.array([57,67,77,87,107,117,127,147,157.])
XLIM=[-.1,.21]
TEXT={
 'EN':dict(figure='Supplementary Figure S11',rowhead='p-tau181\nadjustment',effect='Memory\nassociation',qhead='Difference\nq',
   rowlabels=['None','Linear log₂','Linear,\nuntransformed','Spline log₂','None','Linear log₂','Spline log₂','None','Linear log₂'],
   groups=['No background scores','SomaScan + TMT-MS background scores','All three assay background scores'],
   units='ADNI-MEM/year\nper protein residual SD',q='−log₁₀ q'),
 'ZH':dict(figure='补充图 S11',rowhead='p-tau181\n调整方式',effect='记忆关联',qhead='检测间\n差异q',
   rowlabels=['不调整','对数线性','原值线性','对数样条','不调整','对数线性','对数样条','不调整','对数线性'],
   groups=['无背景分数','SomaScan和TMT-MS背景分数','三种检测的背景分数'],
   units='ADNI-MEM分/年\n每1蛋白残差SD',q='−log₁₀ q')
}
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.linewidth':.7,
    'xtick.major.width':.7,'ytick.major.width':.7,'xtick.major.size':3,'ytick.major.size':0,
    'axes.unicode_minus':False,'savefig.facecolor':'white','figure.facecolor':'white'})

class Figure:
    def __init__(self,lang):
        self.lang=lang;self.t=TEXT[lang]
        self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=160)
        self.axes={};self.primitives={}
    def font(self,size=9,bold=False,math=False):
        filename='segoeui.ttf' if math else (('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN' else ('msyhbd.ttc' if bold else 'msyh.ttc'))
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/filename),size=size,weight='bold' if bold else 'normal')
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
        for assay,x in [('NULISA',27),('SomaScan',76),('TMT-MS',126)]:
            self.fig.lines.append(Line2D([x/W],[1-15.5/H],transform=self.fig.transFigure,linestyle='none',marker=MARKERS[assay],
                markersize=3.7,markerfacecolor=COLORS[assay],markeredgecolor=EDGES[assay],markeredgewidth=.55))
            self.text(x+3.2,15.5,assay,9)
        self.text(2,38,self.t['rowhead'],8.8,True)
        for y,fill,label in zip([47,97,137],['#FAF5E8','#D0E3DB','#E8E2DB'],self.t['groups']):
            self.box(2,y,177,6.3,fill)
            self.text(3.7,y,label,9.2,True)
        for j,y in enumerate(ROWS):
            if j%2==0:self.box(2,y,177,9,'#D0E3DB',.18)
            self.text(3.7,y,self.t['rowlabels'][j],8.8,math=self.lang=='EN' and '₂' in self.t['rowlabels'][j])
        for panel,gene,x,badge in zip('ABC',DATA['proteins'],[33,83,133],['#FAD6C9','#B0C6A8','#92B3CD']):
            g=DATA['candidates'][gene]
            self.box(x,27,46,7.2,'#FAF5E8')
            self.box(x,27,5.5,7.2,badge)
            self.text(x+2.75,27,panel,11.2,True,'center')
            self.text(x+6.8,27,f'{gene}  (n={g["participants"]})',9,True)
            self.text(x+15,38,self.t['effect'],8.8,True,'center')
            self.text(x+39,38,self.t['qhead'],8.8,True,'center')
            ax=self.fig.add_axes([x/W,1-163/H,30/W,111/H],facecolor='none');self.axes[panel]=ax
            ax.spines[['top','right','left']].set_visible(False);ax.spines['bottom'].set_color(INK)
            ax.tick_params(colors=INK,direction='out',pad=3);ax.set_axisbelow(True)
            # Separate the three background groups without joining distinct models.
            for low,high in [(52,92),(102,132),(142,162)]:
                ax.vlines(0,low,high,color='#A6B3C0',linestyle=(0,(3,2)),linewidth=.75,zorder=1)
                for tick in [-.1,.1,.2]:ax.vlines(tick,low,high,color=GRID,linewidth=.5,zorder=0)
            groups=[]
            for k,assay in enumerate(COLORS):
                e=g['assays'][assay];b=np.asarray(e['beta']);lo=np.asarray(e['ci_low']);hi=np.asarray(e['ci_high'])
                yp=ROWS+(k-1)*2
                ax.errorbar(b,yp,xerr=np.vstack([b-lo,hi-b]),fmt=MARKERS[assay],markersize=3.5,
                    mfc=COLORS[assay],mec=EDGES[assay],mew=.55,color=EDGES[assay],elinewidth=.75,capsize=1.4,capthick=.6,zorder=3)
                groups.append(dict(assay=assay,models=DATA['models'],beta=b.tolist(),ci_low=lo.tolist(),ci_high=hi.tolist(),y=yp.tolist()))
            ax.set(xlim=XLIM,ylim=(163,52),xticks=[-.1,0,.1,.2],yticks=[],xlabel=self.t['units'])
            self.style(ax)
            labels=[]
            for q,y in zip(g['q'],ROWS):
                strength=-np.log10(q)
                self.box(x+32,y,14,6.4,CMAP(strength/4))
                if q<=.05:self.box(x+32,y,14,6.4,'none',edge=ROSE,lw=.8)
                label=f'{q:.3g}';labels.append(label)
                self.text(x+39,y,label,8.8,ha='center',math=True)
            self.primitives[panel]=dict(protein=gene,participants=g['participants'],groups=groups,
                q=g['q'],minus_log10_q=(-np.log10(g['q'])).tolist(),displayed_q_text=labels,supported=[q<=.05 for q in g['q']],
                original_family_size=33,x_limits=XLIM,all_models_same_cases_observations_and_protein_SD=True,
                q_colour_range=[0,4],q_colour_truncated=False)
        cax=self.fig.add_axes([124/W,1-188/H,55/W,2.8/H]);self.axes['colourbar']=cax
        cb=self.fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,4),cmap=CMAP),cax=cax,orientation='horizontal',ticks=[0,1,2,3,4])
        cb.solids.set_rasterized(False);cb.solids.set_edgecolor('face')
        cb.outline.set_linewidth(.55);cb.outline.set_edgecolor(INK);cb.set_label(self.t['q'],labelpad=3)
        cax.tick_params(colors=INK,pad=2);self.style(cax)
        cax.xaxis.label.set_fontproperties(self.font(9,math=True))
        self.fig.canvas.draw()
        return self
    def layout(self):
        self.fig.canvas.draw();renderer=self.fig.canvas.get_renderer();records=[];seen=set()
        for artist in self.fig.findobj(Text):
            if not artist.get_visible() or not artist.get_text().strip() or id(artist) in seen:continue
            seen.add(id(artist));box=artist.get_window_extent(renderer)
            mm=np.array([box.x0,box.y0,box.x1,box.y1])/self.fig.dpi*25.4
            assert mm[0]>=-.15 and mm[2]<=W+.15 and mm[1]>=-.15 and mm[3]<=H+.15,('clipped text',self.lang,artist.get_text(),mm.tolist())
            records.append(dict(text=artist.get_text(),point_size=float(artist.get_fontsize()),x0_mm=float(mm[0]),y0_mm=float(H-mm[3]),x1_mm=float(mm[2]),y1_mm=float(H-mm[1])))
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
        title='Supplementary Figure S11: Memory associations under alternative p-tau181 and background adjustment models'
        self.fig.savefig(BASE/f'Figure_S11_{self.lang}.pdf',metadata={'Title':title,'Creator':'Study supplementary figures'})
        self.fig.savefig(BASE/f'editable/Figure_S11_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_S11_{self.lang}.png',dpi=300)
        tif=BASE/f'high_resolution/Figure_S11_{self.lang}.tiff'
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
        palette=PALETTE,assay_colours=COLORS,assay_shapes=MARKERS,
        choice='Continue the apricot, sage and gray-blue assay mapping from Figures S8–S9; rose shows within-model difference q.',
        figure_size_mm=[W,H],panels=3,common_effect_axis=XLIM,model_groups=[4,3,2],original_family_size=33,
        q_colour_range=[0,4],q_colours_truncated=False,colourbar_vector=True,loose_footer_notes=False),ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(results,ensure_ascii=False),flush=True)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');main()
