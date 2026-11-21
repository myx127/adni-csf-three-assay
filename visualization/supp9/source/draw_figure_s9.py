"""Draw four comparable memory-effect panels and evidence with the original correction families."""
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
PALETTE=['#FAD6C9','#B0C6A8','#FAF5E8','#92B3CD','#D99DA3','#D0E3DB','#A6B3C0','#E8E2DB']
COLORS={'NULISA':'#FAD6C9','SomaScan':'#B0C6A8','TMT-MS':'#92B3CD'}
EDGES={'NULISA':'#B5806D','SomaScan':'#6F8C64','TMT-MS':'#547F9E'}
MARKERS={'NULISA':'o','SomaScan':'s','TMT-MS':'D'}
INK='#263A46';GRID='#E9EEEB';ROSE='#A36872'
CMAP=LinearSegmentedColormap.from_list('StudyPeach',['#FAF5E8','#FAD6C9','#D99DA3'])
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.linewidth':.7,
    'xtick.major.width':.7,'ytick.major.width':.7,'xtick.major.size':3,'ytick.major.size':0,
    'axes.unicode_minus':False,'savefig.facecolor':'white','figure.facecolor':'white'})

TEXT={
 'EN':dict(figure='Supplementary Figure S9',A='Primary\nmodel',B='No p-tau181\nadjustment',C='≥2 follow-ups\nspan ≥2 years',D='Equal participant\nweighting',
     E='Cross-assay evidence across sensitivity analyses',protein='Protein',n='n',units='ADNI-MEM/year\nper protein residual SD',
     groups=['Primary model\n68 tests','No p-tau181\n28 tests','≥2 follow-ups; ≥2 years\n28 tests','Equal weighting\n28 tests'],
     qcols=['Shared','Diff.']*4,q='−log₁₀ q'),
 'ZH':dict(figure='补充图 S9',A='主模型',B='不调整\np-tau181',C='≥2次随访\n且跨度≥2年',D='个体\n等权重',
     E='各敏感性分析的跨检测关联证据',protein='蛋白',n='n',units='ADNI-MEM分/年\n每1蛋白残差SD',
     groups=['主模型\n68项检验','不调整p-tau181\n28项检验','≥2次随访且跨度≥2年\n28项检验','个体等权重\n28项检验'],
     qcols=['共同','差异']*4,q='−log₁₀ q')}

class Figure:
    def __init__(self,lang):
        self.lang=lang;self.t=TEXT[lang];self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=160)
        self.axes={};self.primitives={}
    def font(self,size=9,bold=False,math=False):
        if math:return FontProperties(fname='C:/Windows/Fonts/segoeui.ttf',size=size)
        filename=('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN' else ('msyhbd.ttc' if bold else 'msyh.ttc')
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/filename),size=size,weight='bold' if bold else 'normal')
    def text(self,x,y,value,size=9,bold=False,ha='left'):
        return self.fig.text(x/W,1-y/H,value,fontproperties=self.font(size,bold),color=INK,ha=ha,va='center',linespacing=1.05)
    def header(self,key,x,y,width,badge,fill,tall=False):
        height=10 if tall else 6.4
        self.fig.patches.extend([
            Rectangle((x/W,1-(y+height/2)/H),width/W,height/H,transform=self.fig.transFigure,facecolor=fill,edgecolor='none',zorder=-1),
            Rectangle((x/W,1-(y+height/2)/H),5.5/W,height/H,transform=self.fig.transFigure,facecolor=badge,edgecolor='none')])
        self.text(x+2.75,y,key,11.2,True,'center');self.text(x+6.9,y,self.t[key],9.2 if tall else 10,True)
    def axis(self,key,x,y,w,h):
        ax=self.fig.add_axes([x/W,1-(y+h)/H,w/W,h/H]);self.axes[key]=ax
        ax.spines[['top','right','left']].set_visible(False);ax.spines['bottom'].set_color(INK)
        ax.tick_params(colors=INK,direction='out',pad=3);ax.set_axisbelow(True)
        return ax
    def style(self,ax):
        for label in ax.get_xticklabels()+ax.get_yticklabels():label.set_fontproperties(self.font(9));label.set_color(INK)
        for label in [ax.xaxis.label,ax.yaxis.label]:label.set_fontproperties(self.font(9));label.set_color(INK)
    def draw(self):
        self.text(2,5.7,self.t['figure'],13,True)
        for assay,x in [('NULISA',27),('SomaScan',76),('TMT-MS',126)]:
            self.fig.lines.append(Line2D([x/W],[1-15.5/H],transform=self.fig.transFigure,linestyle='none',marker=MARKERS[assay],
                markersize=3.7,markerfacecolor=COLORS[assay],markeredgecolor=EDGES[assay],markeredgewidth=.55))
            self.text(x+3.2,15.5,assay,9)
        configs=[('A',18,'#FAD6C9','#FAF5E8'),('B',58,'#B0C6A8','#D0E3DB'),('C',98,'#92B3CD','#FAF5E8'),('D',138,'#D99DA3','#E8E2DB')]
        yy=np.arange(14)
        self.text(2,37.7,self.t['protein'],9,True)
        for (panel,x,badge,fill),analysis in zip(configs,DATA['analysis_order']):
            self.header(panel,x,27,39,badge,fill,True)
            self.text(x+7,37.7,self.t['n'],9,True,'right')
            ax=self.axis(panel,x+8,43,31,83)
            for i in range(14):
                if i%2==0:ax.axhspan(i-.48,i+.48,color='#D0E3DB',alpha=.2,zorder=0)
            ax.axvline(0,color='#A6B3C0',linestyle=(0,(3,2)),linewidth=.75,zorder=1)
            groups=[];a=DATA['analyses'][analysis]
            for j,assay in enumerate(COLORS):
                e=a['assays'][assay];beta=np.asarray(e['beta']);low=np.asarray(e['ci_low']);high=np.asarray(e['ci_high']);yp=yy+(j-1)*.22
                ax.errorbar(beta,yp,xerr=np.vstack([beta-low,high-beta]),fmt=MARKERS[assay],markersize=3.2,
                    mfc=COLORS[assay],mec=EDGES[assay],mew=.5,color=EDGES[assay],elinewidth=.7,capsize=1.3,capthick=.6,zorder=3)
                groups.append(dict(assay=assay,proteins=DATA['proteins'],beta=beta.tolist(),ci_low=low.tolist(),ci_high=high.tolist(),y=yp.tolist()))
            ax.set(xlim=(-.13,.29),ylim=(13.5,-.5),xticks=[-.1,0,.1,.2],yticks=[],xlabel=self.t['units'])
            ax.xaxis.grid(True,color=GRID,linewidth=.5);self.style(ax)
            for i,n in enumerate(a['participants']):
                y=43+(i+.5)*83/14
                self.text(x+7,y,str(n),8.8,False,'right')
            self.primitives[panel]=dict(analysis=analysis,participants=a['participants'],groups=groups,x_limits=[-.13,.29],
                original_family_size=a['family_size'],adjusts_ptau=a['adjusts_ptau'],same_original_protein_sd=True)
        for i,gene in enumerate(DATA['proteins']):
            self.text(2,43+(i+.5)*83/14,gene,8.8)
        self.header('E',2,151.5,176,'#D99DA3','#E8E2DB')
        columns=[(a,t) for a in DATA['analysis_order'] for t in ['Shared direction','Assay differences']]
        q=np.column_stack([DATA['analyses'][a]['tests'][t]['q'] for a,t in columns]);strength=-np.log10(q)
        ax=self.fig.add_axes([20/W,1-237/H,156/W,56/H]);self.axes['E']=ax
        labels=[]
        for i in range(14):
            row=[]
            for j in range(8):
                ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor=CMAP(strength[i,j]/3),edgecolor='white',linewidth=.7))
                if q[i,j]<=.05:ax.add_patch(Rectangle((j-.47,i-.44),.94,.88,facecolor='none',edgecolor=ROSE,linewidth=.8))
                label=f'{q[i,j]:.3g}';row.append(label)
                ax.text(j,i,label,fontproperties=self.font(8.8,math=True),color=INK,ha='center',va='center',zorder=3)
            labels.append(row)
        ax.set(xlim=(-.5,7.5),ylim=(13.5,-.5),xticks=[],yticks=[])
        for spine in ax.spines.values():spine.set_visible(False)
        for j,label in enumerate(self.t['groups']):self.text(39.5+39*j,162.5,label,8.8,False,'center')
        for j,label in enumerate(self.t['qcols']):self.text(20+19.5*(j+.5),174,label,8.8,False,'center')
        for i,gene in enumerate(DATA['proteins']):self.text(2,181+(i+.5)*4,gene,8.8)
        self.primitives['E']=dict(proteins=DATA['proteins'],columns=[dict(analysis=a,test=t) for a,t in columns],
            q=q.tolist(),minus_log10_q=strength.tolist(),displayed_text=labels,supported=(q<=.05).tolist(),
            original_family_sizes=DATA['family_sizes'],family_count=4,displayed_test_count=112,colour_range=[0,3],colours_truncated=False)
        cax=self.fig.add_axes([120/W,1-245.8/H,56/W,2.8/H]);self.axes['colourbar']=cax
        cb=self.fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,3),cmap=CMAP),cax=cax,orientation='horizontal',ticks=[0,1,2,3])
        cb.solids.set_rasterized(False);cb.solids.set_edgecolor('face');cb.outline.set_linewidth(.55);cb.outline.set_edgecolor(INK)
        cb.set_label(self.t['q'],labelpad=3);cax.tick_params(colors=INK,pad=2);self.style(cax)
        cax.xaxis.label.set_fontproperties(self.font(9,math=True))
        self.fig.canvas.draw();return self
    def layout(self):
        self.fig.canvas.draw();renderer=self.fig.canvas.get_renderer();records=[];seen=set()
        for artist in self.fig.findobj(Text):
            if not artist.get_visible() or not artist.get_text().strip() or id(artist) in seen:continue
            seen.add(id(artist));box=artist.get_window_extent(renderer);mm=np.array([box.x0,box.y0,box.x1,box.y1])/self.fig.dpi*25.4
            assert mm[0]>=-.15 and mm[2]<=W+.15 and mm[1]>=-.15 and mm[3]<=H+.15,('text outside canvas',self.lang,artist.get_text(),mm.tolist())
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
        title='Supplementary Figure S9: Sensitivity of memory change associations among fourteen proteins with shared p-tau181 associations'
        self.fig.savefig(BASE/f'Figure_S9_{self.lang}.pdf',metadata={'Title':title,'Creator':'Study supplementary figures'})
        self.fig.savefig(BASE/f'editable/Figure_S9_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_S9_{self.lang}.png',dpi=300)
        tif=BASE/f'high_resolution/Figure_S9_{self.lang}.tiff';self.fig.savefig(tif,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
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
        choice='Continue Figure S8 apricot, sage and gray-blue assay colors for consistency; muted rose for the q evidence scale.',
        figure_size_mm=[W,H],panels=5,common_effect_axis=[-.13,.29],original_family_sizes=[68,28,28,28],
        colourbar_vector=True,explanatory_footnotes_outside_plot=False),ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(results,ensure_ascii=False),flush=True)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');main()
