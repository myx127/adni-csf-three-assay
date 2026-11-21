"""Draw bilingual primary memory effects and distinct shared/difference evidence."""
from pathlib import Path
import json,logging,sys,warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from PIL import Image

logging.getLogger('fontTools').setLevel(logging.ERROR)
BASE=Path(__file__).resolve().parent.parent
AUDIT=BASE/'audit'
DATA=json.loads((AUDIT/'PLOT_INPUT.json').read_text(encoding='utf-8'))
W,H=180,237
PALETTE=['#DFF0F5','#A3C4DE','#AFD8C2','#7EA98B','#F6EBBA','#E9B0A6','#C0C4C1','#F5F8FA']
COLORS={'NULISA':'#A3C4DE','SomaScan':'#AFD8C2','TMT-MS':'#7EA98B'}
EDGES={'NULISA':'#537F9F','SomaScan':'#649879','TMT-MS':'#486D52'}
MARKERS={'NULISA':'o','SomaScan':'s','TMT-MS':'D'}
INK='#263D54';GRID='#E6ECEF';GRAY='#C0C4C1'
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.linewidth':.7,
    'xtick.major.width':.7,'ytick.major.width':.7,'xtick.major.size':3,'ytick.major.size':0,
    'axes.unicode_minus':False,'savefig.facecolor':'white','figure.facecolor':'white'})

TEXT={
 'EN':dict(figure='Supplementary Figure S6',D='Cross-assay\nevidence',protein='Protein',n='n',
     beta='β (score/year)',q='−log₁₀ q',pathology='Shared p-tau181 association (14)',
     shared='Shared direction',difference='Assay difference',cutoff='q = 0.05',unavailable='—†'),
 'ZH':dict(figure='补充图 S6',D='跨检测\n关联证据',protein='蛋白',n='n',
     beta='β（分/年）',q='−log₁₀ q',pathology='共同p-tau181关联（14个蛋白）',
     shared='共同方向',difference='检测间差异',cutoff='q = 0.05',unavailable='—†')}

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
    def header(self,letter,title,x,width,badge,fill):
        self.fig.patches.extend([
            Rectangle((x/W,1-23/H),width/W,11/H,transform=self.fig.transFigure,facecolor=fill,edgecolor='none',zorder=-1),
            Rectangle((x/W,1-23/H),5.8/W,11/H,transform=self.fig.transFigure,facecolor=badge,edgecolor='none')])
        self.text(x+2.9,17.5,letter,11.5,True,'center');self.text(x+7,17.5,title,9.5,True)
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
        columns=[('A','NULISA',32,32),('B','SomaScan',66.5,32),('C','TMT-MS',101,32)]
        for panel,assay,x,width in columns:self.header(panel,assay,x-1,width+1,COLORS[assay],'#F5F8FA')
        self.header('D',self.t['D'],137,41,'#7EA98B','#F6EBBA')
        self.fig.patches.append(Rectangle((2/W,1-29.6/H),4.8/W,2.8/H,transform=self.fig.transFigure,
            facecolor='#DFF0F5',edgecolor='none'))
        self.text(8,28.2,self.t['pathology'],8.8)
        self.fig.lines.extend([
            Line2D([81/W],[1-28.2/H],transform=self.fig.transFigure,linestyle='none',marker='o',markersize=3.6,
                markerfacecolor='white',markeredgecolor='#537F9F',markeredgewidth=.6),
            Line2D([112/W],[1-28.2/H],transform=self.fig.transFigure,linestyle='none',marker='D',markersize=3.5,
                markerfacecolor='#E9B0A6',markeredgecolor='#A96256',markeredgewidth=.6),
            Line2D([147/W,152/W],[1-28.2/H,1-28.2/H],transform=self.fig.transFigure,
                color='#687784',linestyle=(0,(3,2)),linewidth=.8)])
        self.text(84,28.2,self.t['shared'],8.8);self.text(115,28.2,self.t['difference'],8.8)
        self.text(154,28.2,self.t['cutoff'],8.8)
        self.text(2,35,self.t['protein'],9,True);self.text(29.3,35,self.t['n'],9,True,'right')
        yy=np.arange(34);focus=np.asarray(DATA['pathology_shared'],bool);qualified=np.asarray(DATA['qualified'],bool)
        for panel,assay,x,width in columns:
            ax=self.axis(panel,x,40,width,174)
            e=DATA['assays'][assay];beta=np.asarray(e['beta']);lo=np.asarray(e['ci_low']);hi=np.asarray(e['ci_high'])
            for i,highlight in enumerate(focus):
                if highlight:ax.axhspan(i-.48,i+.48,color='#DFF0F5',alpha=.48,zorder=0)
            ax.axvline(0,color='#8D9AA2',linestyle=(0,(3,2)),linewidth=.75,zorder=1)
            for available in [True,False]:
                use=qualified if available else ~qualified
                ax.errorbar(beta[use],yy[use],xerr=np.vstack([beta[use]-lo[use],hi[use]-beta[use]]),fmt=MARKERS[assay],
                    markersize=3.9,mfc=COLORS[assay] if available else 'white',mec=EDGES[assay] if available else '#9CA59F',
                    mew=.55,color=EDGES[assay] if available else GRAY,elinewidth=.8,capsize=1.9,capthick=.65,zorder=3)
            ax.set(xlim=(-.12,.22),ylim=(33.65,-.65),xticks=[-.1,0,.1,.2],yticks=[],xlabel=self.t['beta'])
            ax.xaxis.grid(True,color=GRID,linewidth=.5)
            self.primitives[panel]=dict(assay=assay,proteins=DATA['proteins'],beta=beta.tolist(),ci_low=lo.tolist(),ci_high=hi.tolist(),
                x_limits=[-.12,.22],y=yy.tolist(),pathology_highlighted=focus.tolist(),qualified=qualified.tolist(),descriptive_style='Open gray symbol and gray CI for CST3')
        anchor=self.axes['A']
        for i,(gene,n) in enumerate(zip(DATA['proteins'],DATA['participants'])):
            ypx=anchor.transData.transform((0,i))[1];y_mm=H-ypx/self.fig.dpi*25.4
            if focus[i]:self.fig.patches.append(Rectangle((1/W,1-(y_mm+2.4)/H),29/W,4.8/H,
                transform=self.fig.transFigure,facecolor='#DFF0F5',alpha=.48,edgecolor='none',zorder=-1))
            self.text(2,y_mm,gene+('†' if not qualified[i] else ''),8.8,bool(focus[i]))
            self.text(29.3,y_mm,str(n),8.8,False,'right')
        ax=self.axis('D',138,40,39,174)
        for i,highlight in enumerate(focus):
            if highlight:ax.axhspan(i-.48,i+.48,color='#DFF0F5',alpha=.48,zorder=0)
        ax.axvline(-np.log10(.05),color='#687784',linestyle=(0,(3,2)),linewidth=.8,zorder=1)
        groups=[]
        for j,(test,face,edge,marker) in enumerate([
            ('Shared direction','#A3C4DE','#537F9F','o'),('Assay differences','#E9B0A6','#A96256','D')]):
            q=np.array([np.nan if v is None else v for v in DATA['tests'][test]['q']]);show=qualified&np.isfinite(q)
            strength=-np.log10(q[show]);yp=yy[show]+(-.15 if j==0 else .15);passed=q[show]<=.05
            for value,y in zip(strength,yp):ax.plot([0,value],[y,y],color=face,linewidth=.65,zorder=2)
            ax.scatter(strength,yp,s=17,marker=marker,facecolors=[face if ok else 'white' for ok in passed],edgecolors=edge,linewidths=.6,zorder=3)
            groups.append(dict(test=test,proteins=[g for g,keep in zip(DATA['proteins'],show) if keep],q=q[show].tolist(),minus_log10_q=strength.tolist(),
                y=yp.tolist(),supported=passed.tolist(),unavailable_proteins=DATA['unavailable_proteins']))
        for i in np.flatnonzero(~qualified):
            ax.text(1.7,i,self.t['unavailable'],fontproperties=self.font(8.8),color='#7D8883',ha='center',va='center')
        ax.set(xlim=(-.15,3.25),ylim=(33.65,-.65),xticks=[0,1,2,3],yticks=[],xlabel=self.t['q'])
        ax.xaxis.grid(True,color=GRID,linewidth=.5)
        self.primitives['D']=dict(groups=groups,cutoff_q=.05,cutoff_minus_log10=float(-np.log10(.05)),display_truncated=False,
            family_size=68,unavailable_proteins=DATA['unavailable_proteins'],unavailable_placeholder=self.t['unavailable'])
        for axis in self.axes.values():self.style(axis)
        self.axes['D'].xaxis.label.set_fontproperties(self.font(9,math=True))
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
        title='Supplementary Figure S6: Associations of CSF proteins with annual memory change across three assays'
        self.fig.savefig(BASE/f'Figure_S6_{self.lang}.pdf',metadata={'Title':title,'Creator':'Study supplementary figures'})
        self.fig.savefig(BASE/f'editable/Figure_S6_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_S6_{self.lang}.png',dpi=300)
        tif=BASE/f'high_resolution/Figure_S6_{self.lang}.tiff';self.fig.savefig(tif,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
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
        user_reference='User-specified clear blue, bamboo green, pale yellow and coral palette for Figure S6.',
        figure_size_mm=[W,H],common_effect_axis=[-.12,.22],q_evidence_axis=[-.15,3.25],explanatory_footnotes_outside_plot=False,
        unqualified_CST3_shown_as_descriptive=True,unqualified_joint_tests_not_plotted_as_q1=True),ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(results,ensure_ascii=False),flush=True)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');main()
