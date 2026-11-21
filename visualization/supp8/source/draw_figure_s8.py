"""Draw two distinct secondary cognitive outcomes and their original shared evidence family."""
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
W,H=180,252
PALETTE=['#FAD6C9','#B0C6A8','#FAF5E8','#92B3CD','#D99DA3','#D0E3DB','#A6B3C0','#E8E2DB']
COLORS={'NULISA':'#FAD6C9','SomaScan':'#B0C6A8','TMT-MS':'#92B3CD'}
EDGES={'NULISA':'#B5806D','SomaScan':'#6F8C64','TMT-MS':'#547F9E'}
MARKERS={'NULISA':'o','SomaScan':'s','TMT-MS':'D'}
INK='#263A46';GRID='#E8ECEA';ROSE='#A36872'
CMAP=LinearSegmentedColormap.from_list('StudyPeach',['#FAF5E8','#FAD6C9','#D99DA3'])
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.linewidth':.7,
    'xtick.major.width':.7,'ytick.major.width':.7,'xtick.major.size':3,'ytick.major.size':0,
    'axes.unicode_minus':False,'savefig.facecolor':'white','figure.facecolor':'white'})

TEXT={
 'EN':dict(figure='Supplementary Figure S8',A='Baseline\nmemory',B='Executive function\nchange',C='Cross-assay\nevidence',
     protein='Protein',n='n',pathology='Shared p-tau181 (14)',baseline_x='ADNI-MEM score\nper protein residual SD',executive_x='ADNI-EF/year\nper protein residual SD',
     groups=['Baseline\nmemory','Executive\nfunction change'],qcols=['Shared','Diff.','Shared','Diff.'],q='−log₁₀ q'),
 'ZH':dict(figure='补充图 S8',A='基线\n记忆',B='执行功能\n年变化',C='跨检测\n关联证据',
     protein='蛋白',n='n',pathology='共同p-tau181（14）',baseline_x='ADNI-MEM分数\n每1蛋白残差SD',executive_x='ADNI-EF分/年\n每1蛋白残差SD',
     groups=['基线\n记忆','执行功能\n年变化'],qcols=['共同','差异','共同','差异'],q='−log₁₀ q')}

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
    def header(self,key,x,width,badge,fill):
        self.fig.patches.extend([
            Rectangle((x/W,1-23/H),width/W,11/H,transform=self.fig.transFigure,facecolor=fill,edgecolor='none',zorder=-1),
            Rectangle((x/W,1-23/H),5.8/W,11/H,transform=self.fig.transFigure,facecolor=badge,edgecolor='none')])
        self.text(x+2.9,17.5,key,11.5,True,'center');self.text(x+7,17.5,self.t[key],9.5,True)
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
        self.header('A',29,46,'#FAD6C9','#FAF5E8')
        self.header('B',80,46,'#92B3CD','#D0E3DB')
        self.header('C',132,46,'#D99DA3','#E8E2DB')
        self.fig.patches.append(Rectangle((2/W,1-31.2/H),4.8/W,2.8/H,transform=self.fig.transFigure,facecolor='#D0E3DB',edgecolor='none'))
        self.text(8,29.8,self.t['pathology'],8.8)
        for assay,x in [('NULISA',45),('SomaScan',70),('TMT-MS',98)]:
            self.fig.lines.append(Line2D([x/W],[1-29.8/H],transform=self.fig.transFigure,linestyle='none',marker=MARKERS[assay],
                markersize=3.5,markerfacecolor=COLORS[assay],markeredgecolor=EDGES[assay],markeredgewidth=.55))
            self.text(x+3,29.8,assay,8.8)
        self.text(2,40,self.t['protein'],9,True);self.text(26.5,40,self.t['n'],9,True,'right')
        focus=np.asarray(DATA['pathology_shared'],bool);yy=np.arange(34)
        specs=[('A','Baseline memory',30,[-.36,.69],[-.3,0,.3,.6],self.t['baseline_x']),
               ('B','Executive function change',81,[-.11,.23],[-.1,0,.1,.2],self.t['executive_x'])]
        for panel,outcome,x,limits,ticks,xlabel in specs:
            ax=self.axis(panel,x,46,44,177)
            for i,highlight in enumerate(focus):
                if highlight:ax.axhspan(i-.48,i+.48,color='#D0E3DB',alpha=.38,zorder=0)
            ax.axvline(0,color='#A6B3C0',linestyle=(0,(3,2)),linewidth=.8,zorder=1)
            groups=[]
            for j,assay in enumerate(COLORS):
                e=DATA['outcomes'][outcome]['assays'][assay]
                beta=np.asarray(e['beta']);low=np.asarray(e['ci_low']);high=np.asarray(e['ci_high']);yp=yy+(j-1)*.22
                ax.errorbar(beta,yp,xerr=np.vstack([beta-low,high-beta]),fmt=MARKERS[assay],markersize=3.1,
                    mfc=COLORS[assay],mec=EDGES[assay],mew=.5,color=EDGES[assay],elinewidth=.7,capsize=1.2,capthick=.6,zorder=3)
                groups.append(dict(assay=assay,proteins=DATA['proteins'],beta=beta.tolist(),ci_low=low.tolist(),ci_high=high.tolist(),y=yp.tolist()))
            ax.set(xlim=limits,ylim=(33.5,-.5),xticks=ticks,yticks=[],xlabel=xlabel)
            ax.xaxis.grid(True,color=GRID,linewidth=.5);self.style(ax)
            self.primitives[panel]=dict(outcome=outcome,groups=groups,x_limits=limits,pathology_highlighted=focus.tolist(),
                units='ADNI-MEM score per fixed protein residual SD' if panel=='A' else 'ADNI-EF score/year per fixed protein residual SD',all_designs_qualified=True)
        for i,(gene,n) in enumerate(zip(DATA['proteins'],DATA['participants'])):
            y=46+(i+.5)*177/34
            if focus[i]:self.fig.patches.append(Rectangle((1/W,1-(y+2.4)/H),27/W,4.8/H,
                transform=self.fig.transFigure,facecolor='#D0E3DB',alpha=.38,edgecolor='none',zorder=-1))
            self.text(2,y,gene,8.8,bool(focus[i]));self.text(26.5,y,str(n),8.8,False,'right')
        # Four columns share one original 136-hypothesis family, rather than four separate families.
        columns=[(o,t) for o in ['Baseline memory','Executive function change'] for t in ['Shared direction','Assay differences']]
        q=np.column_stack([DATA['outcomes'][o]['tests'][t]['q'] for o,t in columns]);strength=-np.log10(q)
        ax=self.fig.add_axes([133/W,1-223/H,44/W,177/H]);self.axes['C']=ax
        labels=[]
        for i in range(34):
            row=[]
            for j in range(4):
                ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor=CMAP(strength[i,j]/2),edgecolor='white',linewidth=.7))
                if q[i,j]<=.05:ax.add_patch(Rectangle((j-.46,i-.45),.92,.9,facecolor='none',edgecolor=ROSE,linewidth=.9))
                label=f'{q[i,j]:.3g}';row.append(label)
                ax.text(j,i,label,fontproperties=self.font(8.8,math=True),color=INK,ha='center',va='center',zorder=3)
            labels.append(row)
        ax.set(xlim=(-.5,3.5),ylim=(33.5,-.5),xticks=[],yticks=[])
        for spine in ax.spines.values():spine.set_visible(False)
        for j,label in enumerate(self.t['groups']):self.text(144+22*j,29.5,label,8.8,False,'center')
        for j,label in enumerate(self.t['qcols']):self.text(133+11*(j+.5),39.5,label,8.8,False,'center')
        self.primitives['C']=dict(proteins=DATA['proteins'],columns=[dict(outcome=o,test=t) for o,t in columns],q=q.tolist(),minus_log10_q=strength.tolist(),
            displayed_text=labels,supported=(q<=.05).tolist(),family_size=136,colour_range=[0,2],colours_truncated=False,all_designs_qualified=True)
        cax=self.fig.add_axes([133/W,1-237.8/H,44/W,2.8/H]);self.axes['colourbar']=cax
        cb=self.fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,2),cmap=CMAP),cax=cax,orientation='horizontal',ticks=[0,1,2])
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
        title='Supplementary Figure S8: CSF protein associations with baseline memory and annual executive function change'
        self.fig.savefig(BASE/f'Figure_S8_{self.lang}.pdf',metadata={'Title':title,'Creator':'Study supplementary figures'})
        self.fig.savefig(BASE/f'editable/Figure_S8_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_S8_{self.lang}.png',dpi=300)
        tif=BASE/f'high_resolution/Figure_S8_{self.lang}.tiff';self.fig.savefig(tif,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
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
        user_reference='User-specified apricot, pale olive, cream, gray blue, muted rose, pale green and warm gray palette for Figure S8.',
        figure_size_mm=[W,H],panels=3,baseline_effect_axis=[-.36,.69],annual_executive_effect_axis=[-.11,.23],
        shared_secondary_family=136,colourbar_vector=True,explanatory_footnotes_outside_plot=False,CST3_qualified_for_both_outcomes=True),ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(results,ensure_ascii=False),flush=True)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');main()
