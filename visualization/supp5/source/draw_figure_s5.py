"""Draw joint-evidence sensitivity and assay-specific coefficient comparisons."""
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
PALETTE=['#35594F','#668578','#A5BBAA','#E0E9DF','#6A627F','#9D93AD','#D4CDDF','#F7F6F9']
INK='#2D493F';PURPLE='#6A627F';GRID='#E8EAE9'
MODELS=['Primary','APOE','APOE and education']
CMAP=LinearSegmentedColormap.from_list('StudyGreen',['#F7F6F9','#E0E9DF','#A5BBAA','#668578','#35594F'])
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
    'axes.linewidth':.7,'xtick.major.width':.7,'ytick.major.width':.7,'xtick.major.size':3,'ytick.major.size':3,
    'axes.unicode_minus':False,'savefig.facecolor':'white','figure.facecolor':'white'})

TEXT={
 'EN':dict(figure='Supplementary Figure S5',A='Shared association evidence',
     headings=['Primary\nmodel','APOE ε4','APOE ε4\n+ education'],threshold='q ≤ 0.05',
     legend=['APOE ε4','APOE ε4 + education'],primary_beta='Primary-model β',adjusted_beta='Adjusted-model β',q='−log₁₀ q'),
 'ZH':dict(figure='补充图 S5',A='共同关联证据的调整敏感性',
     headings=['主模型','APOE ε4','APOE ε4\n＋教育'],threshold='q≤0.05',
     legend=['APOE ε4','APOE ε4＋教育'],primary_beta='主模型β',adjusted_beta='调整后模型β',q='−log₁₀ q')}

def qlabel(value):
    if value==1:return '1'
    if value>=.001:return f'{value:.3g}'
    mantissa,exponent=f'{value:.2e}'.split('e')
    superscript=str(int(exponent)).translate(str.maketrans('-0123456789','⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))
    return mantissa+'×10'+superscript

class Figure:
    def __init__(self,lang):
        self.lang=lang;self.t=TEXT[lang]
        self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=160)
        self.axes={};self.primitives={}
    def font(self,size=9,bold=False,math=False):
        if math:return FontProperties(fname='C:/Windows/Fonts/segoeui.ttf',size=size)
        filename=('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN' else ('msyhbd.ttc' if bold else 'msyh.ttc')
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/filename),size=size,weight='bold' if bold else 'normal')
    def text(self,x,y,value,size=9,bold=False,ha='left'):
        return self.fig.text(x/W,1-y/H,value,fontproperties=self.font(size,bold),color=INK,ha=ha,va='center',linespacing=1.05)
    def header(self,key,title,x,y,width,badge,fill):
        self.fig.patches.extend([
            Rectangle((x/W,1-(y+3.2)/H),width/W,6.4/H,transform=self.fig.transFigure,facecolor=fill,edgecolor='none',zorder=-1),
            Rectangle((x/W,1-(y+3.2)/H),5.8/W,6.4/H,transform=self.fig.transFigure,facecolor=badge,edgecolor='none')])
        self.fig.text((x+2.9)/W,1-y/H,key,fontproperties=self.font(11.5,True),color='white',ha='center',va='center')
        self.text(x+7.2,y,title,9.8,True)
    def axis(self,key,x,y,w,h):
        ax=self.fig.add_axes([x/W,1-(y+h)/H,w/W,h/H]);self.axes[key]=ax
        ax.spines[['top','right']].set_visible(False)
        for spine in ax.spines.values():spine.set_color(INK)
        ax.tick_params(colors=INK,direction='out',pad=3);ax.set_axisbelow(True)
        return ax
    def style(self,ax):
        for label in ax.get_xticklabels()+ax.get_yticklabels():label.set_fontproperties(self.font(9));label.set_color(INK)
        for label in [ax.xaxis.label,ax.yaxis.label]:label.set_fontproperties(self.font(9));label.set_color(INK)
    def draw(self):
        self.text(2,5.7,self.t['figure'],13,True)
        self.header('A',self.t['A'],2,17.8,84,'#35594F','#E0E9DF')
        self.header('B','NULISA',96,17.8,82,'#35594F','#E0E9DF')
        self.header('C','SomaScan',96,99.8,82,'#6A627F','#D4CDDF')
        self.header('D','TMT-MS',96,181.8,82,'#668578','#E0E9DF')
        self.fig.patches.append(Rectangle((2/W,1-31.3/H),3.6/W,3.4/H,transform=self.fig.transFigure,
            facecolor='white',edgecolor=PURPLE,linewidth=.75))
        self.text(6.5,29.5,self.t['threshold'],8.8)
        for j,label in enumerate(self.t['headings']):self.text(23+20.8*(j+.5),29.5,label,8.8,False,'center')
        handles=[Line2D([],[],linestyle='none',marker='o',markersize=5.2,markerfacecolor='white',markeredgecolor='#668578',markeredgewidth=.9),
                 Line2D([],[],linestyle='none',marker='^',markersize=4.7,markerfacecolor='#9D93AD',markeredgecolor=PURPLE,markeredgewidth=.65)]
        legend=self.fig.legend(handles,self.t['legend'],loc='center left',bbox_to_anchor=(97/W,1-28.2/H),frameon=False,ncol=2,
            prop=self.font(8.8),handlelength=.75,handletextpad=.3,columnspacing=.9,borderaxespad=0)
        for item in legend.get_texts():item.set_color(INK)
        # The three q-value columns retain their original, separate 34-protein BY families.
        ax=self.fig.add_axes([23/W,1-225/H,62.4/W,187/H]);self.axes['A']=ax
        q=np.column_stack([DATA['models'][m]['q'] for m in MODELS]);strength=-np.log10(q)
        support=q<=.05
        # Native vector rectangles retain editable colours and text in SVG/PDF.
        cell_text=[]
        for i in range(34):
            line=[]
            for j in range(3):
                face=CMAP(strength[i,j]/20)
                ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor=face,edgecolor='white',linewidth=.75))
                if support[i,j]:ax.add_patch(Rectangle((j-.47,i-.45),.94,.9,facecolor='none',edgecolor=PURPLE,linewidth=.75))
                label=qlabel(q[i,j]);line.append(label)
                ax.text(j,i,label,fontproperties=self.font(8.8,math=True),ha='center',va='center',color='white' if strength[i,j]>=15 else INK,zorder=3)
            cell_text.append(line)
        ax.set(xlim=(-.5,2.5),ylim=(33.5,-.5),xticks=[],yticks=[])
        for spine in ax.spines.values():spine.set_visible(False)
        for i,gene in enumerate(DATA['proteins']):
            y_mm=38+(i+.5)*187/34
            self.text(21.3,y_mm,gene,8.8,gene in DATA['supported_proteins'],'right')
        self.primitives['A']=dict(proteins=DATA['proteins'],models=MODELS,q=q.tolist(),minus_log10_q=strength.tolist(),
            supported=support.tolist(),displayed_text=cell_text,colours_truncated=False,colour_range=[0,20])
        cax=self.fig.add_axes([23/W,1-238.5/H,62.4/W,3/H]);self.axes['colourbar']=cax
        sm=plt.cm.ScalarMappable(norm=Normalize(0,20),cmap=CMAP)
        cb=self.fig.colorbar(sm,cax=cax,orientation='horizontal',ticks=[0,5,10,15,20])
        cb.solids.set_rasterized(False)
        cb.solids.set_edgecolor('face')
        cb.set_label(self.t['q'],labelpad=3);cb.outline.set_linewidth(.55);cb.outline.set_edgecolor(INK)
        cax.tick_params(colors=INK,pad=2);self.style(cax)
        cax.xaxis.label.set_fontproperties(self.font(9,math=True))
        specs=[('B','NULISA',39,[-.15,1.05],[0,.5,1]),
               ('C','SomaScan',121,[-.6,2.1],[-.5,0,1,2]),
               ('D','TMT-MS',203,[-.45,.95],[-.4,0,.4,.8])]
        for panel,assay,top,limits,ticks in specs:
            ax=self.axis(panel,119,top,54,43)
            base=np.asarray(DATA['assays'][assay]['Primary']['beta'])
            ax.plot(limits,limits,linestyle=(0,(4,3)),color='#A5BBAA',linewidth=.9,zorder=1)
            groups=[]
            for model in ['APOE','APOE and education']:
                adjusted=np.asarray(DATA['assays'][assay][model]['beta'])
                is_apoe=model=='APOE'
                ax.scatter(base,adjusted,s=30 if is_apoe else 19,marker='o' if is_apoe else '^',
                    facecolor='white' if is_apoe else '#9D93AD',edgecolors='#668578' if is_apoe else PURPLE,
                    linewidths=.85 if is_apoe else .6,zorder=3 if is_apoe else 4,alpha=1 if is_apoe else .87)
                groups.append(dict(model=model,proteins=DATA['proteins'],x=base.tolist(),y=adjusted.tolist(),jitter_applied=False))
            ax.set(xlim=limits,ylim=limits,xticks=ticks,yticks=ticks,xlabel=self.t['primary_beta'],ylabel=self.t['adjusted_beta'])
            ax.xaxis.grid(True,color=GRID,linewidth=.5);ax.yaxis.grid(True,color=GRID,linewidth=.5)
            self.style(ax)
            self.primitives[panel]=dict(assay=assay,groups=groups,identical_xy_limits=limits,identity_line=True)
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
        title='Supplementary Figure S5: Sensitivity of shared protein-p-tau181 associations to APOE and education adjustment'
        self.fig.savefig(BASE/f'Figure_S5_{self.lang}.pdf',metadata={'Title':title,'Creator':'Study supplementary figures'})
        self.fig.savefig(BASE/f'editable/Figure_S5_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_S5_{self.lang}.png',dpi=300)
        tif=BASE/f'high_resolution/Figure_S5_{self.lang}.tiff';self.fig.savefig(tif,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
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
    (AUDIT/'STYLE_SPECIFICATION.json').write_text(json.dumps(dict(palette=PALETTE,
        user_reference='User-specified cedar-green and smoky-purple palette for Figure S5.',
        model_styles={'APOE':dict(colour='#668578',marker='open circle'),'APOE and education':dict(colour='#6A627F',marker='triangle')},
        figure_size_mm=[W,H],q_colour_range=[0,20],q_colour_truncated=False,scatter_coordinates_jittered=False,
        explanatory_footnotes_outside_plot=False),ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(results,ensure_ascii=False),flush=True)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');main()
