"""Draw bilingual vector figures using aggregate longitudinal follow-up summaries."""
from pathlib import Path
import json, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from PIL import Image

BASE=Path(__file__).resolve().parent.parent
AUDIT=BASE/'audit'
DATA=json.loads((AUDIT/'PLOT_INPUT.json').read_text(encoding='utf-8'))
W,H=180,258
PALETTE=['#A8D6C5','#BDDDEA','#C9DEB5','#F5E3A5','#F3C8AA','#E6B5C6']
COLORS={'Memory':'#A8D6C5','Executive function':'#BDDDEA'}
EDGES={'Memory':'#487969','Executive function':'#4D7590'}
MARKERS={'Memory':'o','Executive function':'D'}
INK='#263C59';GRID='#E8EEEF'
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
    'axes.linewidth':.7,'xtick.major.width':.7,'ytick.major.width':.7,'xtick.major.size':3,
    'ytick.major.size':0,'axes.unicode_minus':False,'savefig.facecolor':'white','figure.facecolor':'white'})

TEXT={
 'EN':dict(figure='Supplementary Figure S3',A='Analysis sample size',B='Observed follow-up span',
     C='Repeated cognitive assessments',D='Longer follow-up coverage',
     people_x='Participants, n',span_x='Years from baseline\nto last assessment',assess_y='Cognitive assessments, n',
     coverage_y='Participants, %',domain_labels=['Memory','Executive function'],
     coverage_labels=['≥2 follow-ups\nover ≥2 years','Last assessment\n≥4.5 years after CSF']),
 'ZH':dict(figure='补充图 S3',A='逐蛋白分析人数',B='实际随访跨度',
     C='重复认知评估的总次数',D='较长随访的覆盖比例',
     people_x='参与者人数',span_x='基线至末次评估的时间（年）',assess_y='认知评估总次数',
     coverage_y='参与者比例（%）',domain_labels=['记忆','执行功能'],
     coverage_labels=['≥2次随访\n且跨度≥2年','末次评估在\nCSF后≥4.5年'])}

class Figure:
    def __init__(self,lang):
        self.lang=lang;self.t=TEXT[lang]
        self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=160)
        self.axes={};self.primitives={}
    def font(self,size=9,bold=False):
        filename=('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN' else ('msyhbd.ttc' if bold else 'msyh.ttc')
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/filename),size=size,weight='bold' if bold else 'normal')
    def text(self,x,y,value,size=9,bold=False,ha='left'):
        return self.fig.text(x/W,1-y/H,value,fontproperties=self.font(size,bold),color=INK,ha=ha,va='center')
    def header(self,key,x,y,badge,fill):
        self.fig.patches.extend([
            Rectangle((x/W,1-(y+3.1)/H),84/W,6.2/H,transform=self.fig.transFigure,facecolor=fill,alpha=.45,edgecolor='none',zorder=-1),
            Rectangle((x/W,1-(y+3.1)/H),5.7/W,6.2/H,transform=self.fig.transFigure,facecolor=badge,edgecolor='none')])
        self.text(x+2.85,y,key,11.5,True,'center');self.text(x+7.2,y,self.t[key],9.8,True)
    def axis(self,key,x,y,w,h):
        ax=self.fig.add_axes([x/W,1-(y+h)/H,w/W,h/H]);self.axes[key]=ax
        ax.spines[['top','right']].set_visible(False)
        for s in ax.spines.values():s.set_color(INK)
        ax.tick_params(colors=INK,direction='out',pad=3)
        ax.set_axisbelow(True)
        return ax
    def domain_legend(self,x,y):
        handles=[Line2D([],[],linestyle='none',marker=MARKERS[k],markersize=4.6,
            markerfacecolor=COLORS[k],markeredgecolor=EDGES[k],markeredgewidth=.6) for k in COLORS]
        legend=self.fig.legend(handles,self.t['domain_labels'],loc='center left',bbox_to_anchor=(x/W,1-y/H),
            ncol=2,frameon=False,prop=self.font(8.8),handlelength=.8,handletextpad=.35,columnspacing=1.05,borderaxespad=0)
        for item in legend.get_texts():item.set_color(INK)
    def style(self,ax):
        for item in ax.get_xticklabels()+ax.get_yticklabels():item.set_fontproperties(self.font(9));item.set_color(INK)
        for item in [ax.xaxis.label,ax.yaxis.label]:item.set_fontproperties(self.font(9.2));item.set_color(INK)
    def box(self,ax,values,pos,domain,width=.32,jitter=.085,seed=1):
        a=np.asarray(values,dtype=float)
        ax.boxplot([a],positions=[pos],widths=width,patch_artist=True,showfliers=False,manage_ticks=False,
            boxprops=dict(facecolor=COLORS[domain],edgecolor=EDGES[domain],linewidth=.8,alpha=.45),
            medianprops=dict(color=INK,linewidth=1.15,zorder=5),
            whiskerprops=dict(color=EDGES[domain],linewidth=.75),capprops=dict(color=EDGES[domain],linewidth=.75))
        rng=np.random.default_rng(seed);xx=pos+rng.uniform(-jitter,jitter,len(a))
        ax.scatter(xx,a,s=13,marker=MARKERS[domain],color=COLORS[domain],edgecolors=EDGES[domain],linewidths=.35,alpha=.83,zorder=3)
        q1,med,q3=np.quantile(a,[.25,.5,.75]);lo=a[a>=q1-1.5*(q3-q1)].min();hi=a[a<=q3+1.5*(q3-q1)].max()
        return dict(domain=domain,x=xx.tolist(),y=a.tolist(),position=pos,q25=float(q1),median=float(med),q75=float(q3),lower_whisker=float(lo),upper_whisker=float(hi))
    def draw(self):
        self.text(2,6.5,self.t['figure'],13,True)
        self.header('A',2,18,'#A8D6C5','#C9DEB5')
        self.header('B',94,18,'#BDDDEA','#BDDDEA')
        self.header('C',2,186,'#C9DEB5','#F5E3A5')
        self.header('D',94,186,'#F3C8AA','#E6B5C6')
        self.domain_legend(112,26.2);self.domain_legend(112,194.4)
        genes=DATA['proteins'];yy=np.arange(len(genes));n=np.asarray(DATA['domains']['Memory']['participants'])
        ax=self.axis('A',19,33,63,130)
        ax.barh(yy,n,height=.57,facecolor='#A8D6C5',edgecolor='#487969',linewidth=.4,alpha=.78)
        for i,value in enumerate(n):
            ax.text(value+10,i,str(value),fontproperties=self.font(8.8),color=INK,ha='left',va='center')
        ax.set(xlim=(0,650),ylim=(33.65,-.65),xticks=[0,200,400,600],yticks=yy,yticklabels=genes,xlabel=self.t['people_x'])
        ax.xaxis.grid(True,color=GRID,linewidth=.5)
        self.primitives['A']=dict(proteins=genes,participants=n.tolist(),identical_between_outcomes=True)
        ax=self.axis('B',111,33,63,130)
        spans=[]
        for j,domain in enumerate(COLORS):
            values=DATA['domains'][domain];med=np.asarray(values['span_median']);q1=np.asarray(values['span_q25']);q3=np.asarray(values['span_q75'])
            yp=yy+(-.145 if j==0 else .145)
            ax.errorbar(med,yp,xerr=np.vstack([med-q1,q3-med]),fmt=MARKERS[domain],markersize=3.6,
                color=EDGES[domain],mfc=COLORS[domain],mec=EDGES[domain],mew=.45,elinewidth=.7,capsize=1.7,capthick=.6,zorder=3)
            spans.append(dict(domain=domain,proteins=genes,median=med.tolist(),q25=q1.tolist(),q75=q3.tolist(),y=yp.tolist()))
        ax.set(xlim=(0,5.3),ylim=(33.65,-.65),xticks=[0,1,2,3,4,5],yticks=yy,yticklabels=[],xlabel=self.t['span_x'])
        ax.xaxis.grid(True,color=GRID,linewidth=.5)
        for i in range(len(genes)):
            if i%2==0:
                for aa in (self.axes['A'],ax):aa.axhspan(i-.48,i+.48,facecolor='#BDDDEA',alpha=.10,zorder=0)
        self.primitives['B']=spans
        ax=self.axis('C',20,205,62,35)
        assays=[]
        for j,domain in enumerate(COLORS):
            assays.append(self.box(ax,DATA['domains'][domain]['assessments'],j,domain,width=.35,jitter=.13,seed=20261101+j))
        ax.set(xlim=(-.5,1.5),ylim=(0,3150),xticks=[0,1],xticklabels=self.t['domain_labels'],yticks=[0,1000,2000,3000],ylabel=self.t['assess_y'])
        ax.yaxis.grid(True,color=GRID,linewidth=.5)
        self.primitives['C']=assays
        ax=self.axis('D',112,205,62,35)
        coverage=[]
        for k,metric in enumerate(['repeated_followup_pct','late_followup_pct']):
            for j,domain in enumerate(COLORS):
                z=self.box(ax,DATA['domains'][domain][metric],k+(-.17 if j==0 else .17),domain,width=.23,jitter=.075,seed=20261104+2*k+j)
                z['metric']=metric;coverage.append(z)
        ax.set(xlim=(-.53,1.53),ylim=(0,100),xticks=[0,1],xticklabels=self.t['coverage_labels'],yticks=[0,25,50,75,100],ylabel=self.t['coverage_y'])
        ax.yaxis.grid(True,color=GRID,linewidth=.5)
        self.primitives['D']=coverage
        for aa in self.axes.values():self.style(aa)
        self.fig.canvas.draw();return self
    def layout(self):
        self.fig.canvas.draw();renderer=self.fig.canvas.get_renderer();records=[]
        artists=[]
        for artist in self.fig.findobj(Text):
            if artist.get_visible() and artist.get_text().strip() and id(artist) not in [id(x) for x in artists]:artists.append(artist)
        for artist in artists:
            box=artist.get_window_extent(renderer);mm=np.array([box.x0,box.y0,box.x1,box.y1])/self.fig.dpi*25.4
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
        title='Supplementary Figure S3: Cognitive follow-up in the analysis samples for individual proteins'
        self.fig.savefig(BASE/f'Figure_S3_{self.lang}.pdf',metadata={'Title':title,'Creator':'Study supplementary figures'})
        self.fig.savefig(BASE/f'editable/Figure_S3_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_S3_{self.lang}.png',dpi=300)
        tif=BASE/f'high_resolution/Figure_S3_{self.lang}.tiff'
        self.fig.savefig(tif,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        with Image.open(tif) as im:
            rgb=im.convert('RGB');rgb.save(tif,compression='tiff_lzw',dpi=(600,600));rgb.close()
        (AUDIT/f'LAYOUT_{self.lang}.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
        (AUDIT/f'PLOT_PRIMITIVES_{self.lang}.json').write_text(json.dumps(self.primitives,ensure_ascii=False,indent=2),encoding='utf-8')
        return dict(language=self.lang,width_mm=W,height_mm=H,minimum_font_pt=min(r['point_size'] for r in records),checked_text_elements=len(records),text_clipped=False,text_overlap=False)

def main():
    results=[]
    for lang in ['EN','ZH']:
        figure=Figure(lang).draw();results.append(figure.save());plt.close(figure.fig)
    (AUDIT/'LAYOUT_VALIDATION.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
    (AUDIT/'STYLE_SPECIFICATION.json').write_text(json.dumps(dict(palette=PALETTE,domain_colours=COLORS,domain_shapes=MARKERS,
        user_reference='Fresh mint, pale blue, sage, pale yellow, peach and rose; specified for Figure S3.',
        explanatory_footnotes_outside_plot=False,figure_size_mm=[W,H]),ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(results,ensure_ascii=False),flush=True)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');main()
