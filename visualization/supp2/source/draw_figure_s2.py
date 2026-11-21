"""Publication-style descriptive measurement panels with the user-supplied apricot/green palette."""
from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.text import Text
from PIL import Image

BASE=Path(__file__).resolve().parent.parent
AUDIT=BASE/'audit'
DATA=json.loads((AUDIT/'PLOT_INPUT.json').read_text(encoding='utf-8'))
W,H=180,188
PALETTE=['#E7AF98','#F2D2B7','#F3E6CE','#82B6AB','#B2D0BF','#ACBFD0']
COLORS={'NULISA':'#E7AF98','SomaScan':'#82B6AB','TMT-MS':'#ACBFD0'}
BOX_COLORS={'NULISA':'#F2D2B7','SomaScan':'#B2D0BF','TMT-MS':'#ACBFD0'}
INK='#33464B';EDGE='#5B6667';GRID='#E6EBE8'
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
   'axes.linewidth':.7,'xtick.major.width':.7,'ytick.major.width':.7,'xtick.major.size':3,
   'ytick.major.size':3,'axes.unicode_minus':False,'savefig.facecolor':'white','figure.facecolor':'white'})

LABELS={
 'EN':dict(figure='Supplementary Figure S2',A='Protein measurement availability',B='NULISA detectability',
   C='Participants per technical batch',D='TMT plex-associated variation',
   availability_y='Participants with a usable\nmeasurement, %',detect_x='Participants with detectable NULISA, %',
   detect_y='Proteins, n',criterion='80% criterion',batch_y='Participants per plate or plex, n',
   batches=['NULISA\n19 plates','SomaScan\n36 plates','TMT-MS\n65 plexes'],
   before_x='Before GIS normalization, R² (%)',after_y='After GIS normalization, R² (%)',median='Median R²: {before:.1f}% → {after:.1f}%'),
 'ZH':dict(figure='补充图 S2',A='蛋白测量的可用率',B='NULISA检出率',C='每个技术批次的参与者数',D='TMT批次相关变异',
   availability_y='具有可用测量值的\n参与者比例（%）',detect_x='达到NULISA检出限的\n参与者比例（%）',
   detect_y='蛋白数',criterion='80%检出率要求',batch_y='每块板或复用批次的参与者数',
   batches=['NULISA\n19块板','SomaScan\n36块板','TMT-MS\n65个批次'],
   before_x='内标归一化前的R²（%）',after_y='内标归一化后的R²（%）',median='R²中位数：{before:.1f}% → {after:.1f}%')}

class Figure:
    def __init__(self,lang):
        self.lang=lang;self.t=LABELS[lang]
        self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=150)
        self.axes={};self.primitives={}
    def font(self,size=9,bold=False):
        name=('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN' else ('msyhbd.ttc' if bold else 'msyh.ttc')
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/name),size=size,weight='bold' if bold else 'normal')
    def figtext(self,x,y,text,size=9,bold=False,ha='left'):
        return self.fig.text(x/W,1-y/H,text,fontproperties=self.font(size,bold),color=INK,ha=ha,va='center')
    def header(self,letter,title,x,y,badge):
        self.fig.patches.extend([
            Rectangle((x/W,1-(y+3.1)/H),84/W,6.2/H,transform=self.fig.transFigure,facecolor='#F3E6CE',alpha=.42,edgecolor='none',zorder=-1),
            Rectangle((x/W,1-(y+3.1)/H),5.7/W,6.2/H,transform=self.fig.transFigure,facecolor=badge,edgecolor='none')])
        self.figtext(x+2.85,y,letter,11.5,True,'center')
        self.figtext(x+7.2,y,title,10.1,True)
    def axis(self,key,x,y,w,h):
        ax=self.fig.add_axes([x/W,1-(y+h)/H,w/W,h/H]);self.axes[key]=ax
        ax.spines[['top','right']].set_visible(False)
        for spine in ax.spines.values():spine.set_color(INK)
        ax.tick_params(colors=INK,direction='out',pad=3)
        ax.set_axisbelow(True)
        return ax
    def style(self,ax):
        for label in ax.get_xticklabels()+ax.get_yticklabels():label.set_fontproperties(self.font(9));label.set_color(INK)
        for label in (ax.xaxis.label,ax.yaxis.label):label.set_fontproperties(self.font(9.2));label.set_color(INK)
    def boxdots(self,ax,values,colors,seed):
        rng=np.random.default_rng(seed)
        stats=[];points=[]
        for i,(name,a) in enumerate(values.items()):
            a=np.asarray(a,dtype=float)
            artists=ax.boxplot([a],positions=[i],widths=.35,patch_artist=True,showfliers=False,manage_ticks=False,
               boxprops={'facecolor':BOX_COLORS[name],'edgecolor':EDGE,'linewidth':.7,'alpha':.65},
               medianprops={'color':INK,'linewidth':1.2,'zorder':5},
               whiskerprops={'color':EDGE,'linewidth':.75},capprops={'color':EDGE,'linewidth':.75})
            xx=i+rng.uniform(-.19,.19,len(a))
            ax.scatter(xx,a,s=15,color=colors[name],alpha=.8,edgecolors=INK,linewidths=.22,zorder=3)
            q1,med,q3=np.quantile(a,[.25,.5,.75]);lo=a[a>=q1-1.5*(q3-q1)].min();hi=a[a<=q3+1.5*(q3-q1)].max()
            stats.append(dict(assay=name,n=len(a),q25=float(q1),median=float(med),q75=float(q3),lower_whisker=float(lo),upper_whisker=float(hi)))
            points.append(dict(assay=name,x=xx.tolist(),y=a.tolist()))
        return stats,points
    def draw(self):
        t=self.t
        self.figtext(2,6.5,t['figure'],13,True)
        for letter,x,y,badge in [('A',2,18,'#E7AF98'),('B',94,18,'#F2D2B7'),('C',2,105,'#82B6AB'),('D',94,105,'#ACBFD0')]:
            self.header(letter,t[letter],x,y,badge)
        # A: the same 43 candidate proteins, with all percentages using 571 participants.
        ax=self.axis('A',20,29,62,54)
        stats,points=self.boxdots(ax,DATA['availability'],COLORS,20261031)
        self.primitives['A']=dict(boxes=stats,points=points)
        ax.set(xlim=(-.5,2.5),ylim=(0,108),xticks=[0,1,2],xticklabels=list(COLORS),yticks=np.arange(0,101,20),ylabel=t['availability_y'])
        ax.yaxis.grid(True,color=GRID,linewidth=.5)
        # B: NPQ at or above the target/plate limit of detection, not mere numerical availability.
        ax=self.axis('B',112,29,62,54)
        counts=np.asarray(DATA['histogram_counts']);edges=np.asarray(DATA['histogram_edges'])
        ax.bar(edges[:-1]+5,counts,width=9.2,color='#F2D2B7',edgecolor='#A4745F',linewidth=.65,zorder=2)
        ax.axvline(80,color='#527F70',linestyle=(0,(4,3)),linewidth=1.0,zorder=3)
        for center,n in zip(edges[:-1]+5,counts):
            if n:ax.text(center,n+.65,str(int(n)),ha='center',va='bottom',fontproperties=self.font(8.8),color=INK)
        ax.set(xlim=(0,100),ylim=(0,40),xticks=np.arange(0,101,20),yticks=[0,10,20,30,40],xlabel=t['detect_x'],ylabel=t['detect_y'])
        ax.yaxis.grid(True,color=GRID,linewidth=.5)
        handle=Line2D([0],[0],color='#527F70',linewidth=1.0,linestyle=(0,(4,3)))
        legend=ax.legend([handle],[t['criterion']],loc='upper left',frameon=False,prop=self.font(8.8),handlelength=1.8,borderaxespad=.3)
        for x in legend.get_texts():x.set_color(INK)
        self.primitives['B']=dict(bin_edges=edges.tolist(),counts=counts.tolist(),criterion=80)
        # C: counts of paired participants, one point per plate or plex.
        ax=self.axis('C',20,117,62,53)
        stats,points=self.boxdots(ax,DATA['batch_sizes'],COLORS,20261101)
        self.primitives['C']=dict(boxes=stats,points=points)
        ax.set(xlim=(-.5,2.5),ylim=(0,60),xticks=[0,1,2],xticklabels=t['batches'],yticks=[0,10,20,30,40,50,60],ylabel=t['batch_y'])
        ax.yaxis.grid(True,color=GRID,linewidth=.5)
        # D: every before/after pair uses identical GIS-available observations.
        ax=self.axis('D',112,117,62,53)
        before=np.asarray(DATA['tmt_r2_before']);after=np.asarray(DATA['tmt_r2_after'])
        ax.plot([0,100],[0,100],linestyle=(0,(4,3)),color='#7C8588',linewidth=.8,zorder=1)
        ax.scatter(before,after,s=23,color='#ACBFD0',edgecolors='#526978',linewidths=.55,alpha=.92,zorder=3)
        summary=t['median'].format(before=DATA['tmt_r2_before_median'],after=DATA['tmt_r2_after_median'])
        ax.text(.035,.965,summary,transform=ax.transAxes,ha='left',va='top',fontproperties=self.font(8.8,True),color=INK)
        j=DATA['proteins'].index('S100A12')
        ax.annotate('S100A12',xy=(before[j],after[j]),xytext=(34,66),textcoords='data',
          fontproperties=self.font(9),color=INK,ha='center',va='center',
          arrowprops={'arrowstyle':'-','color':INK,'linewidth':.6,'shrinkA':2,'shrinkB':4})
        ax.set(xlim=(0,100),ylim=(0,100),xticks=np.arange(0,101,20),yticks=np.arange(0,101,20),xlabel=t['before_x'],ylabel=t['after_y'])
        ax.xaxis.grid(True,color=GRID,linewidth=.45);ax.yaxis.grid(True,color=GRID,linewidth=.45)
        self.primitives['D']=dict(proteins=DATA['proteins'],x=before.tolist(),y=after.tolist(),identity_line=[0,100],labelled_protein='S100A12',median_text=summary)
        for ax in self.axes.values():self.style(ax)
        self.fig.canvas.draw()
        return self
    def check_layout(self):
        self.fig.canvas.draw();renderer=self.fig.canvas.get_renderer();records=[]
        all_text=[]
        for artist in self.fig.findobj(Text):
            if artist.get_visible() and artist.get_text().strip() and id(artist) not in [id(x) for x in all_text]:all_text.append(artist)
        for i,artist in enumerate(all_text):
            box=artist.get_window_extent(renderer)
            mm=np.array([box.x0,box.y0,box.x1,box.y1])/self.fig.dpi*25.4
            check=mm[0]>=-.15 and mm[2]<=W+.15 and mm[1]>=-.15 and mm[3]<=H+.15
            assert check,(self.lang,'text outside canvas',artist.get_text(),mm.tolist())
            records.append(dict(text=artist.get_text(),point_size=float(artist.get_fontsize()),
               x0_mm=float(mm[0]),y0_mm=float(H-mm[3]),x1_mm=float(mm[2]),y1_mm=float(H-mm[1])))
        overlap=[]
        for i,a in enumerate(records):
            for b in records[i+1:]:
                xx=min(a['x1_mm'],b['x1_mm'])-max(a['x0_mm'],b['x0_mm'])
                yy=min(a['y1_mm'],b['y1_mm'])-max(a['y0_mm'],b['y0_mm'])
                if xx>.35 and yy>.35:overlap.append((a['text'],b['text'],round(xx,2),round(yy,2)))
        assert not overlap,(self.lang,'text overlap',overlap)
        return records
    def save(self):
        records=self.check_layout()
        for dirname in ('editable','high_resolution'): (BASE/dirname).mkdir(exist_ok=True)
        title='Supplementary Figure S2: Coverage and technical characteristics of CSF protein measurements across three assays'
        self.fig.savefig(BASE/f'Figure_S2_{self.lang}.pdf',metadata={'Title':title,'Creator':'Study supplementary figures'})
        self.fig.savefig(BASE/f'editable/Figure_S2_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_S2_{self.lang}.png',dpi=300)
        tiff=BASE/f'high_resolution/Figure_S2_{self.lang}.tiff'
        self.fig.savefig(tiff,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        with Image.open(tiff) as im:
            rgb=im.convert('RGB');rgb.save(tiff,compression='tiff_lzw',dpi=(600,600));rgb.close()
        (AUDIT/f'LAYOUT_{self.lang}.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
        (AUDIT/f'PLOT_PRIMITIVES_{self.lang}.json').write_text(json.dumps(self.primitives,ensure_ascii=False,indent=2),encoding='utf-8')
        return dict(language=self.lang,width_mm=W,height_mm=H,minimum_font_pt=min(r['point_size'] for r in records),
                    checked_text_elements=len(records),text_clipped=False,text_overlap=False)

def main():
    results=[]
    for lang in ('EN','ZH'):
        fig=Figure(lang).draw();results.append(fig.save());plt.close(fig.fig)
    (AUDIT/'LAYOUT_VALIDATION.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
    (AUDIT/'STYLE_SPECIFICATION.json').write_text(json.dumps(dict(palette=PALETTE,assay_colours=COLORS,
       user_reference='Apricot, cream, mint green and gray blue; current Figure S2 instruction overrides the Figure S1 palette.',
       explanatory_footnotes_outside_plot=False),ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(results,ensure_ascii=False),flush=True)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');main()
