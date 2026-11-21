"""Draw paired deletion estimates with full-sample anchors and equal-axis identity references."""
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
W,H=180,220
INK='#263A46';GRID='#E9EEEB';IDENTITY='#A6B3C0'
COLORS={'NULISA':'#FAD6C9','SomaScan':'#B0C6A8','TMT-MS':'#92B3CD'}
EDGES={'NULISA':'#B5806D','SomaScan':'#6F8C64','TMT-MS':'#547F9E'}
MARKERS={'NULISA':'o','SomaScan':'s','TMT-MS':'D'}
PALETTE=['#FAD6C9','#B0C6A8','#92B3CD','#FAF5E8','#D0E3DB','#E8E2DB','#A6B3C0']
TEXT={
 'EN':dict(figure='Supplementary Figure S14',legendhead='Omitted\nbatch',full='Full sample',identity='No change with adjustment',
   rowtitles=['SomaScan − NULISA memory association contrast','TMT-MS − NULISA memory association contrast'],
   x='Without p-tau181\n(ADNI-MEM/year per protein residual SD)',y='With log₂ p-tau181'),
 'ZH':dict(figure='补充图 S14',legendhead='删除的\n批次',full='完整样本',identity='调整后对比不变',
   rowtitles=['SomaScan减NULISA的记忆关联对比','TMT-MS减NULISA的记忆关联对比'],
   x='不调整p-tau181\n（ADNI-MEM分/年，每1蛋白残差SD）',y='调整对数p-tau181')
}
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.linewidth':.7,
    'xtick.major.width':.7,'ytick.major.width':.7,'xtick.major.size':3,'ytick.major.size':3,
    'axes.unicode_minus':False,'savefig.facecolor':'white','figure.facecolor':'white'})

def ticklabel(value):
    return '0' if value==0 else f'{value:.2f}'.replace('-','−')

class Figure:
    def __init__(self,lang):
        self.lang=lang;self.t=TEXT[lang]
        self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=160)
        self.axes={};self.primitives={}
    def font(self,size=9,bold=False,math=False):
        name='segoeui.ttf' if math else (('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN' else ('msyhbd.ttc' if bold else 'msyh.ttc'))
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/name),size=size,weight='bold' if bold else 'normal')
    def text(self,x,y,label,size=9,bold=False,ha='left',math=False,rotation=0):
        return self.fig.text(x/W,1-y/H,label,fontproperties=self.font(size,bold,math),color=INK,ha=ha,va='center',
                             linespacing=1.06,rotation=rotation,rotation_mode='anchor')
    def box(self,x,y,width,height,fill):
        self.fig.patches.append(Rectangle((x/W,1-(y+height/2)/H),width/W,height/H,transform=self.fig.transFigure,
            facecolor=fill,edgecolor='none',zorder=-1))
    def style(self,ax):
        for label in ax.get_xticklabels()+ax.get_yticklabels():
            label.set_fontproperties(self.font(8.8));label.set_color(INK)
    def draw(self):
        self.text(2,5.7,self.t['figure'],13,True)
        self.text(2,15.5,self.t['legendhead'],8.8,True)
        for assay,x in [('NULISA',29),('SomaScan',77),('TMT-MS',130)]:
            self.fig.lines.append(Line2D([x/W],[1-15.5/H],transform=self.fig.transFigure,linestyle='none',
                marker=MARKERS[assay],markersize=3.8,markerfacecolor=COLORS[assay],
                markeredgecolor=EDGES[assay],markeredgewidth=.6))
            self.text(x+3,15.5,f'{assay} ({DATA["batch_counts"][assay]})',8.8)
        self.fig.lines.append(Line2D([46/W],[1-24.5/H],transform=self.fig.transFigure,linestyle='none',
            marker='*',markersize=7,markerfacecolor=INK,markeredgecolor='white',markeredgewidth=.45))
        self.text(49,24.5,self.t['full'],8.8)
        self.fig.lines.append(Line2D([102/W,109/W],[1-24.5/H]*2,transform=self.fig.transFigure,
            color=IDENTITY,linewidth=1,linestyle=(0,(4,3))))
        self.text(112,24.5,self.t['identity'],8.8)
        for row_index,contrast in enumerate(DATA['contrasts']):
            band_y=36+96*row_index;head_y=47+96*row_index;top=56+96*row_index
            xlabel_y=110+96*row_index
            fill='#D0E3DB' if row_index==0 else '#E8E2DB'
            self.box(2,band_y,176,6.4,fill)
            self.text(3.5,band_y,self.t['rowtitles'][row_index],9.5,True)
            self.text(3.4,top+21,self.t['y'],9,ha='center',math=self.lang=='EN',rotation=90)
            self.text(99,xlabel_y,self.t['x'],9,ha='center')
            scale=DATA['scales'][contrast];lo,hi=scale['limits'];ticks=scale['ticks']
            for j,(gene,x) in enumerate(zip(DATA['proteins'],[20,77,134])):
                panel='ABCDEF'[row_index*3+j];candidate=DATA['candidates'][gene]
                self.box(x,head_y,42,6.4,'#FAF5E8')
                self.box(x,head_y,5.5,6.4,COLORS[DATA['omitted_assays'][j]])
                self.text(x+2.75,head_y,panel,10.5,True,'center')
                self.text(x+6.8,head_y,f'{gene} (n={candidate["full_participants"]})',8.8,True)
                ax=self.fig.add_axes([x/W,1-(top+42)/H,42/W,42/H])
                self.axes[panel]=ax;ax.spines[['top','right']].set_visible(False)
                ax.spines['left'].set_color(INK);ax.spines['bottom'].set_color(INK)
                ax.tick_params(colors=INK,direction='out',pad=2)
                ax.set_axisbelow(True);ax.grid(color=GRID,linewidth=.5)
                ax.plot([lo,hi],[lo,hi],color=IDENTITY,linestyle=(0,(4,3)),linewidth=.9,zorder=1)
                groups=[]
                # Draw the numerous TMT omissions first, leaving the other groups visible at overlaps.
                for assay in reversed(DATA['omitted_assays']):
                    group=candidate['omissions'][assay];values=group['contrasts'][contrast]
                    ax.plot(values['before'],values['after'],linestyle='none',marker=MARKERS[assay],markersize=3.3,
                        mfc=COLORS[assay],mec=EDGES[assay],mew=.5,alpha=.68,zorder=3)
                for assay in DATA['omitted_assays']:
                    group=candidate['omissions'][assay];values=group['contrasts'][contrast]
                    groups.append(dict(omitted_assay=assay,batch_labels=group['batch_labels'],
                        participants_retained=group['participants_retained'],
                        before=values['before'],after=values['after'],change=values['change']))
                ref=candidate['reference'][contrast]
                ax.plot(ref['before'],ref['after'],linestyle='none',marker='*',markersize=7,
                    mfc=INK,mec='white',mew=.55,zorder=7)
                ax.set(xlim=(lo,hi),ylim=(lo,hi),xticks=ticks,yticks=ticks)
                ax.set_xticklabels([ticklabel(v) for v in ticks]);ax.set_yticklabels([ticklabel(v) for v in ticks])
                ax.set_aspect('equal',adjustable='box');self.style(ax)
                self.primitives[panel]=dict(protein=gene,contrast=contrast,full_sample_participants=candidate['full_participants'],
                    groups=groups,reference=ref,axis_limits=[lo,hi],same_x_y_scale=True,equal_aspect=True,
                    batch_points=sum(len(g['before']) for g in groups),
                    coordinate_jitter=False,confidence_intervals_drawn=False,new_tests=False)
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
        title='Supplementary Figure S14: Sensitivity of candidate memory associations to individual technical batches'
        self.fig.savefig(BASE/f'Figure_S14_{self.lang}.pdf',metadata={'Title':title,'Creator':'Study supplementary figures'})
        self.fig.savefig(BASE/f'editable/Figure_S14_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_S14_{self.lang}.png',dpi=300)
        tif=BASE/f'high_resolution/Figure_S14_{self.lang}.tiff'
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
        fig=Figure(lang).draw();results.append(fig.save());plt.close(fig.fig)
    (AUDIT/'LAYOUT_VALIDATION.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
    (AUDIT/'STYLE_SPECIFICATION.json').write_text(json.dumps(dict(
        palette=PALETTE,omitted_assay_colours=COLORS,omitted_assay_shapes=MARKERS,
        choice='Continue apricot, sage and gray-blue for the assay whose batch is omitted; charcoal stars mark full-sample references.',
        figure_size_mm=[W,H],panels=6,paired_points=720,full_sample_stars=6,
        row_scales=DATA['scales'],equal_x_y_scales_within_each_row=True,no_jitter=True,no_intervals=True,no_new_P_or_q=True,
        loose_footer_notes=False),ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(results,ensure_ascii=False),flush=True)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');main()
