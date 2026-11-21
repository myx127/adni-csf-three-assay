"""Draw nominal calibration, complete tail-assessed designs, the CST3 gate and unresolved resolution."""
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
W,H=180,258
INK='#263A46';GRID='#E9EEEB';BLUE='#527F94';GREEN='#6F8C64';ROSE='#A36872'
PALETTE=['#92B3CD','#B0C6A8','#D99DA3','#FAF5E8','#D0E3DB','#E8E2DB','#527F94']
ROWS=np.arange(8)*12+47.
TEXT={
 'EN':dict(figure='Supplementary Figure S7',analysis='Analysis and\ncalibration',
   A='False-positive rate\nat α=0.05',B='BY-threshold\nassessment',
   C='CST3: calibration of\nmemory association tests',D='Monte Carlo resolution\nfor relative p-tau181 responses',
   rate='Rate',upper='Upper bound',passed='Qualified',failed='Not qualified',
   labels=['Shared p-tau181: age/sex\n35/35 analyses','+ APOE ε4\n35/35 analyses',
           '+ Education/APOE ε4\n35/35 analyses','Memory change\n33/34 analyses',
           'Secondary cognition\n68/68 analyses','14-protein sensitivities\n42/42 analyses',
           'p-tau181 adjustment\n42/42 comparisons','Relative p-tau181 responses\n35/35 nominal only'],
   unitsA='False-positive probability',unitsB='Largest Wilson upper bound\n/ BY first-step threshold',
   unresolved='Tail not resolved',
   shared='Shared direction',diff='Assay difference',
   gate_rows=['Nominal\n(limit 0.06)','BY: 68 tests','BY: 136 tests','BY: 204 tests'],
   unitsC='Upper bound / required limit',
   resolution_rows=['BY first-step threshold\n(34 tests)','Wilson upper bound\n(0/4,000 rejections)'],
   unitsD='Rejection probability'),
 'ZH':dict(figure='补充图 S7',analysis='分析及\n校准情况',
   A='显著性水平α=0.05\n下的假阳性率',B='BY阈值\n校准评价',
   C='CST3记忆关联检验\n的阈值校准表现',D='相对p-tau181响应检验\n的模拟分辨率',
   rate='拒绝率',upper='上界',passed='合格',failed='未合格',
   labels=['共同p-tau181：年龄/性别\n35/35项分析','增加APOE ε4\n35/35项分析',
           '增加教育/APOE ε4\n35/35项分析','记忆年变化\n33/34项分析',
           '次要认知结局\n68/68项分析','14蛋白敏感性分析\n42/42项分析',
           'p-tau181调整比较\n42/42比较','相对p-tau181响应\n35/35仅名义检查'],
   unitsA='假阳性概率',unitsB='最大Wilson上界\n／BY首步阈值',
   unresolved='尾部尚未解析',
   shared='共同方向',diff='检测间差异',
   gate_rows=['名义水平\n（上界限0.06）','BY：68项','BY：136项','BY：204项'],
   unitsC='上界／要求的限值',
   resolution_rows=['BY首步阈值\n（34项）','最低Wilson上界\n（0/4,000次拒绝）'],
   unitsD='拒绝概率')
}
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.linewidth':.7,
    'xtick.major.width':.7,'ytick.major.width':.7,'xtick.major.size':3,'ytick.major.size':0,
    'axes.unicode_minus':False,'savefig.facecolor':'white','figure.facecolor':'white'})

class Figure:
    def __init__(self,lang):
        self.lang=lang;self.t=TEXT[lang]
        self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=160);self.axes={};self.primitives={}
    def font(self,size=9,bold=False,math=False):
        name='segoeui.ttf' if math else (('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN' else ('msyhbd.ttc' if bold else 'msyh.ttc'))
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/name),size=size,weight='bold' if bold else 'normal')
    def text(self,x,y,label,size=9,bold=False,ha='left',math=False):
        return self.fig.text(x/W,1-y/H,label,fontproperties=self.font(size,bold,math),color=INK,ha=ha,va='center',linespacing=1.06)
    def box(self,x,y,width,height,fill):
        self.fig.patches.append(Rectangle((x/W,1-(y+height/2)/H),width/W,height/H,
            transform=self.fig.transFigure,facecolor=fill,edgecolor='none',zorder=-1))
    def header(self,key,x,y,w,fill,badge):
        self.box(x,y,w,10,fill);self.box(x,y,5.5,10,badge)
        self.text(x+2.75,y,key,11.2,True,'center');self.text(x+6.8,y,self.t[key],9.2,True)
    def axis(self,key,x,top,w,h):
        ax=self.fig.add_axes([x/W,1-(top+h)/H,w/W,h/H],facecolor='none');self.axes[key]=ax
        ax.spines[['top','right','left']].set_visible(False);ax.spines['bottom'].set_color(INK)
        ax.tick_params(colors=INK,direction='out',pad=3);ax.set_axisbelow(True)
        return ax
    def style(self,ax):
        for label in ax.get_xticklabels()+ax.get_yticklabels()+[ax.xaxis.label,ax.yaxis.label]:
            label.set_fontproperties(self.font(8.8));label.set_color(INK)
    def marker_legend(self,x,y,label,marker,face,edge):
        self.fig.lines.append(Line2D([x/W],[1-y/H],transform=self.fig.transFigure,linestyle='none',
            marker=marker,markersize=3.8,markerfacecolor=face,markeredgecolor=edge,markeredgewidth=.6))
        self.text(x+3,y,label,8.8)
    def draw(self):
        self.text(2,5.7,self.t['figure'],13,True)
        self.text(2,27,self.t['analysis'],9.2,True)
        self.header('A',58,27,53,'#FAF5E8','#92B3CD')
        self.header('B',121,27,57,'#D0E3DB','#B0C6A8')
        self.marker_legend(61,37,self.t['rate'],'o','#92B3CD',BLUE)
        self.fig.lines.extend([
            Line2D([86/W,90/W],[1-37/H]*2,transform=self.fig.transFigure,color=ROSE,linewidth=.7),
            Line2D([90/W]*2,[1-35.8/H,1-38.2/H],transform=self.fig.transFigure,color=ROSE,linewidth=.8)])
        self.text(93,37,self.t['upper'],8.8)
        self.marker_legend(124,37,self.t['passed'],'o','#B0C6A8',GREEN)
        self.marker_legend(151,37,self.t['failed'],'D','#D99DA3',ROSE)
        for i,(g,y) in enumerate(zip(DATA['groups'],ROWS)):
            if i%2==0:self.box(2,y,176,10,'#E8E2DB')
            self.text(2,y,self.t['labels'][i],8.8)
        a=self.axis('A',58,41,53,96)
        rates=np.array([g['rate'] for g in DATA['groups']]);upper=np.array([g['upper'] for g in DATA['groups']])
        a.hlines(ROWS,rates,upper,color=BLUE,linewidth=.85,zorder=2)
        a.plot(rates,ROWS,'o',ms=3.9,mfc='#92B3CD',mec=BLUE,mew=.6,zorder=3)
        a.plot(upper,ROWS,linestyle='none',marker='|',ms=6,color=ROSE,mew=.85,zorder=4)
        a.set(xlim=(0,.0125),ylim=(137,41),xticks=[0,.005,.01],yticks=[],xlabel=self.t['unitsA'])
        a.set_xticklabels(['0','0.005','0.010']);a.xaxis.grid(True,color=GRID,linewidth=.5);self.style(a)
        self.primitives['A']=dict(group_keys=[g['key'] for g in DATA['groups']],rates=rates.tolist(),upper=upper.tolist(),
            nominal_level=.05,nominal_upper_limit=.06,x_limits=[0,.0125],bounds_cellwise_not_simultaneous=True)
        b=self.axis('B',122,41,56,96)
        b.vlines(1,41,124,color=ROSE,linestyle=(0,(4,3)),lw=.8,zorder=1)
        bpoints=[];failed_label=None
        for g,y in zip(DATA['groups'],ROWS):
            if not g['tail_assessed']:
                self.text(150.5,y,self.t['unresolved'],8.8,ha='center');continue
            offsets=np.linspace(-2.2,2.2,len(g['design_records']))
            for row,off in zip(g['design_records'],offsets):
                ratio=row['tail_ratio'];qualified=row['qualified'];yp=float(y+off)
                b.plot(ratio,yp,linestyle='none',marker='o' if qualified else 'D',ms=3.1 if qualified else 4.5,
                    mfc='#B0C6A8' if qualified else '#D99DA3',mec=GREEN if qualified else ROSE,mew=.5,alpha=.65 if qualified else 1,zorder=4 if qualified else 6)
                bpoints.append(dict(group=g['key'],protein=row['protein'],model=row['model'],ratio=ratio,qualified=qualified,y=yp))
                if not qualified:failed_label=(ratio,yp)
        b.set(xlim=(0,1.55),ylim=(137,41),xticks=[0,.5,1,1.5],yticks=[],xlabel=self.t['unitsB'])
        b.vlines([0,.5,1.5],41,124,color=GRID,linewidth=.5,zorder=0);self.style(b)
        if failed_label:
            b.annotate('CST3',xy=failed_label,xytext=(failed_label[0]-.05,failed_label[1]-5.5),ha='right',
                fontproperties=self.font(8.8),color=ROSE,arrowprops=dict(arrowstyle='-',lw=.7,color=ROSE),zorder=8)
        self.primitives['B']=dict(points=bpoints,ratio_limit=1,x_limits=[0,1.55],unresolved_group='relative',
            unresolved_designs=35,only_categorical_vertical_spread=True,horizontal_values_unchanged=True)
        self.header('C',2,167,84,'#E8E2DB','#D99DA3')
        self.header('D',94,167,84,'#FAF5E8','#B0C6A8')
        self.marker_legend(3.5,181.5,self.t['shared'],'o','#92B3CD',BLUE)
        self.marker_legend(40,181.5,self.t['diff'],'s','#B0C6A8',GREEN)
        cy=np.array([198,210,222,234.])
        c=self.axis('C',32,191,52,50);c.axvline(1,color=ROSE,linestyle=(0,(4,3)),lw=.8,zorder=1)
        for y,label in zip(cy,self.t['gate_rows']):self.text(2,y,label,8.8)
        cpoints=[]
        thresholds=['NOMINAL','BY68_SUFFICIENT','BY136_SUFFICIENT','BY204_SUFFICIENT']
        for j,(test,marker,face,edge) in enumerate([('Shared direction','o','#92B3CD',BLUE),('Assay differences','s','#B0C6A8',GREEN)]):
            for threshold,y in zip(thresholds,cy):
                row=next(r for r in DATA['cst3'] if r['test']==test and r['threshold']==threshold)
                yp=float(y+(j-.5)*3)
                c.plot(row['ratio'],yp,linestyle='none',marker=marker,ms=4.2,
                    mfc=face if row['passed'] else '#D99DA3',mec=edge if row['passed'] else ROSE,mew=.7,zorder=4)
                cpoints.append(dict(**row,y=yp))
        c.set(xlim=(0,1.6),ylim=(241,191),xticks=[0,.5,1,1.5],yticks=[],xlabel=self.t['unitsC'])
        c.xaxis.grid(True,color=GRID,linewidth=.5);self.style(c)
        self.primitives['C']=dict(points=cpoints,ratio_limit=1,x_limits=[0,1.6],nominal_limit=.06,tail_limits_are_original_cutoffs=True)
        r=DATA['resolution'];d=self.axis('D',137,191,40,50)
        values=[r['by34_cutoff'],r['zero_event_upper']];dy=[202.,226.]
        for y,label in zip(dy,self.t['resolution_rows']):self.text(96,y,label,8.8)
        for value,y,marker,face,edge in zip(values,dy,['s','D'],['#B0C6A8','#D99DA3'],[GREEN,ROSE]):
            d.hlines(y,0,value,color=edge,lw=1.2);d.plot(value,y,linestyle='none',marker=marker,ms=4.6,mfc=face,mec=edge,mew=.65)
        self.text(156.5,195.5,'3.57 × 10⁻⁴',8.8,ha='center',math=True)
        self.text(156.5,219.5,'9.59 × 10⁻⁴',8.8,ha='center',math=True)
        d.set(xlim=(0,.00115),ylim=(241,191),xticks=[0,.0005,.001],yticks=[],xlabel=self.t['unitsD'])
        d.set_xticklabels(['0','0.0005','0.001']);d.xaxis.grid(True,color=GRID,linewidth=.5);self.style(d)
        self.primitives['D']=dict(values=values,quantities=['Existing BY-34 first-step cutoff','Saved zero-of-4000 Wilson upper bound'],
            draws=4000,events=0,tail_not_resolved=True,not_an_empirical_BY_tail_estimate=True,x_limits=[0,.00115])
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
        title='Supplementary Figure S7: Calibration of statistical tests under the observed study designs'
        self.fig.savefig(BASE/f'Figure_S7_{self.lang}.pdf',metadata={'Title':title,'Creator':'Study supplementary figures'})
        self.fig.savefig(BASE/f'editable/Figure_S7_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_S7_{self.lang}.png',dpi=300)
        tif=BASE/f'high_resolution/Figure_S7_{self.lang}.tiff'
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
        palette=PALETTE,choice='Soft blue for nominal rejection, sage for qualified designs, muted rose for failed or unresolved limits.',
        figure_size_mm=[W,H],panels=4,nominal_groups=8,tail_design_points=291,cst3_checks=8,
        unresolved_relative_designs=35,power_not_shown=True,only_categorical_vertical_spread=True,
        quantitative_coordinates_not_jittered=True,loose_footer_notes=False),ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(results,ensure_ascii=False),flush=True)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');main()
