"""Model-shift arrows, signed change intervals, correlation triads and paired evidence bars."""
from pathlib import Path
import json,logging,warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch,FancyBboxPatch
from matplotlib.font_manager import FontProperties
from matplotlib.text import Text
from PIL import Image
BASE=Path(__file__).resolve().parents[1]
D=json.loads((BASE/'audit/PLOT_INPUT.json').read_text(encoding='utf-8'))
W,H=180,230
INK='#293C4B';GRID='#E8EDF0';ROSE='#A96789'
FILLS={'NULISA':'#D4977D','SomaScan':'#76ACA0','TMT-MS':'#92A5C3'}
EDGES={'NULISA':'#A3674F','SomaScan':'#3D7F73','TMT-MS':'#5B759C'}
SHAPES={'NULISA':'o','SomaScan':'s','TMT-MS':'D'}
DIFF='#607F96';SHIFT='#BB87A3'
logging.getLogger('fontTools').setLevel(logging.ERROR)
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','font.family':'Arial',
 'figure.facecolor':'white','savefig.facecolor':'white','axes.unicode_minus':False,'axes.linewidth':.65,
 'xtick.major.size':2.6,'ytick.major.size':2.6,'xtick.major.width':.65,'ytick.major.width':.65})
TEXT={
 'EN':dict(A='Memory associations before and after adjustment',B='Direct changes in assay differences',C='Observed residual correlations',D='After background adjustment',
  before='Without p-tau181',after='With p-tau181',memory='Memory association',unit='ADNI-MEM units/year\nper protein residual SD',
  bfamily='Joint change: nine comparisons',changeunit='Change in assay contrast\nADNI-MEM units/year per protein residual SD',
  csign='Solid: positive; dashed: negative',cwidth='Line width follows |r|',slope='Memory slope',
  dsub='Three background scores; 33 comparisons',remaining='Assay difference',additional='p-tau adjustment change',evidence='−log$_{10}$ q'),
 'ZH':dict(A='调整前后的记忆关联',B='跨检测差异的直接变化',C='观察到的残差相关',D='蛋白背景调整后的结果',
  before='不调整p-tau181',after='调整p-tau181',memory='记忆关联系数',unit='ADNI-MEM分/年\n每蛋白残差SD',
  bfamily='联合变化检验：原始9项校正范围',changeunit='跨检测关联对比的变化\nADNI-MEM分/年/蛋白残差SD',
  csign='实线：正相关；虚线：负相关',cwidth='线宽随|r|增加',slope='记忆斜率',
  dsub='纳入三个背景分数；33项校正范围',remaining='跨检测差异',additional='p-tau调整变化',evidence='−log$_{10}$ q')}

def qtext(q):
    if q<.001:
        m,e=f'{q:.2e}'.split('e');return '$'+m+r'\times10^{'+str(int(e))+'}$'
    return f'{q:.3g}'

class Figure:
    def __init__(self,lang):
        self.lang=lang;self.t=TEXT[lang];self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=180);self.axes={};self.records={}
    def font(self,size=8,bold=False):
        name=('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN' else ('msyhbd.ttc' if bold else 'msyh.ttc')
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/name),family='Arial' if self.lang=='EN' else 'Microsoft YaHei',size=size,weight='bold' if bold else 'normal')
    def text(self,x,y,s,size=8,bold=False,ha='left',color=INK,rotation=0):
        return self.fig.text(x/W,1-y/H,s,fontproperties=self.font(size,bold),color=color,ha=ha,va='center',linespacing=1.1,rotation=rotation,zorder=12)
    def panel(self,letter,title,x,y):
        self.text(x,y,letter,11.5,True);self.text(x+7,y,title,8.9,True)
    def axis(self,key,x,y,w,h):
        ax=self.fig.add_axes([x/W,1-(y+h)/H,w/W,h/H]);self.axes[key]=ax
        ax.spines[['top','right']].set_visible(False);ax.spines[['bottom','left']].set_color('#7F929D');ax.tick_params(colors=INK,pad=3,direction='out')
        return ax
    def ticks(self,ax,size=7.4):
        for t in ax.get_xticklabels()+ax.get_yticklabels():t.set_fontproperties(self.font(size));t.set_color(INK)
    def glyph(self,x,y,marker,color,filled=True,size=3.8):
        self.fig.lines.append(Line2D([x/W],[1-y/H],transform=self.fig.transFigure,marker=marker,markersize=size,
             markerfacecolor=color if filled else 'white',markeredgecolor=color,markeredgewidth=.7,linestyle='none',zorder=10))
    def shifts(self):
        self.panel('A',self.t['A'],3,7)
        for x,label,filled in [(8,self.t['before'],False),(48,self.t['after'],True)]:
            self.glyph(x,17,'o','#718795',filled,3.8);self.text(x+3,17,label,7.3)
        self.text(57.5,24.5,self.t['memory'],7.8,True,'center')
        ax=self.axis('A',29,29,57,83);ax.set(xlim=(-.025,.15),ylim=(11.8,-.6),xticks=[-.02,0,.05,.10,.15],yticks=[])
        ax.set_xticklabels(['−0.02','0','0.05','0.10','0.15']);ax.spines['left'].set_visible(False)
        ax.axvline(0,color='#A2AFB7',linestyle=(0,(3,2)),linewidth=.7,zorder=0)
        positions=[]
        for k,g in enumerate(D['genes']):
            header_y=29+(.6+4*k)*83/12.4
            self.text(3,header_y,g,8.2,True)
            n=next(r['Participants'] for r in D['coefficient_shifts'] if r['Protein']==g)
            self.text(18.5,header_y,f'n={n}',7.2)
            if k%2==0:ax.axhspan(4*k+.5,4*k+3.5,color='#F6F8F9',zorder=-1)
            for j,a in enumerate(D['assays']):
                r=next(r for r in D['coefficient_shifts'] if r['Protein']==g and r['Assay']==a)
                yy=4*k+j+1;physical_y=29+(yy+.6)*83/12.4
                self.text(3,physical_y,a,7.5,color=EDGES[a])
                start,end=r['Without p-tau181'],r['With p-tau181']
                ax.plot([start,end],[yy,yy],color=FILLS[a],alpha=.35,linewidth=3,zorder=1)
                ax.add_patch(FancyArrowPatch((start,yy),(end,yy),arrowstyle='-|>',color=EDGES[a],linewidth=.8,mutation_scale=6,shrinkA=2,shrinkB=2,zorder=3))
                ax.plot(start,yy,marker=SHAPES[a],markersize=4.0,markerfacecolor='white',markeredgecolor=EDGES[a],markeredgewidth=.75,linestyle='none',zorder=5)
                ax.plot(end,yy,marker=SHAPES[a],markersize=4.0,markerfacecolor=FILLS[a],markeredgecolor=EDGES[a],markeredgewidth=.65,linestyle='none',zorder=6)
                positions.append({**r,'y':yy})
        self.ticks(ax);self.text(57.5,125.5,self.t['unit'],7.6,False,'center')
        self.records['A']={'pairs':positions,'x_limits':[-.025,.15],'point_estimates_only':True,'same_participants_visits_weights_scales':True}
    def direct_changes(self):
        self.panel('B',self.t['B'],95,7)
        for x,label,assay in [(99,'SomaScan − NULISA','SomaScan'),(144,'TMT-MS − NULISA','TMT-MS')]:
            self.glyph(x,17,SHAPES[assay],EDGES[assay],True,3.7);self.text(x+3,17,label,7.1,color=EDGES[assay])
        self.text(136,27,self.t['bfamily'],7.4,False,'center')
        ax=self.axis('B',99,39,78,73);ax.set(xlim=(-.06,.10),ylim=(2.75,-.75),xticks=[-.05,0,.05,.10],yticks=[])
        ax.set_xticklabels(['−0.05','0','0.05','0.10']);ax.spines['left'].set_visible(False)
        ax.axvline(0,color='#8497A3',linestyle=(0,(3,2)),linewidth=.8,zorder=0)
        points=[]
        for i,g in enumerate(D['genes']):
            py=39+(i+.75)*73/3.5
            self.text(95,py-9,g,8.2,True)
            q=next(r['Joint q'] for r in D['joint_change_tests'] if r['Protein']==g)
            self.text(154,py-9,'q='+qtext(q),7.7,False,'center')
            for offset,assay,contrast in [(-.14,'SomaScan','SomaScan − NULISA'),(.14,'TMT-MS','TMT-MS − NULISA')]:
                r=next(r for r in D['contrast_changes'] if r['Protein']==g and r['Contrast']==contrast)
                b=r['Change'];lo=r['CI lower'];hi=r['CI upper']
                ax.errorbar(b,i+offset,xerr=[[b-lo],[hi-b]],fmt=SHAPES[assay],markersize=4.6,mfc=FILLS[assay],mec=EDGES[assay],
                    mew=.7,color=EDGES[assay],elinewidth=1.15,capsize=2.7,capthick=.8,zorder=4)
                points.append({**r,'y':i+offset})
        self.ticks(ax);self.text(138,125.5,self.t['changeunit'],7.5,False,'center')
        self.records['B']={'intervals':points,'joint_tests':D['joint_change_tests'],'family':9,'x_limits':[-.06,.10],
            'positive_is_increase_in_signed_contrast':True,'intervals_are_pointwise_not_joint':True}
    def network_line(self,pts,color,r,parallel_offset=0):
        p0,p1=np.asarray(pts[0],float),np.asarray(pts[1],float);v=p1-p0;normal=np.array([-v[1],v[0]])/np.linalg.norm(v)
        p0=p0+parallel_offset*normal;p1=p1+parallel_offset*normal
        self.fig.lines.append(Line2D([p0[0]/W,p1[0]/W],[1-p0[1]/H,1-p1[1]/H],transform=self.fig.transFigure,
            color=color,linewidth=.55+2.1*abs(r),linestyle='-' if r>=0 else (0,(3,2)),zorder=2))
    def node(self,x,y,w,h,label,fill):
        self.fig.patches.append(FancyBboxPatch(((x-w/2)/W,1-(y+h/2)/H),w/W,h/H,
            transform=self.fig.transFigure,boxstyle='round,pad=0,rounding_size=.008',facecolor=fill,edgecolor='#D2DCE3',linewidth=.6,zorder=5))
        self.text(x,y,label,7.4,True,'center')
    def correlations(self):
        self.panel('C',self.t['C'],3,142)
        self.text(44,150,self.t['csign'],7.2,False,'center');self.text(44,155.4,self.t['cwidth'],7.2,False,'center')
        network=[]
        for center,g in zip([169,192,215],D['genes']):
            values=[r for r in D['residual_correlations'] if r['Protein']==g]
            for offset,assay in zip([-.9,0,.9],D['assays']):
                r=next(r for r in values if r['Assay']==assay)
                self.network_line([(29,center),(59,center-7)],EDGES[assay],r['Protein–p-tau r'],offset)
                self.network_line([(29,center),(59,center+7)],EDGES[assay],r['Protein–memory slope r'],offset)
                network.append(r)
            rho=values[0]['p-tau–memory slope r']
            self.network_line([(69,center-4),(69,center+4)],ROSE,rho)
            self.node(18,center,22,7,g,'#F2F5F7');self.node(69,center-7,20,6,'p-tau181','#F7EDF3')
            self.node(69,center+7,20,6,self.t['slope'],'#F0EDF6')
            self.text(86,center,'r='+f'{rho:.3f}'.replace('-','−'),7.1,False,'right',ROSE)
        self.records['C']={'correlations':network,'line_width_rule':'0.55 + 2.1 * abs(r) in points',
             'all_edges_undirected':True,'common_covariates_adjusted_only':True,'third_variable_not_conditioned_on':True,
             'no_significance_tests_added':True}
    def background(self):
        self.panel('D',self.t['D'],95,142)
        self.text(137,153.3,self.t['dsub'],7.3,False,'center')
        for x,label,color,hatch in [(99,self.t['remaining'],DIFF,None),(136,self.t['additional'],SHIFT,'///')]:
            self.fig.patches.append(plt.Rectangle((x/W,1-161.2/H),2.8/W,2.5/H,transform=self.fig.transFigure,
                facecolor=color,edgecolor='#875C79' if hatch else color,linewidth=.6,hatch=hatch,zorder=4))
            self.text(x+4,160,label,7.2)
        ax=self.axis('D',107,170,69,46);ax.set(xlim=(-.5,2.5),ylim=(0,2),xticks=[0,1,2],yticks=[0,1,2])
        ax.set_xticklabels(D['genes']);ax.spines['left'].set_color('#80939E')
        ax.axhline(D['cutoff_strength'],color='#9E93A0',linewidth=.8,linestyle=(0,(3,2)),zorder=1)
        for tick in [1,2]:ax.axhline(tick,color=GRID,linewidth=.5,zorder=0)
        bars=[]
        for i,g in enumerate(D['genes']):
            for offset,test,color,hatch in [(-.16,'Remaining assay differences',DIFF,None),(.16,'Additional p-tau181 change',SHIFT,'///')]:
                r=next(r for r in D['background_tests'] if r['Protein']==g and r['Test']==test)
                value=r['Evidence strength'];x=i+offset
                ax.bar(x,value,width=.26,color=color,edgecolor='#875C79' if hatch else color,linewidth=.65,hatch=hatch,zorder=3)
                if value==0:ax.plot(x,0,marker='_',markersize=7,color=color,markeredgewidth=1.3,zorder=5,clip_on=False)
                ax.text(x,value+.065,'q='+f"{r['q']:.3g}",fontproperties=self.font(7.1),color=INK,ha='center',va='bottom')
                bars.append({**r,'x':x,'height':value})
        self.ticks(ax,7.5);self.text(96.5,193,self.t['evidence'],7.6,False,'center',rotation=90)
        self.records['D']={'bars':bars,'family':33,'y_limits':[0,2],'cutoff_strength':D['cutoff_strength'],
             'two_distinct_hypotheses':True,'all_three_background_scores_included':True,'bar_labels_are_q_values':True}
    def layout(self):
        self.fig.canvas.draw();renderer=self.fig.canvas.get_renderer();records=[];seen=set();overlap=[]
        for a in self.fig.findobj(Text):
            if not a.get_visible() or not a.get_text().strip() or id(a) in seen:continue
            seen.add(id(a));b=a.get_window_extent(renderer);v=np.asarray([b.x0,b.y0,b.x1,b.y1])/self.fig.dpi*25.4
            r={'text':a.get_text(),'font_pt':float(a.get_fontsize()),'x0_mm':float(v[0]),'y0_mm':float(H-v[3]),'x1_mm':float(v[2]),'y1_mm':float(H-v[1])}
            assert r['x0_mm']>=.1 and r['x1_mm']<=W-.1 and r['y0_mm']>=.1 and r['y1_mm']<=H-.1,('Clipped label',self.lang,r)
            records.append(r)
        for i,a in enumerate(records):
            for b in records[i+1:]:
                xx=min(a['x1_mm'],b['x1_mm'])-max(a['x0_mm'],b['x0_mm']);yy=min(a['y1_mm'],b['y1_mm'])-max(a['y0_mm'],b['y0_mm'])
                if xx>.22 and yy>.22:overlap.append([a['text'],b['text'],round(xx,2),round(yy,2)])
        assert not overlap,('Text overlap',self.lang,overlap)
        return records
    def save(self):
        self.shifts();self.direct_changes();self.correlations();self.background();labels=self.layout()
        title='Figure 4. p-tau181 adjustment and between-assay differences in memory associations'
        self.fig.savefig(BASE/f'Figure_4_{self.lang}.pdf',metadata={'Title':title,'Creator':'Paired ADNI CSF proteomic study'})
        self.fig.savefig(BASE/f'editable/Figure_4_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_4_{self.lang}.png',dpi=300);self.fig.savefig(BASE/f'Figure_4_{self.lang}_preview.png',dpi=180)
        tif=BASE/f'high_resolution/Figure_4_{self.lang}.tiff';self.fig.savefig(tif,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        with Image.open(tif) as im:
            rgb=im.convert('RGB');rgb.save(tif,compression='tiff_lzw',dpi=(600,600));rgb.close()
        (BASE/f'audit/PLOT_PRIMITIVES_{self.lang}.json').write_text(json.dumps(self.records,ensure_ascii=False,indent=2),encoding='utf-8')
        (BASE/f'audit/LAYOUT_{self.lang}.json').write_text(json.dumps(labels,ensure_ascii=False,indent=2),encoding='utf-8')
        return {'language':self.lang,'size_mm':[W,H],'panels':4,'minimum_font_pt':min(r['font_pt'] for r in labels),
           'text_overlap':False,'text_clipped':False,'model_shift_pairs':9,'contrast_intervals':6,'background_test_bars':6}
def main():
    warnings.filterwarnings('error',message='Glyph.*missing');reports=[]
    for lang in ['EN','ZH']:
        f=Figure(lang);reports.append(f.save());plt.close(f.fig)
    (BASE/'audit/LAYOUT_VALIDATION.json').write_text(json.dumps(reports,ensure_ascii=False,indent=2),encoding='utf-8')
    (BASE/'audit/STYLE_SPECIFICATION.json').write_text(json.dumps({'size_mm':[W,H],'assay_fills':FILLS,'assay_edges':EDGES,
         'assay_shapes':SHAPES,'evidence_bars':{'Remaining assay differences':DIFF,'Additional p-tau181 change':SHIFT},
         'formats':['paired model-shift arrows','signed paired-change intervals','weighted residual correlation triads','two-question evidence bars'],
         'no_causal_arrowheads_in_networks':True,'outside_grey_explanatory_footer':False},ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(reports))
if __name__=='__main__':main()
