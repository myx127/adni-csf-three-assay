"""Original vector illustrations of study design and methods; no observed results."""
from pathlib import Path
import json,logging,warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Ellipse,Polygon,FancyBboxPatch,PathPatch,FancyArrowPatch
from matplotlib.path import Path as MPath
from matplotlib.lines import Line2D
from matplotlib.font_manager import FontProperties
from PIL import Image
BASE=Path(__file__).resolve().parents[1]
DATA=json.loads((BASE/'audit/PLOT_INPUT.json').read_text(encoding='utf-8'))
W,H=180,158
INK='#293C4B';EDGE='#D0D9E0';LINE='#8496A1';DEEP='#3E6178'
FILLS={'NULISA':'#D4977D','SomaScan':'#76ACA0','TMT-MS':'#92A5C3'}
EDGES={'NULISA':'#A3674F','SomaScan':'#3D7F73','TMT-MS':'#5B759C'}
TINTS={'NULISA':'#F7EBE5','SomaScan':'#EAF3EF','TMT-MS':'#EDF1F7'}
ASSAYS=DATA['assays'];COG='#9783B0';PTAU='#B47791'
logging.getLogger('fontTools').setLevel(logging.ERROR)
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white',
 'savefig.facecolor':'white','font.family':'Arial','axes.unicode_minus':False})
WORDS={
 'EN':{
  'head':['ADNI CSF cohort','Three assay technologies','Paired protein comparisons','Shared p-tau181 associations','Cognitive associations','Exploratory follow-up'],
  'cohort':'571 paired participants','repeated':'Repeated cognition','window':'Follow-up up to five years\nafter CSF collection','scores':'ADNI-MEM · ADNI-EF',
  'mechanisms':['Antibody–DNA detection','Aptamer detection','Peptide mass spectrometry'],
  'matched':'Paired by participant and CSF date','proteins':'34 eligible proteins','units':'Comparable effect units','sd':'Protein residual SD',
  'elecsys':'Elecsys p-tau181','direction':'Shared direction?','positive':'Positive in all three','negative':'Negative in all three','or':'or',
  'primary':'Annual memory change','same':'Shared direction?','different':'Different\nmagnitudes?',
  'secondary':['Baseline memory','Executive function change'],
  'adjust':'p-tau181 adjustment','without':'Without p-tau181','with':'With p-tau181',
  'paired':'Direct paired comparison','question':'Does the assay difference change?',
  'background':'Background protein scores','fixed':['Same participants, visits and scales'],'xlab':'p-tau181','ylab':'Protein signal'},
 'ZH':{
  'head':['ADNI脑脊液队列','三种检测技术','逐蛋白配对比较','共同p-tau181关联','认知关联','探索性后续分析'],
  'cohort':'571名配对参与者','repeated':'重复认知评价','window':'CSF采集后最长五年\n认知随访','scores':'ADNI-MEM · ADNI-EF',
  'mechanisms':['抗体识别与DNA读出','适配体识别','肽段质谱定量'],
  'matched':'同一参与者 · 同一采集日期','proteins':'34个符合条件的蛋白','units':'可比较的关联效应单位','sd':'蛋白残差标准差',
  'elecsys':'Elecsys测定p-tau181','direction':'是否具有共同方向？','positive':'三种检测均为正向','negative':'三种检测均为负向','or':'或',
  'primary':'年度记忆变化','same':'方向是否共同？','different':'关联强度\n是否不同？',
  'secondary':['基线记忆','年度执行功能变化'],
  'adjust':'p-tau181调整','without':'不调整p-tau181','with':'调整p-tau181',
  'paired':'直接配对比较','question':'调整是否改变跨检测差异？',
  'background':'蛋白背景分数','fixed':['同一参与者、访视与固定蛋白尺度'],'xlab':'p-tau181','ylab':'蛋白信号'}}
class Figure:
    def __init__(self,lang):
        self.lang=lang;self.t=WORDS[lang]
        self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=180)
        self.ax=self.fig.add_axes([0,0,1,1]);self.ax.set(xlim=(0,W),ylim=(H,0));self.ax.axis('off')
        self.texts=[];self.boxes=[];self.primitives={'stages':DATA['stages'],'connections':DATA['connections'],
            'all_graphics_schematic':True,'observed_associations_plotted':False,'clinical_analysis_gated_on_discoveries':False}
    def font(self,size=8,bold=False):
        filename=('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN' else ('msyhbd.ttc' if bold else 'msyh.ttc')
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/filename),family='Arial' if self.lang=='EN' else 'Microsoft YaHei',size=size,weight='bold' if bold else 'normal')
    def text(self,x,y,s,size=8,bold=False,ha='center',color=INK,box=None,role='label',rotation=0):
        a=self.ax.text(x,y,s,fontproperties=self.font(size,bold),ha=ha,va='center',color=color,linespacing=1.15,zorder=20,rotation=rotation)
        self.texts.append({'artist':a,'box':box,'role':role});return a
    def rect(self,x,y,w,h,fill='white',edge=EDGE,r=.9,lw=.65,z=1):
        self.ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f'round,pad=0,rounding_size={r}',facecolor=fill,edgecolor=edge,linewidth=lw,zorder=z))
        return [x,y,x+w,y+h]
    def line(self,pts,color=LINE,lw=.8,z=3):
        self.ax.add_line(Line2D([p[0] for p in pts],[p[1] for p in pts],color=color,linewidth=lw,solid_capstyle='round',zorder=z))
    def arrow(self,a,b,double=False,color=LINE,lw=.8,scale=6):
        self.ax.add_patch(FancyArrowPatch(a,b,arrowstyle='<->' if double else '-|>',mutation_scale=scale,
            color=color,linewidth=lw,shrinkA=0,shrinkB=0,zorder=4))
    def circle(self,x,y,r,fill,edge='none',lw=.6,z=5):
        self.ax.add_patch(Circle((x,y),r,facecolor=fill,edgecolor=edge,linewidth=lw,zorder=z))
    def person(self,x,y,s=1,color='#7095A6'):
        self.circle(x,y+1.7*s,1.55*s,color)
        self.ax.add_patch(Polygon([(x-2.7*s,y+4.1*s),(x+2.7*s,y+4.1*s),(x+1.9*s,y+9.0*s),(x-1.9*s,y+9.0*s)],
             closed=True,facecolor=color,edgecolor='none',zorder=5))
        for dx in [-1,1]:self.line([(x+dx*s,y+8.7*s),(x+dx*s,y+10.6*s)],color,1.7*s,5)
        self.line([(x-2.2*s,y+5*s),(x-4.1*s,y+9.3*s)],color,1.2*s,5)
        self.line([(x+2.2*s,y+5*s),(x+4.1*s,y+9.3*s)],color,1.2*s,5)
    def tube(self,x,y,w=3.5,h=11,color='#A7CDDF'):
        self.rect(x,y,w,h,'#FCFDFE','#6B8190',w/2,.65,5)
        self.rect(x+.4,y+h*.40,w-.8,h*.52,color,'none',.55,.65,6)
        self.rect(x-.4,y-.8,w+.8,1.7,'#A9BBC8','#6B8190',.3,.6,7)
    def brain(self,x,y,s=1,color=COG):
        for dx,dy,rx,ry in [(-2.8,-.9,2.3,2.7),(-1.9,-3.0,2.3,2.0),(.6,-3.3,2.5,1.9),(3.1,-1.7,2.2,2.4),(2.8,1,2.1,2.3),(.2,2,2.4,2),(-2.7,1.2,2,2.2)]:
            self.ax.add_patch(Ellipse((x+dx*s,y+dy*s),2*rx*s,2*ry*s,facecolor=color,edgecolor='none',zorder=5))
        pale='#F9F6FC'
        self.line([(x-.3*s,y-4.2*s),(x-.7*s,y-1.5*s),(x+.4*s,y+.1*s),(x+.1*s,y+3.2*s)],pale,.65,6)
        self.line([(x-3.8*s,y-1.5*s),(x-2.1*s,y-.2*s),(x-2.8*s,y+1.5*s)],pale,.55,6)
        self.line([(x+2.1*s,y-3.0*s),(x+1.4*s,y-1.1*s),(x+3.4*s,y-.4*s)],pale,.55,6)
    def form(self,x,y,w=5.5,h=5):
        self.rect(x,y,w,h,'#F7F4FA','#AB9ABB',.4,.6,5)
        for d in [1.4,2.6,3.8]:self.line([(x+1,y+d),(x+w-1,y+d)],'#9680AD',.5,6)
    def protein(self,x,y,s=1,color='#B4A4C5',edge='none'):
        pts=np.asarray([[-2.4,-.5],[-1.7,-1.9],[-.3,-2.2],[.7,-1.5],[2.2,-1.1],[2.5,.2],[1.8,1.5],[.4,1.3],[-.8,2.0],[-2,1.2]])*s+[x,y]
        self.ax.add_patch(Polygon(pts,closed=True,facecolor=color,edgecolor=edge,linewidth=.55,zorder=7))
    def antibody(self,x,y,s=1,color=EDGES['NULISA']):
        self.line([(x,y+3.0*s),(x,y),(x-2.2*s,y-2.1*s)],color,1.2,6)
        self.line([(x,y),(x+2.2*s,y-2.1*s)],color,1.2,6)
        self.line([(x-.9*s,y+2.6*s),(x-.9*s,y+.2*s),(x-2.9*s,y-1.6*s)],color,.65,6)
        self.line([(x+.9*s,y+2.6*s),(x+.9*s,y+.2*s),(x+2.9*s,y-1.6*s)],color,.65,6)
    def dna(self,x,y,w=6,h=3,color=EDGES['NULISA']):
        t=np.linspace(0,1,50);xx=x+t*w;yy=y+np.sin(t*3*np.pi)*h/2
        self.line(list(zip(xx,yy)),color,.65,6);self.line(list(zip(xx,2*y-yy)),color,.65,6)
        for t in np.linspace(.08,.92,6):
            q=np.sin(t*3*np.pi)*h/2;self.line([(x+t*w,y+q),(x+t*w,y-q)],color,.45,6)
    def aptamer(self,x,y,color=EDGES['SomaScan']):
        path=MPath([(x-3,y+4),(x-2.3,y+1),(x-3.4,y-1.3),(x-1.8,y-4.2),(x+.4,y-4),(x+1.7,y-2.7),(x+.8,y-.8),(x-1.1,y-.4),(x+.2,y+2.3),(x+3,y+3.9)],
            [MPath.MOVETO]+[MPath.CURVE3]*8+[MPath.LINETO])
        self.ax.add_patch(PathPatch(path,facecolor='none',edgecolor=color,linewidth=1.1,zorder=6))
        self.protein(x+3.2,y-2.1,.8,'#C5A7BD')
        for xx,yy in [(x-2.9,y+3),(x-2.9,y-.9),(x-.4,y-3.9),(x+.9,y-1.6),(x,y+2.3)]:self.circle(xx,yy,.38,'#C6E0D7',color,.35,7)
    def mass_spec(self,x,y):
        self.rect(x-5,y-5,10,10.8,'#E7ECF1','#7A8D9C',1,.65,5)
        self.rect(x-3.6,y-3.9,7.2,4.9,'#2A4357','none',.5,.65,6)
        self.line([(x-2.7,y-.1),(x-2,y-.4),(x-1.6,y-2.8),(x-1,y-.3),(x-.3,y-1.5),(x+.5,y-.3),(x+1.5,y-2.3),(x+2.3,y-.2)],'#B7D4E9',.55,7)
        self.circle(x-1.8,y+3.2,.65,'#8194A4');self.circle(x+1.7,y+3.2,.65,'#8194A4')
    def tiny_trend(self,x,y,w,h,color,up=True):
        self.line([(x,y),(x,y+h),(x+w,y+h)],'#B2BEC6',.55,4)
        yy=[y+h*.8,y+h*.57,y+h*.22] if up else [y+h*.2,y+h*.45,y+h*.8]
        self.line([(x+w*.18,yy[0]),(x+w*.52,yy[1]),(x+w*.87,yy[2])],color,1.1,6)
    def header(self,id,x,y,w,h):
        box=self.rect(x,y,w,h,'white',EDGE,.7,.65)
        self.rect(x+.2,y+.2,w-.4,10.2,'#F4F7F9','none',.6,0,2)
        self.circle(x+4.5,y+5.1,2.6,DEEP,z=6)
        self.text(x+4.5,y+5.05,str(id),9.1,True,color='white',role='stage number')
        self.text(x+8.5,y+5.1,self.t['head'][id-1],8.35,True,'left',box=box,role='stage title')
        self.boxes.append(box);return box
    def cohort(self,b):
        for x,c in zip([10,18,26,34],['#C6D5DF','#8CB0C1','#4D859D','#8EB6C5']):self.person(x,16,.84,c)
        self.tube(46,17,3.4,10)
        self.text(31,31,self.t['cohort'],8.6,True,box=b)
        self.text(31,38.2,self.t['repeated'],8.0,True,box=b)
        self.line([(12,47),(49,47)],'#9BAAB4',.7)
        for x in [13,29,47]:
            self.circle(x,47,.7,'#8799A7');self.brain(x,44,.45);self.form(x-2.6,49.6,5.2,4.5)
        self.text(31,59,self.t['scores'],7.8,box=b)
        self.text(31,65,self.t['window'],7.9,True,box=b)
    def platforms(self,b):
        for i,assay in enumerate(ASSAYS):
            y=16.5+i*16.2
            self.rect(65,y,50,14.6,TINTS[assay],'none',.9)
            if i==0:
                self.antibody(71,y+5,.8);self.dna(67.5,y+11,7,2)
            elif i==1:self.aptamer(71,y+7)
            else:self.mass_spec(71,y+7)
            self.text(80,y+4.2,assay,8.3,True,'left',EDGES[assay],box=b)
            self.text(80,y+10,self.t['mechanisms'][i],7.45,False,'left',box=b)
    def pairing(self,b):
        self.text(149,18.3,self.t['matched'],7.7,box=b)
        self.text(149,25.3,self.t['proteins'],8.8,True,box=b)
        for y in [33,39,45]:
            self.line([(130,y),(168,y)],'#D7E1E6',.65)
            for j,assay in enumerate(ASSAYS):self.protein(134+j*15,y,.75,FILLS[assay],EDGES[assay])
        self.text(149,54.8,self.t['units'],8.0,True,box=b)
        for j,assay in enumerate(ASSAYS):
            x=130+j*15;self.line([(x,61),(x+8,61)],EDGES[assay],.85)
            for dx in [0,2,4,6,8]:self.line([(x+dx,60.4),(x+dx,61.6)],EDGES[assay],.6)
        self.text(149,66.5,self.t['sd'],7.7,box=b)
    def pathology(self,b):
        self.tube(8.5,93,2.9,8);self.protein(12.7,94.2,.65,PTAU)
        self.text(35,96,self.t['elecsys'],7.8,True,box=b)
        self.text(31,104,self.t['direction'],8.1,True,box=b)
        for i,assay in enumerate(ASSAYS):self.tiny_trend(9+i*15,110,10,7,EDGES[assay],True)
        self.text(31,122,self.t['positive'],7.8,box=b)
        self.text(31,128,self.t['or'],7.6,box=b)
        for i,assay in enumerate(ASSAYS):self.tiny_trend(9+i*15,133,10,7,EDGES[assay],False)
        self.text(31,146,self.t['negative'],7.8,box=b)
        self.text(31,152,self.t['xlab'],7.3,box=b)
        self.text(6,128,self.t['ylab'],7.3,box=b,rotation=90)
    def cognition(self,b):
        self.text(90,96,self.t['primary'],8.1,True,box=b)
        self.brain(78,105,.85)
        self.circle(103,104.4,4.2,'#F3F6F8','#98ACB8',.65)
        self.line([(103,101.4),(103,104.4),(105.5,105.7)],'#708A9B',.8,6)
        self.circle(103,104.4,.3,'#708A9B',z=7)
        self.line([(85,105),(97.7,105)],'#B0BDC7',.65)
        self.rect(66,114,48,12,'#EDF5F1','#C7DCD3',1)
        for i,assay in enumerate(ASSAYS):self.arrow((70+i*4.2,122.7),(70+i*4.2,117),color=EDGES[assay],lw=.95,scale=5.0)
        self.text(97,120,self.t['same'],7.8,True,box=b)
        self.rect(66,129,48,13,'#F0F3F8','#CDD6E5',1)
        for i,(assay,length) in enumerate(zip(ASSAYS,[7.4,4.6,9.7])):
            self.line([(69,132+i*2.9),(69+length,132+i*2.9)],EDGES[assay],1.7)
        self.text(97,135.5,self.t['different'],7.8,True,box=b)
        self.text(90,147.1,self.t['secondary'][0],7.6,box=b)
        self.text(90,152,self.t['secondary'][1],7.6,box=b)
    def model_glyph(self,x,y,adjusted=False):
        for i,assay in enumerate(ASSAYS):self.protein(x+i*4.2,y,.50,FILLS[assay],EDGES[assay])
        self.line([(x+10.5,y),(x+24.5,y)],'#9EADB7',.65)
        self.brain(x+28,y,.58)
        if adjusted:
            self.circle(x+7,y-3,1.3,'#EFDCE7',PTAU,.6)
            self.text(x+7,y-2.95,'p',7.1,True,color='#8C5873',role='p-tau symbol')
    def conditioning(self,b):
        self.text(149,95,self.t['adjust'],8.3,True,box=b)
        self.rect(125,100,48,16,'#F6F8FA','#D5DEE4',1)
        self.text(149,104.3,self.t['without'],7.8,True,box=b)
        self.model_glyph(131.5,111)
        self.arrow((127.8,116.6),(127.8,121.7),double=True,lw=.8,scale=5)
        self.text(151,119.5,self.t['paired'],7.3,True,box=b)
        self.rect(125,123,48,16,'#F8F0F5','#E0CBD6',1)
        self.text(149,127.2,self.t['with'],7.8,True,box=b)
        self.model_glyph(131.5,134,True)
        self.text(149,144,self.t['question'],7.65,True,box=b)
        for x,y in [(127.1,149.5),(130.0,148.7),(132.4,150)]:self.protein(x,y,.35,'#BCCDBF')
        self.text(152,149.3,self.t['background'],7.4,box=b)
        self.text(149,153,self.t['fixed'][0],7.35,box=b)
    def draw(self):
        specs=[(1,3,3,56,67),(2,62,3,56,67),(3,121,3,56,67),
               (4,3,79,56,76),(5,62,79,56,76),(6,121,79,56,76)]
        panels=[self.header(*s) for s in specs]
        for fn,b in zip([self.cohort,self.platforms,self.pairing,self.pathology,self.cognition,self.conditioning],panels):fn(b)
        self.arrow((59.2,36),(61.8,36),scale=5.5)
        self.arrow((118.2,36),(120.8,36),scale=5.5)
        self.line([(149,70.2),(149,74.8),(31,74.8)],'#7B939F',.8)
        self.arrow((31,74.8),(31,78.8),scale=5.5)
        self.arrow((90,74.8),(90,78.8),scale=5.5)
        self.arrow((118.2,116),(120.8,116),scale=5.5)
        self.primitives.update({'figure_size_mm':[W,H],'diagram_counts':{'participants':571,'proteins':34,'window_years':5},
            'shared_direction_depicts_both_possible_signs':True,'same_participants_visits_scales_for_adjustment':True,
            'p_tau_symbol_is_a_covariate_not_causal_path':True,'assay_colour_map':FILLS})
        return self
    def layout(self):
        self.fig.canvas.draw();renderer=self.fig.canvas.get_renderer();records=[];overlaps=[]
        for spec in self.texts:
            a=spec['artist'];b=a.get_window_extent(renderer)
            x0,y0,x1,y1=np.asarray([b.x0,b.y0,b.x1,b.y1])/self.fig.dpi*25.4
            top,bottom=H-y1,H-y0
            assert x0>=.35 and x1<=W-.35 and top>=.35 and bottom<=H-.35,('Canvas clipping',a.get_text(),[x0,top,x1,bottom])
            if spec['box']:
                l,t,r,bb=spec['box']
                assert x0>=l+.5 and x1<=r-.5 and top>=t+.35 and bottom<=bb-.35,('Panel clipping',a.get_text(),[x0,top,x1,bottom],spec['box'])
            records.append({'text':a.get_text(),'role':spec['role'],'font_pt':float(a.get_fontsize()),'x0_mm':float(x0),'y0_mm':float(top),'x1_mm':float(x1),'y1_mm':float(bottom)})
        for i,a in enumerate(records):
            for b in records[i+1:]:
                xx=min(a['x1_mm'],b['x1_mm'])-max(a['x0_mm'],b['x0_mm'])
                yy=min(a['y1_mm'],b['y1_mm'])-max(a['y0_mm'],b['y0_mm'])
                if xx>.25 and yy>.25:overlaps.append([a['text'],b['text'],round(xx,2),round(yy,2)])
        assert not overlaps,('Text overlap',self.lang,overlaps)
        return records
    def save(self):
        self.draw();records=self.layout()
        title='Figure 1. Study design and analytical framework'
        self.fig.savefig(BASE/f'Figure_1_{self.lang}.pdf',metadata={'Title':title,'Creator':'Paired ADNI CSF study; original vector illustrations'})
        self.fig.savefig(BASE/f'editable/Figure_1_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_1_{self.lang}.png',dpi=300)
        self.fig.savefig(BASE/f'Figure_1_{self.lang}_preview.png',dpi=180)
        tif=BASE/f'high_resolution/Figure_1_{self.lang}.tiff'
        self.fig.savefig(tif,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        with Image.open(tif) as im:
            rgb=im.convert('RGB');rgb.save(tif,compression='tiff_lzw',dpi=(600,600));rgb.close()
        (BASE/f'audit/LAYOUT_{self.lang}.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
        (BASE/f'audit/PLOT_PRIMITIVES_{self.lang}.json').write_text(json.dumps(self.primitives,ensure_ascii=False,indent=2),encoding='utf-8')
        return {'language':self.lang,'figure_size_mm':[W,H],'minimum_font_pt':min(r['font_pt'] for r in records),'stages':6,
                'text_overlap':False,'text_clipped':False,'text_within_panels':True,'observed_results_plotted':False}
def main():
    warnings.filterwarnings('error',message='Glyph.*missing')
    outputs=[]
    for lang in ['EN','ZH']:
        f=Figure(lang);outputs.append(f.save());plt.close(f.fig)
    (BASE/'audit/LAYOUT_VALIDATION.json').write_text(json.dumps(outputs,ensure_ascii=False,indent=2),encoding='utf-8')
    (BASE/'audit/STYLE_SPECIFICATION.json').write_text(json.dumps({'version':2,'reference':'User-specified JTM Figure_1.pdf',
      'structure':'Six illustrated stages, with pathology and cognition in parallel','original_vector_icons':True,
      'assay_fills':FILLS,'assay_edges':EDGES,'assay_tints':TINTS,'background':'white','width_mm':W,'height_mm':H,
      'result_estimates_or_discovery_counts':False,'external_small_grey_footnotes':False},ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(outputs))
if __name__=='__main__':main()
