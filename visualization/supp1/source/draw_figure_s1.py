"""Draw an editable, data-linked medical study flow in English and Chinese."""
from pathlib import Path
import csv,json,math,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
from matplotlib.font_manager import FontProperties
from PIL import Image

BASE=Path(__file__).resolve().parent.parent
AUDIT=BASE/'audit'
EDITABLE=BASE/'editable'
HIGH=BASE/'high_resolution'
DATA=json.loads((AUDIT/'STUDY_FLOW_VALUES.json').read_text(encoding='utf-8'))
W,H=180,244
MM=1/25.4
P2MM=25.4/72
# User-specified scientific palette. Text remains dark for comfortable reading.
PALETTE=['#527F94','#9DBCC7','#D7E5E6','#C289A0','#E2B9C7','#A5B59A']
C=dict(ink='#304754',muted='#304754',rule='#D7E5E6',blue='#527F94',blue_fill='#D7E5E6',
       teal='#C289A0',teal_fill='#F5E9EF',grey='#304754',grey_fill='#F5F8F8',grey_edge='#9DBCC7',
       purple='#A5B59A',purple_fill='#EAF0E6',nulisa='#527F94',soma='#C289A0',tmt='#A5B59A')
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
                     'savefig.facecolor':'white','figure.facecolor':'white','axes.unicode_minus':False})

LABELS={
'EN':dict(
 figure='Supplementary Figure S1',
 A='Paired CSF measurements\nand participant selection',B='Protein mapping\nand selection',C='Pathology and cognitive association analyses',
 assay1='Antibody assay',assay2='Aptamer assay',assay3='Mass spectrometry',
 paired_title='Paired CSF measurements',paired_detail='Matched by participant and date',
 ptau_title='Same-day p-tau181',ptau_detail='Elecsys measurement available',
 path_title='Clinical information available',path_detail='p-tau181, age and sex recorded',
 cog_title='Memory follow-up available',cog_detail='Baseline and ≥1 follow-up\nwithin five years',
 excluded='excluded',ptau_missing='Same-day result\nunavailable',age_missing='Age unavailable',memory_missing='No eligible\nmemory data',
 nominal_title='Candidate protein targets',nominal_detail='Measured by all three assays',
 single_title='Proteins with a unique match',single_detail='One matched SomaScan target',
 eligible_title='Proteins meeting study criteria',eligible_detail='Adequate detection and data',
 formal_title='Proteins in the main analyses',formal_detail='Pathology and cognitive associations',
 mapping_missing='26 with multiple\nSomaScan matches\n1 combined protein\nannotation',
 design_missing='Low detection or\ninsufficient data',mapt='Analysed separately\n1 protein',
 formal_parent='All 34 proteins included in both analyses',selected_parent='Exploratory analyses',
 pathology='Shared p-tau181 associations',pathology_test='Associations in the same\ndirection across all three assays',
 interpretation='Highlighted for clinical\ninterpretation',
 cognition='Cognitive associations',primary='Primary:\nAnnual change in ADNI-MEM',
 secondary='Secondary:\nBaseline ADNI-MEM and\nannual ADNI-EF change',
 cog_tests='Shared associations and\ndifferences between assays',
 adjustment='p-tau181 adjustment',selected='Selected because memory\nassociations differed by assay',
 models='Nine adjustment models',paired_changes='Paired comparison of\ndifferences between assays\nbefore and after adjustment'),
'ZH':dict(
 figure='补充图 S1',A='配对CSF测量\n与参与者筛选',B='蛋白映射\n与筛选',C='病理与认知关联的分析框架',
 assay1='抗体检测',assay2='适配体检测',assay3='质谱检测',
 paired_title='配对CSF测量',paired_detail='按参与者与采样日期匹配',
 ptau_title='同日p-tau181检测',ptau_detail='具有可用Elecsys检测结果',
 path_title='临床资料完整的参与者',path_detail='具有p-tau181、年龄及性别资料',
 cog_title='具有记忆随访的参与者',cog_detail='基线及至少1次随访\nCSF采集后五年内',
 excluded='排除',ptau_missing='无合格同日\np-tau181',age_missing='年龄不可用',memory_missing='无合格记忆\n评估资料',
 nominal_title='三检测共同蛋白靶标',nominal_detail='依据检测靶标注释确定',
 single_title='唯一匹配的蛋白',single_detail='SomaScan靶标唯一对应',
 eligible_title='满足分析要求的蛋白',eligible_detail='检出率与样本条件满足',
 formal_title='纳入主要分析的蛋白',formal_detail='病理与认知关联分析',
 mapping_missing='26个SomaScan多重匹配\n1个复合注释',
 design_missing='检出率偏低或\n资料不足',mapt='单独分析\n1个蛋白',
 formal_parent='全部34个蛋白均纳入以下两类分析',selected_parent='探索性分析',
 pathology='共同p-tau181关联',pathology_test='三种检测中\n关联方向一致',interpretation='临床结果解释中\n重点关注的蛋白',
 cognition='认知关联',primary='主要结局：\nADNI-MEM年变化',secondary='次要结局：\n基线ADNI-MEM及\nADNI-EF年变化',
 cog_tests='共同关联与\n跨检测关联差异',
 adjustment='p-tau181调整',selected='因记忆关联在不同检测间\n存在差异而选定',models='9种调整模型',
 paired_changes='配对比较调整前后\n检测间关联差异的变化')}

class Diagram:
    def __init__(self,lang):
        self.lang=lang;self.labels=LABELS[lang]
        self.fig=plt.figure(figsize=(W*MM,H*MM),dpi=150)
        self.ax=self.fig.add_axes([0,0,1,1]);self.ax.set(xlim=(0,W),ylim=(H,0));self.ax.axis('off')
        self.fig.canvas.draw();self.renderer=self.fig.canvas.get_renderer();self.texts=[];self.boxes=[]
    def font(self,size,bold=False):
        name=('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN' else ('msyhbd.ttc' if bold else 'msyh.ttc')
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/name),size=size,weight='bold' if bold else 'normal')
    def width(self,text,size,bold=False):
        return self.renderer.get_text_width_height_descent(text,self.font(size,bold),False)[0]/self.fig.dpi*25.4
    def wrap(self,text,width,size,bold=False):
        out=[]
        for raw in str(text).split('\n'):
            if not raw:out.append('');continue
            tokens=raw.split(' ') if self.lang=='EN' else list(raw)
            join=' ' if self.lang=='EN' else ''
            line=''
            for token in tokens:
                trial=token if not line else line+join+token
                if line and self.width(trial,size,bold)>width:
                    out.append(line);line=token
                else:line=trial
            out.append(line)
        return out
    def text(self,x,y,text,size=9,bold=False,color=None,ha='left',width=None,region=None,tag=''):
        lines=self.wrap(text,width,size,bold) if width else str(text).split('\n')
        gap=size*P2MM*1.16
        for i,line in enumerate(lines):
            obj=self.ax.text(x,y+i*gap,line,fontproperties=self.font(size,bold),ha=ha,va='center',
                             color=color or C['ink'],clip_on=False)
            self.texts.append(dict(obj=obj,tag=tag,region=region,text=line,point_size=size))
        return len(lines)*gap
    def rect(self,x,y,w,h,edge,fill='white',lw=.8,rounding=.6,tag=''):
        box=FancyBboxPatch((x,y),w,h,boxstyle=f'round,pad=0,rounding_size={rounding}',
                          linewidth=lw,edgecolor=edge,facecolor=fill,zorder=1)
        self.ax.add_patch(box);self.boxes.append(dict(tag=tag,x=x,y=y,w=w,h=h))
    def node(self,x,y,w,h,parts,edge,fill='white',tag='',pad=2):
        self.rect(x,y,w,h,edge,fill,tag=tag)
        prepared=[]
        for text,size,bold,color in parts:
            lines=self.wrap(text,w-2*pad,size,bold)
            prepared.append((lines,size,bold,color,len(lines)*size*P2MM*1.16))
        total=sum(r[-1] for r in prepared)+.65*(len(parts)-1)
        assert total<=h-1.6,(self.lang,tag,'text block too tall',total,h,[(r[0],r[1]) for r in prepared])
        pos=y+(h-total)/2
        for lines,size,bold,color,hh in prepared:
            self.text(x+w/2,pos+size*P2MM*1.16/2,'\n'.join(lines),size,bold,color,ha='center',
                      region=(x+.5,y+.5,x+w-.5,y+h-.5),tag=tag)
            pos+=hh+.65
    def arrow(self,start,end,color=None,lw=.8):
        self.ax.add_patch(FancyArrowPatch(start,end,arrowstyle='-|>',mutation_scale=6.6,
                            linewidth=lw,color=color or C['gray'],shrinkA=1,shrinkB=1,zorder=0))
    def line(self,points,color=None,lw=.65):
        self.ax.plot([p[0] for p in points],[p[1] for p in points],color=color or C['grey_edge'],linewidth=lw,zorder=0)
    def header(self,letter,title,x,y,width):
        badge={'A':'#9DBCC7','B':'#E2B9C7','C':'#A5B59A'}[letter]
        self.rect(x,y-2.5,5,5,badge,badge,lw=0,rounding=.25,tag='panel badge '+letter)
        self.text(x+2.5,y,letter,11.5,True,ha='center',tag='panel '+letter)
        self.text(x+7,y,title,10.3,True,width=width-7,tag='panel '+letter+' title')
    def exclusion(self,x,y,w,h,count,description,tag,reference=False):
        title='MAPT' if reference else (f'{count} excluded' if self.lang=='EN' else f'排除{count}'+('人' if tag.startswith('A') else '个'))
        self.node(x,y,w,h,[(title,9.3,True,C['ink']),
                           (description,8.8,False,C['ink'])],
                  C['teal'] if reference else C['grey_edge'],C['teal_fill'] if reference else C['grey_fill'],tag,pad=1.5)
    def draw(self):
        t=self.labels;v=DATA
        self.text(2,6.5,t['figure'],13,True,tag='figure label')
        self.header('A',t['A'],2,17.5,84)
        self.header('B',t['B'],95,17.5,83)
        self.line([(2,27),(86,27)],C['rule']);self.line([(95,27),(178,27)],C['rule'])
        # A: three released proteomic datasets converge on a matched index visit.
        for x,name,detail,color in [(2,'NULISA',t['assay1'],C['nulisa']),(31,'SomaScan',t['assay2'],C['soma']),(60,'TMT-MS',t['assay3'],C['tmt'])]:
            self.node(x,32,25,13.5,[(name,9.5,True,C['ink']),(detail,8.8,False,C['ink'])],color,tag='A assay '+name,pad=1)
            self.line([(x+12.5,45.5),(x+12.5,47.5)],C['grey_edge'])
        self.line([(14.5,47.5),(72.5,47.5)],C['grey_edge'])
        self.arrow((28.5,47.5),(28.5,52),C['blue'])
        patient_specs=[(52,18,'paired_title','paired_participants','paired_detail'),
                       (80,18,'ptau_title','valid_same_day_ptau','ptau_detail'),
                       (108,19,'path_title','pathology_reference','path_detail'),
                       (137,20.5,'cog_title','cognitive_reference','cog_detail')]
        for y,h,title,quantity,detail in patient_specs:
            count=f'{v[quantity]} participants' if self.lang=='EN' else f'{v[quantity]}人'
            self.node(2,y,53,h,[(t[title],9.3,True,C['ink']),(count,12.2,True,C['ink']),
                              (t[detail],8.8,False,C['ink'])],C['blue'],C['blue_fill'],tag='A '+quantity)
        for start,end,count,desc,h in [(70,80,v['no_valid_same_day_ptau'],'ptau_missing',17),
                                      (98,108,v['age_unavailable'],'age_missing',16),
                                      (127,137,v['no_eligible_memory'],'memory_missing',20)]:
            mid=(start+end)/2
            self.arrow((28.5,start),(28.5,end),C['blue'])
            self.arrow((28.5,mid),(60,mid),C['grey_edge'],.6)
            self.exclusion(60,mid-h/2,25,h,count,t[desc],'A '+desc)
        # B: candidate identities, measurement eligibility and the separate MAPT reference.
        protein_specs=[(32,18,'nominal_title','nominal_targets','nominal_detail'),
                       (72,20,'single_title','single_mapping','single_detail'),
                       (107,20,'eligible_title','eligible_proteins','eligible_detail'),
                       (137,20,'formal_title','formal_proteins','formal_detail')]
        for y,h,title,quantity,detail in protein_specs:
            suffix=' targets' if quantity=='nominal_targets' else ' proteins'
            count=str(v[quantity])+suffix if self.lang=='EN' else f'{v[quantity]}个'
            self.node(95,y,53,h,[(t[title],9.3,True,C['ink']),(count,12.2,True,C['ink']),
                               (t[detail],8.8,False,C['ink'])],C['teal'],C['teal_fill'],tag='B '+quantity)
        for start,end,count,desc,h,ref in [(50,72,v['mapping_exclusions'],'mapping_missing',29,False),
                                          (92,107,v['design_exclusions'],'design_missing',18,False),
                                          (127,137,1,'mapt',17,True)]:
            mid=(start+end)/2
            self.arrow((121.5,start),(121.5,end),C['teal'])
            self.arrow((121.5,mid),(153,mid),C['teal'] if ref else C['grey_edge'],.6)
            self.exclusion(153,mid-h/2,25,h,count,t[desc],'B '+desc,ref)
        # C: the full family feeds both formal branches; candidate selection follows cognition.
        self.line([(2,165),(178,165)],C['rule'])
        self.header('C',t['C'],2,171,176)
        self.node(2,180,115,7.5,[(t['formal_parent'],9.6,True,C['ink'])],C['blue'],'#9DBCC7',tag='C full family',pad=1)
        self.node(124,180,54,7.5,[(t['selected_parent'],9.6,True,C['ink'])],C['purple'],'#A5B59A',tag='C exploratory label',pad=1)
        for cx in (29,90):self.arrow((cx,187.5),(cx,195),C['blue'])
        count=f"{v['shared_ptau_proteins']} proteins" if self.lang=='EN' else f"{v['shared_ptau_proteins']}个蛋白"
        self.node(2,195,54,41.5,[(t['pathology'],9.5,True,C['ink']),(t['pathology_test'],9.0,False,C['ink']),
                   (count,11.5,True,C['ink']),(t['interpretation'],9.0,False,C['ink'])],C['blue'],tag='C shared pathology')
        self.node(63,195,54,41.5,[(t['cognition'],9.5,True,C['ink']),(t['primary'],9.0,False,C['ink']),
                   (t['secondary'],9.0,False,C['ink']),(t['cog_tests'],9.0,False,C['ink'])],C['blue'],tag='C full family cognition')
        self.node(124,195,54,41.5,[(t['adjustment'],9.5,True,C['ink']),('NPTX2 · TAFA5 · TIMP3',9.1,True,C['ink']),
                   (t['selected'],9.0,False,C['ink']),(t['models'],9.0,False,C['ink']),
                   (t['paired_changes'],9.0,False,C['ink'])],C['purple'],tag='C selected candidates')
        self.arrow((117,215.75),(124,215.75),C['purple'])
        return self
    def check_layout(self):
        self.fig.canvas.draw();renderer=self.fig.canvas.get_renderer();records=[]
        for entry in self.texts:
            box=entry['obj'].get_window_extent(renderer)
            xy=self.ax.transData.inverted().transform([[box.x0,box.y0],[box.x1,box.y1]])
            x0,x1=sorted(xy[:,0]);y0,y1=sorted(xy[:,1])
            assert x0>=-.15 and x1<=W+.15 and y0>=-.15 and y1<=H+.15,(self.lang,entry['tag'],'outside canvas',entry['text'],(x0,y0,x1,y1))
            if entry['region']:
                a,b,c,d=entry['region']
                assert x0>=a-.25 and x1<=c+.25 and y0>=b-.4 and y1<=d+.4,(self.lang,entry['tag'],'outside text container',entry['text'],(x0,y0,x1,y1),entry['region'])
            records.append(dict(tag=entry['tag'],text=entry['text'],size_pt=entry['point_size'],region=entry['region'],x0=float(x0),y0=float(y0),x1=float(x1),y1=float(y1)))
        overlaps=[]
        for i,a in enumerate(records):
            for b in records[i+1:]:
                x=min(a['x1'],b['x1'])-max(a['x0'],b['x0']);y=min(a['y1'],b['y1'])-max(a['y0'],b['y0'])
                if x>.15 and y>.15:overlaps.append((a['tag'],a['text'],b['tag'],b['text'],x,y))
        assert not overlaps,(self.lang,'text overlaps',overlaps)
        return records
    def save(self):
        records=self.check_layout()
        EDITABLE.mkdir(exist_ok=True);HIGH.mkdir(exist_ok=True)
        self.fig.savefig(BASE/f'Figure_S1_{self.lang}.pdf',format='pdf',metadata={'Title':'Supplementary Figure S1: Study design, participant selection, and analytical framework','Creator':'Study figure authoring'})
        self.fig.savefig(EDITABLE/f'Figure_S1_{self.lang}.svg',format='svg',metadata={'Title':'Study design, participant selection, and analytical framework'})
        self.fig.savefig(BASE/f'Figure_S1_{self.lang}.png',dpi=300)
        target=HIGH/f'Figure_S1_{self.lang}.tiff'
        self.fig.savefig(target,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        with Image.open(target) as im:
            rgb=im.convert('RGB');rgb.save(target,compression='tiff_lzw',dpi=(600,600));rgb.close()
        (AUDIT/f'LAYOUT_{self.lang}.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
        return dict(language=self.lang,width_mm=W,height_mm=H,minimum_text_pt=min(r['size_pt'] for r in records),
                    checked_text_lines=len(records),text_overlap=False,text_clipped=False,svg_text_editable=True,tiff_dpi=600)

def main():
    results=[]
    for lang in ('EN','ZH'):
        drawing=Diagram(lang).draw();results.append(drawing.save());plt.close(drawing.fig)
    (AUDIT/'FIGURE_LAYOUT_VALIDATION.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(results,ensure_ascii=False),flush=True)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');main()
