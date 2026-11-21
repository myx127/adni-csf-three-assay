"""Nested retention rings, an error-scenario matrix, peptide maps and vertical interval glyphs."""
from pathlib import Path
import json,logging,warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Wedge,Rectangle
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.font_manager import FontProperties
from matplotlib.text import Text
from PIL import Image
BASE=Path(__file__).resolve().parents[1]
D=json.loads((BASE/'audit/PLOT_INPUT.json').read_text(encoding='utf-8'))
W,H=180,250
INK='#293C4B';GRID='#E7ECF0';PEPTIDE='#B46F8D'
FILLS={'NULISA':'#D4977D','SomaScan':'#76ACA0','TMT-MS':'#92A5C3'}
EDGES={'NULISA':'#A3674F','SomaScan':'#3D7F73','TMT-MS':'#5B759C'}
PROTEIN={'NPTX1':'#568F88','KLK6':'#9980AF'}
DOMAIN={'NPTX1':'#BBD9D0','KLK6':'#D5C8E2'}
CMAP=LinearSegmentedColormap.from_list('signed_contrast',['#7194B4','#FCFCFC','#BA829C'])
logging.getLogger('fontTools').setLevel(logging.ERROR)
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','font.family':'Arial',
 'figure.facecolor':'white','savefig.facecolor':'white','axes.unicode_minus':False,'axes.linewidth':.65,
 'xtick.major.size':2.6,'ytick.major.size':2.6,'xtick.major.width':.65,'ytick.major.width':.65})
TEXT={
 'EN':dict(A='Protein variation retained after adjustment',B='One signal, unequal measurement error',
    C='Molecular coverage of targeted peptides',D='Peptide associations in additional ADNI participants',
    outer='Outer: p-tau181 only',inner='Inner: + three background scores',
    hypothetical='Illustrative model with assigned parameters',reference='Reference error variance',
    comparison='Comparison error variance',contrast='Contrast\nchange',
    residue='Amino-acid position',targeted='Targeted mass spectrometry',specific='Five gene-specific peptides',
    PTX='Pentraxin',S1='Peptidase S1',cohort='n=95; no overlap with the three-assay sample',
    unit='Peptide-signal change\nper SD log$_2$ p-tau181'),
 'ZH':dict(A='调整后保留的蛋白变异',B='同一信号与不等的测量误差',
    C='靶向肽段的分子覆盖',D='额外ADNI参与者中的肽段关联',
    outer='外环：仅调整p-tau181',inner='内环：再纳入三个背景分数',
    hypothetical='假设模型；参数未依据ADNI估计',reference='参照检测的误差方差',
    comparison='比较检测的误差方差',contrast='关联对比\n的变化',
    residue='氨基酸位置',targeted='靶向质谱',specific='五条基因特异性肽段',
    PTX='Pentraxin结构域',S1='Peptidase S1结构域',cohort='n=95；与三检测主样本无参与者重叠',
    unit='肽段信号变化\n每SD log$_2$ p-tau181')}
class Figure:
    def __init__(self,lang):
        self.lang=lang;self.t=TEXT[lang];self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=180);self.records={}
    def font(self,size=8,bold=False):
        name=('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN' else ('msyhbd.ttc' if bold else 'msyh.ttc')
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/name),family='Arial' if self.lang=='EN' else 'Microsoft YaHei',size=size,weight='bold' if bold else 'normal')
    def text(self,x,y,s,size=8,bold=False,ha='left',color=INK,rotation=0):
        return self.fig.text(x/W,1-y/H,s,fontproperties=self.font(size,bold),color=color,ha=ha,va='center',linespacing=1.08,rotation=rotation,zorder=12)
    def panel(self,letter,title,x,y):
        self.text(x,y,letter,11.5,True);self.text(x+7,y,title,8.7,True)
    def axis(self,x,y,w,h):
        ax=self.fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
        ax.spines[['top','right']].set_visible(False);ax.spines[['bottom','left']].set_color('#8294A0');ax.tick_params(colors=INK,pad=3,direction='out')
        return ax
    def ticks(self,ax,size=7.4):
        for t in ax.get_xticklabels()+ax.get_yticklabels():t.set_fontproperties(self.font(size));t.set_color(INK)
    def ring(self,x,y,r,width,fraction,color):
        ax=self.fig.add_axes([(x-r)/W,1-(y+r)/H,2*r/W,2*r/H]);ax.set(xlim=(-r,r),ylim=(-r,r));ax.set_aspect('equal');ax.axis('off')
        ax.add_patch(Wedge((0,0),r,0,360,width=width,facecolor=GRID,edgecolor='none'))
        ax.add_patch(Wedge((0,0),r,90-360*fraction,90,width=width,facecolor=color,edgecolor='none'))
    def retention(self):
        self.panel('A',self.t['A'],3,7)
        self.text(23,17,self.t['outer'],7.5);self.text(23,22,self.t['inner'],7.5)
        self.ring(19,17,1.8,.45,.8,'#7F929E');self.ring(19,22,1.2,.4,.4,'#7F929E')
        rec=[]
        for x,a in zip([29,52,75],D['assays']):self.text(x,30,a,7.8,True,'center',EDGES[a])
        for y,g in zip([45,68,91],D['genes']):
            rows=[r for r in D['retention'] if r['Protein']==g]
            self.text(3,y-1.8,g,8.1,True);self.text(3,y+2.8,'n='+str(rows[0]['Participants']),7.1)
            for x,a in zip([29,52,75],D['assays']):
                r=next(r for r in rows if r['Assay']==a)
                self.ring(x,y,7.7,1.35,r['p-tau alone'],FILLS[a])
                self.ring(x,y,5.4,1.25,r['p-tau and three backgrounds'],EDGES[a])
                label=f"{100*r['p-tau and three backgrounds']:.0f}%"
                self.text(x,y,label,7.5,True,'center')
                rec.append({**r,'outer_fraction':r['p-tau alone'],'inner_fraction':r['p-tau and three backgrounds'],'displayed_inner_percentage':label,'outer_angle':360*r['p-tau alone'],'inner_angle':360*r['p-tau and three backgrounds']})
        self.records['A']={'ring_pairs':rec,'scale':[0,1],'full_circle_is_reference_residual_sum_of_squares':True,'not_an_assay_reliability_estimate':True}
    def error_grid(self):
        self.panel('B',self.t['B'],95,7)
        self.text(137,18,self.t['hypothetical'],7.5,False,'center')
        levels=[.1,.5,1,2,4];values=np.empty((5,5));cells=[]
        for r in D['error_grid']:
            i=levels.index(r['Reference-assay error variance']);j=levels.index(r['Comparison-assay error variance'])
            values[i,j]=r['Adjusted minus unadjusted contrast']
            cells.append({**r,'row':i,'column':j,'displayed':f"{r['Adjusted minus unadjusted contrast']:.3f}"})
        ax=self.axis(114,30,57,57);ax.spines[:].set_visible(False)
        norm=Normalize(vmin=-.06,vmax=.06)
        im=ScalarMappable(norm=norm,cmap=CMAP)
        for i in range(5):
            for j in range(5):
                ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor=CMAP(norm(values[i,j])),edgecolor='white',linewidth=.8))
        ax.set(xlim=(-.5,4.5),ylim=(4.5,-.5),aspect='equal')
        ax.set_xticks(range(5),[str(v) for v in levels]);ax.set_yticks(range(5),[str(v) for v in levels])
        ax.set_xticks(np.arange(-.5,5,1),minor=True);ax.set_yticks(np.arange(-.5,5,1),minor=True)
        ax.grid(which='minor',color='white',linewidth=1);ax.tick_params(which='both',length=0);self.ticks(ax,7.5)
        for c in cells:
            label=c['displayed'].replace('-','−')
            if abs(float(c['displayed']))==0:label='0.000'
            ax.text(c['column'],c['row'],label,fontproperties=self.font(7.3),color=INK,ha='center',va='center')
        self.text(101,58,self.t['reference'],7.5,False,'center',rotation=90);self.text(142.5,96,self.t['comparison'],7.5,False,'center')
        cax=self.fig.add_axes([124/W,1-105/H,45/W,2.5/H])
        cb=self.fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-.06,0,.06]);cb.outline.set_visible(False)
        if cb.solids is not None:cb.solids.set_rasterized(False)
        cax.set_xticklabels(['−0.06','0','0.06']);cax.tick_params(length=2,pad=2,width=.6);self.ticks(cax,7.2)
        self.text(113,104,self.t['contrast'],7.1,False,'center')
        self.records['B']={'cells':cells,'reference_levels':levels,'comparison_levels':levels,
            'colour_limits':[-.06,.06],'is_illustrative_population_model':True,'not_fitted_to_ADNI':True,
            'single_underlying_signal_and_equal_loadings':True,'signed_change_not_absolute_change':True}
    def line(self,points,color=INK,lw=.7,style='-'):
        self.fig.lines.append(Line2D([p[0]/W for p in points],[1-p[1]/H for p in points],transform=self.fig.transFigure,color=color,linewidth=lw,linestyle=style,zorder=3))
    def rectangle(self,x,y,w,h,fill,edge='none',lw=.4,z=4):
        self.fig.patches.append(Rectangle((x/W,1-(y+h/2)/H),w/W,h/H,transform=self.fig.transFigure,facecolor=fill,edgecolor=edge,linewidth=lw,zorder=z))
    def molecular_map(self):
        self.panel('C',self.t['C'],3,120)
        sx=lambda pos:31+pos/450*145
        # Both canonical proteins share this residue scale.
        for pos in [0,100,200,300,400]:
            self.text(sx(pos),129,str(pos),7.2,False,'center')
        self.text(175,124.8,self.t['residue'],7.2,False,'right')
        self.line([(sx(0),132),(sx(450),132)],'#AFBBC3',.5)
        for pos in [0,100,200,300,400]:self.line([(sx(pos),131.4),(sx(pos),132.6)],'#AFBBC3',.5)
        rows=[];tracks=[]
        for gene,y in [('NPTX1',143),('KLK6',165)]:
            domain=next(r for r in D['domains'] if r['Protein']==gene);length=domain['Canonical length']
            self.text(3,y-1,gene,8.3,True);self.text(3,y+3.4,domain['UniProt accession'],7.1)
            self.rectangle(sx(0),y,sx(length)-sx(0),2.7,'#E4E9ED')
            start,end=domain['Domain start'],domain['Domain end']
            # Inclusive residue intervals extend half a residue at either side.
            self.rectangle(sx(start-.5),y,sx(end+.5)-sx(start-.5),4.1,DOMAIN[gene])
            self.text(sx(length)+1.4,y,str(length),7.2)
            for r in [r for r in D['peptides'] if r['Protein']==gene]:
                self.rectangle(sx(r['Start']-.5),y,sx(r['End']+.5)-sx(r['Start']-.5),5.0,PEPTIDE,z=6)
                center=sx((r['Start']+r['End'])/2)
                placements={'FQLTFPLR':(center,137),'LENLEQYSR':(center,137),
                    'ESSQEQSSVVR':(53,172),'LSELIQPLPLER':(73,158),'YTNWIQK':(108,158)}
                lx,ly=placements[r['Peptide']]
                edge=y-2.7 if ly<y else y+2.7
                self.line([(center,edge),(lx,ly+(1.8 if ly<y else -1.8))],PEPTIDE,.55)
                self.text(lx,ly,r['Peptide'],7.4,False,'center')
                rows.append({**r,'track_y_mm':y,'segment_start_mm':sx(r['Start']-.5),'segment_end_mm':sx(r['End']+.5)})
            if gene=='NPTX1':self.text((sx(start)+sx(end))/2,150,self.t['PTX'],7.3,False,'center',PROTEIN[gene])
            else:self.text(85,172,self.t['S1'],7.3,False,'center',PROTEIN[gene])
            tracks.append({**domain,'start_mm':sx(0),'end_mm':sx(length)})
        self.text(148,160,self.t['targeted'],7.8,True,'center')
        self.text(148,167,self.t['specific'],7.3,False,'center')
        self.records['C']={'peptides':rows,'proteins':tracks,'common_amino_acid_scale':[0,450],
          'coordinates_one_based_and_inclusive':True,'NPTX1_is_distinct_from_NPTX2':True,'mapping_does_not_identify_proteoforms':True}
    def peptide_associations(self):
        self.panel('D',self.t['D'],3,181)
        self.text(101,188,self.t['cohort'],7.5,False,'center')
        ax=self.axis(27,200,149,32)
        ax.set(xlim=(-.5,4.5),ylim=(-.06,.40),xticks=range(5),yticks=[-.05,0,.2,.4])
        ax.set_yticklabels(['−0.05','0','0.2','0.4'])
        ax.set_xticklabels(D['peptide_order'])
        ax.axvspan(1.5,4.5,color='#F7F5FA',zorder=-2)
        ax.axhline(0,color='#8999A5',linestyle=(0,(3,2)),linewidth=.7,zorder=-1)
        for y in [.2,.4]:ax.axhline(y,color=GRID,linewidth=.5,zorder=-1)
        rows=[]
        for i,r in enumerate(D['peptides']):
            b,lo,hi=r['Coefficient'],r['CI lower'],r['CI upper'];colour=PROTEIN[r['Protein']]
            ax.vlines(i,lo,hi,colors=colour,linewidth=2.0,zorder=3)
            ax.hlines([lo,hi],i-.065,i+.065,colors=colour,linewidth=.9,zorder=3)
            ax.plot(i,b,marker='o',mfc='white',mec=colour,mew=1.4,markersize=5.5,zorder=4)
            px=27+(i+.5)/5*149
            self.text(px,245,'q='+f"{r['q']:.3g}",7.5,False,'center')
            rows.append({**r,'x':i,'y':b,'interval_low':lo,'interval_high':hi})
        self.ticks(ax,7.4);ax.tick_params(axis='x',length=0,pad=8)
        self.text(56.8,195,'NPTX1',8.1,True,'center',PROTEIN['NPTX1'])
        self.text(131.3,195,'KLK6',8.1,True,'center',PROTEIN['KLK6'])
        self.text(8,218,self.t['unit'],7.4,False,'center',rotation=90)
        self.records['D']={'intervals':rows,'family':5,'n':95,'y_limits':[-.06,.4],
            'q_one_sided_positive_tests':True,'CI_two_sided_95_percent':True,'t_degrees_of_freedom':4,
            'not_external_cohort_validation':True,'signal_scale_is_released_normalized_log_peptide_not_protein_SD':True}
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
        self.retention();self.error_grid();self.molecular_map();self.peptide_associations();labels=self.layout()
        title='Figure 5. Measurement context and peptide-level support for assay interpretation'
        self.fig.savefig(BASE/f'Figure_5_{self.lang}.pdf',metadata={'Title':title,'Creator':'Paired ADNI CSF proteomic study'})
        self.fig.savefig(BASE/f'editable/Figure_5_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_5_{self.lang}.png',dpi=300);self.fig.savefig(BASE/f'Figure_5_{self.lang}_preview.png',dpi=180)
        tif=BASE/f'high_resolution/Figure_5_{self.lang}.tiff';self.fig.savefig(tif,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        with Image.open(tif) as im:
            rgb=im.convert('RGB');rgb.save(tif,compression='tiff_lzw',dpi=(600,600));rgb.close()
        (BASE/f'audit/PLOT_PRIMITIVES_{self.lang}.json').write_text(json.dumps(self.records,ensure_ascii=False,indent=2),encoding='utf-8')
        (BASE/f'audit/LAYOUT_{self.lang}.json').write_text(json.dumps(labels,ensure_ascii=False,indent=2),encoding='utf-8')
        return {'language':self.lang,'size_mm':[W,H],'panels':4,'minimum_font_pt':min(r['font_pt'] for r in labels),
            'text_overlap':False,'text_clipped':False,'retention_ring_pairs':9,'hypothetical_cells':25,'mapped_peptides':5,'association_intervals':5}
def main():
    warnings.filterwarnings('error',message='Glyph.*missing');reports=[]
    for lang in ['EN','ZH']:
        f=Figure(lang);reports.append(f.save());plt.close(f.fig)
    (BASE/'audit/LAYOUT_VALIDATION.json').write_text(json.dumps(reports,ensure_ascii=False,indent=2),encoding='utf-8')
    (BASE/'audit/STYLE_SPECIFICATION.json').write_text(json.dumps({'size_mm':[W,H],
       'assay_fills':FILLS,'assay_edges':EDGES,'protein_colours':PROTEIN,'peptide_segment_colour':PEPTIDE,
       'formats':['nested retained-variation rings','signed error-scenario matrix','common-scale protein architecture tracks','vertical peptide association intervals'],
       'simulation_distinguished_from_observed_data':True,'outside_grey_explanatory_footer':False,
       'no_assay_accuracy_ranking':True},ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(reports))
if __name__=='__main__':main()
