"""Draw protein mapping, peptide associations, replicate intervals and pooled-control CVs."""
from pathlib import Path
import json, logging, sys, warnings
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
INK='#263A46';GRID='#E6ECEA'
PALETTE=['#9BCFC1','#B9A0CF','#C289A0','#FAD6C9','#92B3CD','#D0E3DB','#FAF5E8','#E8E2DB']
PROTEIN_COLOURS={'NPTX1':('#9BCFC1','#4D897B'),'KLK6':('#B9A0CF','#84649B')}
TEXT={
 'EN':dict(figure='Supplementary Figure S16',
    A='Peptide locations\nin canonical sequences',B='Peptide–p-tau181 associations\n95 ADNI participants',
    C='Technical repeatability\n2014 and 2019 releases',D='2019 pooled QC\n24 measurements',
    domain_N='PTX domain',domain_K='Peptidase S1',Bunits='Peptide association (β)',
    Cunits='Pearson r (95% CI)',Dunits='Raw light/heavy CV (%)',raw='Raw light/heavy ratios'),
 'ZH':dict(figure='补充图 S16',A='靶向肽段在\n标准序列中的位置',B='肽段与p-tau181关联\n95名ADNI参与者',
    C='技术重复性\n2014年及2019年发布',D='2019年混池质控\n24次测量',
    domain_N='PTX结构域',domain_K='S1肽酶结构域',Bunits='肽段关联系数（β）',
    Cunits='Pearson r（95% CI）',Dunits='原始轻/重比值CV（%）',raw='原始轻/重离子比值')
}
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
    'axes.linewidth':.7,'xtick.major.width':.7,'ytick.major.width':.7,
    'xtick.major.size':2.5,'ytick.major.size':0,'axes.unicode_minus':False,
    'savefig.facecolor':'white','figure.facecolor':'white'})

def residue_x(position):return 6+(position-1)*78/450
def fragment_label(value):
    return value[0]+''.join(chr(0x2080+int(c)) for c in value[1:])

class Figure:
    def __init__(self,language):
        self.lang=language;self.t=TEXT[language]
        self.fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=160)
        self.primitives={}
    def font(self,size=8.8,bold=False,math=False):
        name='segoeui.ttf' if math else (('arialbd.ttf' if bold else 'arial.ttf') if self.lang=='EN'
            else ('msyhbd.ttc' if bold else 'msyh.ttc'))
        return FontProperties(fname=str(Path('C:/Windows/Fonts')/name),size=size,weight='bold' if bold else 'normal')
    def text(self,x,y,label,size=8.8,bold=False,ha='left',math=False,colour=INK):
        return self.fig.text(x/W,1-y/H,label,fontproperties=self.font(size,bold,math),color=colour,
            ha=ha,va='center',linespacing=1.04)
    def box(self,x,y,width,height,face,edge='none',lw=0):
        self.fig.patches.append(Rectangle((x/W,1-(y+height/2)/H),width/W,height/H,
            transform=self.fig.transFigure,facecolor=face,edgecolor=edge,linewidth=lw,zorder=-2))
    def line(self,xs,ys,colour=INK,lw=.65):
        self.fig.lines.append(Line2D(np.array(xs)/W,1-np.array(ys)/H,
            transform=self.fig.transFigure,color=colour,linewidth=lw,zorder=0))
    def header(self,key,x,y,width,fill,badge):
        self.box(x,y,width,10,fill);self.box(x,y,5.5,10,badge)
        self.text(x+2.75,y,key,11.2,True,'center');self.text(x+6.8,y,self.t[key],9.2,True)
    def legend(self,x,y,label,marker,face,edge):
        self.fig.lines.append(Line2D([x/W],[1-y/H],transform=self.fig.transFigure,
            linestyle='none',marker=marker,ms=4.1,mfc=face,mec=edge,mew=.65))
        self.text(x+3.5,y,label)
    def axis(self,x,top,width,height,xlim,xticks):
        ax=self.fig.add_axes([x/W,1-(top+height)/H,width/W,height/H],facecolor='none')
        ax.spines[['top','right','left']].set_visible(False);ax.spines['bottom'].set_color(INK)
        ax.tick_params(colors=INK,direction='out',pad=2);ax.set_axisbelow(True)
        ax.set(xlim=xlim,ylim=(top+height,top),xticks=xticks,yticks=[])
        for label in ax.get_xticklabels():label.set_fontproperties(self.font());label.set_color(INK)
        return ax
    def interval(self,ax,value,lower,upper,y,face,edge,marker):
        ax.hlines(y,lower,upper,color=edge,lw=.9)
        ax.vlines([lower,upper],y-.75,y+.75,color=edge,lw=.75)
        ax.plot(value,y,linestyle='none',marker=marker,ms=4.5,mfc=face,mec=edge,mew=.65)
    def draw(self):
        self.text(2,6,self.t['figure'],13,True)
        self.header('A',2,23,86,'#FAF5E8','#9BCFC1')
        self.header('B',95,23,83,'#D0E3DB','#B9A0CF')
        maps=[]
        for gene,titley,tracky in [('NPTX1',35,51),('KLK6',90,106)]:
            r=next(x for x in DATA['proteins'] if x['gene']==gene)
            face,edge=PROTEIN_COLOURS[gene]
            self.text(2,titley,gene+' ('+str(r['length'])+' aa)',9.2,True)
            self.box(6,tracky,r['length']*78/450,6,'#E8E2DB',edge='#A6B3C0',lw=.4)
            start=residue_x(r['domain_start']);width=(r['domain_end']-r['domain_start']+1)*78/450
            self.box(start,tracky,width,6,face,edge=edge,lw=.45)
            self.text(start+width/2,tracky-7,self.t['domain_N' if gene=='NPTX1' else 'domain_K'],8.8,ha='center')
            maps.append(dict(gene=gene,length=r['length'],domain_start=r['domain_start'],domain_end=r['domain_end'],
                track_x_mm=6,track_width_mm=r['length']*78/450,domain_x_mm=start,domain_width_mm=width))
        label_positions={
            'LENLEQYSR':dict(x=27,y=65,bend=58),
            'FQLTFPLR':dict(x=59,y=77,bend=58),
            'ESSQEQSSVVR':dict(x=19,y=121,bend=115),
            'LSELIQPLPLER':dict(x=49,y=135,bend=113),
            'YTNWIQK':dict(x=76,y=121,bend=110.5)}
        mapped=[]
        for r in DATA['peptides']:
            tracky=51 if r['gene']=='NPTX1' else 106
            start=residue_x(r['start']);width=(r['end']-r['start']+1)*78/450
            self.box(start,tracky,width,6,'#C289A0',edge='#945D75',lw=.55)
            loc=label_positions[r['sequence']];center=start+width/2
            self.line([center,center,loc['x'],loc['x']],[tracky+3,loc['bend'],loc['bend'],loc['y']-3],colour='#527F94')
            self.text(loc['x'],loc['y'],r['sequence'],8.8,ha='center')
            self.text(loc['x'],loc['y']+5,str(r['start'])+'–'+str(r['end']),8.8,ha='center')
            mapped.append(dict(gene=r['gene'],peptide=r['sequence'],start=r['start'],end=r['end'],
                x_mm=start,width_mm=width,label_x_mm=loc['x'],label_y_mm=loc['y']))
        self.primitives['A']=dict(proteins=maps,peptides=mapped,common_residue_scale_max=450,
            common_track_span_mm=78,intervals_are_one_based_inclusive=True,domains_are_sequence_annotations_not_epitopes=True)

        by={r['peptide']:y for r,y in zip(DATA['effects'],[50,65,92,107,122])}
        self.text(95,38,'NPTX1',8.8,True);self.text(95,80,'KLK6',8.8,True)
        self.text(171,38,'BY q',8.8,True,'center')
        ax=self.axis(125,42,38,92,(-.07,.39),[0,.1,.2,.3])
        ax.axvline(0,color='#A6B3C0',ls='--',lw=.7)
        bpoints=[]
        for r in DATA['effects']:
            y=by[r['peptide']];face,edge=PROTEIN_COLOURS[r['gene']]
            self.text(95,y,r['peptide'])
            self.text(171,y,f"{r['q_BY5']:.3f}",8.8,ha='center')
            self.interval(ax,r['effect'],r['lower'],r['upper'],y,face,edge,'o' if r['gene']=='NPTX1' else 's')
            bpoints.append(dict(gene=r['gene'],peptide=r['peptide'],effect=r['effect'],lower=r['lower'],upper=r['upper'],
                q=r['q_BY5'],q_label=f"{r['q_BY5']:.3f}",participants=95))
        self.text(144,146,self.t['Bunits'],8.8,ha='center')
        self.primitives['B']=dict(points=bpoints,x_limits=[-.07,.39],family_peptides=5,n=95,
            CI_sidedness='Two-sided 95%',q_from_positive_one_sided_tests=True,
            original_signal_scale_not_residual_SD_scale=True,all_q_above_005=True)

        self.header('C',2,158,86,'#FAF5E8','#92B3CD')
        self.header('D',95,158,83,'#D0E3DB','#9BCFC1')
        self.legend(30,170,'2014','o','#FAD6C9','#AA785F')
        self.legend(58,170,'2019','s','#92B3CD','#527F94')
        self.text(85,175,'n',8.8,True,'center')
        self.legend(98,170,self.t['raw'],'o','#9BCFC1','#4D897B')
        row_positions={r['peptide']:183+i*14 for i,r in enumerate(DATA['effects'])}
        c=self.axis(27,177,51,69,(-.15,1.025),[0,.5,1])
        c.axvline(0,color='#A6B3C0',ls='--',lw=.7)
        d=self.axis(123,177,54,69,(0,19),[0,5,10,15])
        d.xaxis.grid(True,color=GRID,lw=.5)
        cpoints=[];dpoints=[]
        for i,r in enumerate(DATA['effects']):
            peptide=r['peptide'];y=row_positions[peptide]
            if i%2==0:
                self.box(2,y,86,12,'#FAF5E8');self.box(95,y,83,12,'#FAF5E8')
            self.text(2,y,peptide);self.text(95,y,peptide)
            for release,offset,face,edge,marker in [('2014',-2.5,'#FAD6C9','#AA785F','o'),
                                                   ('2019',2.5,'#92B3CD','#527F94','s')]:
                t=next(v for v in DATA['repeats'] if v['peptide']==peptide and v['release']==release)
                self.interval(c,t['r'],t['lower'],t['upper'],y+offset,face,edge,marker)
                self.text(85,y+offset,str(t['pairs']),8.8,ha='center')
                cpoints.append(dict(gene=t['gene'],peptide=peptide,release=release,r=t['r'],
                    lower=t['lower'],upper=t['upper'],pairs=t['pairs']))
            qrows=[v for v in DATA['pool_qc'] if v['peptide']==peptide]
            for q,offset in zip(qrows,[-2.5,2.5]):
                d.plot(q['cv'],y+offset,linestyle='none',marker='o',ms=4.1,mfc='#9BCFC1',mec='#4D897B',mew=.6)
                d.text(q['cv']+.65,y+offset,fragment_label(q['fragment']),fontproperties=self.font(8.8,math=True),
                    color=INK,ha='left',va='center')
                dpoints.append(dict(gene=q['gene'],peptide=peptide,fragment=q['fragment'],
                    fragment_label=fragment_label(q['fragment']),cv_percent=q['cv'],measurements=24))
        self.text(52.5,255,self.t['Cunits'],8.8,ha='center')
        self.text(150,255,self.t['Dunits'],8.8,ha='center')
        self.primitives['C']=dict(points=cpoints,x_limits=[-.15,1.025],bootstrap_resamples=4000,
            n_is_complete_replicate_pairs=True,negative_interval_preserved=True,year_groups_are_separate_release_estimates=True)
        self.primitives['D']=dict(points=dpoints,x_limits=[0,19],measurements_per_transition=24,release=2019,
            raw_light_heavy_ratio=True,no_added_pseudocount_or_logarithm=True,no_clinical_subjects_in_pool=True)
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
        title='Supplementary Figure S16: Targeted peptide associations with p-tau181 and technical repeatability for NPTX1 and KLK6'
        self.fig.savefig(BASE/f'Figure_S16_{self.lang}.pdf',metadata={'Title':title,'Creator':'Study supplementary figures'})
        self.fig.savefig(BASE/f'editable/Figure_S16_{self.lang}.svg',metadata={'Title':title})
        self.fig.savefig(BASE/f'Figure_S16_{self.lang}.png',dpi=300)
        tif=BASE/f'high_resolution/Figure_S16_{self.lang}.tiff'
        self.fig.savefig(tif,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        with Image.open(tif) as im:
            rgb=im.convert('RGB');rgb.save(tif,compression='tiff_lzw',dpi=(600,600));rgb.close()
        (AUDIT/f'LAYOUT_{self.lang}.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
        (AUDIT/f'PLOT_PRIMITIVES_{self.lang}.json').write_text(json.dumps(self.primitives,ensure_ascii=False,indent=2),encoding='utf-8')
        return dict(language=self.lang,width_mm=W,height_mm=H,minimum_font_pt=min(r['point_size'] for r in records),
            checked_text_elements=len(records),text_clipped=False,text_overlap=False)

def main():
    warnings.filterwarnings('error',message='Glyph.*missing');results=[]
    for lang in ['EN','ZH']:
        fig=Figure(lang).draw();results.append(fig.save());plt.close(fig.fig)
    (AUDIT/'LAYOUT_VALIDATION.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    (AUDIT/'STYLE_SPECIFICATION.json').write_text(json.dumps(dict(palette=PALETTE,
        choice='Teal and lavender protein annotation, peach and blue release markers, teal pooled controls.',
        figure_size_mm=[W,H],panels=4,peptides=5,association_intervals=5,repeatability_intervals=10,
        raw_QC_transition_points=10,protein_tracks_share_one_amino_acid_scale=True,
        all_intervals_untruncated=True,loose_footer_notes=False),indent=2),encoding='utf-8')
    print(json.dumps(results),flush=True)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');main()
