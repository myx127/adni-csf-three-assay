from config import *
import pickle

def run():
 frozen=json.loads((RESULTS/'STAGE8_PLAN_FROZEN.json').read_text(encoding='utf-8'));assert frozen['protocol_sha256']==sha(ROOT/'STAGE8_PROTOCOL_v1.0.md')
 designs=pickle.load((S6/'00_restricted/ATTRIBUTION_MODELS_RESTRICTED.pkl').open('rb'))
 # Existing released second moments only; do not calculate new real shared P here.
 cols=['accession','gene','model','N_residual_SD','S_residual_SD','T_residual_SD','c_NS','c_NT','c_ST']
 existing=pd.read_csv(S6/'results/COVARIATE_FOUR_MODEL_RESULTS.csv',usecols=cols);templates={};rows=[]
 for r in existing.loc[existing.model.isin(MODELS)].itertuples():
  d=designs[(r.accession,r.model)];sigma=np.array([r.N_residual_SD,r.S_residual_SD,r.T_residual_SD])
  C=np.diag(sigma**2);C[0,1]=C[1,0]=r.c_NS;C[0,2]=C[2,0]=r.c_NT;C[1,2]=C[2,1]=r.c_ST
  assert np.linalg.eigvalsh(C).min()>0
  key=r.model+'__'+r.accession
  templates[key]={'accession':r.accession,'gene':r.gene,'model':r.model,'d':d,'C':C,'sigma':sigma,'L':np.linalg.cholesky(C)}
  rows.append({'design_id':key,'accession':r.accession,'gene':r.gene,'model':r.model,'in_family34':r.accession in FAMILY,'n':d['n'],'rank_X':d['rank'],'residual_df':d['df'],'N_G':len(d['onehots'][0]),'S_G':len(d['onehots'][1]),'T_G':len(d['onehots'][2]),'min_cov_eigenvalue':np.linalg.eigvalsh(C).min(),'single_factor_not_required':True})
 assert len(templates)==105 and len({t['accession'] for t in templates.values()})==35
 for acc in {t['accession'] for t in templates.values()}:
  base=templates['M0__'+acc]['d']
  for model in MODELS:assert np.array_equal(base['indices'],templates[model+'__'+acc]['d']['indices'])
 path=PRIVATE/'SIMULATION_TEMPLATES_RESTRICTED.pkl'
 if path.exists():
  previous=pickle.load(path.open('rb'));assert list(previous)==list(templates)
 else:
  with path.open('wb') as f:pickle.dump(templates,f,pickle.HIGHEST_PROTOCOL)
 csvout(RESULTS/'TEMPLATE_DESIGN_AUDIT.csv',rows);csvout(RESULTS/'EFFECT_PATTERN_GRID.csv',[{**r,'effect':json.dumps(r['effect'])} for r in PATTERNS])
 log('105 observed-covariance templates prepared; all model comparisons use identical subjects; no latent-factor restriction')

if __name__=='__main__':run()
