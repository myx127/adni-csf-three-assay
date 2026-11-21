from __future__ import annotations
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from audit_common import *
from run_coverage import match_records

def design(part):
 if len(part)==0:return {'technical_design_rank':0,'n_minus_technical_rank':0}
 blocks=[np.ones((len(part),1))]
 cols=['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex']
 for c in cols:blocks.append(pd.get_dummies(part[c].astype(str),drop_first=True,dtype=float).to_numpy())
 x=np.column_stack(blocks);rank=int(np.linalg.matrix_rank(x))
 out={'technical_design_rank':rank,'n_minus_technical_rank':len(part)-rank,
      'technical_nominal_columns':x.shape[1],'not_a_final_covariate_design':True}
 for col,prefix in zip(cols,['N','S','T']):
  n=part[col].value_counts().to_numpy()
  out.update({prefix+'_clusters':len(n),prefix+'_cluster_min_n':int(n.min()),prefix+'_cluster_max_n':int(n.max()),
      prefix+'_singleton_clusters':int((n==1).sum()),prefix+'_effective_clusters':float(n.sum()**2/(n*n).sum())})
 return out

def run():
 clinical=pd.read_csv(RESTRICTED/'INDEX_CLINICAL_COVERAGE_FLAGS_RESTRICTED.csv').set_index('key')
 masks=pd.read_parquet(STAGE1/'00_restricted/TRIPLE_CANDIDATE_AVAILABILITY_RESTRICTED.parquet')
 meta=pd.read_csv(STAGE1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
 routes={'GIS_NUMERIC':'complete_GIS_S_all','GIS_LOD':'complete_GIS_detect_all_N_S_all',
         'RAW_NUMERIC':'complete_raw_S_all','GIS_UNAMBIGUOUS':'complete_GIS_unambiguous_S_all'}
 rows_out=[];tech=[];at=[]
 for accession,part in masks.groupby('accession'):
  part=part.set_index('key').join(clinical,how='left',validate='one_to_one')
  assert len(part)==571 and part.RID.notna().all()
  for name,col in routes.items():
   available=part[col];p=available & part.ELECSYS_PTAU_W0
   pmin=p & part.COV_MINIMAL;pfull=p & part.COV_COGNITIVE
   out={'accession':accession,'gene':meta.loc[accession,'nulisa_target'],'annotation_single_feature':bool(meta.loc[accession,'annotation_strict_candidate']),
    'route':name,'n_three_assays':int(available.sum()),'n_same_day_ptau':int(p.sum()),
    'n_same_day_ptau_age_sex':int(pmin.sum()),'n_same_day_ptau_fullcov':int(pfull.sum()),
    'n_p_fullcov_cog_minimum':int((pfull & part.COG_MINIMUM_5Y).sum()),
    'n_p_fullcov_cog_highinfo':int((pfull & part.COG_HIGH_INFO_5Y).sum()),
    'n_amy_pet_90':int((available & part.AMY_PET_CENTILOIDS_W90).sum()),
    'n_amy_pet_365':int((available & part.AMY_PET_CENTILOIDS_W365).sum()),
    'n_p_fullcov_amy_pet_90':int((pfull & part.AMY_PET_CENTILOIDS_W90).sum()),
    'not_final_eligibility':True}
   rows_out.append(out)
   if name=='GIS_NUMERIC':
    subset=part.loc[pmin]
    tech.append({'accession':accession,'gene':out['gene'],'annotation_single_feature':out['annotation_single_feature'],
                 'n':len(subset),**design(subset)})
    for cat,n in subset.PTAU_ANCHOR_AT_W0.value_counts().items():
     at.append({'accession':accession,'gene':out['gene'],'category':cat,'n':int(n),'operational_label_not_tauPET':True})
 data=pd.DataFrame(rows_out);write_csv(RESULTS/'PER_PROTEIN_JOINT_COVERAGE.csv',data)
 write_csv(RESULTS/'PER_PROTEIN_TECHNICAL_DESIGN.csv',tech);write_csv(RESULTS/'PER_PROTEIN_OPERATIONAL_AT_COUNTS.csv',at)
 summaries=[]
 for group,groupdata in [('ALL70',data),('SINGLE_ANNOTATION43',data.loc[data.annotation_single_feature])]:
  for route,part in groupdata.groupby('route'):
   for col in [c for c in data if c.startswith('n_')]:
    values=part[col]
    summaries.append({'group':group,'route':route,'metric':col,'n_proteins':len(part),'min':int(values.min()),
       'median':float(values.median()),'max':int(values.max()),**{f'n_proteins_ge_{n}':int((values>=n).sum()) for n in [100,200,300,400,457]},
       'thresholds_descriptive_only':True})
 write_csv(RESULTS/'PER_PROTEIN_COVERAGE_SUMMARY.csv',summaries)
 # Explain absent contemporaneous tau PET without changing its analysis eligibility.
 index=clinical.reset_index();index['EXAMDATE']=index.EXAMDATE.astype(str)
 tau=pd.read_csv(RESTRICTED/'TAU_PET_FLAGS_RESTRICTED.csv');tau['date']=tau.date.map(parse_date)
 extended=[]
 for tracer in ['ALL','FTP','MK6240','PI2620']:
  frame=tau if tracer=='ALL' else tau.loc[tau.tracer.eq(tracer)]
  for window in [730,1096,1826,3653,100000]:
   matches=match_records(index,frame,'META_TEMPORAL_SUVR_EXACT',window)
   ds=[r['delta_days'] for r in matches.values()]
   extended.append({'tracer':tracer,'window_days':window,'n':len(matches),'n_before':sum(d<0 for d in ds),
       'n_after':sum(d>0 for d in ds),'minimum_absolute_days':min(map(abs,ds)) if ds else None,
       'median_absolute_days':float(np.median(np.abs(ds))) if ds else None,
       'purpose':'Time-gap explanation only; not contemporaneous validation'})
 write_csv(RESULTS/'TAU_PET_EXTENDED_TIMING.csv',extended)
 cov=pd.read_csv(RESTRICTED/'COVARIATE_FLAGS_RESTRICTED.csv')
 flags=[c for c in cov if 'CONFLICT' in c or c.startswith('AGE_SOURCE')]
 write_csv(RESULTS/'COVARIATE_CONFLICT_COUNTS.csv',[{'flag':c,'n':int(cov[c].sum())} for c in flags])
 chosen=clinical.loc[clinical.ELECSYS_PTAU_W0].copy()
 batch=chosen.PTAU_ANCHOR_BATCH_W0.value_counts().rename_axis('batch').reset_index(name='n')
 write_csv(RESULTS/'SELECTED_PTAU_ASSAY_BATCH_COUNTS.csv',batch)
 # Descriptive category occupancy in technical batches. No association test.
 occupancy=[]
 for c in ['NULISA_PlateID','SOMASCAN_PlateId','TMT_plex']:
  for label,part in chosen.groupby(c):
   cats=part.PTAU_ANCHOR_AT_W0
   row={'technical_factor':c,'batch_label':str(label),'n':len(part),'n_operational_AT_categories':int(cats[cats.ne('UNKNOWN')].nunique())}
   row.update({cat:int(cats.eq(cat).sum()) for cat in ['A-T-','A-T+','A+T-','A+T+','UNKNOWN']});occupancy.append(row)
 write_csv(RESULTS/'TECHNICAL_BATCH_CATEGORY_OCCUPANCY.csv',occupancy)
 write_csv(RESULTS/'COHORT_TECHNICAL_DESIGN.csv',[
  {'cohort':'TRIPLE571','n':len(clinical),**design(clinical)},
  {'cohort':'SAME_DAY_PTAU_AGE_SEX','n':int((clinical.ELECSYS_PTAU_W0&clinical.COV_MINIMAL).sum()),
   **design(clinical.loc[clinical.ELECSYS_PTAU_W0&clinical.COV_MINIMAL])}])
 save_json(RESULTS/'PROTEIN_COVERAGE_SUMMARY.json',{'index_n':571,'annotation_candidates':43,'extended_candidates':70,
   'protein_value_columns_read':0,'technical_design_includes_actual_age_values':False,
   'index_date_min':str(index.EXAMDATE.min()),'index_date_max':str(index.EXAMDATE.max()),
   'per_protein_rows':len(data),'main_min_median_max':data.loc[data.annotation_single_feature&data.route.eq('GIS_NUMERIC'),'n_same_day_ptau_age_sex'].agg(['min','median','max']).to_dict()})
 log('Per-protein intersections and technical design audit complete; boolean masks only')
 print((RESULTS/'PROTEIN_COVERAGE_SUMMARY.json').read_text(encoding='utf-8'))

if __name__=='__main__':run()
