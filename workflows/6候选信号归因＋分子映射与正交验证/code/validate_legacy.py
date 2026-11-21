from config import *
from scipy.linalg import lstsq
from scipy.stats import t
from statsmodels.stats.multitest import multipletests
import zipfile,pickle,re,gzip

def validate(check):
 for r in pd.read_csv(RESULTS/'LEGACY_INPUT_MANIFEST.csv').itertuples():check('legacy_source_hash_'+r.input_id,sha(r.path)==r.sha256)
 ref=json.loads((RESULTS/'HUMAN_REFERENCE_PROTEOME_MANIFEST.json').read_text(encoding='utf-8'));check('reference_proteome_hash',sha(ROOT/ref['path'])==ref['sha256'])
 first=json.loads((RESULTS/'FIRST_LEGACY_QUANTITATIVE_ACCESS.json').read_text(encoding='utf-8'));check('legacy_SAP_linked_before_association',first['plan_sha256_before_association']==sha(RESULTS/'LEGACY_VERIFICATION_SAP_FROZEN.json'))
 reg=pd.read_csv(RESULTS/'LEGACY_CANDIDATE_TRANSITION_REGISTRY.csv');manifest=pd.read_csv(RESULTS/'LEGACY_INPUT_MANIFEST.csv').set_index('input_id')
 with zipfile.ZipFile(manifest.loc['L04','path']) as z:
  with z.open('csfproteomics/adni_csfproteomics.csv') as f:raw=pd.read_csv(f,low_memory=False)
 check('legacy_raw750_rows',len(raw)==750)
 date=pd.to_datetime(raw.EXAMDATE,errors='coerce');valid=date.notna()&pd.to_numeric(raw.RID,errors='coerce').gt(0)
 keys=raw.RID.astype('Int64').astype(str)+'|'+date.dt.strftime('%Y-%m-%d').fillna('')
 check('legacy_unique729_dates_198_RID',keys[valid].nunique()==729 and raw.loc[valid,'RID'].nunique()==198)
 frame=pd.read_csv(RESTRICTED/'LEGACY_L04_CLINICAL_DESIGN_RESTRICTED.csv');cache=np.load(RESTRICTED/'LEGACY_PEPTIDE_CACHE_RESTRICTED.npz')
 check('legacy_cache_keys729',len(frame)==729 and frame.key.is_unique and np.array_equal(frame.key,cache['keys']))
 selected=raw.iloc[frame.source_row.to_numpy(int)-2]
 check('legacy_raw_date_RID_alignment',np.array_equal(keys.iloc[frame.source_row.to_numpy(int)-2],frame.key))
 for j,pep in enumerate(cache['peptides']):
  columns=reg.loc[reg.peptide.eq(pep),'transition_field'].tolist();a=selected[columns].apply(pd.to_numeric,errors='coerce').to_numpy(float)
  independent=(a[:,0]+a[:,1])/2;check('legacy_log_mean_no_second_log_'+pep,np.allclose(cache['Y'][:,j],independent,equal_nan=True,atol=1e-12))
 # Re-select blind duplicates solely using reported technical metadata.
 selected_meta=raw.loc[valid,['RUN.NUMBER','SUBRUN','RUN.INJECTION.NUMBER']].copy();selected_meta['key']=keys[valid];selected_meta['line']=selected_meta.index+2
 earliest=selected_meta.sort_values(['key','RUN.NUMBER','SUBRUN','RUN.INJECTION.NUMBER','line']).drop_duplicates('key').set_index('key')
 check('legacy_duplicate_selection_metadata_only',np.array_equal(frame.source_row,earliest.reindex(frame.key).line))
 old=pd.read_csv(S5/'00_restricted/CLINICAL_DESIGN_RESTRICTED.csv').set_index('key');part=frame.loc[frame.key.isin(old.index)].set_index('key')
 for c in ['P_LOG','age_month_precision','sex_code2']:check('legacy_same_day_clinical_matches_S5_'+c,np.allclose(part[c],old.reindex(part.index)[c],equal_nan=True,atol=1e-10))
 ds=pickle.load((RESTRICTED/'LEGACY_CLINICAL_DESIGNS_RESTRICTED.pkl').open('rb'));assoc=pd.read_csv(RESULTS/'LEGACY_NEW_RID_PEPTIDE_PTAU_ASSOCIATION.csv').set_index('design_id')
 numerics=[]
 for ident,d in ds.items():
  row=assoc.loc[ident];j=list(cache['peptides']).index(row.peptide);y=cache['Y'][d['indices'],j];X=np.c_[d['Z'],d['P']]
  coef,_,rank,_=lstsq(X,y);E=y-X@coef;gcoef,_,_,_=lstsq(d['Z'],d['P']);rP=d['P']-d['Z']@gcoef;u=E*rP/(rP@rP);V=0.
  for axis in range(2):
   lab=d['groups'][:,axis];sums=np.array([u[lab==g].sum() for g in np.unique(lab)]);sums-=sums.mean();V+=sums@sums*len(sums)/(len(sums)-1)
  se=np.sqrt(V);tdf=min(len(set(d['groups'][:,0])),len(set(d['groups'][:,1])))-1;p=t.sf(coef[-1]/se,tdf)
  check('legacy_independent_OLS_SUM2_t_'+ident,np.allclose([coef[-1],se,p],[row.b_P,row.se_SUM2_CR1,row.p_positive_direction],rtol=1e-8,atol=1e-12))
  subjects=frame.loc[d['indices'],'RID'];check('legacy_95_subject_disjoint_'+ident,len(subjects)==95 and subjects.is_unique and not set(subjects)&set(old.RID))
  dates=frame.loc[d['indices'],['RID','date']].set_index('RID').date;earliestdate=frame.groupby('RID').date.min();check('legacy_earliest_date_before_clinical_selection_'+ident,(dates==earliestdate.reindex(dates.index)).all())
  numerics.append({'design_id':ident,'b_independent':coef[-1],'se_independent':se,'t_df_independent':tdf,'p_independent':p})
 csvout(RESULTS/'INDEPENDENT_LEGACY_NUMERIC_VALIDATION.csv',numerics)
 for cohort,part in assoc.reset_index().groupby('cohort'):
  check('legacy_BY5_reference_'+cohort,len(part)==5 and np.allclose(part.q_BY5_single_assay,multipletests(part.p_for_family,method='fdr_by')[1]))
 relative=pd.read_csv(RESULTS/'LEGACY_SAME_DATE_RELATIVE_DESCRIPTIVE.csv');check('legacy_no_underpowered_relative_inference',relative.formal_relative_p_value.isna().all() and relative.residual_df.lt(40).all())
 oldy=np.load(S5/'00_restricted/PROTEIN_TRANSFORMED_RESTRICTED.npz');meta=pd.read_csv(S1/'results/EXTENDED_CORE_MAPPING_DRAFT.csv').set_index('accession')
 from designs import make
 masks=pickle.load((RESTRICTED/'LEGACY_CASE_MASKS_RESTRICTED.pkl').open('rb'))
 for r in relative.itertuples():
  mask=masks[r.design_id];d=make(frame,mask,tech=('NULISA_PlateID','SOMASCAN_PlateId','MRM_PLATE','MRM_RUN'),additional_factors=('ptau_batch','MRM_INSTRUMENT'))
  idx=frame.loc[mask,'S5_index'].to_numpy(int);j=list(oldy['accessions']).index(r.accession);sj=list(oldy['soma_features']).index(meta.loc[r.accession,'soma_features']);k=list(cache['peptides']).index(r.peptide)
  Y=np.c_[oldy['N'][idx,j],oldy['S'][idx,sj],cache['Y'][mask,k]];X=np.c_[d['Z'],d['P']]
  # Nested technical indicators are exactly dependent. Use a dimension-aware SVD
  # threshold, independently of the production QR basis (default eps admits a
  # spurious ~9e-15 singular value and unstable coefficients on this small set).
  coef,_,rank,_=lstsq(X,Y,cond=max(X.shape)*np.finfo(float).eps);E=Y-X@coef;R=np.corrcoef(E,rowvar=False)
  check('legacy_independent_SVD_rank_'+r.design_id,rank==d['rank'])
  check('legacy_descriptive_correlations_'+r.design_id,np.allclose([R[0,2],R[1,2]],[r.corr_residual_N_targetMS,r.corr_residual_S_targetMS],atol=1e-8))
 # Independent combined-regex scan of the full FASTA, using I/L equivalence.
 expected=pd.read_csv(RESULTS/'HUMAN_REFERENCE_PEPTIDE_ALL_MATCHES.csv');spec=pd.read_csv(RESULTS/'HUMAN_REFERENCE_PEPTIDE_SPECIFICITY.csv')
 peps=sorted(set(spec.peptide));lookup={p.replace('I','L'):p for p in peps};pattern=re.compile('(?=('+'|'.join(map(re.escape,lookup))+'))');hits=set();n=0
 def examine(header,sequence):
  accession=header.split('|')[1];m=re.search(r'\bGN=(\S+)',header);gene=m[1] if m else '';s=sequence.replace('I','L')
  for p in {lookup[m[1]] for m in pattern.finditer(s)}:hits.add((p,accession,gene))
 header=None;seq=[]
 with gzip.open(ROOT/ref['path'],'rt',encoding='utf-8') as f:
  for line in f:
   if line.startswith('>'):
    if header is not None:examine(header,''.join(seq));n+=1
    header=line[1:].strip();seq=[]
   else:seq.append(line.strip())
 if header is not None:examine(header,''.join(seq));n+=1
 comparison=set(expected.loc[expected.match_mode.eq('I_L_EQUIVALENT'),['peptide','matched_accession','matched_gene']].itertuples(index=False,name=None))
 check('independent_reference_proteome_IL_scan_all_hits',hits==comparison,f'{n} sequences; {len(hits)} hits')
 check('reference_sequence_count169651',n==169651)
 for r in spec.itertuples():
  obj=json.loads((EVIDENCE/(r.accession+'_UNIPROT.json')).read_text(encoding='utf-8'));sequence=obj['sequence']['value'];start=sequence.find(r.peptide)+1
  check('independent_canonical_position_'+r.accession+'_'+r.peptide,int(r.canonical_start_1based)==start and start>0)
 status=json.loads((RESULTS/'DIRECT_DISCOVERY_VALIDATION_STATUS.json').read_text(encoding='utf-8'))
 check('final_direct_coverage_two_not_three',status['direct_discovery_accessions_covered_count']==2 and status['NPTX1']['direct_2019_Caprion_peptides']==2 and status['KLK6']['direct_2019_Caprion_peptides']==3 and status['PDGFRB']['direct_2019_Caprion_peptides']==0)
 check('no_overclaimed_relative_or_external_replication',status['relative_disagreement_replications_completed']==0 and not status['external_cohort_validation_completed'])
 check('L01_L02_not_counted_as_distinct_experiments','same 2014 MRM release family' in json.loads((RESULTS/'LEGACY_VERIFICATION_SAP_FROZEN.json').read_text(encoding='utf-8'))['L01_L02'])
 return {'direct_discovery_genes_with_legacy_peptides':2,'direct_legacy_peptides':5,'new_subject_disjoint_RID':95,'reference_proteome_sequences':n}
