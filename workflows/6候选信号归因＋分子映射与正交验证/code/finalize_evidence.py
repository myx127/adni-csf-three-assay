"""Combine all audited sources without turning coverage or association into replication."""
from config import *

def run():
 initial=RESULTS/'EMORY_PANEL_DIRECT_VALIDATION_STATUS.json';old=RESULTS/'DIRECT_DISCOVERY_VALIDATION_STATUS.json'
 if not initial.exists():
  obj=json.loads(old.read_text(encoding='utf-8'));assert obj['NPTX1']=='NO_DIRECT_TARGETED_ASSAY_COVERAGE'
  obj['scope']='Emory 2022/2023 V01/V02 only; initial snapshot before later legacy source extension';save_json(initial,obj)
 associations=pd.read_csv(RESULTS/'LEGACY_NEW_RID_PEPTIDE_PTAU_ASSOCIATION.csv');screen=pd.read_csv(RESULTS/'LEGACY_DIRECT_DESIGN_SCREEN.csv')
 relative=pd.read_csv(RESULTS/'LEGACY_SAME_DATE_RELATIVE_DESCRIPTIVE.csv');specificity=pd.read_csv(RESULTS/'HUMAN_REFERENCE_PEPTIDE_SPECIFICITY.csv')
 regions=pd.read_csv(RESULTS/'UNIPROT_REGIONS_AND_PROCESSING.csv');mapping=[];status={};rows=[]
 for acc,gene in CANDIDATES.items():
  part=specificity.loc[specificity.accession.eq(acc)&specificity.panel.eq('CAPRION_2019')]
  a=associations.loc[associations.accession.eq(acc)&associations.cohort.eq('NEW_RID_CLINICAL')]
  r=relative.loc[relative.accession.eq(acc)&relative.cohort.eq('S5_SAME_DATE')]
  state={'accession':acc,'gene':gene,'direct_Emory_panel_peptides':0,'direct_2019_Caprion_peptides':len(part),
    'direct_2014_release_peptides':len(part),'2014_vs_2019_distinct_experiments_not_release_aliases':True,
    '2014_exact_date_and_technical_crosswalk_resolved':False,'2014_candidate_associations_fitted':False,
    'strict_relative_disagreement_validation_completed':False,'external_cohort_validation_completed':False,
    '2019_same_day_n_min':int(r.n.min()) if len(r) else 0,'2019_same_day_n_max':int(r.n.max()) if len(r) else 0,
    '2019_same_day_residual_df_min':int(r.residual_df.min()) if len(r) else 0,'2019_same_day_residual_df_max':int(r.residual_df.max()) if len(r) else 0,
    'new_RID_single_assay_n':int(a.n.min()) if len(a) else 0,'new_RID_peptides_nominal_positive_p05':int(a.p_positive_direction.le(.05).sum()) if len(a) else 0,
    'new_RID_peptides_BY_positive_q05':int(a.positive_association_BY.sum()) if len(a) else 0,
    'new_RID_best_q_BY5':float(a.q_BY5_single_assay.min()) if len(a) else None,
    'available_conclusion':'DIRECT_PEPTIDE_COVERAGE_WITH_EXPLORATORY_PATHOLOGY_ASSOCIATION_ONLY' if len(part) else 'NO_DIRECT_TARGETED_RECEPTOR_COVERAGE',
    'no_validated_mechanism_or_binding_epitope':True}
  status[gene]=state;rows.append(state)
  for p in part.to_dict('records'):
   start=int(p['canonical_start_1based']);end=int(p['canonical_end_1based']);reg=regions.loc[regions.accession.eq(acc)&regions.start.le(end)&regions.end.ge(start)]
   mapping.append({**p,'overlapping_annotated_regions':'; '.join(f'{r.feature_type}:{int(r.start)}-{int(r.end)}:{r.description if isinstance(r.description,str) else ""}' for r in reg.itertuples()),
    'released_two_transitions_are_same_peptide':True,'does_not_establish_cleaved_active_or_intact_species':True})
 save_json(old,{**status,'direct_discovery_accessions_covered_count':2,'direct_three_discovery_validation_completed':False,
   'relative_disagreement_replications_completed':0,'external_cohort_validation_completed':False,
   'status':'AVAILABLE_SOURCE_ANALYSES_COMPLETE; direct molecular forms and relative-pattern replication unresolved',
   'reason':'Two discoveries have direct legacy peptides; exact-date common multi-assay design fails frozen df40 gate. New RID peptide-pTau association tests a different question; no BY5 positive associations.',
   'Emory_initial_scope_record':'EMORY_PANEL_DIRECT_VALIDATION_STATUS.json'})
 csvout(RESULTS/'THREE_DISCOVERY_FINAL_EVIDENCE_STATUS.csv',rows);csvout(RESULTS/'LEGACY_PEPTIDE_CANONICAL_MAPPING.csv',mapping)
 # Explicitly register the nine distinct evidence uses. No repeated release counts as replication.
 evidence=[
  ('S5_FROZEN','34 protein family; N/S/TMT and Elecsys pTau','Discovery','571 indexed; 542-555 for three discoveries','Original three discoveries under frozen conditional mean/covariance assumptions','No arbitrary-error-covariance identification'),
  ('ATTRIBUTION','same cases; four covariate models','Postdiscovery explanation','35 proteins x 4 models','Exact algebraic changes in H,V,Q and BY; not causal mediation','No mechanistic percentage explained'),
  ('EMORY_V01_V02','two Emory targeted releases, duplicate handling fixed','Additional measurement','S5 overlap 538-560; new RID84-87','4 other genes /5 peptides, not NPTX1/PDGFRB/KLK6','Same cohort; new relative-design df insufficient'),
  ('FOUR_ASSAY','N/S/original TMT/Emory targetedMS','Posthoc covariance-direction exploration','5 peptides of4 other genes','Cross-block rank1 under zero cross-block error covariance','Does not test original three discoveries or prove pTau-response equivalence'),
  ('LEGACY_2014','L01 and L02 are one MRM release family','Coverage and source audit','303 released rows;287 RID','Direct NPTX1 two and KLK6 three peptide labels','No verified exact date/technical crosswalk; associations not fitted; not two independent assays'),
  ('CAPRION_2019_SAME_DATE','L04 separate longitudinal MRM experiment','Direct same-person measurement','90-94 matched cases; residual df20-24','Direct peptide descriptive correlations and slopes','Frozen df40 gate fails; no relative-constraint P values'),
  ('CAPRION_2019_NEW_RID','earliest reported date; exclude all571 S5 RID','Subject-disjoint same-ADNI peptide-pTau association','95 RID;5 peptides','Separate fixed five-peptide BY family; calibrated SUM2/t4','No BY rejection; different hypothesis; no external cohort'),
  ('MOLECULAR_MAPPING','official UniProt 2026_03 reference and12 entry JSONs','Sequence interpretation','169651 canonical/isoform sequences;16 peptide targets','Gene-specificity within defined reference; canonical regions','Sequence identity does not prove analytical specificity, binding epitope or active/mature form'),
  ('EQUIVALENCE','S5 frozen35 designs; working margins .2/.35/.5SD','Exploratory positive-consensus test','34 family;11 eligible;0 BY equivalent','Fieller90 joint interval/IUT with finite boundary calibration','Postdiscovery working margins; absence of rejection is not consensus or global disagreement')]
 csvout(RESULTS/'EVIDENCE_USE_AND_INDEPENDENCE_REGISTRY.csv',[dict(zip(['evidence_id','source','role','scope','permitted_interpretation','limit'],e)) for e in evidence])
 log('Final source/evidence matrix reconciled: two directly covered genes; zero formal relative-pattern replications and zero external-cohort validations')

if __name__=='__main__':run()
