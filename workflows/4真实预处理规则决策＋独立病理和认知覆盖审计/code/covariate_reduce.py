from __future__ import annotations
from collections import defaultdict,Counter
import re
import pandas as pd
from audit_common import *
from clinical_reduce import stream,STRUCTURE,consolidate

def covariates(index):
 targets={int(r.RID):parse_date(r.EXAMDATE) for r in index.itertuples()};demo=defaultdict(list);counts=Counter()
 fields=['RID','VISDATE','VISCODE2','PTDOB','PTDOBYY','PTGENDER','PTEDUCAT','HAS_QC_ERROR']
 for line,row in stream('I13',fields):
  counts['rows']+=1;person=rid(row.get('RID'))
  if person not in targets:continue
  if text(row.get('HAS_QC_ERROR'))=='1':counts['explicit_QC_error']+=1;continue
  day=parse_date(row.get('VISDATE'));dob=parse_date(row.get('PTDOB'));year=numeric(row.get('PTDOBYY'))
  # This release stores PTDOB as MM/YYYY and PTDOBYY as YYYY-MM-DD.
  # The latter is treated as year information, never as proof of an exact birth day.
  month_year=re.fullmatch(r'(\d{1,2})/(\d{4})',text(row.get('PTDOB')))
  month_year_valid=bool(month_year and 1<=int(month_year[1])<=12 and 1800<int(month_year[2])<2027)
  year_date=parse_date(row.get('PTDOBYY'))
  year_value=int(year.low) if year.exact and year.low.is_integer() and 1800<year.low<2027 else year_date.year if year_date else None
  birth_year=dob.year if dob else int(month_year[2]) if month_year_valid else year_value
  birth_month=dob.month if dob else int(month_year[1]) if month_year_valid else None
  dob_conflict=bool(birth_year and year_value and birth_year!=year_value)
  counts['birth_year_field_disagreement']+=int(dob_conflict)
  approximate=(targets[person].year-birth_year) if birth_year else None
  ageok=approximate is not None and 0<=approximate<=120 and not dob_conflict
  sex=text(row.get('PTGENDER'));sex=sex if sex in ['1','2'] else ''
  education=numeric(row.get('PTEDUCAT'));eduok=education.exact and 0<=education.low<=30
  demo[person].append({'date':day,'age_valid':ageok,'birth_signature':hashlib.sha256(str(birth_year).encode()).hexdigest() if ageok else '',
    'birth_month':birth_month,'birth_field_conflict':dob_conflict,
    'age_source':'FULL_DATE' if dob else 'MONTH_YEAR' if month_year_valid else 'YEAR_ONLY' if ageok else 'UNAVAILABLE','gender_category':sex,
    'education_available':eduok,'source_row':line,'education_signature':hashlib.sha256(str(education.low).encode()).hexdigest() if eduok else ''})
 STRUCTURE.append({'source':'I13',**counts,'index_RID_with_any_record':len(demo)})
 apo=defaultdict(list);counts=Counter()
 for line,row in stream('I14',['RID','GENOTYPE','APUSABLE']):
  counts['rows']+=1;person=rid(row.get('RID'))
  if person not in targets:continue
  if text(row.get('APUSABLE'))=='0':counts['explicit_unusable']+=1;continue
  raw=text(row.get('GENOTYPE'));alleles=re.findall(r'[234]',raw)
  if len(alleles)!=2:apo[person].append(None);continue
  apo[person].append(sum(a=='4' for a in alleles))
 STRUCTURE.append({'source':'I14',**counts,'index_RID_with_any_record':len(apo)})
 outputs=[];category_counts=Counter()
 for person,day in targets.items():
  records=demo.get(person,[]);birth={r['birth_signature'] for r in records if r['age_valid']};sex={r['gender_category'] for r in records if r['gender_category']}
  dated=[r for r in records if r['date'] and r['education_available'] and abs((r['date']-day).days)<=365]
  dated.sort(key=lambda r:(abs((r['date']-day).days),r['date']>day,r['source_row']))
  education=bool(dated);same_date_conflict=False
  if dated:
   first=dated[0];same=[r for r in dated if r['date']==first['date']]
   same_date_conflict=len({r['education_signature'] for r in same})>1
   education=not same_date_conflict
  copies={v for v in apo.get(person,[]) if v is not None};carrier={v>0 for v in copies}
  months={r['birth_month'] for r in records if r['age_valid'] and r['birth_month'] is not None}
  birth_conflict=len(birth)>1 or len(months)>1 or any(r['birth_field_conflict'] for r in records)
  age=len(birth)==1 and not birth_conflict;gender=len(sex)==1;e4=len(carrier)==1
  if gender:category_counts['SEX_'+next(iter(sex))]+=1
  if e4:category_counts['APOE4_CARRIER_'+str(int(next(iter(carrier))))]+=1
  if len(copies)==1:category_counts['APOE4_COPIES_'+str(next(iter(copies)))]+=1
  out={'RID':person,'AGE_AVAILABLE':age,'AGE_BIRTH_INFO_CONFLICT':birth_conflict,
       'SEX_AVAILABLE':gender,'SEX_CONFLICT':len(sex)>1,'EDUCATION_365_AVAILABLE':education,
       'EDUCATION_SAME_DATE_CONFLICT':same_date_conflict,'APOE4_AVAILABLE':e4,'APOE4_DOSAGE_AVAILABLE':len(copies)==1,
       'APOE4_CARRIER_CONFLICT':len(carrier)>1,'EDUCATION_ABS_DAYS':abs((dated[0]['date']-day).days) if dated else None,
       'AGE_SOURCE_YEAR_ONLY':age and all(r['age_source']=='YEAR_ONLY' for r in records if r['age_valid']),
       'AGE_SOURCE_MONTH_YEAR':age and any(r['age_source']=='MONTH_YEAR' for r in records if r['age_valid'])}
  out['COV_MINIMAL']=age and gender
  out['COV_COGNITIVE']=age and gender and education and e4
  outputs.append(out)
 write_csv(RESULTS/'COVARIATE_CATEGORY_COUNTS.csv',[{'category':k,'n_participants':v} for k,v in sorted(category_counts.items())])
 data=pd.DataFrame(outputs)
 write_csv(RESULTS/'COVARIATE_AVAILABILITY.csv',[{'covariate_flag':c,'n_available':int(data[c].sum()),'denominator':len(data)} for c in data if c.endswith('_AVAILABLE') or c.startswith('COV_')])
 return data

def load_diagnosis():
 records=[];counts=Counter();people=set()
 for line,row in stream('I15',['RID','EXAMDATE','VISCODE2','DIAGNOSIS','HAS_QC_ERROR']):
  counts['rows']+=1;person=rid(row.get('RID'));day=parse_date(row.get('EXAMDATE'))
  if not person or not day:counts['invalid_ID_or_date']+=1;continue
  people.add(person);qc=text(row.get('HAS_QC_ERROR'))!='1';code=text(row.get('DIAGNOSIS'))
  value={'1':'CN','2':'MCI','3':'DEMENTIA','10':'TEAM_NODX'}.get(code,'UNKNOWN') if qc else 'UNKNOWN'
  records.append({'RID':person,'date':day,'source_row':line,'DX_CAT':value,'DX_VALID':value in ['CN','MCI','DEMENTIA'],
                   'VISCODE':text(row.get('VISCODE2')),'_sig_DX':value if value!='UNKNOWN' else ''})
 STRUCTURE.append({'source':'I15',**counts,'unique_RID':len(people)})
 result=consolidate(records,['RID','date'],'I15')
 result['DX_VALID']=result.DX_CAT.isin(['CN','MCI','DEMENTIA'])
 return result
