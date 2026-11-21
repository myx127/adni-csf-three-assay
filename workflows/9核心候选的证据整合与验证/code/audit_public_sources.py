from config import *
import requests,re
from pypdf import PdfReader

def run():
    method_rows=[]
    for name in ['ADNI_EMORY_CSF_TMT_MS_Methods_20251007.pdf','ADNI_CruchagaLab_Methods_CSF_SomaScan7K_Proteomic_Data_20230817.pdf','ADNI_BSHRI_PLA_CSF_NULISA_CNS_METHODS_20260116.pdf']:
        p=RAW/'CSF'/name
        reader=PdfReader(p,strict=False);text='\n'.join(f'[PAGE {i+1}]\n'+(page.extract_text() or '') for i,page in enumerate(reader.pages))
        out=EVIDENCE/(p.stem+'_EXTRACTED.txt');out.write_text(text,encoding='utf-8')
        method_rows.append({'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size,'evidence_extract':str(out.relative_to(ROOT)),'synapse_ids':'|'.join(sorted(set(re.findall(r'syn\d+',text))))})
    csvout(RESULTS/'ORIGINAL_ASSAY_METHOD_SOURCES.csv',method_rows)
    frozen=EVIDENCE/'PUBLIC_REPOSITORY_AUDIT_SNAPSHOT.json'
    if frozen.exists() and '--refresh' not in sys.argv:
        snapshot=json.loads(frozen.read_text(encoding='utf-8'))
        for key,name in [('requests','PUBLIC_REPOSITORY_REQUEST_LOG.csv'),('entities','PUBLIC_REPOSITORY_ENTITY_AUDIT.csv'),('children','PUBLIC_REPOSITORY_CHILDREN_METADATA.csv')]:csvout(RESULTS/name,snapshot[key])
        external=pd.read_csv(S7/'results/EXTERNAL_VALIDATION_OPPORTUNITIES.csv')
        external['stage9_individual_validation_completed']=False
        external['stage9_action']='Prior source applicability reviewed; no account access or author contact'
        external.loc[external.source_id.eq('X01'),'stage9_action']='Frozen public Synapse entity/children metadata and403 file-access snapshot; no individual-level replication'
        csvout(RESULTS/'EXTERNAL_VALIDATION_READINESS.csv',external)
        log('Original methods extracted; frozen public repository metadata retained; no network/account action')
        return
    requests_log=[];entities=[];children=[]
    session=requests.Session()
    for ident in ['syn30549757','syn59804727']:
        url=f'https://repo-prod.prod.sagebase.org/repo/v1/entity/{ident}'
        try:
            response=session.get(url,timeout=15)
            record={'resource':ident,'url':url,'method':'GET','status_code':response.status_code,'time':now(),'used_authentication':False}
            if response.ok:
                obj=response.json();save_json(EVIDENCE/(ident+'_entity.json'),obj);entities.append({'id':ident,'name':obj.get('name'),'type':obj.get('concreteType'),'access':'PUBLIC_METADATA_ONLY'})
            else:record['reason']=response.text[:300]
            requests_log.append(record)
            if response.ok and obj.get('concreteType','').endswith(('Project','Folder')):
                cu='https://repo-prod.prod.sagebase.org/repo/v1/entity/children'
                body={'parentId':ident,'includeTypes':['folder','file'],'sortBy':'NAME','sortDirection':'ASC','includeTotalChildCount':True}
                cr=session.post(cu,json=body,timeout=15)
                rec={'resource':ident+'_children','url':cu,'method':'POST_READ_ONLY_LIST','status_code':cr.status_code,'time':now(),'used_authentication':False}
                if cr.ok:
                    co=cr.json();save_json(EVIDENCE/(ident+'_children.json'),co)
                    for child in co.get('page',[]):
                        children.append({'parent':ident,'id':child.get('id'),'name':child.get('name'),'type':child.get('type'),'benefactorId':child.get('benefactorId'),'listed_metadata_only':True})
                    rec['totalChildCount']=co.get('totalChildCount')
                else:rec['reason']=cr.text[:300]
                requests_log.append(rec)
        except Exception as e:requests_log.append({'resource':ident,'url':url,'method':'GET_OR_READ_ONLY_LIST','time':now(),'error':str(e),'used_authentication':False})
    csvout(RESULTS/'PUBLIC_REPOSITORY_REQUEST_LOG.csv',requests_log)
    csvout(RESULTS/'PUBLIC_REPOSITORY_ENTITY_AUDIT.csv',entities)
    csvout(RESULTS/'PUBLIC_REPOSITORY_CHILDREN_METADATA.csv',children)
    external=pd.read_csv(S7/'results/EXTERNAL_VALIDATION_OPPORTUNITIES.csv')
    external['stage9_individual_validation_completed']=False
    external['stage9_action']='Prior source applicability reviewed; no account access or author contact'
    external.loc[external.source_id.eq('X01'),'stage9_action']='Public Synapse entity/children metadata requested; individual-level replication not performed'
    csvout(RESULTS/'EXTERNAL_VALIDATION_READINESS.csv',external)
    log(f'Original assay methods extracted;{len(requests_log)} public read-only repository requests audited; no account action')
    print(pd.DataFrame(entities).to_string(index=False))
    print(pd.DataFrame(children).to_string(index=False))
    print(pd.DataFrame(requests_log)[[c for c in ['resource','status_code','totalChildCount','error'] if c in pd.DataFrame(requests_log)]].to_string(index=False))

if __name__=='__main__':
    assert_unsealed()
    run()
