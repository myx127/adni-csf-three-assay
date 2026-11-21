from config import *
import requests

def run():
    target=EVIDENCE/'PUBLIC_DATA_DOWNLOAD_AUDIT.json'
    if target.exists() and '--refresh' not in sys.argv:
        log('Frozen external file-access snapshot retained; no network download attempted');return
    session=requests.Session();records=[]
    for ident in ['syn30549894','syn30550193','syn50877009']:
        url='https://repo-prod.prod.sagebase.org/repo/v1/entity/'+ident
        try:
            response=session.get(url,timeout=15);rec={'entity':ident,'metadata_status':response.status_code,'authenticated':False}
            if response.ok:
                obj=response.json();save_json(EVIDENCE/(ident+'_entity.json'),obj);rec['name']=obj.get('name')
                file_response=session.get(url+'/file',params={'redirect':'false'},timeout=15,allow_redirects=False)
                rec['file_access_status']=file_response.status_code
                if not file_response.ok:rec['file_access_response']=file_response.text[:300]
                # Only the verified public permission response is recorded.
                # Newly accessible data require a new frozen analysis branch, rather than silently entering old analyses.
                if file_response.ok:rec['next_step']='File now accessible; register data and protocol in a new branch before ingestion'
            records.append(rec)
        except Exception as exc:records.append({'entity':ident,'error':str(exc),'authenticated':False})
    save_json(target,records)
    log('Unauthenticated file-access checks archived; no signed URLs or credentials saved')

if __name__=='__main__':
    assert_unsealed()
    run()
