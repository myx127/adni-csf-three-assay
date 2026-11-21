from __future__ import annotations

import hashlib
import json
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
STAGE1 = ROOT.parent / '1候选逐蛋白可评价性审计'
RESULTS, LOGS, FIGURES = ROOT/'results', ROOT/'logs', ROOT/'figures'
for directory in (RESULTS, LOGS, FIGURES):
    directory.mkdir(parents=True, exist_ok=True)
SEED = 20261027
METHODS = ['PCA_IN_SAMPLE', 'PCA_CROSSFIT_2FOLD', 'COMMON_LOADING_MEAN_WALD',
           'INVARIANCE_JOINT_IID', 'INVARIANCE_MEAN_IID', 'INVARIANCE_COV_IID',
           'INVARIANCE_JOINT_TMT_CLUSTER', 'INVARIANCE_JOINT_3WAY_CLUSTER', 'ORACLE_COMPLETE_TRUTH']
COMPATIBLE = ['N01','N02','N03','N04','N05','N06','N11']

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def save_json(path, value):
    def clean(x):
        if isinstance(x,dict): return {str(k):clean(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)): return [clean(y) for y in x]
        if isinstance(x,np.ndarray): return clean(x.tolist())
        if isinstance(x,np.generic): return clean(x.item())
        if isinstance(x,float) and not np.isfinite(x): return None
        return x
    Path(path).write_text(json.dumps(clean(value),ensure_ascii=False,indent=2),encoding='utf-8')

def csv(path, frame):
    if not isinstance(frame,pd.DataFrame): frame=pd.DataFrame(frame)
    frame.to_csv(path,index=False,encoding='utf-8-sig')

def log(message):
    line=datetime.now().isoformat(timespec='seconds')+' '+str(message)
    print(line,flush=True)
    with (LOGS/'SIMULATION_RUN.log').open('a',encoding='utf-8') as f: f.write(line+'\n')

def scenarios():
    rows=[
        ('N01','无病理作用','NO_EXTRA_RELATIVE_EFFECT','无病理均值作用，等精度',{'beta':0.,'sd':[.8,.8,.8]}),
        ('N02','共同均值_等精度','NO_EXTRA_RELATIVE_EFFECT','仅共同均值变化，等精度',{'sd':[.8,.8,.8]}),
        ('N03','共同均值_不等精度','NO_EXTRA_RELATIVE_EFFECT','仅共同均值变化，不等精度',{}),
        ('N04','仅共同方差变化','NO_EXTRA_RELATIVE_EFFECT','均值不变，潜在方差随组变化',{'beta':0.,'uvar':[.5,1.,1.8]}),
        ('N05','共同均值及方差变化','NO_EXTRA_RELATIVE_EFFECT','共同均值和方差均变化',{'uvar':[.6,1.,1.6]}),
        ('N06','误差方差随病理变化','NO_EXTRA_RELATIVE_EFFECT','M3噪声SD乘exp(0.35P)',{'hetero':.35}),
        ('N07','固定相关误差','NO_EXTRA_RELATIVE_EFFECT','M1/M2误差相关0.35',{'rho':.35}),
        ('N08','变化的相关误差','NO_EXTRA_RELATIVE_EFFECT','M1/M2误差相关按组为0.05/0.25/0.50',{'rho':[.05,.25,.50]}),
        ('N09','固定二次响应','NO_EXTRA_RELATIVE_EFFECT','M3=T+0.35*T²，响应不随病理改变',{'response':'quadratic'}),
        ('N10','固定饱和响应','NO_EXTRA_RELATIVE_EFFECT','M3=1.8*tanh(T/1.8)，响应不随病理改变',{'response':'saturation'}),
        ('N11','MCAR20百分比','NO_EXTRA_RELATIVE_EFFECT','每平台独立20%缺失，完整记录分析',{'missing':'mcar'}),
        ('N12','LOD30删除','NO_EXTRA_RELATIVE_EFFECT','M1低于理论边际30%阈值设缺失',{'missing':'lod30'}),
        ('N13','LOD70删除','NO_EXTRA_RELATIVE_EFFECT','M1低于理论边际70%阈值设缺失',{'missing':'lod70'}),
        ('N14','LOD30地板替代','NO_EXTRA_RELATIVE_EFFECT','M1低于阈值以阈值替代，示范错误处理风险',{'missing':'floor30'}),
        ('N15','平衡GIS共享噪声','NO_EXTRA_RELATIVE_EFFECT','P与plex无关，同plex GIS噪声SD0.6被共享',{'batch':'gis_balanced'}),
        ('N16','未校正plex混杂','NO_EXTRA_RELATIVE_EFFECT','P与plex偏移相关，未用GIS去除偏移',{'batch':'raw_confounded'}),
        ('N17','GIS后plex混杂设计','NO_EXTRA_RELATIVE_EFFECT','P与plex偏移相关；GIS去除偏移但保留共享参照噪声',{'batch':'gis_confounded'}),
        ('N18','选择性缺失','NO_EXTRA_RELATIVE_EFFECT','M3是否可用依赖其读数与P',{'missing':'selective'}),
        ('N19','弱共同信号','NO_EXTRA_RELATIVE_EFFECT','共同残差方差0.05，beta0.10，SD0.8/1/1.5',{'beta':.10,'uvar':[.05]*3,'sd':[.8,1.,1.5]}),
        ('N20','交叉批次共享误差','NO_EXTRA_RELATIVE_EFFECT','三平台分别共享19/36/65批次扰动SD0.45',{'batch':'crossed'}),
        ('A01','真实相对截距偏移','EXTRA_RELATIVE_EFFECT','M2额外0.35*P',{'delta':[0.,.35,0.]}),
        ('A02','真实相对载荷变化','EXTRA_RELATIVE_EFFECT','M3额外0.40*P*T',{'loading_delta':[0.,0.,.40]}),
        ('A03','截距与载荷组合偏离','EXTRA_RELATIVE_EFFECT','M2额外0.25P；M3额外0.30P*T',{'delta':[0.,.25,0.],'loading_delta':[0.,0.,.30]}),
        ('U01','平行共同方向直接项','OBSERVATIONALLY_UNIDENTIFIABLE','三个平台同时额外0.30P；可被共同信号均值吸收',{'delta':[.30]*3}),
    ]
    result=[]
    for index,(sid,label,truth,description,settings) in enumerate(rows):
        result.append({'scenario_id':sid,'scenario_index':index,'label_zh':label,'truth_class':truth,
                       'description':description,'settings':settings,'compatible_linear_independent_null':sid in COMPATIBLE})
    return result

def design_grid(registry):
    cells=[]
    def add(phase, row, n, reps, continuous=False, effect=None):
        # Seed components are independent of execution order and process scheduling.
        key=[SEED, {'MAIN':1,'N_GRID':2,'CONTINUOUS':3,'EFFECT_CURVE':4}[phase],row['scenario_index'],n,int(continuous),0 if effect is None else int(round(effect*1000))+1]
        seed=int(np.random.SeedSequence(key).generate_state(1,dtype=np.uint64)[0])
        cells.append({'cell_id':f'{phase}_{row["scenario_id"]}_n{n}'+('' if effect is None else f'_e{effect:.2f}'),
                      'phase':phase,'scenario_id':row['scenario_id'],'n':n,'repetitions':reps,
                      'continuous_P':continuous,'offset_effect_override':effect,'seed':str(seed)})
    lookup={r['scenario_id']:r for r in registry}
    for row in registry: add('MAIN',row,571,2000)
    for sid in ['N02','N03','N07','N09','N19','A01']:
        for n in [200,400]: add('N_GRID',lookup[sid],n,500)
    for sid in ['N02','N03','N05','N07','N09','N12','A01','A02']: add('CONTINUOUS',lookup[sid],571,1000,True)
    for effect in [0.,.10,.20,.35,.50]: add('EFFECT_CURVE',lookup['A01'],571,1000,effect=effect)
    assert sum(r['repetitions'] for r in cells)==67000
    return cells
