from config import *
from openpyxl.styles import Font,PatternFill,Alignment
import platform,importlib.metadata

WORKBOOK='候选信号归因_分子映射与正交验证总表.xlsx'
PRIORITY=['THREE_DISCOVERY_FINAL_EVIDENCE_STATUS','COVARIATE_FOUR_MODEL_RESULTS','EDUCATION_APOE_SHAPLEY_Q',
 'H_V_AND_B_C_DECOMPOSITION','BY_Q_OWN_VS_FAMILY_DECOMPOSITION','THREE_DISCOVERY_MOLECULAR_EVIDENCE_MATRIX',
 'LEGACY_PEPTIDE_CANONICAL_MAPPING','HUMAN_REFERENCE_PEPTIDE_SPECIFICITY','LEGACY_DIRECT_DESIGN_SCREEN',
 'LEGACY_NEW_RID_PEPTIDE_PTAU_ASSOCIATION','LEGACY_SAME_DATE_RELATIVE_DESCRIPTIVE','ORTHOGONAL_MEASUREMENT_SUPPORT_RESULTS',
 'FOUR_ASSAY_CROSSBLOCK_RESULTS','FOUR_ASSAY_RESIDUAL_CORRELATIONS','EQUIVALENCE_WORKING_MARGIN_SUMMARY',
 'OPERATIONAL_EQUIVALENCE_ALL35','RELATIVE_RESPONSE_FIELLER_INTERVALS','EVIDENCE_USE_AND_INDEPENDENCE_REGISTRY',
 'FIGURE_INDEX','INDEPENDENT_QA_CHECKS']

def run():
 environment={'python':sys.version,'platform':platform.platform(),'packages':{p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','statsmodels','matplotlib','openpyxl','pypdf','pyarrow']},
  'main_variance_version':0,'primary_family_unchanged':True,'simulation_seed_names':SEEDS,
  'other_simulation_seed_roots':{'four_assay':2026102906,'four_direction_boundary':2026102907,'legacy_single_assay':2026102908},
  'seeds_are_arbitrary_integer_labels_not_execution_dates':True}
 save_json(RESULTS/'SOFTWARE_ENVIRONMENT.json',environment)
 paths=sorted(RESULTS.glob('*.csv'),key=lambda p:(PRIORITY.index(p.stem) if p.stem in PRIORITY else len(PRIORITY),p.name))
 # Workbook/delivery QA is created after this file. It remains a separate CSV to avoid a self-referential workbook.
 paths=[p for p in paths if p.name not in ['DELIVERY_QA_CHECKS.csv','WORKBOOK_SHEET_INDEX.csv','FILE_MANIFEST_VERIFICATION.csv']]
 records=[];data=[]
 for i,p in enumerate(paths,1):
  f=pd.read_csv(p);sheet=f'S{i:02d}_'+p.stem[:26];assert len(sheet)<=31
  records.append({'sheet':sheet,'source_csv':'results/'+p.name,'rows':len(f),'columns':len(f.columns),'description':p.stem})
  data.append((sheet,f))
 info=pd.DataFrame([
  ['研究入口','先读STAGE6_MAIN_REPORT_ZH.md；简版见DECISION_MEMO_ZH.md'],
  ['主发现','Stage5的PDGFRB/NPTX1/KLK6与34家族不改；本轮已知发现后的解释及有限验证'],
  ['直接覆盖','旧版MRM覆盖NPTX1两肽段、KLK6三肽段；PDGFRB无直接覆盖'],
  ['新增人群','95名不在Stage5索引571人中的ADNI受试者；仍是同队列，不是外部验证'],
  ['验证结论','五肽段新增人群关联均BY q>0.05；同日相对约束设计df不足，未作正式复现'],
  ['共识结论','34家族仅11个等效可评价；三档工作容差均0 BY等效；非显著不判共识'],
  ['独立性','2014 CSV和整合XLSX是一个实验家族；2019为另一实验；Emory两发布重叠不重复计数'],
  ['方法边界','任意误差相关下真实共同丰度与检测特异病理机制未识别；表位/原TMT肽段不可得'],
  ['模拟','1010个独立单元；3160000次独立模拟，P值版本不重复计数'],
  ['信息边界','无个体记录进入本总表；00_restricted留作本地复现，目录名不代表OS访问控制'],
  ['读表','保留英文机器字段名；INDEX链接全部结果表；字段释义见FIELD_GUIDE_ZH.md'],
  ['复核','INDEPENDENT_QA_SUMMARY.json；交付检查及散列另存文件，不循环嵌入工作簿']
 ],columns=['事项','说明'])
 with pd.ExcelWriter(ROOT/WORKBOOK,engine='openpyxl') as writer:
  info.to_excel(writer,sheet_name='README',index=False);pd.DataFrame(records).to_excel(writer,sheet_name='INDEX',index=False)
  for sheet,f in data:f.to_excel(writer,sheet_name=sheet,index=False)
  book=writer.book
  for ws in book:
   ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions;ws.sheet_view.showGridLines=False
   for cell in ws[1]:cell.fill=PatternFill('solid',fgColor='254A68');cell.font=Font(color='FFFFFF',bold=True);cell.alignment=Alignment(wrap_text=True,vertical='center')
   ws.row_dimensions[1].height=34
   for column in ws.iter_cols():
    width=min(58,max(11,max(len(str(c.value or '')) for c in list(column)[:200])+2));ws.column_dimensions[column[0].column_letter].width=width
    for cell in column[1:]:
     if isinstance(cell.value,str) and cell.value.startswith('='):cell.data_type='s'
     if isinstance(cell.value,float):cell.number_format='0.000000E+00' if 'p_' in str(column[0].value).lower() else '0.000000'
  ws=book['INDEX']
  for i,r in enumerate(records,2):ws.cell(i,1).hyperlink=f"#'{r['sheet']}'!A1";ws.cell(i,1).font=Font(color='1565A3',underline='single')
  book['README'].column_dimensions['A'].width=18;book['README'].column_dimensions['B'].width=120
  for row in book['README'].iter_rows(min_row=2):
   for cell in row:cell.alignment=Alignment(wrap_text=True,vertical='top')
   book['README'].row_dimensions[row[0].row].height=34
 csvout(RESULTS/'WORKBOOK_SHEET_INDEX.csv',records)
 save_json(RESULTS/'WORKBOOK_BUILD.json',{'filename':WORKBOOK,'result_sheets':len(records),'total_sheets':len(records)+2,'total_table_rows':sum(r['rows'] for r in records),'includes_no_individual_data':True})
 log(f'Aggregate workbook built: {len(records)+2} sheets; all sheets indexed to source CSV')

if __name__=='__main__':run()
