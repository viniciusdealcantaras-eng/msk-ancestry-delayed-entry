"""Build the integrated MSK-IMPACT 50K + MSK-CHORD patient-level dataset.

Input : raw cBioPortal files in data/msk_impact_50k_2026/ and data/msk_chord_2024/
Output: out/banco_ancestralidade_MSK.csv (local only, git-ignored, never redistributed:
        both studies are licensed CC BY-NC-ND 4.0)
"""
import os, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
R50,RC,OUT,MO='data/msk_impact_50k_2026/','data/msk_chord_2024/','out/',30.4375
os.makedirs(OUT,exist_ok=True)

# ---------- MSK-IMPACT 50K ----------
p=pd.read_csv(R50+'data_clinical_patient.txt',sep='\t',comment='#')
s=pd.read_csv(R50+'data_clinical_sample_50k.tsv',sep='\t')
idx=s.sort_values('SAMPLE_ID').groupby('PATIENT_ID').first().reset_index()
nsamp=s.groupby('PATIENT_ID').size()

f=pd.DataFrame({'patient_id':p.PATIENT_ID})
f['ancestry']=p.ANCESTRY_LABEL.values
f['sex']=p.SEX.values
f['age_dx']=pd.to_numeric(p.AGE_AT_DX,errors='coerce').values
f['n_other_cancer_types']=pd.to_numeric(p.N_OTHER_CANCER_TYPES,errors='coerce').values
f['os_status']=p.OS_STATUS.map({'DECEASED':1,'LIVING':0}).values
f['os_months_seq']=pd.to_numeric(p.OS_MONTHS,errors='coerce').values
f['n_samples']=f.patient_id.map(nsamp).fillna(0).astype(int)
ix=idx.set_index('PATIENT_ID')
cols={'SAMPLE_ID':'index_sample_id','CANCER_TYPE':'cancer_type','CANCER_TYPE_DETAILED':'cancer_type_detailed',
 'ONCOTREE_CODE':'oncotree_code','SAMPLE_TYPE':'sample_type','DISEASE_STATUS':'disease_status','PRIMARY_SITE':'primary_site',
 'METASTATIC_SITE':'metastatic_site','GENE_PANEL':'gene_panel','TMB_SCORE':'tmb','MUTATION_COUNT':'mutation_count',
 'MSI_SCORE':'msi_score','MSI_TYPE':'msi_type','FRACTION_GENOME_ALTERED':'fga','FACETS_WGD':'wgd','FACETS_QC':'facets_qc',
 'FACETS_PURITY':'facets_purity','FACETS_PLOIDY':'facets_ploidy','TUMOR_PURITY':'tumor_purity_path','SAMPLE_COVERAGE':'sample_coverage'}
for a,b in cols.items(): f[b]=f.patient_id.map(ix[a])
f['wgd']=f.wgd.map({'TRUE':1,'FALSE':0,True:1,False:0})
f['msi_high']=np.where(f.msi_type.eq('Instable'),1,np.where(f.msi_type.eq('Stable'),0,np.nan))
f['tmb_high10']=np.where(f.tmb.notna(),(f.tmb>=10).astype(float),np.nan)

# ---------- alteracoes condutoras na amostra indice ----------
isamp=set(f.index_sample_id.dropna())
m=pd.read_csv(R50+'mutations_driver_genes_50k.tsv',sep='\t',dtype=str)
m=m[m.sampleId.isin(isamp)]
keep=~m.mutationType.isin(["5'Flank",'Splice_Region'])|((m.hugo=='TERT')&(m.mutationType=="5'Flank"))
m=m[keep]
genes=sorted(m.hugo.unique())
mut=pd.crosstab(m.sampleId,m.hugo).clip(upper=1)
cna_amp=['ERBB2','EGFR','MET','MYC','CCND1','CDK4','MDM2','FGFR2','KRAS','AR']
cna_del=['CDKN2A','PTEN','RB1','SMAD4','BRCA2','BRCA1','ATM']
cn=pd.read_csv(R50+'data_cna.txt',sep='\t',index_col=0)
cn=cn.loc[cn.index.isin(cna_amp+cna_del)]
cn=cn[[c for c in cn.columns if c in isamp]].T
sv=pd.read_csv(R50+'structural_variants_50k.tsv',sep='\t',dtype=str)
sv=sv[sv.sampleId.isin(isamp)]
fus=['ALK','ROS1','RET','NTRK1','NTRK2','NTRK3','FGFR2','FGFR3','ERG','BRAF','MET','NRG1','EGFR']
svl=pd.concat([sv[['sampleId','site1HugoSymbol']].rename(columns={'site1HugoSymbol':'g'}),sv[['sampleId','site2HugoSymbol']].rename(columns={'site2HugoSymbol':'g'})])
svl=svl[svl.g.isin(fus)].drop_duplicates()
svx=pd.crosstab(svl.sampleId,svl.g).clip(upper=1)
sid=f.index_sample_id
for g in genes: f['mut_'+g]=sid.map(mut[g]).fillna(0).astype(int).where(sid.notna())
for g in cna_amp:
    if g in cn: f['amp_'+g]=np.where(sid.isin(cn.index),(sid.map(cn[g])==2).astype(float),np.nan)
for g in cna_del:
    if g in cn: f['del_'+g]=np.where(sid.isin(cn.index),(sid.map(cn[g])==-2).astype(float),np.nan)
for g in fus: f['sv_'+g]=sid.map(svx[g]).fillna(0).astype(int).where(sid.notna()) if g in svx else np.where(sid.notna(),0,np.nan)

# ---------- MSK-CHORD ----------
pc=pd.read_csv(RC+'data_clinical_patient.txt',sep='\t',comment='#')
sc=pd.read_csv(RC+'data_clinical_sample.txt',sep='\t',comment='#')
dx=pd.read_csv(RC+'data_timeline_diagnosis.txt',sep='\t',dtype=str)
dx['START_DATE']=pd.to_numeric(dx.START_DATE,errors='coerce')
first=dx.sort_values(['PATIENT_ID','START_DATE']).groupby('PATIENT_ID').first()
ctc=sc.sort_values('SAMPLE_ID').groupby('PATIENT_ID')[['SAMPLE_ID','CANCER_TYPE']].first()
stage={'Stage 1-3':'I-III','Stage 4':'IV'}
c=pd.DataFrame({'patient_id':pc.PATIENT_ID})
c['chord_race']=pc.RACE.values; c['chord_ethnicity']=pc.ETHNICITY.values
c['chord_index_sample_id']=c.patient_id.map(ctc.SAMPLE_ID); c['chord_cancer_type']=c.patient_id.map(ctc.CANCER_TYPE)
c['chord_stage_dx']=c.patient_id.map(first.STAGE_CDM_DERIVED).map(stage)
c['chord_stage_highest']=pc.STAGE_HIGHEST_RECORDED.map(stage).values
c['chord_os_status']=pc.OS_STATUS.astype(str).str.startswith('1').astype(int).values
c['chord_os_months_seq']=pd.to_numeric(pc.OS_MONTHS,errors='coerce').values
off=c.patient_id.map(first.START_DATE)
c['chord_dx_valid']=(off.notna()&(off<=0)).astype(int)
c['dx_to_seq_months']=np.where(c.chord_dx_valid==1,-off/MO,np.nan)
c['os_months_dx']=c.dx_to_seq_months+c.chord_os_months_seq

f=f.merge(c,on='patient_id',how='left')
f['in_chord']=f.chord_race.notna().astype(int)|f.chord_os_status.notna().astype(int)
f['in_chord']=f.patient_id.isin(set(c.patient_id)).astype(int)
f['same_index_sample']=np.where(f.in_chord==1,(f.index_sample_id==f.chord_index_sample_id).astype(float),np.nan)
order={'nonASJ-EUR':1,'ASJ-EUR':2,'ADMIX_OTHER':3,'AFR':4,'EAS':5,'SAS':6,'NAM':7}
f['ancestry_cat']=f.ancestry.map(order)
f.to_csv(OUT+'banco_ancestralidade_MSK.csv',index=False)
print(f.shape, f.in_chord.sum(), int(((f.in_chord==1)&(f.chord_dx_valid==1)).sum()))
print(f.ancestry.value_counts(dropna=False).to_dict())
