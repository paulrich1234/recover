import xlrd
import pandas as pd
import numpy as np
import copy
import seaborn as sns  ; sns.set()
import matplotlib.pyplot as plt



excel_name1 = 'mmc7.xlsx'
excel_name2 = 'mmc5.xlsx'
excel_name3 = 'mmc4.xlsx'
wb = xlrd.open_workbook(excel_name3)



df3 = pd.read_excel(excel_name1, index=False, encoding='utf8')
df3=df3.fillna(0)

# pathway_average3=df3.sum(axis=0)/df3.shape[0]
# print(pathway_average3)

df1 = pd.read_excel(excel_name1, index=False, encoding='utf8')
df1=df1.fillna(0)

pathway_average=df1.sum(axis=0)/df1.shape[0]
# print(pathway_average)
df2 = pd.read_excel(excel_name2, index=True, skiprows=2,iencoding='utf8')
Subtype_Sample = {}
for row in df2.iterrows():
    if row[1].SUBTYPE == 'Not_Applicable':
        if row[1].DISEASE not in Subtype_Sample.keys():
            Subtype_Sample.setdefault(row[1].DISEASE,[]).append(row[1].SAMPLE_BARCODE)
        else:
            Subtype_Sample[row[1].DISEASE].append(row[1].SAMPLE_BARCODE)
    else:
        SUBTYPE=str(row[1].DISEASE)+' '+str(row[1].SUBTYPE)
        if SUBTYPE not in Subtype_Sample.keys():
            Subtype_Sample.setdefault(SUBTYPE,[]).append(row[1].SAMPLE_BARCODE)
        else:
            Subtype_Sample[SUBTYPE].append(row[1].SAMPLE_BARCODE)


STES_CIN=[]
STES_GS=[]
STES_EBV=[]
STES_MSI_POLE=[]
STES_Squamous=[]
CRC_MSI_POLE=[]
CRC_GS=[]
CRC_CIN=[]
GBM = []
# for key ,value in Subtype_Sample.items():
#         print(key)
print('8'*100)


new_dict=copy.deepcopy(Subtype_Sample)
for key ,value in Subtype_Sample.items():
    if key=='ESCA CIN' or key=='STAD CIN':
        STES_CIN+=new_dict.pop(key)
    elif key=='ESCA GS' or key=='STAD GS':
        STES_GS+=new_dict.pop(key)
    elif 'STAD EBV' == key:
        STES_EBV+=new_dict.pop(key)
    elif 'STAD MSI' == key or 'ESCA POLE' == key or 'ESCA MSI' == key or 'STAD POLE' == key:
        STES_MSI_POLE+=new_dict.pop(key)
    elif 'ESCA ESCC' == key:
        STES_Squamous += new_dict.pop(key)
    else:
        pass

for key ,value in Subtype_Sample.items():
    if key=='COAD GS' or key=='READ GS':
        CRC_GS+=new_dict.pop(key)
    elif key=='READ CIN' or key=='COAD CIN':
        CRC_CIN+=new_dict.pop(key)
    elif 'READ MSI' == key or 'COAD MSI'==key or 'READ POLE' == key or 'COAD POLE' == key:
        CRC_MSI_POLE +=new_dict.pop(key)
    elif 'GBM nan'==key or 'GBM IDHwt' ==key or 'GBM IDHmut-non-codel'==key:
        GBM+=new_dict.pop(key)
    else:
        pass

new_dict['STES_CIN']=STES_CIN
new_dict['STES_GS']=STES_GS
new_dict['STES_EBV']=STES_EBV
new_dict['STES_MSI_POLE']=STES_MSI_POLE
new_dict['STES_Squamous']=STES_Squamous
new_dict['CRC_MSI_POLE']=CRC_MSI_POLE
new_dict['CRC_GS']=CRC_GS
new_dict['CRC_CIN']=CRC_CIN
new_dict['GBM']=GBM




for pathway in df1.columns:
    print(pathway)
Disease_Pathway2 =pd.DataFrame(index=new_dict.keys(),columns=list(df1.columns),data=np.zeros(shape=(len(new_dict.keys()),len(list(df1.columns)))))
for key ,values in new_dict.items():
    for value in values:
        for pathway in df1.columns:
            if df1.loc[value,pathway]==1:
                Disease_Pathway2.loc[key,pathway]+=1

Disease_Pathway2=Disease_Pathway2.loc[['GBM','LGG IDHwt','LGG IDHmut-non-codel','LGG IDHmut-codel','UVM','HNSC HPV+','HNSC HPV-','THCA','ACC','PCPG','THYM',
                                      'LUAD','MESO','LUSC','BRCA LumA','BRCA LumB','BRCA Her2','BRCA Basal','BRCA Normal','STES_Squamous',
                                     'STES_CIN','STES_EBV','STES_GS','STES_MSI_POLE',
                                     'CRC_MSI_POLE','CRC_GS','CRC_CIN',
                                     'LIHC','CHOL','PAAD','KIRC','KIRP','KICH',
                                     'BLCA','PRAD','TGCT seminoma','TGCT non-seminoma','OV',
                                     'UCEC CN_HIGH','UCEC CN_LOW','UCEC POLE','UCS','CESC AdenoCarcinoma','CESC SquamousCarcinoma',
                                     'SKCM','SARC DDLPS','SARC LMS','SARC MFS/UPS','SARC Other','DLBC','LAML']]


# down=Disease_Pathway2.sum(axis=0)/9125
# print(down)



Disease_Pathway =pd.DataFrame(index=new_dict.keys(),columns=list(df1.columns),data=np.zeros(shape=(len(new_dict.keys()),len(list(df1.columns)))))
for key ,values in new_dict.items():
    for value in values:
        for pathway in df1.columns:
            if df1.loc[value,pathway]==1:
                Disease_Pathway.loc[key,pathway]+=1
    Disease_Pathway.loc[key,:]=Disease_Pathway.loc[key,:]/len(values)

Disease_Pathway=Disease_Pathway*100+0.5
Disease_Pathway=Disease_Pathway.astype('int')
reset_columns = ['RTK RAS','Cell Cycle','PI3K','TP53','NOTCH','WNT','MYC','HIPPO','TGF-Beta','NRF2']
Disease_Pathway=Disease_Pathway[reset_columns] # STES=ESCA+STAD  CRC = COAD +READ
print(Disease_Pathway)
Disease_Pathway.to_csv('result22.csv')
# STES_Squamous = Disease_Pathway.loc[''].sum(axis=0)


Disease_Pathway=Disease_Pathway.loc[['GBM','LGG IDHwt','LGG IDHmut-non-codel','LGG IDHmut-codel','UVM','HNSC HPV+','HNSC HPV-','THCA','ACC','PCPG','THYM',
                                      'LUAD','MESO','LUSC','BRCA LumA','BRCA LumB','BRCA Her2','BRCA Basal','BRCA Normal','STES_Squamous',
                                     'STES_CIN','STES_EBV','STES_GS','STES_MSI_POLE',
                                     'CRC_MSI_POLE','CRC_GS','CRC_CIN',
                                     'LIHC','CHOL','PAAD','KIRC','KIRP','KICH',
                                     'BLCA','PRAD','TGCT seminoma','TGCT non-seminoma','OV',
                                     'UCEC CN_HIGH','UCEC CN_LOW','UCEC POLE','UCS','CESC AdenoCarcinoma','CESC SquamousCarcinoma',
                                     'SKCM','SARC DDLPS','SARC LMS','SARC MFS/UPS','SARC Other','DLBC','LAML']]




# Disease_Pathway =Disease_Pathway.drop(['READ POLE','COAD POLE','READ CIN','PRAD','COAD GS','ESCA GS','STAD GS','ESCA CIN'],axis=0)


cmap = sns.cubehelix_palette(start = 1.5, rot = 1, gamma=0.8, as_cmap = True)
print(Disease_Pathway)
ax = sns.heatmap(Disease_Pathway,annot=True, fmt="d",cbar=0,annot_kws={'size':8,'weight':'bold', 'color':'black'},linewidths=0.1,cmap="YlGnBu",xticklabels=True, yticklabels=True)
ax.set_title('Alteration frequencies',fontsize=18,loc='left')
plt.show()
# Disease_Pathway.to_csv('result33.csv')










