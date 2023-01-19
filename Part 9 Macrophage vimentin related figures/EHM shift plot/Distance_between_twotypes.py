# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 13:53:40 2022

@author:    QiuXinyao & LiShuai
"""
#Python 3.9.12 
import numpy as np   #1.22.4
import pandas as pd   #1.4.2
from numpy import *   #1.22.4
###Input data
table=pd.read_csv(r'B:\codex_density\All_CN_TP_0424 metastasis.csv')
a=table.Allsubtypes.value_counts()  
cn=a.index
final1=[]
for i in cn:
    if "Tumor" in i:
        cc=1
        final1.append(i)
name=['CD4T_FOXP3+','Macrophages_Vimentin+']

###The slices' name of different group
P_classname=table[(table['Metastasis'])=='P']['Class'].drop_duplicates()
P_LM_classname=table[(table['Metastasis'])=='P_LM']['Class'].drop_duplicates()

###filter the slices which do not have two celltypes in one slice.
s=list()
for i in P_classname:
    select_table=table[(table['Class']==i)]
    select_cell1 = select_table[(select_table['Allsubtypes']==name[0])]
    select_cell2 = select_table[(select_table['Allsubtypes']==name[1])]
    a=select_cell1.empty
    b=select_cell2.empty
    if a==True or b==True:
        c=True
    else:
        c=False
    s.append(c)

opposite=list()
for i in s:
    if i == False:
        a=0
        opposite.append(True)
    if i == True:
        opposite.append(False)
P_classname_select=P_classname[opposite].reset_index()

###Caculate the distance bewteen two celltypes
oneslice_final=list()
for i in range(0,P_classname_select.shape[0]):
    select_table=table[(table['Class']==P_classname_select['Class'][i])]
    select_cell = select_table[(select_table['Allsubtypes']==name[0])|(select_table['Allsubtypes']==name[1])]

    select_cell1= select_cell[(select_cell['Allsubtypes']==name[0])][['XMin','YMin']]
    select_cell2= select_cell[(select_cell['Allsubtypes']==name[1])][['XMin','YMin']]

    select_cell1=select_cell1.reset_index(drop=True)
    select_cell2=select_cell2.reset_index(drop=True)
    for i in range(0,select_cell1.shape[0]):
        core_cell=select_cell1.loc[i,].T
        core_cell_distance=tile(core_cell,(select_cell2.shape[0],1))-select_cell2
        sqDiffMat=core_cell_distance**2
        sqDistances=sqDiffMat.sum(axis=1)
        #sqDistances=sqDiffMat**0.5
        finalDistances= sqDistances.min()
        oneslice_final.append(finalDistances)
    
oneslice_final= np.asarray(oneslice_final)
oneslice_final=oneslice_final**0.5

###Repeat to caculate another group

s=list()
for i in P_LM_classname:
    select_table=table[(table['Class']==i)]
    select_cell1 = select_table[(select_table['Allsubtypes']==name[0])]
    select_cell2 = select_table[(select_table['Allsubtypes']==name[1])]
    a=select_cell1.empty
    b=select_cell2.empty
    if a==True or b==True:
        c=True
    else:
        c=False
    s.append(c)

opposite=list()
for i in s:
    if i == False:
        a=0
        opposite.append(True)
    if i == True:
        opposite.append(False)
P_LM_classname_select=P_LM_classname[opposite].reset_index()


oneslice_LM_final=list()
for i in range(0,P_LM_classname_select.shape[0]):
    select_table=table[(table['Class']==P_LM_classname_select['Class'][i])]
    select_cell = select_table[(select_table['Allsubtypes']==name[0])|(select_table['Allsubtypes']==name[1])]

    select_cell1= select_cell[(select_cell['Allsubtypes']==name[0])][['XMin','YMin']]
    select_cell2= select_cell[(select_cell['Allsubtypes']==name[1])][['XMin','YMin']]

    select_cell1=select_cell1.reset_index(drop=True)
    select_cell2=select_cell2.reset_index(drop=True)
    for i in range(0,select_cell1.shape[0]):
        core_cell=select_cell1.loc[i,].T
        core_cell_distance=tile(core_cell,(select_cell2.shape[0],1))-select_cell2
        sqDiffMat=core_cell_distance**2
        sqDistances=sqDiffMat.sum(axis=1)
        #sqDistances=sqDiffMat**0.5
        finalDistances= sqDistances.min()
        oneslice_LM_final.append(finalDistances)

oneslice_LM_final= np.asarray(oneslice_LM_final)
oneslice_LM_final=oneslice_LM_final**0.5

path1='B:\\codex_density\\oneslice_LM_final_'
cell_path=name[1]
path_total1=path1+cell_path+'.csv'

path2='B:\codex_density\\oneslice_final_'
path_total2=path2+cell_path+'.csv'

pd.DataFrame(oneslice_LM_final).to_csv(path_total1,index=False)
pd.DataFrame(oneslice_final).to_csv(path_total2,index=False)







