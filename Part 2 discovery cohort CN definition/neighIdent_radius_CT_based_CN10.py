# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 18:05:03 2021

@author: qiuxi
"""


###Python 3.8.8

import pandas as pd  ###v1.1.4
import numpy as np    ###v1.22.4
from sklearn.neighbors import NearestNeighbors   ###v0.22
import time
import sys
import matplotlib.pyplot as plt    ###v3.5.0
from sklearn.cluster import MiniBatchKMeans   ###v0.22
import seaborn as sns    ###v0.11.2
import random     

random.seed(123)

def get_windows(job,radius):
    '''
    For each region and each individual cell in dataset, return the indices of the nearest neighbors.
    
    'job:  meta data containing the start time,index of region, region name, indices of region in original dataframe
    n_neighbors:  the number of neighbors to find for each cell
    '''
    start_time,idx,tissue_name,indices = job
    job_start = time.time()
    
    print ("Starting:", str(idx+1)+'/'+str(len(exps)),': ' + exps[idx])

    tissue = tissue_group.get_group(tissue_name)
    to_fit = tissue.loc[indices][[X,Y]].values


#     fit = NearestNeighbors(n_neighbors=n_neighbors+1).fit(tissue[[X,Y]].values)
    fit = NearestNeighbors(radius=radius).fit(tissue[[X,Y]].values)
    m = fit.radius_neighbors(to_fit)
#     m = m[0][:,1:], m[1][:,1:]
    m = m[0], m[1]
    a= pd.DataFrame(m[1])
    a.columns = ["X"]
    a['N'] = a["X"].apply(lambda x : len(x))
    f=pd.DataFrame()
    nn=len(a)
    for i in range(nn):
        b = a["X"][i]
        b = tissue.iloc[b,:].index
        c= values[b]
        d=c.sum(axis=0)
        e=pd.DataFrame(d/d.sum(axis=0))#reg/cell finished
        aa = pd.DataFrame(m[0])
        aa.columns = ["CellID"]
        bb = aa["CellID"][i]
        cc = np.where(bb == 0)
        dd = b[cc]
        e.columns = dd
        e=pd.DataFrame(e.values.T, index=e.columns, columns=e.index)
        f=pd.concat([f,e])
    end_time = time.time()
   
    print ("Finishing:", str(idx+1)+"/"+str(len(exps)),": "+ exps[idx],end_time-job_start,end_time-start_time)
    return f

def get_windows_CN(job,radius):
    '''
    For each region and each individual cell in dataset, return the indices of the nearest neighbors.
    
    'job:  meta data containing the start time,index of region, region name, indices of region in original dataframe
    n_neighbors:  the number of neighbors to find for each cell
    '''
    start_time,idx,tissue_name,indices = job
    job_start = time.time()
    
    print ("Starting:", str(idx+1)+'/'+str(len(exps)),': ' + exps[idx])

    tissue = tissue_group.get_group(tissue_name)
    to_fit = tissue.loc[indices][[X,Y]].values


#     fit = NearestNeighbors(n_neighbors=n_neighbors+1).fit(tissue[[X,Y]].values)
    fit = NearestNeighbors(radius=radius).fit(tissue[[X,Y]].values)
    m = fit.radius_neighbors(to_fit)
#     m = m[0][:,1:], m[1][:,1:]
    m = m[0], m[1]
    a= pd.DataFrame(m[1])
    a.columns = ["X"]
    a['N'] = a["X"].apply(lambda x : len(x))
    nn=len(a)
    aa = pd.DataFrame(m[0])
    aa.columns = ["CellID"]
    dd = []
    for i in range(nn):
        if a["N"][i] > 5:
            b = a["X"][i]
            b = tissue.iloc[b,:].index
            bb = aa["CellID"][i]
            cc = np.where(bb == 0)
            dd += list(b[cc])
        else: continue
        
    end_time = time.time()
   
    print ("Finishing:", str(idx+1)+"/"+str(len(exps)),": "+ exps[idx],end_time-job_start,end_time-start_time)
    return dd


#ks = [5,10,20] # k=5 means it collects 5 nearest neighbors for each center cell
radius = 50
path_to_data = r'E:\myPythonProject\TrainingData\20220423 celltype CN\testing_allReg_seurat_._analysis_20220315_AllSubtypes_withoutartifact_celltype.csv'
X = 'XMin'
Y = 'YMin'
reg = 'Class'
file_type = 'csv'

cluster_col = 'celltype' #CODEX Clusters Nucleus Intensity
keep_cols = [X,Y,reg,cluster_col]
save_path = ''


#read in data and do some quick data rearrangement
#n_neighbors = max(ks)
assert (file_type=='csv' or file_type =='pickle') #


if file_type == 'pickle':
    cells = pd.read_pickle(path_to_data)
if file_type == 'csv':
    cells = pd.read_csv(path_to_data)

cells = pd.concat([cells,pd.get_dummies(cells[cluster_col])],1)  #将每个clster在后面单独作为一列


#cells = cells.reset_index() #Uncomment this line if you do any subsetting of dataframe such as removing dirt etc or will throw error at end of next next code block (cell 6)

sum_cols = cells[cluster_col].unique()  #计算多少个cluster
values = cells[sum_cols].values  #将每个细胞的clusterID放入相应的列中

#find windows for each cell in each tissue region
tissue_group = cells[[X,Y,reg]].groupby(reg)
exps = list(cells[reg].unique())
tissue_chunks = [(time.time(),exps.index(t),t,a) for t,indices in tissue_group.groups.items() for a in np.array_split(indices,1)] 
tissues = [get_windows(job,radius) for job in tissue_chunks]
tissues_a = [get_windows_CN(job,radius) for job in tissue_chunks]

tissues_aa = []
for i in range(288):
    tissues_aa += tissues_a[i]

windows2 = pd.concat(tissues[0:288])
windows2.columns = sum_cols
windows2 = pd.concat([cells[keep_cols],windows2],1)
windows2 = windows2.loc[tissues_aa,:]
cells = cells.loc[tissues_aa,:]


# k = 10
# n_neighborhoods =20#cluster
# neighborhood_name = "neighborhood"+str(k)
# k_centroids = {}

# km = MiniBatchKMeans(n_clusters = n_neighborhoods,random_state=0)

# labelskm = km.fit_predict(windows2[sum_cols].values)

# k_centroids[k] = km.cluster_centers_
# cells['neighborhood10'] = labelskm
# cells[neighborhood_name] = cells[neighborhood_name].astype('category')
# cells.to_csv("E:\\myPythonProject\\TrainingData\\20220308\\cells_r=50_CN=10.csv")

#a.to_csv("/Users/taozhou/Desktop/00_SMMU/Projects/Codex/Python_CN/a_1_reg001.csv")
#['reg064_A','reg066_A','reg018_B','reg023_A']

#cell_order = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
#              16,17,18,19,20,21,22,23]


# this plot shows the types of cells (ClusterIDs) in the different niches (0-7)
# k_to_plot = 10
# niche_clusters = (k_centroids[k_to_plot])
# tissue_avgs = values.mean(axis = 0)
# fc = np.log2(((niche_clusters+tissue_avgs)/(niche_clusters+tissue_avgs).sum(axis = 1, keepdims = True))/tissue_avgs)
# fc = pd.DataFrame(fc,columns = sum_cols)
# s=sns.clustermap(fc, vmin =-3,vmax = 3,cmap = 'bwr',row_cluster = False, figsize=(20,10))
# s.savefig("E:\\myPythonProject\\TrainingData\\20220308\\r=50_cn=20.pdf")
# s1=sns.clustermap(fc, vmin =-3,vmax = 3,cmap = 'bwr',row_cluster = True, figsize=(20,15))
# s1.savefig("E:\\myPythonProject\\TrainingData\\20220308\\r=50_cn=20_clustering.pdf")


# cells['neighborhood10'] = cells['neighborhood10'].astype('category')
# sns.lmplot(data = cells[cells['groups']==1],x = 'XMin',y='YMin',hue = 'neighborhood10',palette = 'bright',height = 8,col = reg,col_wrap = 10,fit_reg = False)

# cells['neighborhood10'] = cells['neighborhood10'].astype('category')
# sns.lmplot(data = cells[cells['groups']==2],x = 'X:X',y='Y:Y',hue = 'neighborhood10',palette = 'bright',height = 8,col = reg,col_wrap = 10,fit_reg = False)

k = 10
for i in range(5, 16, 1):
    n_neighborhoods = i #cluster
    neighborhood_name = "neighborhood"+str(k)
    k_centroids = {}
    km = MiniBatchKMeans(n_clusters = n_neighborhoods,random_state=0)
    labelskm = km.fit_predict(windows2[sum_cols].values)
    k_centroids[k] = km.cluster_centers_
    cells['neighborhood10'] = labelskm
    cells[neighborhood_name] = cells[neighborhood_name].astype('category')
    cells.to_csv("E:\\myPythonProject\\TrainingData\\20220423 celltype CN\\cells_r=50_CN="+str(i)+".csv")
    k_to_plot = 10
    niche_clusters = (k_centroids[k_to_plot])
    tissue_avgs = values.mean(axis = 0)
    fc = np.log2(((niche_clusters+tissue_avgs)/(niche_clusters+tissue_avgs).sum(axis = 1, keepdims = True))/tissue_avgs)
    fc = pd.DataFrame(fc,columns = sum_cols)
    s=sns.clustermap(fc, vmin =-3,vmax = 3,cmap = 'bwr',row_cluster = False, figsize=(8,8))
    s.savefig("E:\\myPythonProject\\TrainingData\\20220423 celltype CN\\r=50_CN="+str(i)+".pdf")
    s1=sns.clustermap(fc, vmin =-3,vmax = 3,cmap = 'bwr',row_cluster = True, figsize=(10,8))
    s1.savefig("E:\\myPythonProject\\TrainingData\\20220423 celltype CN\\r=50_CN="+str(i)+"clustering.pdf")



