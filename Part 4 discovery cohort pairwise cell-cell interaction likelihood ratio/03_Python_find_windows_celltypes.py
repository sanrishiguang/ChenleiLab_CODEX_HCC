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


def get_windows(job,n_neighbors):
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
    fit = NearestNeighbors(n_neighbors=n_neighbors).fit(tissue[[X,Y]].values)
    m = fit.kneighbors(to_fit)
#     m = m[0][:,1:], m[1][:,1:]
    m = m[0], m[1]
    

    #sort_neighbors
    args = m[0].argsort(axis = 1)
    add = np.arange(m[1].shape[0])*m[1].shape[1]
    sorted_indices = m[1].flatten()[args+add[:,None]]

    neighbors = tissue.index.values[sorted_indices]
   
    end_time = time.time()
   
    print ("Finishing:", str(idx+1)+"/"+str(len(exps)),": "+ exps[idx],end_time-job_start,end_time-start_time)
    return neighbors.astype(np.int32)


#######################################################################################
#####################################TP1##################################
#####################################################################################
#ks = [5,10,20] # k=5 means it collects 5 nearest neighbors for each center cell
ks = [5,10,15] # k=5 means it collects 5 nearest neighbors for each center cell
path_to_data = r'E:\myPythonProject\TrainingData\AllCN_r50_CN40_0331\spatialanalysis\data_SP_LII.csv'
#path_to_data = r'E:\myPythonProject\CODEX TLS 0830\TLS_seurat_python.csv'
X = 'XMin'   #对照表格修改
Y = 'YMin'    #对照表格修改
reg = 'Class'    #对照表格修改
file_type = 'csv'   #对照表格修改
cluster_col = 'celltype'      #对照表格修改
keep_cols = [X,Y,reg,cluster_col]
save_path = ''
#read in data and do some quick data rearrangement
n_neighbors = max(ks)
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
tissues = [get_windows(job,n_neighbors) for job in tissue_chunks]

#for each cell and its nearest neighbors, reshape and count the number of each cell type in those neighbors.
out_dict = {}
for k in ks:
    for neighbors,job in zip(tissues,tissue_chunks):

        chunk = np.arange(len(neighbors))#indices
        tissue_name = job[2]
        indices = job[3]
        window = values[neighbors[chunk,:k].flatten()].reshape(len(chunk),k,len(sum_cols)).sum(axis = 1)
        out_dict[(tissue_name,k)] = (window.astype(np.float16),indices)
        
#concatenate the summed windows and combine into one dataframe for each window size tested.
windows = {}
for k in ks:
   
    window = pd.concat([pd.DataFrame(out_dict[(exp,k)][0],index = out_dict[(exp,k)][1].astype(int),columns = sum_cols) for exp in exps],0)
    window = window.loc[cells.index.values]
    window = pd.concat([cells[keep_cols],window],1)
    windows[k] = window
#windows 作为字典，包含K=5， 10， 20 的数据，每一个都是一个dataframe
windows2 = windows[5]   #提取k=10的数据作为一个单独的数据框
# windows2[cluster_col] = cells[cluster_col]
##把细胞提出来，按细胞类型计算interaction的细胞个数： interaction count
arr = windows2.loc[:,"Class":]
group = arr.groupby(["Class","celltype"]).agg('sum')
group.shape
##根据interaction count计算logs Odds ratio, 并可视化
group.to_csv(r"E:\11. CODEX\TrainingData\20220308\AllCN_r50_CN40_0331\spatial analysis\spatialanalysis_output_celltype_SP_LII_class.csv")



#######################################################################################
#####################################TP2##################################
#####################################################################################
#ks = [5,10,20] # k=5 means it collects 5 nearest neighbors for each center cell
ks = [5,10,15] # k=5 means it collects 5 nearest neighbors for each center cell
path_to_data = r'E:\myPythonProject\TrainingData\AllCN_r50_CN40_0331\spatialanalysis\data_SP_HII.csv'
#path_to_data = r'E:\myPythonProject\CODEX TLS 0830\TLS_seurat_python.csv'
X = 'XMin'   #对照表格修改
Y = 'YMin'    #对照表格修改
reg = 'Class'    #对照表格修改
file_type = 'csv'   #对照表格修改
cluster_col = 'celltype'      #对照表格修改
keep_cols = [X,Y,reg,cluster_col]
save_path = ''
#read in data and do some quick data rearrangement
n_neighbors = max(ks)
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
tissues = [get_windows(job,n_neighbors) for job in tissue_chunks]

#for each cell and its nearest neighbors, reshape and count the number of each cell type in those neighbors.
out_dict = {}
for k in ks:
    for neighbors,job in zip(tissues,tissue_chunks):

        chunk = np.arange(len(neighbors))#indices
        tissue_name = job[2]
        indices = job[3]
        window = values[neighbors[chunk,:k].flatten()].reshape(len(chunk),k,len(sum_cols)).sum(axis = 1)
        out_dict[(tissue_name,k)] = (window.astype(np.float16),indices)
        
#concatenate the summed windows and combine into one dataframe for each window size tested.
windows = {}
for k in ks:
   
    window = pd.concat([pd.DataFrame(out_dict[(exp,k)][0],index = out_dict[(exp,k)][1].astype(int),columns = sum_cols) for exp in exps],0)
    window = window.loc[cells.index.values]
    window = pd.concat([cells[keep_cols],window],1)
    windows[k] = window
#windows 作为字典，包含K=5， 10， 20 的数据，每一个都是一个dataframe
windows2 = windows[5]   #提取k=10的数据作为一个单独的数据框
# windows2[cluster_col] = cells[cluster_col]
##把细胞提出来，按细胞类型计算interaction的细胞个数： interaction count
arr = windows2.loc[:,"Class":]
group = arr.groupby(["Class","celltype"]).agg('sum')
group.shape
##根据interaction count计算logs Odds ratio, 并可视化
group.to_csv(r"E:\11. CODEX\TrainingData\20220308\AllCN_r50_CN40_0331\spatial analysis\spatialanalysis_output_celltype_SP_HII_class.csv")



#######################################################################################
#####################################TP3##################################
#####################################################################################
#ks = [5,10,20] # k=5 means it collects 5 nearest neighbors for each center cell
ks = [5,10,15] # k=5 means it collects 5 nearest neighbors for each center cell
path_to_data = r'E:\myPythonProject\TrainingData\AllCN_r50_CN40_0331\spatialanalysis\data_SP_Pf.csv'
#path_to_data = r'E:\myPythonProject\CODEX TLS 0830\TLS_seurat_python.csv'
X = 'XMin'   #对照表格修改
Y = 'YMin'    #对照表格修改
reg = 'Class'    #对照表格修改
file_type = 'csv'   #对照表格修改
cluster_col = 'celltype'      #对照表格修改
keep_cols = [X,Y,reg,cluster_col]
save_path = ''
#read in data and do some quick data rearrangement
n_neighbors = max(ks)
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
tissues = [get_windows(job,n_neighbors) for job in tissue_chunks]

#for each cell and its nearest neighbors, reshape and count the number of each cell type in those neighbors.
out_dict = {}
for k in ks:
    for neighbors,job in zip(tissues,tissue_chunks):

        chunk = np.arange(len(neighbors))#indices
        tissue_name = job[2]
        indices = job[3]
        window = values[neighbors[chunk,:k].flatten()].reshape(len(chunk),k,len(sum_cols)).sum(axis = 1)
        out_dict[(tissue_name,k)] = (window.astype(np.float16),indices)
        
#concatenate the summed windows and combine into one dataframe for each window size tested.
windows = {}
for k in ks:
   
    window = pd.concat([pd.DataFrame(out_dict[(exp,k)][0],index = out_dict[(exp,k)][1].astype(int),columns = sum_cols) for exp in exps],0)
    window = window.loc[cells.index.values]
    window = pd.concat([cells[keep_cols],window],1)
    windows[k] = window
#windows 作为字典，包含K=5， 10， 20 的数据，每一个都是一个dataframe
windows2 = windows[5]   #提取k=10的数据作为一个单独的数据框
# windows2[cluster_col] = cells[cluster_col]
##把细胞提出来，按细胞类型计算interaction的细胞个数： interaction count
arr = windows2.loc[:,"Class":]
group = arr.groupby(["Class","celltype"]).agg('sum')
group.shape
##根据interaction count计算logs Odds ratio, 并可视化
group.to_csv(r"E:\11. CODEX\TrainingData\20220308\AllCN_r50_CN40_0331\spatial analysis\spatialanalysis_output_celltype_SP_Pf_class.csv")


#######################################################################################
#####################################TP4##################################
#####################################################################################
#ks = [5,10,20] # k=5 means it collects 5 nearest neighbors for each center cell
ks = [5,10,15] # k=5 means it collects 5 nearest neighbors for each center cell
path_to_data = r'E:\myPythonProject\TrainingData\AllCN_r50_CN40_0331\spatialanalysis\data_SP_MII_Tumor.csv'
#path_to_data = r'E:\myPythonProject\CODEX TLS 0830\TLS_seurat_python.csv'
X = 'XMin'   #对照表格修改
Y = 'YMin'    #对照表格修改
reg = 'Class'    #对照表格修改
file_type = 'csv'   #对照表格修改
cluster_col = 'celltype'      #对照表格修改
keep_cols = [X,Y,reg,cluster_col]
save_path = ''
#read in data and do some quick data rearrangement
n_neighbors = max(ks)
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
tissues = [get_windows(job,n_neighbors) for job in tissue_chunks]

#for each cell and its nearest neighbors, reshape and count the number of each cell type in those neighbors.
out_dict = {}
for k in ks:
    for neighbors,job in zip(tissues,tissue_chunks):

        chunk = np.arange(len(neighbors))#indices
        tissue_name = job[2]
        indices = job[3]
        window = values[neighbors[chunk,:k].flatten()].reshape(len(chunk),k,len(sum_cols)).sum(axis = 1)
        out_dict[(tissue_name,k)] = (window.astype(np.float16),indices)
        
#concatenate the summed windows and combine into one dataframe for each window size tested.
windows = {}
for k in ks:
   
    window = pd.concat([pd.DataFrame(out_dict[(exp,k)][0],index = out_dict[(exp,k)][1].astype(int),columns = sum_cols) for exp in exps],0)
    window = window.loc[cells.index.values]
    window = pd.concat([cells[keep_cols],window],1)
    windows[k] = window
#windows 作为字典，包含K=5， 10， 20 的数据，每一个都是一个dataframe
windows2 = windows[5]   #提取k=10的数据作为一个单独的数据框
# windows2[cluster_col] = cells[cluster_col]
##把细胞提出来，按细胞类型计算interaction的细胞个数： interaction count
arr = windows2.loc[:,"Class":]
group = arr.groupby(["Class","celltype"]).agg('sum')
group.shape
##根据interaction count计算logs Odds ratio, 并可视化
group.to_csv(r"E:\11. CODEX\TrainingData\20220308\AllCN_r50_CN40_0331\spatial analysis\spatialanalysis_output_celltype_SP_MII_Tumor_class.csv")



#######################################################################################
#####################################TP5##################################
#####################################################################################
#ks = [5,10,20] # k=5 means it collects 5 nearest neighbors for each center cell
ks = [5,10,15] # k=5 means it collects 5 nearest neighbors for each center cell
path_to_data = r'E:\myPythonProject\TrainingData\AllCN_r50_CN40_0331\spatialanalysis\data_SP_MII_Fibro.csv'
#path_to_data = r'E:\myPythonProject\CODEX TLS 0830\TLS_seurat_python.csv'
X = 'XMin'   #对照表格修改
Y = 'YMin'    #对照表格修改
reg = 'Class'    #对照表格修改
file_type = 'csv'   #对照表格修改
cluster_col = 'celltype'      #对照表格修改
keep_cols = [X,Y,reg,cluster_col]
save_path = ''
#read in data and do some quick data rearrangement
n_neighbors = max(ks)
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
tissues = [get_windows(job,n_neighbors) for job in tissue_chunks]

#for each cell and its nearest neighbors, reshape and count the number of each cell type in those neighbors.
out_dict = {}
for k in ks:
    for neighbors,job in zip(tissues,tissue_chunks):

        chunk = np.arange(len(neighbors))#indices
        tissue_name = job[2]
        indices = job[3]
        window = values[neighbors[chunk,:k].flatten()].reshape(len(chunk),k,len(sum_cols)).sum(axis = 1)
        out_dict[(tissue_name,k)] = (window.astype(np.float16),indices)
        
#concatenate the summed windows and combine into one dataframe for each window size tested.
windows = {}
for k in ks:
   
    window = pd.concat([pd.DataFrame(out_dict[(exp,k)][0],index = out_dict[(exp,k)][1].astype(int),columns = sum_cols) for exp in exps],0)
    window = window.loc[cells.index.values]
    window = pd.concat([cells[keep_cols],window],1)
    windows[k] = window
#windows 作为字典，包含K=5， 10， 20 的数据，每一个都是一个dataframe
windows2 = windows[5]   #提取k=10的数据作为一个单独的数据框
# windows2[cluster_col] = cells[cluster_col]
##把细胞提出来，按细胞类型计算interaction的细胞个数： interaction count
arr = windows2.loc[:,"Class":]
group = arr.groupby(["Class","celltype"]).agg('sum')
group.shape
##根据interaction count计算logs Odds ratio, 并可视化
group.to_csv(r"E:\11. CODEX\TrainingData\20220308\AllCN_r50_CN40_0331\spatial analysis\spatialanalysis_output_celltype_SP_MII_Fibro_class.csv")




######总的########
ks = [5,10,15] # k=5 means it collects 5 nearest neighbors for each center cell
path_to_data = r'E:\11. CODEX\TrainingData\20220308\AllCN_r50_CN40_0331\All_CN_TP_0424.csv'
#path_to_data = r'E:\myPythonProject\CODEX TLS 0830\TLS_seurat_python.csv'
X = 'XMin'   #对照表格修改
Y = 'YMin'    #对照表格修改
reg = 'Class'    #对照表格修改
file_type = 'csv'   #对照表格修改
cluster_col = 'celltype'      #对照表格修改
keep_cols = [X,Y,reg,cluster_col]
save_path = ''
#read in data and do some quick data rearrangement
n_neighbors = max(ks)
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
tissues = [get_windows(job,n_neighbors) for job in tissue_chunks]

#for each cell and its nearest neighbors, reshape and count the number of each cell type in those neighbors.
out_dict = {}
for k in ks:
    for neighbors,job in zip(tissues,tissue_chunks):

        chunk = np.arange(len(neighbors))#indices
        tissue_name = job[2]
        indices = job[3]
        window = values[neighbors[chunk,:k].flatten()].reshape(len(chunk),k,len(sum_cols)).sum(axis = 1)
        out_dict[(tissue_name,k)] = (window.astype(np.float16),indices)
        
#concatenate the summed windows and combine into one dataframe for each window size tested.
windows = {}
for k in ks:
   
    window = pd.concat([pd.DataFrame(out_dict[(exp,k)][0],index = out_dict[(exp,k)][1].astype(int),columns = sum_cols) for exp in exps],0)
    window = window.loc[cells.index.values]
    window = pd.concat([cells[keep_cols],window],1)
    windows[k] = window
#windows 作为字典，包含K=5， 10， 20 的数据，每一个都是一个dataframe
windows2 = windows[5]   #提取k=10的数据作为一个单独的数据框
# windows2[cluster_col] = cells[cluster_col]
##把细胞提出来，按细胞类型计算interaction的细胞个数： interaction count
arr = windows2.loc[:,"Class":]
group = arr.groupby(["Class","celltype"]).agg('sum')
group.shape
##根据interaction count计算logs Odds ratio, 并可视化
group.to_csv(r"E:\11. CODEX\TrainingData\20220308\AllCN_r50_CN40_0331\spatial analysis\spatialanalysis_output_celltype_class.csv")





