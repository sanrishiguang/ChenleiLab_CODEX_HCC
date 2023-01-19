

###R 4.1.0
rm(list=ls())
gc()

library(tidyverse)  #‘1.3.2’
library(tidyr) #‘1.2.0’
library(ggpubr) # ‘0.4.0’
library(ggthemes) #‘4.2.4’
library(Rtsne)  # ‘0.16’
library(Rcpp)  #‘1.0.8.3’
library(igraph)  #‘1.3.2’
library(ggthemes)  # ‘4.2.4’
library(survival)   #‘3.3.1’
library(survminer)   # ‘0.4.9’
library(pheatmap)  # ‘1.0.12’
library(ggsci)  # ‘2.9’
library(reshape2) #‘1.4.4’

##### 1: data after definition 
library(tidyverse)
library(stringr)
df <- read.csv(file="testing_allReg_seurat_._analysis_20220411_AllSubtypes.csv",header = T,row.names = 1,stringsAsFactors = F)

names(df)
unique(df$Class) 
unique(df$TumorSubtype)
unique(df$StromalSubtype)

#remove aritifact
df_filter <- filter(df, StromalSubtype != 'CD4T_artifact')

#create a column: Allsubtypes
df_filter$Allsubtypes <- df_filter$StromalSubtype
df_filter$Allsubtypes[df_filter$Allsubtypes=='Tumor'] <- df_filter$TumorSubtype[df_filter$Allsubtypes=='Tumor']

unique(df_filter$Allsubtypes)

#create a column: celltype
df_filter$celltype <- df_filter$StromalSubtype
df_filter$celltype[df_filter$StromalSubtype=='Tumor'] <- 'Tumor'

sort(unique(df_filter$StromalSubtype))
df_filter$celltype[df_filter$StromalSubtype=='CD4T_FOXP3+'] <- 'Treg'
df_filter$celltype[df_filter$StromalSubtype=='Macrophages_CD11c+_HLA-DR+'|df_filter$StromalSubtype=='Macrophages_HLA-DR+'] <- 'APC'
sort(unique(df_filter$celltype))

df_filter$celltype[str_detect(df_filter$celltype,'^CD4T')] <- "CD4T"
df_filter$celltype[str_detect(df_filter$celltype,'^CD8T')] <- "CD8T"
df_filter$celltype[str_detect(df_filter$celltype,'^Macrophage')] <- "Macrophages"

##check
sort(unique(df_filter$Allsubtypes))
sort(unique(df_filter$celltype))

write.csv(df_filter,file='testing_allReg_seurat_._analysis_20220412_AllSubtypes_celltype.csv')  ##as input for CN definition in python
