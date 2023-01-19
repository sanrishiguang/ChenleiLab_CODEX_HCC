
###R 4.1.0

# generate celltype column for defining CT-based-CN in Python -------------------------------------------------


library(tidyverse)   #‘1.3.2’
library(stringr) # ‘1.4.0’
df <- read.csv(file="E:\\11. CODEX\\TrainingData\\2022\\testing_allReg_seurat_._analysis_20220315_AllSubtypes_withoutartifact.csv",header = T,row.names = 1,stringsAsFactors = F)

names(df)
unique(df$Class) ##100 samples
unique(df$TumorSubtype)
unique(df$StromalSubtype)

unique(df$Allsubtypes)
df_filter <- df


df_filter$celltype <- df_filter$StromalSubtype
df_filter$celltype[df_filter$StromalSubtype=='Tumor'] <- 'Tumor'

sort(unique(df_filter$StromalSubtype))
df_filter$celltype[df_filter$StromalSubtype=='CD4T_FOXP3+'] <- 'Treg'
df_filter$celltype[df_filter$StromalSubtype=='Macrophages_CD11c+_HLA-DR+'|df_filter$StromalSubtype=='Macrophages_HLA-DR+'] <- 'APC'
sort(unique(df_filter$celltype))

df_filter$celltype[str_detect(df_filter$celltype,'^CD4T')] <- "CD4T"
df_filter$celltype[str_detect(df_filter$celltype,'^CD8T')] <- "CD8T"
df_filter$celltype[str_detect(df_filter$celltype,'^Macrophage')] <- "Macrophages"

unique(df_filter$celltype)


write.csv(df_filter,file='E:\\myPythonProject\\TrainingData\\20220423 celltype CN\\testing_allReg_seurat_._analysis_20220315_AllSubtypes_withoutartifact_celltype.csv')

