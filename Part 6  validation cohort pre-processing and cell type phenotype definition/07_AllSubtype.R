# R Version: 3.6.3
# rm(list=ls())

library(Seurat) # Version: 3.2.3
library(tidyverse)
library(patchwork)
library(ggplot2)
library(harmony)
library(dplyr)
library(circlize)
library(RColorBrewer)
library(pheatmap)
library(data.table)
library(reshape)
library(tidyr)
library(stringr)
library(Rphenograph)
library(ggpubr)
library(flowCore)
library(Rcpp)
library(cytofkit)
library(igraph)
library(ggthemes)
library(Rtsne)
library(cytofexplorer)
library(survival)
library(survminer)
library(pvclust)
library(ComplexHeatmap)
library(ggbreak)
library(sva)

#### load all subtypes ####
load("~/OneDrive/01_Project/01_Codex/Bioinfo_Analysis/RealRun/Analysis_07v_Validation/03_DefineTumorSubtype/Outputs_v1/20220411_TumorSubtype.RData")
load("~/OneDrive/01_Project/01_Codex/Bioinfo_Analysis/RealRun/Analysis_07v_Validation/04_DefineMacrophageSubtype/Outputs_v1/20220411_MacrophageSubtype.RData")
load("~/OneDrive/01_Project/01_Codex/Bioinfo_Analysis/RealRun/Analysis_07v_Validation/05_DefineCD4TSubtype/Outputs_v1/20220411_CD4TSubtype.RData")
load("~/OneDrive/01_Project/01_Codex/Bioinfo_Analysis/RealRun/Analysis_07v_Validation/06_DefineCD8TSubtype/Outputs_v1/20220411_CD8TSubtype.RData")

identical(MacrophageSubtype$cellid,CD8TSubtype$cellid)
table(MacrophageSubtype$MacrophagesSubtype)
table(CD8TSubtype$CD8TSubtype)

identical(MacrophageSubtype$cellid,CD4TSubtype$cellid)
identical(CD8TSubtype$cellid,CD4TSubtype$cellid)
identical(TumorSubtype$cellid,CD4TSubtype$cellid)
identical(TumorSubtype$cellid,CD8TSubtype$cellid)
identical(TumorSubtype$cellid,MacrophageSubtype$cellid)
identical(TumorSubtype$Class,MacrophageSubtype$Class)
identical(TumorSubtype$`actin Nucleus Intensity`,MacrophageSubtype$`actin Nucleus Intensity`)
identical(TumorSubtype$YMax,MacrophageSubtype$YMax)
identical(rownames(TumorSubtype),rownames(MacrophageSubtype))

Allsubtypes = cbind(TumorSubtype,CD4TSubtype[,42:43],CD8TSubtype[,42:43],MacrophageSubtype[,42:43])

load("~/Library/CloudStorage/OneDrive-个人/01_Project/01_Codex/Bioinfo_Analysis/RealRun/Analysis_07v_Validation/02_DefineCelltype/Outputs/codex_F_com_F.RData")
codex_F_com_F@meta.data$cellid = as.numeric(rownames(codex_F_com_F@meta.data))

EndothelialCelltype <- subset(codex_F_com_F,idents="Endothelial cells")
EndothelialCelltype = subset(EndothelialCelltype@meta.data,select = c("celltype","cellid"))

BCelltype <- subset(codex_F_com_F,idents="B cells")
BCelltype = subset(BCelltype@meta.data,select = c("celltype","cellid"))

FibroblastCelltype <- subset(codex_F_com_F,idents="Fibroblasts")
FibroblastCelltype = subset(FibroblastCelltype@meta.data,select = c("celltype","cellid"))

LymEndoCelltype <- subset(codex_F_com_F,idents="Lymphatic endothelial cells")
LymEndoCelltype = subset(LymEndoCelltype@meta.data,select = c("celltype","cellid"))

BiTractCelltype <- subset(codex_F_com_F,idents="Biliary tract cells")
BiTractCelltype = subset(BiTractCelltype@meta.data,select = c("celltype","cellid"))
  
FibroEndoCelltype <- subset(codex_F_com_F,idents="Fibroblasts/Endothelial cells")
FibroEndoCelltype = subset(FibroEndoCelltype@meta.data,select = c("celltype","cellid"))

TumorSubtype_subset = subset(TumorSubtype,select = c("TumorSubtype","cellid"))
TumorSubtype_subset = subset(TumorSubtype_subset,TumorSubtype != "Stromal")
colnames(TumorSubtype_subset)[1] = "celltype"

CD4TSubtype_subset = subset(CD4TSubtype,select = c("CD4TSubtype","cellid"))
CD4TSubtype_subset = subset(CD4TSubtype_subset,CD4TSubtype != "cell_others")
colnames(CD4TSubtype_subset)[1] = "celltype"

CD8TSubtype_subset = subset(CD8TSubtype,select = c("CD8TSubtype","cellid"))
CD8TSubtype_subset = subset(CD8TSubtype_subset,CD8TSubtype != "cell_others")
colnames(CD8TSubtype_subset)[1] = "celltype"

MacrophageSubtype_subset = subset(MacrophageSubtype,select = c("MacrophagesSubtype","cellid"))
MacrophageSubtype_subset = subset(MacrophageSubtype_subset,MacrophagesSubtype != "cell_others")
colnames(MacrophageSubtype_subset)[1] = "celltype"

temp = rbind(TumorSubtype_subset,EndothelialCelltype,CD8TSubtype_subset,BCelltype,FibroEndoCelltype,MacrophageSubtype_subset,CD4TSubtype_subset,FibroblastCelltype,LymEndoCelltype,BiTractCelltype)
temp$cellid = as.numeric(temp$cellid)
temp = temp[order(temp$cellid),]

identical(Allsubtypes$cellid,as.character(temp$cellid))

Allsubtypes$Allsubtypes = temp$celltype

TumorSubtype_subset$celltype = "Tumor"

temp = rbind(TumorSubtype_subset,EndothelialCelltype,CD8TSubtype_subset,BCelltype,FibroEndoCelltype,MacrophageSubtype_subset,CD4TSubtype_subset,FibroblastCelltype,LymEndoCelltype,BiTractCelltype)
temp$cellid = as.numeric(temp$cellid)
temp = temp[order(temp$cellid),]

identical(Allsubtypes$cellid,as.character(temp$cellid))

Allsubtypes$StromalSubtype = temp$celltype

write.table(Allsubtypes, file = "testing_allReg_seurat_._analysis_20220411_AllSubtypes.csv", sep=",", col.names = NA, row.names = TRUE)

