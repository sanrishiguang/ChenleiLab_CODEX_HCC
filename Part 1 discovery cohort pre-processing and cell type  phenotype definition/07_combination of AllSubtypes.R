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
load("~/OneDrive/01_Project/01_Codex/Bioinfo_Analysis/RealRun/Analysis_07v/03_DefineTumorSubtype/Outputs_v2/20220315_TumorSubtype.RData")
load("~/OneDrive/01_Project/01_Codex/Bioinfo_Analysis/RealRun/Analysis_07v/04_DefineMacrophageSubtype/Outputs_v3/20220315_MacrophageSubtype.RData")
load("~/OneDrive/01_Project/01_Codex/Bioinfo_Analysis/RealRun/Analysis_07v/05_DefineCD4TSubtype/Outputs_v4/20220315_CD4TSubtype.RData")
load("~/OneDrive/01_Project/01_Codex/Bioinfo_Analysis/RealRun/Analysis_07v/06_DefineCD8TSubtype/Outputs_v3/20220315_CD8TSubtype.RData")

identical(MacrophageSubtype$cellid,CD8TSubtype$cellid)
table(MacrophageSubtype$MacrophagesSubtype)
table(CD8TSubtype$CD8TSubtype)

CD8TSubtype = CD8TSubtype[CD8TSubtype$cellid %in% CD4TSubtype$cellid,]
MacrophageSubtype = MacrophageSubtype[MacrophageSubtype$cellid %in% CD4TSubtype$cellid,]
identical(MacrophageSubtype$cellid,CD4TSubtype$cellid)
identical(CD8TSubtype$cellid,CD4TSubtype$cellid)

CD4TSubtype = CD4TSubtype[CD4TSubtype$cellid %in% TumorSubtype$cellid,]
CD8TSubtype = CD8TSubtype[CD8TSubtype$cellid %in% TumorSubtype$cellid,]
MacrophageSubtype = MacrophageSubtype[MacrophageSubtype$cellid %in% TumorSubtype$cellid,]
TumorSubtype = TumorSubtype[TumorSubtype$cellid %in% CD4TSubtype$cellid,]

identical(TumorSubtype$cellid,CD4TSubtype$cellid)
identical(TumorSubtype$cellid,CD8TSubtype$cellid)
identical(TumorSubtype$cellid,MacrophageSubtype$cellid)
identical(TumorSubtype$Class,MacrophageSubtype$Class)
identical(TumorSubtype$`actin Nucleus Intensity`,MacrophageSubtype$`actin Nucleus Intensity`)
identical(TumorSubtype$YMax,MacrophageSubtype$YMax)
identical(rownames(TumorSubtype),rownames(MacrophageSubtype))

Allsubtypes = cbind(TumorSubtype,CD4TSubtype[,42:43],CD8TSubtype[,42:43],MacrophageSubtype[,42:43])

load("~/Library/CloudStorage/OneDrive-个人/01_Project/01_Codex/Bioinfo_Analysis/RealRun/Analysis_07v/02_DefineCelltype/Outputs/codex_F_com_F.RData")
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

temp = rbind(TumorSubtype_subset,EndothelialCelltype,CD8TSubtype_subset,BCelltype,MacrophageSubtype_subset,CD4TSubtype_subset,FibroblastCelltype,LymEndoCelltype,BiTractCelltype)
temp$cellid = as.numeric(temp$cellid)
temp = temp[order(temp$cellid),]

identical(Allsubtypes$cellid,as.character(temp$cellid))

Allsubtypes$Allsubtypes = temp$celltype

TumorSubtype_subset$celltype = "Tumor"

temp = rbind(TumorSubtype_subset,EndothelialCelltype,CD8TSubtype_subset,BCelltype,MacrophageSubtype_subset,CD4TSubtype_subset,FibroblastCelltype,LymEndoCelltype,BiTractCelltype)
temp$cellid = as.numeric(temp$cellid)
temp = temp[order(temp$cellid),]

identical(Allsubtypes$cellid,as.character(temp$cellid))

Allsubtypes$StromalSubtype = temp$celltype

write.table(Allsubtypes, file = "testing_allReg_seurat_._analysis_20220315_AllSubtypes_withartifact.csv", sep=",", col.names = NA, row.names = TRUE)

Allsubtypes_withoutartifact = Allsubtypes[-grep("artifact",Allsubtypes$Allsubtypes),]
write.table(Allsubtypes_withoutartifact, file = "testing_allReg_seurat_._analysis_20220315_AllSubtypes_withoutartifact.csv", sep=",", col.names = NA, row.names = TRUE)

#### Survival: HIF1a+ Macrophages ####
PatientID = unique(MacrophageSubtype$Class)
NoMaPercent = c()
PosPatient = c()
PosPercent = c()
LowMaPercent = c()
HighMaPercent = c()

for (i in 1:length(PatientID)) {
  temp = subset(MacrophageSubtype,Class == PatientID[i])
  temp1 = as.data.frame(round((table(temp$MacrophagesSubtype)/nrow(temp))*100,4))
  if(length(grep("HIF1a",temp1$Var1))==0){
    NoMaPercent = c(NoMaPercent,PatientID[i])
  }
  else{
    temp2 = round((nrow(temp[temp$MacrophagesSubtype %in% c("Macrophages_HIF1a+","Macrophages_HIF1a+_c-Myc+","Macrophages_HIF1a+_Vimentin+"),])/nrow(temp))*100,2)
    PosPatient = c(PosPatient,PatientID[i])
    PosPercent = c(PosPercent,temp2)
  }
}

HighMaPercent = subset(data.frame(PosPatient,PosPercent),PosPercent>0.01)$PosPatient
LowMaPercent = subset(data.frame(PosPatient,PosPercent),PosPercent<=0.01)$PosPatient

HighMaPercent = data.frame(HighMaPercent,Group=rep("High",length(HighMaPercent)))
colnames(HighMaPercent)[1] = "Class" 
LowMaPercent = data.frame(LowMaPercent,Group=rep("No",length(LowMaPercent)))
colnames(LowMaPercent)[1] = "Class" 
NoMaPercent = data.frame(NoMaPercent,Group=rep("No",length(NoMaPercent)))
colnames(NoMaPercent)[1] = "Class" 

MaGroup = rbind(HighMaPercent,LowMaPercent,NoMaPercent)

temp = str_split(MaGroup$Class,"_",simplify = T)
temp = paste0(temp[,1],"_",temp[,3])
MaGroup$Class = temp
MaGroup = MaGroup[order(MaGroup$Class),]

load("~/OneDrive/01_Project/01_Codex/Bioinfo_Analysis/RealRun/Analysis_07v/01_GenerateSeuratInput/Outputs/TrainingData.RData")
Surv_data = TrainingData.df
Surv_data$TMA_Reg = paste0("TMA",Surv_data$TMA,"_","reg0", sprintf("%02d",Surv_data$Reg))
Surv_data = Surv_data[Surv_data$TMA_Reg %in% temp,]
Surv_data = Surv_data[order(Surv_data$TMA_Reg),]

identical(MaGroup$Class,Surv_data$TMA_Reg)

Surv_data$MaGroup = MaGroup$Group

res <- pairwise_survdiff(Surv(RFSday, RFS01) ~ MaGroup,
                         data = Surv_data, p.adjust.method = "none")
View(res$p.value)

#### Survival: FOXP3+ CD4T ####
PatientID = unique(CD4TSubtype$Class)
NoMaPercent = c()
PosPatient = c()
PosPercent = c()
LowMaPercent = c()
HighMaPercent = c()

for (i in 1:length(PatientID)) {
  temp = subset(CD4TSubtype,Class == PatientID[i])
  temp1 = as.data.frame(round((table(temp$CD4TSubtype)/nrow(temp))*100,4))
  if(length(grep("KI67",temp1$Var1))==0 & length(grep("FOXP3",temp1$Var1))==0){
    NoMaPercent = c(NoMaPercent,PatientID[i])
  }
  else{
    temp2 = round((nrow(temp[temp$CD4TSubtype %in% c("CD4T_KI67+","CD4T_FOXP3+","CD4T_FOXP3+_KI67+"),])/nrow(temp))*100,2)
    PosPatient = c(PosPatient,PatientID[i])
    PosPercent = c(PosPercent,temp2)
  }
}

HighMaPercent = subset(data.frame(PosPatient,PosPercent),PosPercent>1)$PosPatient
LowMaPercent = subset(data.frame(PosPatient,PosPercent),PosPercent<=1)$PosPatient

HighMaPercent = data.frame(HighMaPercent,Group=rep("High",length(HighMaPercent)))
colnames(HighMaPercent)[1] = "Class" 
LowMaPercent = data.frame(LowMaPercent,Group=rep("Low",length(LowMaPercent)))
colnames(LowMaPercent)[1] = "Class" 
NoMaPercent = data.frame(NoMaPercent,Group=rep("No",length(NoMaPercent)))
colnames(NoMaPercent)[1] = "Class" 

MaGroup = rbind(HighMaPercent,LowMaPercent,NoMaPercent)

temp = str_split(MaGroup$Class,"_",simplify = T)
temp = paste0(temp[,1],"_",temp[,3])
MaGroup$Class = temp
MaGroup = MaGroup[order(MaGroup$Class),]

load("~/OneDrive/01_Project/01_Codex/Bioinfo_Analysis/RealRun/Analysis_07v/01_GenerateSeuratInput/Outputs/TrainingData.RData")
Surv_data = TrainingData.df
Surv_data$TMA_Reg = paste0("TMA",Surv_data$TMA,"_","reg0", sprintf("%02d",Surv_data$Reg))
Surv_data = Surv_data[Surv_data$TMA_Reg %in% temp,]
Surv_data = Surv_data[order(Surv_data$TMA_Reg),]

identical(MaGroup$Class,Surv_data$TMA_Reg)

Surv_data$MaGroup = MaGroup$Group

res <- pairwise_survdiff(Surv(RFSday, RFS01) ~ MaGroup,
                         data = Surv_data, p.adjust.method = "none")
View(res$p.value)

## Add celltype
temp = codex_F_com_F@meta.data[as.numeric(codex_F_com_F@meta.data$cellid) %in% as.numeric(Allsubtypes_withoutartifact$cellid),]
identical(temp$cellid,Allsubtypes_withoutartifact$cellid) # F "e+05"
Allsubtypes_withoutartifact$celltype = temp$celltype

write.table(Allsubtypes_withoutartifact, file = "testing_allReg_seurat_._analysis_20220324_AllSubtypes_withoutartifact_addct.csv", sep=",", col.names = NA, row.names = TRUE)

