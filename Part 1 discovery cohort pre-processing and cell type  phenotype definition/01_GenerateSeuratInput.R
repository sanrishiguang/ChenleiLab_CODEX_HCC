# R Version: 4.2.0
# rm(list=ls())

library(Seurat) # v4.1.1
library(tidyverse) # v1.3.2
library(patchwork) # v1.1.2
library(ggplot2) # v3.3.6
library(dplyr) # v1.0.10
library(circlize) # v0.4.16
library(RColorBrewer) # v1.1-3
library(pheatmap) # v1.0.12
library(data.table) # v1.14.2
library(reshape) # v0.8.9
library(tidyr) # v1.2.1
library(stringr) # v1.4.1
library(Rphenograph) # v0.99.1
library(ggpubr) # v0.4.0
library(flowCore) # v2.8.0
library(Rcpp) # v1.0.9
library(cytofkit) # v0.99.0
library(igraph) # v1.3.5
library(ggthemes) # v4.2.4
library(Rtsne) # v0.16
library(cytofexplorer) # v2.05
library(survival) # v3.4-0
library(survminer) # v0.4.9

#### Extract TMA and reg info of the training data ####
load("/Users/taozhou/Desktop/Project/Codex/Analysis_07v/01_GenerateSeuratInput/HCC401.Rdata")
TrainingData.df = df_training[order(df_training$TMA),]
TMA_reg_info = list()
for (i in 1:length(unique(TrainingData.df$TMA))) {
  TMA_reg_info[[i]] = paste0("reg0",sprintf("%02d",TrainingData.df[which(TrainingData.df$TMA==i),]$Reg),"_qupath_cell_seg_data.txt")
}
rm(df_fix,df_training,df_validation)

#### Generate Seurat input ####
base_path <- '/Volumes/Tao的扩展盘/备份/Codex/Bioinfo_Analysis/RealRun/00_SegData_QupathOK'
dirnames <- stringr::str_sort(list.files(base_path),numeric = TRUE)
Panel <- c("actin","aSMA","c_Myc","Caspase_3","CD107A","CD11C","CD163","CD20","CD21","CD3",
           "CD31","CD4","CD44","CD45","CD45RO","CD68","CD8","DAPI","E_Cadherin","FOXP3",
           "Glypican3","Hepar1","HIF1a","HistoneH3","HLA_DR","Keratin","KI67","p_AMPK","p_mTOR","p_S6",
           "p53","Pan_CK","PD_L1","PD1","Podoplanin","Twist1","Vimentin")
Panel_36 <- c("actin","aSMA","c_Myc","Caspase_3","CD107A","CD11C","CD163","CD20","CD21","CD3",
           "CD31","CD4","CD44","CD45","CD45RO","CD68","CD8","E_Cadherin","FOXP3",
           "Glypican3","Hepar1","HIF1a","HistoneH3","HLA_DR","Keratin","KI67","p_AMPK","p_mTOR","p_S6",
           "p53","Pan_CK","PD_L1","PD1","Podoplanin","Twist1","Vimentin")

## New way to avoid for loop
empty_list <- vector(mode = "list", length = length(dirnames))
i=1
while(i <= length(dirnames)){
  base_path_TMA = file.path(base_path,dirnames[i])
  empty_list[i]=base_path_TMA
  i=i+1
}

list_of_files = list()
i=1
while(i <= length(empty_list)){
  file_path = paste(empty_list[i],unlist(TMA_reg_info[i]),sep = "/")
  list_of_files[i]=list(file_path)
  i=i+1
}

inputdata=lapply(list_of_files,function(list_of_files){
  list_of_files = read_tsv(list_of_files,id="File_Path")
}) ############################################################################## Id parameter needs MacOS system to run !!!

inputdata_csd = lapply(inputdata,function(inputdata){
  inputdata = inputdata %>% select(-contains("Phenotype")) %>% 
    select(contains("Cell"),-contains("DAPI"),-contains('Blank'),-contains('Empty')) %>%
    select(contains("Mean")) %>%
    rename_with(~str_remove_all(.x, ': Cell:')) %>%
    rename_with(~str_remove_all(.x, '\\(.*\\) ')) %>%
    rename_with(~str_remove_all(.x, ' Mean')) %>% 
    na.omit()
})
csd = bind_rows(inputdata_csd)

inputdata_raw = lapply(inputdata,function(inputdata){
  inputdata = inputdata %>% select(-contains("Phenotype")) %>% 
    select(-contains(Panel))
})

inputdata_raw[[7]]$`Nucleus: Circularity` = as.double(inputdata_raw[[7]]$`Nucleus: Circularity`)
inputdata_raw[[7]]$`Cell: Circularity` = as.double(inputdata_raw[[7]]$`Cell: Circularity`)

csd_raw = bind_rows(inputdata_raw)

csd_raw$TMA = str_split(csd_raw$File_Path, "/", simplify = TRUE)[,10]
csd_raw$reg = str_split(csd_raw$File_Path, "/", simplify = TRUE)[,11]
csd_raw$reg = str_split(csd_raw$reg, "_", simplify = TRUE)[,1]
csd_raw$Class = paste(csd_raw$TMA,csd_raw$reg,sep = "_")

save(csd,file="csd.RData")
save(csd_raw,file="csd_raw.RData")
save(TrainingData.df, file = "TrainingData.RData")

#### Filter Data ####
identical(rownames(csd),rownames(csd_raw))
DataFilter = csd
DataFilter$TMA = csd_raw$TMA
DataFilter$Cell.ID = rownames(DataFilter)
DataFilter$rowsum = rowSums(DataFilter[,1:36])

FilterID = DataFilter[84145:3378225,] # filter TMA1_5

pdf("all_cells.pdf",width = 4,height = 4)
plot(sort(FilterID$rowsum))
dev.off()

pdf("50000cells.pdf",width = 4,height = 4)
plot(sort(FilterID$rowsum)[1:50000])
abline(h = 700, col="blue")

pdf("rev50000cells.pdf",width = 4,height = 4)
plot(sort(FilterID$rowsum)[3244082:3294081])
abline(h = 200000, col="red")
dev.off()

FilterID = subset(FilterID, rowsum > 700 & rowsum < 200000) 

limit_Exprs = c(5000,40000,20000,2000,45000,15000,
                25000,15000,20000,10000,10000,20000,
                40000,20000,35000,40000,40000,35000,
                40000,15000,45000,15000,max(FilterID$HistoneH3),35000,25000,max(FilterID$KI67),
                2000,10000,20000,10000,30000,20000,
                20000,max(FilterID$Podoplanin),4000,35000)
for (i in 1:36) {
  pdf(paste0(colnames(FilterID)[i],".pdf"),width = 4,height = 4)
  plot(sort(FilterID[[i]])[3233469:3283468],ylab=paste0(colnames(FilterID)[i],"_[3233469:3283468]"))
  abline(h = limit_Exprs[i], col="red")
  dev.off()
}

FilterID = subset(FilterID, actin < 5000 & aSMA < 40000 & c_Myc < 20000 & Caspase_3 < 2000 & CD107A < 45000 & CD11C < 15000 &
                    CD163 < 25000 & CD20 < 15000 & CD21 < 20000 & CD3 < 10000 & CD31 < 10000 & CD4 < 20000 &
                    CD44 < 40000 & CD45 < 20000 & CD45RO < 35000 & CD68 < 40000 & CD8 < 40000 & E_Cadherin < 35000 &
                    FOXP3 < 40000 & Glypican3 < 15000 & Hepar1 < 45000 & HIF1a < 15000 & HLA_DR < 35000 & Keratin < 25000 &
                    p_AMPK < 2000 & p_mTOR < 10000 & p_S6 < 20000 & p53 < 10000 & Pan_CK < 30000 & PD_L1 < 20000 &
                    PD1 < 20000 & Twist1 < 4000 & Vimentin < 35000)

FilterID = FilterID$Cell.ID

rm(DataFilter)
save(FilterID,file = "FilterID.RData")


