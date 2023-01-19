
#R version 4.1.0

#######  load packages#######
rm(list=ls())library(tidyverse)  #‘1.3.2’
library(patchwork)  # ‘1.1.1’
library(plotly)   #‘4.10.0’
# library(SingleR)
library(ggplot2)  #‘3.3.6’
library(clusterProfiler)  #‘4.5.2’
library(org.Hs.eg.db)  #‘3.14.0’
library(enrichplot)  #‘1.14.2’
library(ggpubr)  # ‘0.4.0’
library(Seurat)  #‘4.1.1’
library(fgsea)   # ‘1.20.0’
library(msigdbr)  #‘7.5.1’
library(survminer)  # ‘0.4.9’
library(survival)   #‘3.3.1’
library(data.table)   #‘1.14.2’
gc()



library(RColorBrewer)
cols=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)

mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), 
                 legend.key = element_rect(fill = "white", colour = "white"), 
                 legend.background = (element_rect(colour= "white", fill = "white")))

`%notin%` <- Negate(`%in%`)



####Part 1: Create the object ####
input = 'HCC_log_tpm_expression_matrix.txt'
if(file_test("-f",input)){
  project.data <- fread(input,sep="\t", header=TRUE,data.table=F)
  rownames(project.data) <- project.data[,1]
  project.data[,1] <- NULL
}else{
  project.data <- Read10X(data.dir=input)
}

Mymetadata <- read.table('HCC_cell_metadata.txt', sep = "\t",header = T,row.names = 1)
Mymetadata <- Mymetadata[-1,]
project <- CreateSeuratObject(counts = project.data, project = "project", meta.data=Mymetadata, min.cells = 0, min.features = 0)
sort(unique(project$cell_type))



Tcell_project <- subset(project, cell_type %in% c('C0_Tcell', "C1_Tcell", "C3_Tcell", "C5_Tcell", "C19_Tcell"))
Mye_project <- subset(project, cell_type %in% c("C2_Mye.", "C8_Mye.", "C11_Mye.", "C15_Mey.", "C20_Mye"))


#######~~~ save project, file = 'originalprojectect.Rdata'############
save(project, file = 'originalprojectect.Rdata')
saveRDS(Tcell_project,file='Tcell_project.Rds')
saveRDS(Mye_project,file='Mye_project.Rds')


####Part 2: Clustering Myeloid Celltypes####

######1. Basic Seurat Pipelines ######

Mye_project <- readRDS('Mye_project.Rds')

Mye_project_processed <- NormalizeData(Mye_project, normalization.method = "LogNormalize", scale.factor = 10000) # nolint # nolint
Mye_project_processed <- FindVariableFeatures(Mye_project_processed, selection.method = "vst", nfeatures = 2000)
Mye_project_processed <- ScaleData(Mye_project_processed)
Mye_project_processed <- RunPCA(Mye_project_processed, features = VariableFeatures(object = Mye_project_processed))

#
# Mye_project_processed <- JackStraw(Mye_project_processed, num.replicate = 100)
# Mye_project_processed <- ScoreJackStraw(Mye_project_processed, dims = 1:30)

#
pcs_number=30

pdf(file="DimHeatmap_PC.pdf",width=15,height=45)
DimHeatmap(object = Mye_project_processed, dims = 1:pcs_number, cells = 500, balanced = TRUE) # nolint
dev.off()

# pdf(file="JackStrawPlot_PC.pdf",width=15)
# p <- JackStrawPlot(object = Mye_project_processed, dims = 1:10)
# print(p)
# dev.off()

pdf(file="ElbowPlot_PC.pdf",width=15)
p <- ElbowPlot(Mye_project_processed, ndims = pcs_number, reduction = "pca")
print(p)
dev.off()

######~~~~ saveRDS:Mye_project_processed_before_cluster.Rds#######
saveRDS(Mye_project_processed,file = "Mye_project_processed_before_cluster.Rds")


##PCA Dims 1：20 for downstreaming analysis
#
Mye_project_processed <- readRDS("Mye_project_processed_before_cluster.Rds")

Mye_project_processed <- FindNeighbors(Mye_project_processed, dims = 1:20)
Mye_project_processed <- RunUMAP(Mye_project_processed, dims = 1:20)
Mye_project_processed <- FindClusters(Mye_project_processed, resolution = 0.5)



##Resolution 0.3 for doswnstreaming analysis
DimPlot(Mye_project_processed, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(Mye_project_processed, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "RNA_snn_res.0.15",ncol = 6) + NoLegend()


FeaturePlot(Mye_project_processed, features = c( "CD68","CD163", "CD14","FCGR3A", # Mono/Macro
                                             "CD1C", # cDC
                                             "APOA2","SPP1",'BAG3', 'FCN1','THBS1', 'IL32','PRF1','CD1E','MAT2A', 'IDO1', 'CLEC9A', 'LAMP3','TOP2A','VIM'
                                             ))
DotPlot(Mye_project_processed, features = c( "CD68","CD163", "CD14","FCGR3A", # Mono/Macro
                                                 "CD1C", # cDC
                                                 "APOA2","SPP1",'BAG3', 'FCN1','THBS1', 'IL32','PRF1','CD1E','MAT2A', 'IDO1', 'CLEC9A', 'LAMP3','TOP2A','VIM'
)) + rotate_x_text()

FeaturePlot(Mye_project_processed, features = c('CD68','CD163'), slot = 'data',min.cutoff = 1.5 )
DotPlot(Mye_project_processed, features = c('CD68','CD163'))

VlnPlot(Mye_project_processed,features = c('VIM', 'SPP1')  )
FeaturePlot(Mye_project_processed, features = c('VIM'), slot = 'data',min.cutoff = 2 )
DotPlot(Mye_project_processed, features = 'VIM', scale.min = 0)


####### 2: SingleR -----------------------------------------------------------------

# library(SingleR)
# library(celldex)
# # hpca.se=HumanPrimaryCellAtlasData() 
# # Blue.se=BlueprintEncodeData() 
# # Immune.se=DatabaseImmuneCellExpressionData()
# # Nover.se=NovershternHematopoieticData()
# # MonacoIm.se=MonacoImmuneData()
# # ImmGen.se=ImmGenData() #(鼠)
# # Mouse.se=MouseRNAseqData() #(鼠)
# ###
# 
# load('../ref_Hematopoietic.RData')
# load('../ref_Monaco_114s.RData')
# load('../ref_Human_all.RData')
# 
# 
# project_for_SingleR <- GetAssayData(Mye_project_processed, slot="data")
# 
# ###ref_Monaco
# pbmc.hesc <- SingleR(test = project_for_SingleR, ref = ref_Human_all, 
#                      labels = ref_Human_all$label.fine) 
# ##
# table(pbmc.hesc$labels,Mye_project_processed$seurat_clusters) 
# # pdf("plotScoreHeatmap.pdf")
# print(plotScoreHeatmap(pbmc.hesc))
# # dev.off()
# Mye_project_processed@meta.data$labels <-pbmc.hesc$labels
# # pdf("Umap.pdf",height=5,width=10)
# print(DimPlot(Mye_project_processed, group.by = c("seurat_clusters", "labels"),reduction = "umap", label = T))
# # dev.off()
# library(reshape2)
# aa=table(pbmc.hesc$labels,Mye_project_processed$seurat_clusters)
# aa= apply(aa,2,function(x) x/sum(x))
# df=as.data.frame(melt(aa))
# df$Var2=as.factor(df$Var2)
# g <- ggplot(df, aes(Var2, Var1)) + geom_point(aes(size = value), colour = "darkred") + theme_bw() 
# # pdf("singleR_match_seurat.pdf",height=5,width=10)
# print(g)
# # dev.off()
# library(pheatmap)
# # pdf(paste(i,"heatmap.pdf",sep ="_"),height=5,width=10)
# pheatmap(log2(aa+10), color=colorRampPalette(c("white", "blue"))(101))
# # dev.off()
# 
# 
# 
# ##
# FeaturePlot(project, features = c('CD163', 'CD3D',
#                                              'TRBC1','TRBC2',
#                                              'LYZ','CD68','FCGR3A','CD163',
#                                              'ALB','CPS1',
#                                              'CD31','CLDN5',
#                                              'COL1A1',
#                                              'CD79A',
#                                              'EPCAM'
#                                              ))
# 
# 


###### 3. Celltype Definition ####

## Definition
new.cluster.ids = c(
  "MC0_Monocyte", # 0 
  "MC1_Macrophage", # 1   #
  "MC2_DC2", # 2
  "MC3_Macrophage", # 3
  "MC4_Macrophage", # 4
  "MC5_Macrophage", # 5
  "MC6_Macrophage", # 6
  "MC7_DC2", # 7
  'MC8_cycling',
  'MC9_DC2',
  'MC10_DC1',
  'MC11_DC3'
)

names(new.cluster.ids) <- levels(Mye_project_processed)
Mye_project_processed <- RenameIdents(Mye_project_processed, new.cluster.ids)
Mye_project_processed@meta.data$Myeloid_subtype = Idents(Mye_project_processed)

DimPlot(Mye_project_processed, reduction = "umap", label = T)
DimPlot(Mye_project_processed, reduction = "umap", label = TRUE, repel = T,split.by = "Myeloid_subtype") + NoLegend()

###~~~ saveRDS ：Mye_project_processed_after_definition.Rds######
saveRDS(Mye_project_processed,file = "Mye_project_processed_after_definition.Rds")


####Figures
FeaturePlot(Mye_project_processed, features = 'CD68', slot = 'data',min.cutoff = 1.5 )  ##macrophage
FeaturePlot(Mye_project_processed, features = 'CD163', slot = 'data',min.cutoff = 1.0 ) ##macrophage
FeaturePlot(Mye_project_processed, features = 'FOLR2', slot = 'data',min.cutoff = 1.0 ) ##macrophage
FeaturePlot(Mye_project_processed, features = 'SELENOP', slot = 'data',min.cutoff = 1.0 ) ##macrophage
# FeaturePlot(Mye_project_processed, features = 'SPP1', slot = 'data',min.cutoff = 1.0 ) ##macrophage
# FeaturePlot(Mye_project_processed, features = 'BAG3', slot = 'data',min.cutoff = 1.0 ) ##macrophage


FeaturePlot(Mye_project_processed, features ='FCN1' , slot = 'data',min.cutoff = 0 )  ##monocyte
FeaturePlot(Mye_project_processed, features ='GZMB' , slot = 'data',min.cutoff = 0.4 )  ##monocyte
# FeaturePlot(Mye_project_processed, features ='PRF1' , slot = 'data',min.cutoff = 0)  ##monocyte

FeaturePlot(Mye_project_processed, features = c( "CD1C",'CD1E') , slot = 'data',min.cutoff = 0 )  #DC2
FeaturePlot(Mye_project_processed, features = 'IDO1' , slot = 'data',min.cutoff = 1.0)  ##DC1
FeaturePlot(Mye_project_processed, features = c('CLEC9A'), slot = 'data',min.cutoff = 0 )  ##DC1
FeaturePlot(Mye_project_processed, features = c( 'LAMP3'), slot = 'data',min.cutoff = 0 )  ##DC3

FeaturePlot(Mye_project_processed, features = c( 'TOP2A'), slot = 'data',min.cutoff = 0 )  ##CYCLING



######4. Macrophage_VIM definition according to VIM threshold###############


Mye_project_processed <- readRDS('Mye_project_processed_after_definition.Rds')
DimPlot(Mye_project_processed, reduction = 'umap', label = T)
Macro_project <- subset(Mye_project_processed, Myeloid_subtype %in% c('MC1_Macrophage', 'MC3_Macrophage',
                                                                      'MC4_Macrophage', 'MC5_Macrophage',
                                                                      'MC6_Macrophage'))

FeaturePlot(Macro_project, features = 'VIM', min.cutoff=0) 

dat_VIM <- FetchData(object = Macro_project, vars = 'VIM')
hist(dat_VIM$VIM)
quantile(dat_VIM$VIM,probs = seq(0, 1, 0.1))

FeaturePlot(Macro_project, features = 'VIM', min.cutoff=2)  
FeaturePlot(Macro_project, features = 'VIM', min.cutoff=2)  
FeaturePlot(Macro_project, features = 'VIM', min.cutoff=2.2) 
# FeaturePlot(Macro_project, features = 'HLA-A', min.cutoff=1) 
# FeaturePlot(Macro_project, features = 'CD4', min.cutoff=1) 
# FeaturePlot(Macro_project, features = 'MKI67', min.cutoff=1) 
# FeaturePlot(Macro_project, features = 'ITGAX', min.cutoff=1) 
# FeaturePlot(Macro_project, features = 'HIF1A', min.cutoff=1) 


Macro_project_processed <- NormalizeData(Macro_project, normalization.method = "LogNormalize", scale.factor = 10000)
Macro_project_processed <- FindVariableFeatures(Macro_project_processed, selection.method = "vst", nfeatures = 2000)
Macro_project_processed <- ScaleData(Macro_project_processed)
Macro_project_processed <- RunPCA(Macro_project_processed, features = VariableFeatures(object = Macro_project_processed))

#
# Macro_project_processed <- JackStraw(Macro_project_processed, num.replicate = 100)
# Macro_project_processed <- ScoreJackStraw(Macro_project_processed, dims = 1:30)

#
pcs_number=30
pdf(file="DimHeatmap_PC_macrophage.pdf",width=15,height=45)
DimHeatmap(object = Macro_project_processed, dims = 1:pcs_number, cells = 500, balanced = TRUE)
dev.off()

# pdf(file="JackStrawPlot_PC.pdf",width=15)
# p <- JackStrawPlot(object = Macro_project_processed, dims = 1:10)
# print(p)
# dev.off()

ElbowPlot(Macro_project_processed, ndims = pcs_number, reduction = "pca")
saveRDS(Macro_project_processed,file = "Macro_project_processed_before_cluster.Rds")


##PCA Dims 1：20 for downstreaming analysis
#
Macro_project_processed <- readRDS("Macro_project_processed_before_cluster.Rds")

Macro_project_processed <- FindNeighbors(Macro_project_processed, dims = 1:12)
Macro_project_processed <- RunUMAP(Macro_project_processed, dims = 1:12)

Macro_project_processed <- FindClusters(Macro_project_processed, resolution = 1)

DimPlot(Macro_project_processed, reduction = 'umap', label = T)
FeaturePlot(Macro_project_processed, features = c('VIM','SPP1'), min.cutoff=2) 
DotPlot(Macro_project_processed, features = 'VIM', scale.min = 0)
DotPlot(Macro_project_processed, features = c('VIM','SPP1'), scale.min = 0)
VlnPlot(Macro_project_processed, features = c('VIM','SPP1'))

diffs <- FindAllMarkers(object = Macro_project_processed, logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.1,assay="RNA")
diffs$cluster <-as.factor(as.numeric(as.vector(diffs$cluster)))

# write.table(diffs, sep="\t", file="project.FindAllMarkers.xls",quote=FALSE, row.names=FALSE)
idents<-as.factor(as.numeric(as.vector(Macro_project_processed@active.ident)))
names(idents) <- rownames(Macro_project_processed@meta.data)
Macro_project_processed@active.ident <- idents

sum(duplicated(diffs$gene))
diffs2<-subset(diffs,subset=!duplicated(gene))

library(dplyr)
diffs2  %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(Macro_project_processed, features = top10$gene,slot = "data") + NoLegend()
DoHeatmap(Macro_project_processed, features = top10$gene,slot = "scale.data") + NoLegend()
VlnPlot(Macro_project_processed, features = 'VIM')
table(Idents(Macro_project_processed)) / length(Idents(Macro_project_processed))

######~~~ write.csv ： diffs macro clusters#######
write.csv(diffs2, file='./Macro_project_processed_diffs2.csv')


dat_VIM <- FetchData(object = Macro_project_processed, vars = 'VIM')
VIM_cut_value <- quantile(dat_VIM$VIM,probs = seq(0, 1, 0.1)) [[10]]

Macro_project_processed@meta.data$VIM_group <- 'Macrophage_VIM_low'
Macro_project_processed@meta.data$VIM_group[dat_VIM > VIM_cut_value] <- 'Macrophage_VIM_hi'

with(Macro_project_processed@meta.data, table(VIM_group, RNA_snn_res.1))
ggplot(data=dat_VIM, mapping = aes(x=VIM, y=..density..)) + 
  geom_histogram(color='black', fill='grey')+
  geom_density(color='blue')

Macro_project_processed@meta.data$VIM_group_cluster <- paste0(Macro_project_processed@meta.data$VIM_group, '_C',Macro_project_processed@meta.data$RNA_snn_res.1)
table(Macro_project_processed@meta.data$VIM_group_cluster)
Idents(Macro_project_processed) <- Macro_project_processed$VIM_group
DimPlot(Macro_project_processed, reduction = 'umap')
with(Macro_project_processed@meta.data, table(VIM_group, tissue_source))


######~~~ saveRDS  Macro_project_processed_after_definition_VIM.Rds#######
saveRDS(Macro_project_processed,file = "Macro_project_processed_after_definition_VIM.Rds")


# Part 3 : T-cell procession ----------------------------------------------

######1. regular processing ######
T_cell_processed <-  readRDS('Tcell_project.Rds')
T_cell_processed <- NormalizeData(T_cell_processed, normalization.method = "LogNormalize", scale.factor = 10000)
T_cell_processed <- FindVariableFeatures(T_cell_processed, selection.method = "vst", nfeatures = 2000)
T_cell_processed <- ScaleData(T_cell_processed)
T_cell_processed <- RunPCA(T_cell_processed, features = VariableFeatures(object = T_cell_processed))

pcs_number=30
pdf(file="T_cell_DimHeatmap_PC.pdf",width=15,height=45)
DimHeatmap(object = T_cell_processed, dims = 1:pcs_number, cells = 500, balanced = TRUE)
dev.off()

ElbowPlot(T_cell_processed, ndims = pcs_number, reduction = "pca")


T_cell_processed_1 <- FindNeighbors(T_cell_processed, dims = 1:15)
T_cell_processed_1 <- RunUMAP(T_cell_processed_1, dims = 1:15)
T_cell_processed_1 <- FindClusters(T_cell_processed_1, resolution = 0.5)


##Resolution 0.5 for doswnstreaming analysis
DimPlot(T_cell_processed_1, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(T_cell_processed_1, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "RNA_snn_res.0.5",ncol = 6) + NoLegend()

######~~~ save "T_cell_processing.Rdata"#####
save(T_cell_processed,T_cell_processed_1,file = "T_cell_processing.Rdata")


# FeaturePlot(T_cell_processed_1, features = c('FOXP3','CD4','CD8A','LAG3','MKI67','FGFBP2','IL2RA'))
# DotPlot(T_cell_processed_1, features = c('FOXP3','CD4','CD8A','LAG3','MKI67','FGFBP2','IL2RA'))

FeaturePlot(T_cell_processed_1, features = c('FOXP3','CD4','CD8A','IL2RA'))
DotPlot(T_cell_processed_1, features = c('FOXP3','CD4','IL2RA'))


# # ######2. SingleR #####
# # 
# library(SingleR)
# library(celldex)
# 
# load('../ref_Hematopoietic.RData')
# load('../ref_Monaco_114s.RData')
# load('../ref_Human_all.RData')
# 
# project <- T_cell_processed_1
# project_for_SingleR <- GetAssayData(project, slot="data")
# 
# 
# ###ref_Human_all
# pbmc.hesc <- SingleR(test = project_for_SingleR, ref = ref_Human_all,
#                        labels = ref_Human_all$label.fine)
# ##
# table(pbmc.hesc$labels,project$seurat_clusters)
# # pdf("plotScoreHeatmap.pdf")
# print(plotScoreHeatmap(pbmc.hesc))
# # dev.off()
# project@meta.data$labels <-pbmc.hesc$labels
# # pdf("Umap.pdf",height=5,width=10)
# print(DimPlot(project, group.by = c("seurat_clusters", "labels"),reduction = "umap", label = T))
# # dev.off()
# library(reshape2)
# aa=table(pbmc.hesc$labels,project$seurat_clusters)
# aa= apply(aa,2,function(x) x/sum(x))
# df=as.data.frame(melt(aa))
# df$Var2=as.factor(df$Var2)
# g <- ggplot(df, aes(Var2, Var1)) + geom_point(aes(size = value), colour = "darkred") + theme_bw()
# # pdf("singleR_match_seurat.pdf",height=5,width=10)
# print(g)
# # dev.off()
# library(pheatmap)
# # pdf(paste(i,"heatmap.pdf",sep ="_"),height=5,width=10)
# pheatmap(log2(aa+10), color=colorRampPalette(c("white", "blue"))(101))
# # dev.off()


######3. CD4/CD8 definition#####
load('T_cell_processing.Rdata')
project <- T_cell_processed_1

## Definition
new.cluster.ids = c(
  "TC0_", # 0
  "TC1_Treg", # 1
  "TC2", # 2
  "TC3", # 3
  "TC4", # 4
  "TC5", # 5
  "TC6",
  'TC7',
  'TC8',
  'TC9',
  'TC10',
  'TC11'
)

names(new.cluster.ids) <- levels(project)
project <- RenameIdents(project, new.cluster.ids)
project@meta.data$T_subtype = Idents(project)

DimPlot(project, reduction = "umap", label = T, repel = T)
DimPlot(project, reduction = "umap", label = TRUE, repel = T,split.by = "T_subtype") + NoLegend()

######~~~ saveRDS "T_cell_subtype_definition.Rds"######
saveRDS(project,file = "T_cell_subtype_definition.Rds")
##C1_Treg



# Part 4 : Re-Annotation the total seurat object -----------------------------------------
library(dplyr)
library(purrr)
load('originalSeuratObject.Rdata')
project <- project0
dim(project)  ##16498
T_cell_project <- readRDS('T_cell_subtype_definition.Rds')   ###7294 cells
Macro_project <- readRDS('Macro_project_processed_after_definition_VIM.Rds')  ###1635 cells

table(T_cell_project$T_subtype)
table(Macro_project$VIM_group_cluster)
table(Macro_project$VIM_group)

View(Macro_project@meta.data)
temp0 <- project@meta.data
tempT <- T_cell_project@meta.data %>% dplyr::select(T_subtype)
tempM <- Macro_project@meta.data %>% dplyr::select(Myeloid_subtype, Macrophage_subtype = RNA_snn_res.1, VIM_group,VIM_group_cluster)

temp0T <- merge(temp0, tempT, by='row.names', all=T)
row.names(temp0T) <- temp0T$Row.names
temp0T <- temp0T[,-1]
temp0TM <- merge(temp0T, tempM, by='row.names', all=T)
row.names(temp0TM) <- temp0TM$Row.names
temp0TM <- temp0TM[,-1]

#新建一列包含Treg，Macrophage-VIM-low/hi
temp0TM$MyIdents <- temp0TM$cell_type
temp0TM$MyIdents [temp0TM$T_subtype == 'TC1_Treg'] <- 'Treg'
temp0TM$MyIdents [temp0TM$VIM_group == "Macrophage_VIM_low"] <- "Macrophage_VIM_low"
temp0TM$MyIdents [temp0TM$VIM_group == "Macrophage_VIM_hi"] <- "Macrophage_VIM_hi"

project@meta.data <-  temp0TM

#########~~~ saveRdata: identifiedSeuratObject#######
save(project, file='identifiedSeuratObject.Rdata')



# Part 5:  Macrophage_vimentin_low vs Macrophage_vimentin_hi----------------------------------------------------------------


library(nichenetr)
library(Seurat) # please update to Seurat V4
library(ggplot2)
library(dplyr)
Macrophage_obj <- readRDS('Macro_project_processed_after_definition_VIM.Rds')

Idents(Macrophage_obj) <- Macrophage_obj@meta.data$VIM_group
table(Idents(Macrophage_obj))

DEGs <- FindMarkers(Macrophage_obj, ident.1 = "Macrophage_VIM_hi", ident.2 = "Macrophage_VIM_low", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
deg <- DEGs
###vocalo
View(deg)
library("EnhancedVolcano")
colnames(deg)

EnhancedVolcano(deg,lab=rownames(deg),
                x="avg_log2FC",y="p_val_adj",
                pCutoff = 10e-2,FCcutoff = 0.585,
                xlim = c(-1,1.5),
                labSize = 2 )



library(clusterProfiler)
library(org.Hs.eg.db)

deg$GeneID <- row.names(deg)
gene_up=rownames(deg[deg$avg_log2FC > 0 & deg$p_val_adj < 0.05, ])
gene_down=rownames(deg[deg$avg_log2FC < 0 & deg$p_val_adj < 0.05, ])
gene_up=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                   keys = gene_up,
                                                   columns = 'ENTREZID',
                                                   keytype = 'SYMBOL')[,2]))
gene_down=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                     keys = gene_down,
                                                     columns = 'ENTREZID',
                                                     keytype = 'SYMBOL')[,2]))
gene_up <- unique(gene_up)
gene_down <- unique(gene_down)

go.up <- enrichGO(gene = gene_up,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP" ,
                  pAdjustMethod = "BH",
                  readabl = TRUE)
go.up@result <- go.up@result[order(go.up@result$Count, decreasing = T),]
dotplot(go.up, title = paste0('Up-regulation genes : ', length(gene_up)),showCategory = 10)

go.down <- enrichGO(gene = gene_down,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP" ,
                    pAdjustMethod = "BH",
                    readabl = TRUE)
go.down@result <- go.down@result[order(go.down@result$Count, decreasing = T),]
dotplot(go.down,title = paste0('Down-regulation genes : ', length(gene_down)),showCategory = 10 )

# KEGG.up <- enrichKEGG(gene = gene_up)
# KEGG.down <- enrichKEGG(gene = gene_down)
# dotplot(KEGG.up, title='Up-regulated genes in Macrophage_VIM: KEGG', showCategory = 20)
# dotplot(KEGG.down, title='Down-regulated genes in Macrophage_VIM: KEGG', showCategory = 20)

#######~~~ saveRdata : Macrophage_vimentin_low vs Macrophage_vimentin_hi######
save(list=ls(), file='Macrophage_vimentin_low vs Macrophage_vimentin_hi.Rdata')

######1. GO.UP#######
# #positive regulation of cytokine production in GO.UP
# Mygenelist.up <- go.up@result$geneID[go.up@result$Description == 'positive regulation of cytokine production']
# Mygenelist.up <- strsplit(Mygenelist.up,'/')[[1]]
# Mygenelist.up
# Mygenelist.up_ <- as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
#                                                    keys = Mygenelist.up,
#                                                    columns = 'ENTREZID',
#                                                    keytype = 'SYMBOL')[,2]))
# 
# 
# go.Mygenelist.up <- enrichGO(gene = Mygenelist.up_,
#                   OrgDb = org.Hs.eg.db,
#                   ont = "BP" ,
#                   pAdjustMethod = "BH",
#                   readabl = TRUE)
# View(go.Mygenelist.up@result)
# 
# #negative regulation of immune system process
# Mygenelist.up <- go.up@result$geneID[go.up@result$Description == 'negative regulation of immune system process']
# Mygenelist.up <- strsplit(Mygenelist.up,'/')[[1]]
# Mygenelist.up
# 
# ######2. GO.DOWN######
# 
# #immune response-regulating signaling pathway
# #positive regulation of cytokine production
# #macroautophagy
# Mygenelist <- go.down@result$geneID[go.down@result$Description == 'macroautophagy']
# Mygenelist <- strsplit(Mygenelist,'/')[[1]]
# Mygenelist
# Mygenelist_ <- as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
#                                                              keys = Mygenelist,
#                                                              columns = 'ENTREZID',
#                                                              keytype = 'SYMBOL')[,2]))
# go.Mygenelist <- enrichGO(gene = Mygenelist_,
#                              OrgDb = org.Hs.eg.db,
#                              ont = "BP" ,
#                              pAdjustMethod = "BH",
#                              readabl = TRUE)
# View(go.Mygenelist@result)


######3. gene features to plot######
VlnPlot(Macrophage_obj, features = c( 'IL1B', 'IL10','TGFB1'))


#######4. Macrophage_VIM vs others: geneset from Zhou Tao#######

deg <- DEGs
deg$GeneID <- row.names(deg)
# deg <- deg %>% filter(avg_log2FC > 0 & p_val_adj < 0.05)


####看geneset
GOBP_PHAGOCYTOSIS = clusterProfiler::read.gmt("../gmt/GOBP_PHAGOCYTOSIS.v7.5.1.gmt")
GOBP_PHAGOCYTOSIS = deg$GeneID[deg$GeneID %in% GOBP_PHAGOCYTOSIS$gene]
GOBP_PHAGOCYTOSIS = list(GOBP_PHAGOCYTOSIS) 

M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN = clusterProfiler::read.gmt("../gmt/GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN.v7.5.1.gmt")
M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN = deg$GeneID[deg$GeneID %in% M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN$gene]
M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN = list(M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN) 

M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP = clusterProfiler::read.gmt("../gmt/GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP.v7.5.1.gmt")
M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP = deg$GeneID[deg$GeneID %in% M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP$gene]
M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP = list(M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP) 

HALLMARK_ANGIOGENESIS = clusterProfiler::read.gmt("../gmt/HALLMARK_ANGIOGENESIS.v7.5.1.gmt")
HALLMARK_ANGIOGENESIS = deg$GeneID[deg$GeneID %in% HALLMARK_ANGIOGENESIS$gene]
HALLMARK_ANGIOGENESIS = list(HALLMARK_ANGIOGENESIS) 

WP_ANGIOGENESIS = clusterProfiler::read.gmt("../gmt/WP_ANGIOGENESIS.v7.5.1.gmt")
WP_ANGIOGENESIS = deg$GeneID[deg$GeneID %in% WP_ANGIOGENESIS$gene]
WP_ANGIOGENESIS = list(WP_ANGIOGENESIS) 

###Addmodule
Macrophage_obj <- AddModuleScore(Macrophage_obj, features = GOBP_PHAGOCYTOSIS, name = 'GOBP_PHAGOCYTOSIS')
Macrophage_obj <- AddModuleScore(Macrophage_obj, features = M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN, name = 'M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN')
Macrophage_obj <- AddModuleScore(Macrophage_obj, features = M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP, name = 'M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP')
Macrophage_obj <- AddModuleScore(Macrophage_obj, features = HALLMARK_ANGIOGENESIS, name = 'HALLMARK_ANGIOGENESIS')
Macrophage_obj <- AddModuleScore(Macrophage_obj, features = WP_ANGIOGENESIS, name = 'WP_ANGIOGENESIS')

####
Macrophage_obj_VIM = subset(Macrophage_obj@meta.data, VIM_group == 'Macrophage_VIM_hi')
Macrophage_obj_VIM = subset(Macrophage_obj_VIM, select = c('GOBP_PHAGOCYTOSIS1',
                                                         'M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN1',
                                                         "M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP1",
                                                         'HALLMARK_ANGIOGENESIS1',
                                                         'WP_ANGIOGENESIS1'))

Macrophage_obj_others = subset(Macrophage_obj@meta.data, VIM_group == 'Macrophage_VIM_low')
Macrophage_obj_others = subset(Macrophage_obj_others, select = c('GOBP_PHAGOCYTOSIS1',
                                                           'M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN1',
                                                           "M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP1",
                                                           'HALLMARK_ANGIOGENESIS1',
                                                           'WP_ANGIOGENESIS1'))

GS_plot = rbind(colMeans(Macrophage_obj_VIM),
                colMeans(Macrophage_obj_others))

rownames(GS_plot) = c("Macrophage_VIM_hi","Macrophage_VIM_low")
GS_plot = t(as.matrix(GS_plot))

pheatmap::pheatmap(GS_plot, cluster_rows = F, cluster_cols = F, angle_col = 45, scale = "row")
pheatmap::pheatmap(GS_plot, cluster_rows = F, cluster_cols = F, angle_col = 45)


p <- pheatmap::pheatmap(GS_plot, cluster_rows = F, cluster_cols = F, angle_col = 45)

pdf(file='E:\\11. CODEX\\2022年7月 整图\\Figure 7\\M1M2.pdf', height=3, width=5)
print(p)
dev.off()



##Macropahge_VIM 

####MC1, MC2 etc
table(Macrophage_obj$VIM_group_cluster)
Macrophage_obj_MC1 = subset(Macrophage_obj@meta.data, VIM_group_cluster == 'Macrophage_VIM_hi_C1')
Macrophage_obj_MC1 = subset(Macrophage_obj_VIM, select = c('GOBP_PHAGOCYTOSIS1',
                                                           'M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN1',
                                                           "M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP1",
                                                           'HALLMARK_ANGIOGENESIS1',
                                                           'WP_ANGIOGENESIS1'))

Macrophage_obj_MC3 = subset(Macrophage_obj@meta.data, VIM_group_cluster == 'Macrophage_VIM_low_C1')
Macrophage_obj_MC3 = subset(Macrophage_obj_MC3, select = c('GOBP_PHAGOCYTOSIS1',
                                                           'M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN1',
                                                           "M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP1",
                                                           'HALLMARK_ANGIOGENESIS1',
                                                           'WP_ANGIOGENESIS1'))

Macrophage_obj_MC4 = subset(Macrophage_obj@meta.data, VIM_group_cluster == 'Macrophage_VIM_hi_C6')
Macrophage_obj_MC4 = subset(Macrophage_obj_MC4, select = c('GOBP_PHAGOCYTOSIS1',
                                                           'M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN1',
                                                           "M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP1",
                                                           'HALLMARK_ANGIOGENESIS1',
                                                           'WP_ANGIOGENESIS1'))

Macrophage_obj_MC5 = subset(Macrophage_obj@meta.data, VIM_group_cluster == 'Macrophage_VIM_low_C6')
Macrophage_obj_MC5 = subset(Macrophage_obj_MC5, select = c('GOBP_PHAGOCYTOSIS1',
                                                           'M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN1',
                                                           "M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP1",
                                                           'HALLMARK_ANGIOGENESIS1',
                                                           'WP_ANGIOGENESIS1'))



GS_plot = rbind(colMeans(Macrophage_obj_MC1),
                colMeans(Macrophage_obj_MC3),
                colMeans(Macrophage_obj_MC4),
                colMeans(Macrophage_obj_MC5))

rownames(GS_plot) = c("Macrophage_VIM_hi_C1","Macrophage_VIM_low_C1","Macrophage_VIM_hi_C6","Macrophage_VIM_low_C6")
GS_plot = t(as.matrix(GS_plot))

pheatmap::pheatmap(GS_plot, cluster_rows = F, cluster_cols = F, angle_col = 45, scale = "row")
pheatmap::pheatmap(GS_plot, cluster_rows = F, cluster_cols = F, angle_col = 45, scale = "none")






temp <- Macrophage_obj@meta.data

unique(temp$VIM_group_cluster)
table((temp$VIM_group_cluster))
temp <- Macrophage_obj@meta.data %>% filter(VIM_group_cluster %in% c('Macrophage_VIM_hi_C0',
                                                                     'Macrophage_VIM_low_C0',
                                                                     'Macrophage_VIM_hi_C1',
                                                                     'Macrophage_VIM_low_C1',
                                                                     'Macrophage_VIM_hi_C4',
                                                                     'Macrophage_VIM_low_C4',
                                                                     'Macrophage_VIM_hi_C5',
                                                                     'Macrophage_VIM_low_C5',
                                                                     'Macrophage_VIM_hi_C6',
                                                                     'Macrophage_VIM_low_C6',
                                                                     'Macrophage_VIM_hi_C8',
                                                                     'Macrophage_VIM_low_C8'))

temp$VIM_group_cluster <- factor(temp$VIM_group_cluster, levels=c('Macrophage_VIM_hi_C0',
                                                                     'Macrophage_VIM_low_C0',
                                                                     'Macrophage_VIM_hi_C1',
                                                                     'Macrophage_VIM_low_C1',
                                                                     'Macrophage_VIM_hi_C4',
                                                                     'Macrophage_VIM_low_C4',
                                                                     'Macrophage_VIM_hi_C5',
                                                                     'Macrophage_VIM_low_C5',
                                                                     'Macrophage_VIM_hi_C6',
                                                                     'Macrophage_VIM_low_C6',
                                                                     'Macrophage_VIM_hi_C8',
                                                                     'Macrophage_VIM_low_C8'
                                                                     ))

ggplot(data=temp, mapping = aes(x=VIM_group_cluster, y=WP_ANGIOGENESIS1, color=VIM_group_cluster)) + 
  geom_violin(aes(fill=VIM_group_cluster)) +
  geom_boxplot()+
  coord_flip()+
  guides(fill=F, color=F)
ggplot(data=temp, mapping = aes(x=VIM_group_cluster, y=HALLMARK_ANGIOGENESIS1, color=VIM_group_cluster)) + 
  geom_violin(aes(fill=VIM_group_cluster)) +
  geom_boxplot()+
  coord_flip()+
  guides(fill=F, color=F)




library(ggsignif)
library(ggpubr)
library(ggsci)
mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), 
                 legend.key = element_rect(fill = "white", colour = "white"), 
                 legend.background = (element_rect(colour= "white", fill = "white")))

ggplot(temp,aes(x=VIM_group, y=WP_ANGIOGENESIS1))+
  geom_boxplot(aes(color = VIM_group),outlier.colour=NA)+
  geom_point(aes(color = VIM_group),position = 'jitter',size=0.5,alpha=0.5)+
  stat_compare_means(method="wilcox.test",label.x=1,label.y.npc='top',size=3,vjust=1)+
  # ggsignif::geom_signif(comparisons = , 
  #             test = "wilcox.test",
  #             map_signif_level=T,
  #             step_increase = 0.1,
  #             size = 0.5,
  #             textsize = 2)+
  labs(title='WP_ANGIOGENESIS1')+
  mytheme+
  scale_color_manual(values=c(pal_npg('nrc',alpha = 1.0)(5)))+
  theme(plot.title = element_text(size=8))+ 
  theme(axis.text.x = element_blank())


ggplot(temp,aes(x=VIM_group, y=HALLMARK_ANGIOGENESIS1))+
  geom_boxplot(aes(color = VIM_group),outlier.colour=NA)+
  geom_point(aes(color = VIM_group),position = 'jitter',size=0.5,alpha=0.5)+
  stat_compare_means(method="wilcox.test",label.x=1,label.y.npc='top',size=3,vjust=1)+
  # ggsignif::geom_signif(comparisons = , 
  #             test = "wilcox.test",
  #             map_signif_level=T,
  #             step_increase = 0.1,
  #             size = 0.5,
  #             textsize = 2)+
  labs(title='HALLMARK_ANGIOGENESIS1')+
  mytheme+
  scale_color_manual(values=c(pal_npg('nrc',alpha = 1.0)(5)))+
  theme(plot.title = element_text(size=8))+ 
  theme(axis.text.x = element_blank())
     

ggplot(temp,aes(x=VIM_group, y=M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN1))+
  geom_boxplot(aes(color = VIM_group),outlier.colour=NA)+
  geom_point(aes(color = VIM_group),position = 'jitter',size=0.5,alpha=0.5)+
  stat_compare_means(method="wilcox.test",label.x=1,label.y.npc='top',size=3,vjust=1)+
  # ggsignif::geom_signif(comparisons = , 
  #             test = "wilcox.test",
  #             map_signif_level=T,
  #             step_increase = 0.1,
  #             size = 0.5,
  #             textsize = 2)+
  labs(title='M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN1')+
  mytheme+
  scale_color_manual(values=c(pal_npg('nrc',alpha = 1.0)(5)))+
  theme(plot.title = element_text(size=8))+ 
  theme(axis.text.x = element_blank())


ggplot(temp,aes(x=VIM_group, y=M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP1))+
  geom_boxplot(aes(color = VIM_group),outlier.colour=NA)+
  geom_point(aes(color = VIM_group),position = 'jitter',size=0.5,alpha=0.5)+
  stat_compare_means(method="wilcox.test",label.x=1,label.y.npc='top',size=3,vjust=1)+
  # ggsignif::geom_signif(comparisons = , 
  #             test = "wilcox.test",
  #             map_signif_level=T,
  #             step_increase = 0.1,
  #             size = 0.5,
  #             textsize = 2)+
  labs(title='M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP1')+
  mytheme+
  scale_color_manual(values=c(pal_npg('nrc',alpha = 1.0)(5)))+
  theme(plot.title = element_text(size=8))+ 
  theme(axis.text.x = element_blank())


ggplot(temp,aes(x=VIM_group, y=GOBP_PHAGOCYTOSIS1))+
  geom_boxplot(aes(color = VIM_group),outlier.colour=NA)+
  geom_point(aes(color = VIM_group),position = 'jitter',size=0.5,alpha=0.5)+
  stat_compare_means(method="wilcox.test",label.x=1,label.y.npc='top',size=3,vjust=1)+
  # ggsignif::geom_signif(comparisons = , 
  #             test = "wilcox.test",
  #             map_signif_level=T,
  #             step_increase = 0.1,
  #             size = 0.5,
  #             textsize = 2)+
  labs(title='GOBP_PHAGOCYTOSIS1')+
  mytheme+
  scale_color_manual(values=c(pal_npg('nrc',alpha = 1.0)(5)))+
  theme(plot.title = element_text(size=8))+ 
  theme(axis.text.x = element_blank())






###### 5. Macropahge_VIM_hi vs low DEG: hallmark enrichment######

load('Macrophage_vimentin_low vs Macrophage_vimentin_hi.Rdata')

library(msigdbr)
library(enrichplot)
msigdbr_show_species()
head(Hum_msigdbr)
library(DOSE)
library(forcats)
##
#C5： GO
#C7: immunologic signature gene sets
#C8: celltype signature gene sets
#H： hallmark gene sets


gene_up=rownames(deg[deg$avg_log2FC > 0 & deg$p_val_adj < 0.05, ])
gene_down=rownames(deg[deg$avg_log2FC < 0 & deg$p_val_adj < 0.05, ])

Hum_hallmark <- msigdbr(species = 'Homo sapiens',category="H") %>%
  dplyr::select(gs_name, gene_symbol)


#### UP
hallmark_up <- enricher(gene_up,TERM2GENE=Hum_hallmark)
barplot(hallmark_up)
##
ego3 <- mutate(hallmark_up, richFactor = Count / as.numeric
               (sub("/\\d+", "", BgRatio)))
ggplot(ego3, showCategory = 10,
       aes(richFactor,
           fct_reorder(Description, richFactor))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        trans = "log10",
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_dose(12) +
  xlab("Rich Factor") +
  ylab(NULL) +
ggtitle("Biological Processes")


cnetplot(hallmark_up,showCategory=20)
heatplot(hallmark_up)


#####DOWN

hallmark_down <- enricher(gene_down,TERM2GENE=Hum_hallmark)
barplot(hallmark_down)
##美化
ego3 <- mutate(hallmark_down, richFactor = Count / as.numeric
               (sub("/\\d+", "", BgRatio)))
ggplot(ego3, showCategory = 10,
       aes(richFactor,
           fct_reorder(Description, richFactor))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        trans = "log10",
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_dose(12) +
  xlab("Rich Factor") +
  ylab(NULL)+
  ggtitle("Biological Processes")


cnetplot(hallmark_down,showCategory=20)
heatplot(hallmark_down)




# ###### 6. Macropahge_VIM_hi vs low DEG: C2######
# 
# # load('Macrophage_vimentin_low vs Macrophage_vimentin_hi.Rdata')
# 
# library(msigdbr)
# library(enrichplot)
# msigdbr_show_species()
# 
# Hum_msigdbr <- msigdbr(species = 'Homo sapiens') 
# Hum_msigdbr %>% filter(gs_cat == 'C2')  %>% with(., table(gs_subcat))
# 
# ##
# #C5： GO
# #C7: immunologic signature gene sets
# #C8: celltype signature gene sets
# #H： hallmark gene sets
# 
# 
# gene_up=rownames(deg[deg$avg_log2FC > 0 & deg$p_val_adj < 0.05, ])
# # gene_down=rownames(deg[deg$avg_log2FC < 0 & deg$p_val_adj < 0.05, ])
# 
# Hum_hallmark <- msigdbr(species = 'Homo sapiens',category="C2") %>% filter(gs_subcat == 'CP:REACTOME') %>% 
#   dplyr::select(gs_name, gene_symbol)
# 
# 
# #### UP
# hallmark_up <- enricher(gene_up,TERM2GENE=Hum_hallmark)
# summarise(hallmark_up)
# 
# 
# 
# 
# barplot(hallmark_up)
# ##美化
# ego3 <- mutate(hallmark_up, richFactor = Count / as.numeric
#                (sub("/\\d+", "", BgRatio)))
# ggplot(ego3, showCategory = 50,
#        aes(richFactor,
#            fct_reorder(Description, richFactor))) +
#   geom_segment(aes(xend=0, yend = Description)) +
#   geom_point(aes(color=p.adjust, size = Count)) +
#   scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
#                         trans = "log10",
#                         guide=guide_colorbar(reverse=TRUE, order=1)) +
#   scale_size_continuous(range=c(2, 10)) +
#   theme_dose(12) +
#   xlab("Rich Factor") +
#   ylab(NULL) +
#   ggtitle("CP:REACTOME")
# 
# 
# cnetplot(hallmark_up,showCategory=20)
# heatplot(hallmark_up)
# 
# 





# Part 6 : Nichenet: Treg -------------------------------------------------------

######1. sender:Mac-Vimentin, revceiver: Treg, condition: Tumor vs adjacent liverr######

library(nichenetr)
library(Seurat) # please update to Seurat V4
library(ggplot2)
library(tibble)
library(dplyr)
library(RColorBrewer)


load('E:/11.R files/Nichenet/NicheNet ligand-target model.Rdata')
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
head(lr_network)
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

load('identifiedSeuratObject.Rdata')

Idents(project) <- project@meta.data$MyIdents

#nichenet
nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = project, 
                                               receiver = "Treg", 
                                               condition_colname = "tissue_source", 
                                               condition_oi = "Tumor", 
                                               condition_reference = "Adjacent liver", 
                                               sender = c('Macrophage_VIM_hi'), 
                                               ligand_target_matrix = ligand_target_matrix, 
                                               lr_network = lr_network, 
                                               weighted_networks = weighted_networks, 
                                               organism = "human")

######~~~ save Rdata: Mac-vim vs Treg Nichenet ######
save(nichenet_output, file='NicheNet_Analysis_0601_part6-1.Rdata')

load('NicheNet_Analysis_0601_part6-1.Rdata')
load('Macrophage_vimentin_low vs Macrophage_vimentin_hi.Rdata')

##ligand_target 总图
nichenet_output$ligand_activity_target_heatmap

##ligand
nichenet_output$ligand_expression_dotplot
nichenet_output$ligand_differential_expression_heatmap
View(nichenet_output$ligand_activities)

temp=nichenet_output$ligand_activities[1:20,]
temp$test_ligand <- factor(temp$test_ligand, levels = temp$test_ligand[order(temp$pearson, decreasing = F)])
ggplot(data=temp, mapping = aes(y=test_ligand))+
  geom_bar(aes(fill = pearson) )+ 
  mytheme +
  scale_fill_gradient(low = 'white', high = '#e65100') +
  theme(axis.ticks.x  = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) +
  labs(y = 'Top Active Ligands')


top_ligands <- nichenet_output$top_ligands
# nichenet_output$ligand_target_matrix %>% .[1:10,1:6]
##top_ligands
DotPlot(Macrophage_obj, features = top_ligands)
VlnPlot(Macrophage_obj, features = top_ligands, ncol=10)
##target
nichenet_output$ligand_target_heatmap



###Protein-Protein Interaction
nichenet_output$ligand_receptor_heatmap
###ligand-receptor interaction
nichenet_output$ligand_receptor_heatmap_bonafide
nichenet_output$ligand_receptor_df
nichenet_output$top_receptors


###Differentail Receptor expression in two groups
Idents(project) <-  project@meta.data$MyIdents
DotPlot(project %>% subset(idents = "Treg"), features = nichenet_output$top_receptors %>% rev() %>% .[1:20], split.by = "tissue_source") + RotatedAxis()
VlnPlot(project %>% subset(idents = "Treg"), features = c('IL1R1','IL1R2','IL1RAP')) + RotatedAxis()

FeaturePlot(Macrophage_obj, features = c('VIM', 'SPP1'),min.cutoff = 2 )
FeaturePlot(Macrophage_obj, features = c('MMP12', 'S100A10'),min.cutoff = 2 )
FeaturePlot(Macrophage_obj, features = c('CD74'),min.cutoff = 2)
FeaturePlot(Macrophage_obj, features = c('MIF'),min.cutoff = 2 )


######2. circos plot： ligand-receptor#######

library(magrittr)
library(tidyverse)
library(circlize)
library(tidyr)
library(readr)
library(purrr)
library(tibble)
library(stringr)
#
# avg_expression_ligands = AverageExpression(project %>% subset(subset = aggregate == "LCMV"),features = nichenet_output$top_ligands) # if want to look specifically in LCMV-only cells
avg_expression_ligands = AverageExpression(Macrophage_obj, features = nichenet_output$top_ligands)
#
sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean()) #+ ligand_expression %>% sd())
}) %>% t()

sender_ligand_assignment[1:4,1:4]
sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)

(sender_ligand_assignment)


all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = nichenet_output$top_ligands %>% setdiff(unique_ligands)

VIMhi_specific_ligands = sender_ligand_assignment$Macrophage_VIM_hi %>% names() #%>% setdiff(general_ligands)
VIMlow_specific_ligands = sender_ligand_assignment$Macrophage_VIM_low %>% names() #%>% setdiff(general_ligands)
# Mono_specific_ligands = sender_ligand_assignment$Mono %>% names() %>% setdiff(general_ligands)
# DC_specific_ligands = sender_ligand_assignment$DC %>% names() %>% setdiff(general_ligands)

ligand_type_indication_df = tibble(
  ligand_type = c(rep("Macrophage-VIM-hi-specific", times = VIMhi_specific_ligands %>% length()),
                  rep("Macrophage-VIM-low-specific", times = VIMlow_specific_ligands %>% length())),
  ligand = c(VIMhi_specific_ligands, VIMlow_specific_ligands))

ligand_type_indication_df %>% head


###
active_ligand_target_links_df = nichenet_output$ligand_receptor_df_bonafide %>% 
  mutate(target_type = "Tumor-DE") %>% 
  inner_join(ligand_type_indication_df) 
# if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.40)
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)
ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$receptor %>% unique(), active_ligand_target_links_df_circos$receptor %>% unique())
circos_links = active_ligand_target_links_df %>% filter(!receptor %in% targets_to_remove &!ligand %in% ligands_to_remove)
circos_links
unique(circos_links$ligand_type)


#
grid_col_ligand =c("Macrophage-VIM-hi-specific" = "royalblue",
                   "Macrophage-VIM-low-specific" = "lawngreen")
grid_col_target =c(
  "Tumor-DE" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% dplyr::select(ligand,receptor, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(receptor,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$receptor)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 


#
target_order = circos_links$receptor %>% unique()
ligand_order = c(VIMhi_specific_ligands, VIMlow_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)


#
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5


##gaps
# gaps = c(
#   # width_ligand_target,
#   rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Macrophage-VIM-hi-specific") %>% distinct(ligand) %>% nrow() -1)),
#   width_different_cell,
#   rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Macrophage-VIM-low-specific") %>% distinct(ligand) %>% nrow() -1))
#   # width_ligand_target
# )
# length(gaps)





# circos.par(gap.degree = gaps)
chordDiagram(links_circle, 
             directional = 1,
             order=order,
             link.sort = TRUE, 
             link.decreasing = FALSE, 
             grid.col = grid_col,
             transparency = transparency, 
             diffHeight = 0.005, 
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", 
             link.visible = links_circle$weight >= cutoff_include_all_ligands,
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
circos.clear()


######3.circos plot： ligand-target#######

library(magrittr)
library(tidyverse)
library(circlize)
library(tidyr)
library(readr)
library(purrr)
library(tibble)
library(stringr)
# avg_expression_ligands = AverageExpression(project %>% subset(subset = aggregate == "LCMV"),features = nichenet_output$top_ligands) # if want to look specifically in LCMV-only cells
avg_expression_ligands = AverageExpression(Macrophage_obj, features = nichenet_output$top_ligands)
sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean()) #+ ligand_expression %>% sd())
}) %>% t()

sender_ligand_assignment[1:4,1:4]
sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)

(sender_ligand_assignment)


all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = nichenet_output$top_ligands %>% setdiff(unique_ligands)

VIMhi_specific_ligands = sender_ligand_assignment$Macrophage_VIM_hi %>% names() #%>% setdiff(general_ligands)
VIMlow_specific_ligands = sender_ligand_assignment$Macrophage_VIM_low %>% names() #%>% setdiff(general_ligands)
# Mono_specific_ligands = sender_ligand_assignment$Mono %>% names() %>% setdiff(general_ligands)
# DC_specific_ligands = sender_ligand_assignment$DC %>% names() %>% setdiff(general_ligands)

ligand_type_indication_df = tibble(
  ligand_type = c(rep("Macrophage-VIM-hi-specific", times = VIMhi_specific_ligands %>% length()),
                  rep("Macrophage-VIM-low-specific", times = VIMlow_specific_ligands %>% length())),
  ligand = c(VIMhi_specific_ligands, VIMlow_specific_ligands))

ligand_type_indication_df %>% head


active_ligand_target_links_df = nichenet_output$ligand_target_df %>% 
  mutate(target_type = "Tumor-DE") %>% 
  inner_join(ligand_type_indication_df) 
# if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.40)
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)
ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())
circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)
circos_links
unique(circos_links$ligand_type)


grid_col_ligand =c("Macrophage-VIM-hi-specific" = "lawngreen",
                   "Macrophage-VIM-low-specific" = "royalblue")
grid_col_target =c(
  "Tumor-DE" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% dplyr::select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 


target_order = circos_links$target %>% unique()
ligand_order = c(VIMhi_specific_ligands, VIMlow_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)


width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

# gaps = c(
#   # width_ligand_target,
#   rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Macrophage-VIM-hi-specific") %>% distinct(ligand) %>% nrow() -1)),
#   width_different_cell,
#   rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Macrophage-VIM-low-specific") %>% distinct(ligand) %>% nrow() -1))
#   # width_ligand_target
# )
# length(gaps)





# circos.par(gap.degree = gaps)
chordDiagram(links_circle, 
             directional = 1,
             order=order,
             link.sort = TRUE, 
             link.decreasing = FALSE, 
             grid.col = grid_col,
             transparency = transparency, 
             diffHeight = 0.005, 
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", 
             link.visible = links_circle$weight >= cutoff_include_all_ligands,
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
circos.clear()


#####4. Enrichment analysis of top target genes#####
temp <- nichenet_output$top_targets

gene_up=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                   keys = temp,
                                                   columns = 'ENTREZID',
                                                   keytype = 'SYMBOL')[,2]))

gene_up <- unique(gene_up)

go.up <- enrichGO(gene = gene_up,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP" ,
                  pAdjustMethod = "BH",
                  readabl = TRUE)
go.up@result <- go.up@result[order(go.up@result$Count, decreasing = T),]
dotplot(go.up, title = paste0('Top target genes : ', length(gene_up)),showCategory = 10)

go.up.filter <- filter(go.up, Description %in% c('positive regulation of cytokine production',
                                                 'response to hypoxia',
                                                 'regulation of angiogenesis',
                                                 'myeloid leukocyte migration',
                                                 'cellular response to interleukin-1',
                                                 'response to interferon-gamma',
                                                 'interleukin-10 production',
                                                 'response to transforming growth factor beta',
                                                 'regulation of transforming growth factor beta production',
                                                 'positive regulation of regulatory T cell differentiation'))
dotplot(go.up.filter, title = paste0('Top target genes : ', length(gene_up)),showCategory = 15)


ego3 <- mutate(go.up.filter, richFactor = Count / as.numeric
               (sub("/\\d+", "", BgRatio)))
ggplot(ego3, showCategory = 10,
       aes(richFactor,
           fct_reorder(Description, richFactor))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        trans = "log10",
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_dose(12) +
  xlab("Rich Factor") +
  ylab(NULL)+
  ggtitle("Biological Processes")


cnetplot(go.up.filter,showCategory=20)
heatplot(go.up.filter)




# Part 7  signature definition---------------------------------

library(msigdbr)
library(enrichplot)
library(DOSE)

##FCtop50

load('Macrophage_vimentin_low vs Macrophage_vimentin_hi.Rdata')

######1. Macrophage-vimentin hi signature#######

Mac_VIM_hi_genesets <- deg %>% .[order(deg$p_val_adj),] %>% filter(p_val_adj <0.01) %>%  filter(avg_log2FC > 0.5)
Mac_VIM_hi_genesets <- Mac_VIM_hi_genesets[!grepl("^MT-",Mac_VIM_hi_genesets$GeneID ), ]
Mac_VIM_hi_genesets <- Mac_VIM_hi_genesets[!grepl("^RP[L|S]",Mac_VIM_hi_genesets$GeneID ), ]
Mac_VIM_hi_genesets <- Mac_VIM_hi_genesets %>% filter(GeneID %notin% c('AP000350.10','ACTB', 'AC006386.1','AC012005.1',
                                                                       'AC090498.1'))
Mac_VIM_hi_genesets_R <- Mac_VIM_hi_genesets[order(Mac_VIM_hi_genesets$avg_log2FC, decreasing = T),]

saveRDS(Mac_VIM_hi_genesets_R, file='Mac_VIM_hi_genesets_R.RDS')


# Mac_VIM_hi_genesets_R <- readRDS(file='./MyGeneSignature/Mac_VIM_hi_genesets_R.RDS') 
Mac_VIM_hi_genesets_R_order <- Mac_VIM_hi_genesets_R[order(Mac_VIM_hi_genesets_R$avg_log2FC, decreasing = T),] 

Mac_VIM_hi_genesets_R_FCtop50 <- Mac_VIM_hi_genesets_R_order$GeneID [1:50]
saveRDS(Mac_VIM_hi_genesets_R_FCtop50, file= 'Mac_VIM_hi_genesets_R_FCtop50.RDS')

Mac_VIM_hi_genesets_R_FCtop50_df <- Mac_VIM_hi_genesets_R_order[1:50, ]
write.csv(Mac_VIM_hi_genesets_R_FCtop50_df, file = 'Table S7.csv')

######2. IL-1 production signature ##################

IL1_Production_genesets <- go.up@result %>% dplyr::filter(Description == 'interleukin-1 beta production') %>% dplyr::select(geneID) %>% 
  apply(., 2, function(x) {str_split(x, '/')}) %>% .[[1]]
saveRDS(IL1_Production_genesets[[1]], file = 'IL1_Production_genesets.RDS')





# Part 8 : HCC bulk-seq data to validate the findings --------------------
library(IOBR)
library(EPIC)
library(estimate) 
library(tidyverse)
library(tidyHeatmap)
library(maftools)
library(ggpubr)
library(ggplot2)
library(survival)


Mac_VIM_hi_FCtop50 <- readRDS(file='Mac_VIM_hi_genesets_R_FCtop50.RDS')
IL1_Production_genesets <-  readRDS(file = 'IL1_Production_genesets.RDS')


####signatures from handsome
sign <- read.delim("./immune_sign.txt")
names(sign)
sign_select <- sign[,c("Treg","CD8_Activation",'Antinflammatory',"Anergy","Cytolytic")] %>% as.list()

for (i in seq_along(c("Treg","CD8_Activation",'Antinflammatory',"Anergy","Cytolytic"))) {
  sign_select[[i]] <- sign_select[[i]][sign_select[[i]] != ""]
}



mysig <- list(Mac_VIM_hi_FCtop50 = Mac_VIM_hi_FCtop50,
              IL1_Production_genesets =IL1_Production_genesets)

names(signature_tme)[1:20]
VIM_signature_tme <- c(mysig,  signature_tme[c("T_cell_regulatory_Peng_et_al", 'Treg_Rooney_et_al')], sign_select)

# selected_gene <- c('WNT5A','IL1B','TGFB2','EFNA4','HAS2','DKK1')

save(list=ls(), file='signatures for Bulk-seq validation.Rdata')


###### 12.1 HCC 239 samples RuiRu ########
load(file='signatures for Bulk-seq validation.Rdata')
HCC_239_dat <- read.csv(file='E:\\11. CODEX\\LIHC bulk seq\\239 HCC bulk-seq from RuiRu TPM\\RNA.exp.csv', header=T)
dim(HCC_239_dat) #colname is sample
HCC_239_dat <- HCC_239_dat[ !duplicated(HCC_239_dat$锘縢ene), ]
row.names(HCC_239_dat) <-  HCC_239_dat$锘縢ene
HCC_239_dat <- HCC_239_dat[, -1]
temp <- HCC_239_dat %>% as.matrix()
temp <- log2(temp + 1)
HCC_239_dat_1 <- as.data.frame(temp)

##
HCC_239_dat_1_t <- HCC_239_dat_1 %>% t() %>% as.data.frame()
HCC_239_dat_1_t$ID <- row.names(HCC_239_dat_1_t) 

##
sig_tme_ruiru <- calculate_sig_score(pdata           = NULL,
                                     eset            = HCC_239_dat_1,
                                     signature       = VIM_signature_tme,
                                     method          = "zscore",
                                     mini_gene_count = 3,
                                     adjust_eset = TRUE)




###clinical information
library(readxl)
df_clinical_1st <- read_excel("E:\\11. CODEX\\LIHC bulk seq\\239 HCC bulk-seq from RuiRu TPM\\临床信息.xlsx")
df_clinical_2st <- read_excel("E:\\11. CODEX\\LIHC bulk seq\\239 HCC bulk-seq from RuiRu TPM\\第二次随访数据.xlsx")


temp1 <- df_clinical_1st %>% dplyr::select( "Sample ID" ,  "BCLC")
temp2 <- df_clinical_2st %>% dplyr::select("Sample ID" , 
                                           RFS01 = "第二次随访是否复发（是／否，1／0）",
                                           RFSday =  "手术至复发间隔时间（d）"
)
df_clinical <- merge(temp1, temp2, by="Sample ID")

sig_tme_ruiru_clinical <- merge(sig_tme_ruiru, df_clinical, by.x='ID', by.y= 'Sample ID')
str(sig_tme_ruiru_clinical)
sig_tme_ruiru_clinical$RFS01 <- as.numeric(sig_tme_ruiru_clinical$RFS01)
sig_tme_ruiru_clinical$RFSday <- as.numeric(sig_tme_ruiru_clinical$RFSday)

save(list=ls(), file='HCC-ruiru-239.Rdata')

#######~~~~Mac-Vim: poor prognosis######
library(survival)
library(survminer)
cycle_variable <- names(sig_tme_ruiru_clinical)[3:8]

for (i in 1:length(cycle_variable)) {
  sur.cut <- surv_cutpoint(sig_tme_ruiru_clinical, time= 'RFSday',event = 'RFS01' , variables = cycle_variable[i], minprop = 0.1)
  #summary(sur.cut)
  sur.cat <- surv_categorize(sur.cut)
  #head(sur.cat)
  #table(sur.cat$AllCN32)
  names(sur.cat) <- c("RFSday",'RFS01','group')
  fit <- survfit(Surv(RFSday, RFS01) ~ group, data = sur.cat)
  filename <- paste0('./HCC bulk data Survival/Ruiru//Optimal_RFS',"_",cycle_variable[i],".png")
  png(file=filename,width=1500,height=1500,res=300)
  p <- ggsurvplot(fit,palette = "npg",
                  risk.table = TRUE, pval = TRUE,
                  conf.int = FALSE, xlab="Time in Days",
                  ggtheme = theme_classic(),
                  title = cycle_variable[i])
  print(p)
  dev.off()
  print(cycle_variable[i])
}


hist(sig_tme_ruiru_clinical$Mac_VIM_hi_FCtop50)



library(psych)
library(corrplot)
names(sig_tme_ruiru_clinical)
temp <- sig_tme_ruiru_clinical[3:11]
##基于Frequency
MM<-corr.test(temp, method="spearman")
# pdf(file=paste0('1_CN_basedon_Freq.pdf'), width = 15, height = 15)
corrplot(MM$r,method = "square",order = 'hclust',
         p.mat = MM$p, sig.level = 0.001, insig = "label_sig")
# dev.off()




library(ggpubr)

subytp_A <-  names(sig_tme_ruiru_clinical)[4:11]
subytp_B <-  "Mac_VIM_hi_FCtop50"   

library(ggpubr)
#pearson correlation
for (i in seq_along(subytp_A)) {
  p <- ggscatter(sig_tme_ruiru_clinical, x = subytp_A[i], y = subytp_B,
                 shape = 21, fill ='#8491B4FF',color='darkgrey',size = 2, # Points color, shape and size
                 add = "reg.line",  # Add regressin line
                 add.params = list(color = "#7E6148FF", fill = "#B09C85FF"), # Customize reg. line
                 conf.int = TRUE, # Add confidence interval
                 cor.coef = T, # Add correlation coefficient. see ?stat_cor
                 cor.coeff.args = list(method = "pearson", label.x.npc = "middle", label.y.npc = "top",label.sep = "\n"),
                 title = 'Pearson correlation')
  ggsave(p,file = paste0("./HCC bulk correlation/Ruiru/Pearson_",cycle_variable[i],".png"),width = 3.5,height = 3)
}





######## 12.5 HCCDB6 --------------------------------------------------------------
load(file='signatures for Bulk-seq validation.Rdata')


HCC_HCCDB6_dat <- read.table(file='E:\\11. CODEX\\LIHC bulk seq\\HCCDB6\\GSE14520-GPL3921.gene.txt', sep='\t',header=T, row.names = 1)
dim(HCC_HCCDB6_dat) #colname is sample
HCC_HCCDB6_dat <- HCC_HCCDB6_dat[ !duplicated(HCC_HCCDB6_dat$Symbol), ]
row.names(HCC_HCCDB6_dat) <-  HCC_HCCDB6_dat$Symbol
HCC_HCCDB6_dat <- HCC_HCCDB6_dat[, -1]
HCC_HCCDB6_dat_1 <- HCC_HCCDB6_dat
sum(is.na(HCC_HCCDB6_dat_1 ))  ##有NA

#NA fixation
library(randomForest)
HCC_HCCDB6_dat_1 <- na.roughfix(HCC_HCCDB6_dat_1)

#rowname fixation
HCC_HCCDB6_dat_1 <- HCC_HCCDB6_dat_1[!grepl('///', rownames(HCC_HCCDB6_dat_1)),]


sig_tme_HCCDB6 <- calculate_sig_score(pdata           = NULL,
                                      eset            = HCC_HCCDB6_dat_1,  ##行是gene， 列是sample
                                      signature       = VIM_signature_tme,
                                      method          = "zscore",
                                      mini_gene_count = 3,
                                      adjust_eset = TRUE)
dim(sig_tme_HCCDB6)

###clinical information
df_HCCDB6_sample <- read.table(file='E:\\11. CODEX\\LIHC bulk seq\\HCCDB6\\HCCDB6.sample.txt', sep='\t', row.names = 1) %>%  t() %>% as.data.frame()
df_HCCDB6_patient <- read.table(file='E:\\11. CODEX\\LIHC bulk seq\\HCCDB6\\HCCDB6.patient.txt', sep='\t', row.names = 1) %>%  t() %>% as.data.frame()


df_HCCDB6_SP <- left_join(df_HCCDB6_sample, df_HCCDB6_patient, by='PATIENT')
df_HCCDB6_SP$ID <- gsub('-','.',df_HCCDB6_SP$SAMPLE_ID)


sig_tme_HCCDB6_SP <- left_join(sig_tme_HCCDB6,df_HCCDB6_SP, by='ID' )


ggboxplot(sig_tme_HCCDB6_SP, x='TYPE', y = 'Mac_VIM_hi_FCtop50',fill='TYPE',palette = "jco",add="jitter", add.params = list(fill="white"),trim = F, size = 0.5) +
  stat_compare_means(method = "wilcox.test")


######~~~~MacVim poor prognosis#######
sig_tme_HCCDB6_SP$SURVIVAL_TIME <- as.numeric(sig_tme_HCCDB6_SP$SURVIVAL_TIME)
sig_tme_HCCDB6_SP$STATUS [sig_tme_HCCDB6_SP$STATUS == 'Alive'] <-  0
sig_tme_HCCDB6_SP$STATUS [sig_tme_HCCDB6_SP$STATUS =='Dead'] <- 1
sig_tme_HCCDB6_SP$STATUS <- as.numeric(sig_tme_HCCDB6_SP$STATUS )
str(sig_tme_HCCDB6_SP)

save(list = ls(), file = 'HCCDB6.Rdada')

cycle_variable <- names(sig_tme_HCCDB6)[3:11]

library(survminer)
library(survival)

for (i in 1:length(cycle_variable)) {
  sur.cut <- surv_cutpoint(sig_tme_HCCDB6_SP, time= 'SURVIVAL_TIME',event = 'STATUS' , variables = cycle_variable[i], minprop = 0.1)
  #summary(sur.cut)
  sur.cat <- surv_categorize(sur.cut)
  #head(sur.cat)
  #table(sur.cat$AllCN32)
  names(sur.cat) <- c("OSday",'OS01','group')
  fit <- survfit(Surv(OSday, OS01) ~ group, data = sur.cat)
  filename <- paste0('./HCC bulk data Survival/HCCDB6//Optimal_RFS',"_",cycle_variable[i],".png")
  png(file=filename,width=1500,height=1500,res=300)
  p <- ggsurvplot(fit,palette = "npg",
                  risk.table = TRUE, pval = TRUE,
                  conf.int = FALSE, xlab="Time in Days",
                  ggtheme = theme_classic(),
                  title = cycle_variable[i])
  print(p)
  dev.off()
  print(cycle_variable[i])
}



##### correlation between signatures

library(ggpubr)

subytp_A <-  names(sig_tme_HCCDB6)[4:11]
subytp_B <-  "Mac_VIM_hi_FCtop50"   

#pearson correlation
for (i in seq_along(subytp_A)) {
  p <- ggscatter(sig_tme_HCCDB6_SP, x = subytp_A[i], y = subytp_B,
                 shape = 21, fill ='#8491B4FF',color='darkgrey',size = 2, # Points color, shape and size
                 add = "reg.line",  # Add regressin line
                 add.params = list(color = "#7E6148FF", fill = "#B09C85FF"), # Customize reg. line
                 conf.int = TRUE, # Add confidence interval
                 cor.coef = T, # Add correlation coefficient. see ?stat_cor
                 cor.coeff.args = list(method = "pearson", label.x.npc = "middle", label.y.npc = "top",label.sep = "\n"),
                 title = 'Pearson correlation')
  ggsave(p,file = paste0("./HCC bulk correlation/HCCDB6/Pearson_",subytp_A[i],".png"),width = 3.5,height = 3)
}




###### 12.2 HCCDB15 : also TCGA dataset ########
load(file='signatures for Bulk-seq validation.Rdata')


HCC_HCCDB15_dat <- read.table(file='E:\\11. CODEX\\LIHC bulk seq\\HCCDB15\\TCGA-HCC.log2_normalized.txt', sep='\t',header=T, row.names = 1)
dim(HCC_HCCDB15_dat) #colname is sample
HCC_HCCDB15_dat <- HCC_HCCDB15_dat[ !duplicated(HCC_HCCDB15_dat$Symbol), ]
row.names(HCC_HCCDB15_dat) <-  HCC_HCCDB15_dat$Symbol
HCC_HCCDB15_dat <- HCC_HCCDB15_dat[, -1]
# temp <- HCC_HCCDB15_dat %>% as.matrix()
# temp <- log2(temp + 1)
HCC_HCCDB15_dat_1 <- HCC_HCCDB15_dat

HCC_HCCDB15_dat_1_t <- HCC_HCCDB15_dat_1 %>% t() %>% as.data.frame()
HCC_HCCDB15_dat_1_t$ID <- row.names(HCC_HCCDB15_dat_1_t) 

sig_tme_HCCDB15 <- calculate_sig_score(pdata           = NULL,
                                       eset            = HCC_HCCDB15_dat_1,
                                       signature       = VIM_signature_tme,
                                       method          = "zscore",
                                       mini_gene_count = 3,
                                       adjust_eset = TRUE)
dim(sig_tme_HCCDB15)

###clinical information
df_HCCDB15_sample <- read.table(file='E:\\11. CODEX\\LIHC bulk seq\\HCCDB15\\HCCDB15.sample.txt', sep='\t') %>%  t() %>% as.data.frame()
df_HCCDB15_patient <- read.table(file='E:\\11. CODEX\\LIHC bulk seq\\HCCDB15\\HCCDB15.patient.txt', sep='\t') %>%  t() %>% as.data.frame()

colnames(df_HCCDB15_sample) <- df_HCCDB15_sample[1,]
df_HCCDB15_sample <- df_HCCDB15_sample[-1, ]
df_HCCDB15_sample$SAMPLE_ID <- gsub('-', '.', df_HCCDB15_sample$SAMPLE_ID )


colnames(df_HCCDB15_patient) <- df_HCCDB15_patient[1,]
df_HCCDB15_patient <- df_HCCDB15_patient[-1,]
df_HCCDB15_patient$PATIENT_ID <- gsub('-', '.', df_HCCDB15_patient$PATIENT_ID )


dim(df_HCCDB15_patient) #356 
dim(df_HCCDB15_sample) #405
dim(sig_tme_HCCDB15)   #400


df_sig_sample <- merge(sig_tme_HCCDB15, df_HCCDB15_sample, by.x='ID' ,by.y= 'SAMPLE_ID' )  ##400
df_sig_sample_tumor <- df_sig_sample %>% filter(TYPE == 'HCC')  ###351



###prognosis
dim(df_sig_sample_tumor) #351
df_sig_sample_tumor_patient <- left_join(df_sig_sample_tumor,df_HCCDB15_patient, by='PATIENT' )  ##351
dim(df_sig_sample_tumor_patient) #351

str(df_sig_sample_tumor_patient)
df_sig_sample_tumor_patient$SURVIVAL_TIME <- as.numeric(df_sig_sample_tumor_patient$SURVIVAL_TIME)
df_sig_sample_tumor_patient$STATUS [df_sig_sample_tumor_patient$STATUS == 'Alive'] <-  0
df_sig_sample_tumor_patient$STATUS [df_sig_sample_tumor_patient$STATUS =='Dead'] <- 1
df_sig_sample_tumor_patient$STATUS <- as.numeric(df_sig_sample_tumor_patient$STATUS )


save(list = ls(), file = 'HCCDB15.Rdata')


#######~~~~MacVim: poor prognosis#######
cycle_variable <- names(sig_tme_HCCDB15)[3:11]

library(survminer)
library(survival)

for (i in 1:length(cycle_variable)) {
  sur.cut <- surv_cutpoint(df_sig_sample_tumor_patient, time= 'SURVIVAL_TIME',event = 'STATUS' , variables = cycle_variable[i], minprop = 0.1)
  #summary(sur.cut)
  sur.cat <- surv_categorize(sur.cut)
  #head(sur.cat)
  #table(sur.cat$AllCN32)
  names(sur.cat) <- c("RFSday",'RFS01','group')
  fit <- survfit(Surv(RFSday, RFS01) ~ group, data = sur.cat)
  filename <- paste0('./HCC bulk data Survival/HCCDB15//Optimal_RFS',"_",cycle_variable[i],".png")
  png(file=filename,width=1500,height=1500,res=300)
  p <- ggsurvplot(fit,palette = "npg",
                  risk.table = TRUE, pval = TRUE,
                  conf.int = FALSE, xlab="Time in Days",
                  ggtheme = theme_classic(),
                  title = cycle_variable[i])
  print(p)
  dev.off()
  print(cycle_variable[i])
}



##### correlation between signatures

library(ggpubr)

subytp_A <-  names(sig_tme_HCCDB15)[4:11]
subytp_B <-  "Mac_VIM_hi_FCtop50"   

#pearson correlation
for (i in seq_along(subytp_A)) {
  p <- ggscatter(df_sig_sample_tumor_patient, x = subytp_A[i], y = subytp_B,
                 shape = 21, fill ='#8491B4FF',color='darkgrey',size = 2, # Points color, shape and size
                 add = "reg.line",  # Add regressin line
                 add.params = list(color = "#7E6148FF", fill = "#B09C85FF"), # Customize reg. line
                 conf.int = TRUE, # Add confidence interval
                 cor.coef = T, # Add correlation coefficient. see ?stat_cor
                 cor.coeff.args = list(method = "pearson", label.x.npc = "middle", label.y.npc = "top",label.sep = "\n"),
                 title = 'Pearson correlation')
  ggsave(p,file = paste0("./HCC bulk correlation/HCCDB15/Pearson_",subytp_A[i],".png"),width = 3.5,height = 3)
}












