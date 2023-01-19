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
library(pvclust) # v2.2-0
library(ComplexHeatmap) # v2.12.1
library(ggbreak) # v0.1.0
library(sva) # v3.44.0

#### V1: subtype: Macrophages ####
load("~/Desktop/Project/Codex/Analysis_07v/02_DefineCelltype/Outputs/codex_F_com_F.RData")
codex_F_com_F_macrophages <- subset(codex_F_com_F,idents="Macrophages") %>%
  Seurat::NormalizeData(verbose = FALSE,normalization.method = NULL) %>%
  ScaleData(verbose = FALSE)

codex_F_com_F_macrophages <- RunPCA(codex_F_com_F_macrophages, features=c("c-Myc","CD107A","CD4","CD44", 
                                                                          "HIF1a","HLA-DR","KI67","p-AMPK","p-mTOR", 
                                                                          "p-S6","p53","PD-L1"))
codex_F_com_F_macrophages <- RunUMAP(codex_F_com_F_macrophages, method = "umap-learn", metric = "correlation", dims = 1:length(codex_F_com_F_macrophages@reductions$pca))

codex_F_com_F_macrophages_phenograph <- subset(codex_F_com_F_macrophages,features=c("c-Myc","CD107A","CD4","CD44", 
                                                                              "HIF1a","HLA-DR","KI67","p-AMPK","p-mTOR", 
                                                                              "p-S6","p53","PD-L1"))
codex_F_com_F_macrophages_phenograph = ScaleData(codex_F_com_F_macrophages_phenograph)
codex_F_com_F_macrophages_phenograph <- codex_F_com_F_macrophages_phenograph@assays$CODEX@scale.data
# PG_elbow(t(codex_F_com_F_macrophages_phenograph), k_from = 10, k_to = 30, k_step=5) # cluster row
PhenoGraph_macrophages_result_k25 <-as.numeric(membership(Rphenograph(t(codex_F_com_F_macrophages_phenograph),k=25)))

save(PhenoGraph_macrophages_result_k25, file = "PhenoGraph_macrophages_result_k25.RData")

Idents(codex_F_com_F_macrophages) = PhenoGraph_macrophages_result_k25
levels(codex_F_com_F_macrophages) = as.character(1:25)

save(codex_F_com_F_macrophages,file = "codex_F_com_F_macrophages_PhenoGraph_PCA_UMAP.RData")

# Heatmap
codex_F_com_F_macrophages_SelectedFeatures <- subset(codex_F_com_F_macrophages,features=c("c-Myc","CD107A","CD4","CD44", 
                                                                              "HIF1a","HLA-DR","KI67","p-AMPK","p-mTOR", 
                                                                              "p-S6","p53","PD-L1"))
codex_F_com_F_macrophages_SelectedFeatures = ScaleData(codex_F_com_F_macrophages_SelectedFeatures)

p = DoHeatmap(subset(codex_F_com_F_macrophages_SelectedFeatures, downsample=100), slot = "scale.data", assay = "CODEX",
              features = rownames(codex_F_com_F_macrophages_SelectedFeatures), angle = 90, size = 4,disp.min = -2, disp.max = 2) + NoLegend()
pdf("20220208_Subtype_Macrophages_Heatmap_100.pdf", width = 12, height = 3)
p
dev.off()

p = DoHeatmap(subset(codex_F_com_F_macrophages_SelectedFeatures, downsample=5000), slot = "scale.data", assay = "CODEX",
              features = rownames(codex_F_com_F_macrophages_SelectedFeatures), angle = 90, size = 4,disp.min = -2, disp.max = 2) + NoLegend()
pdf("20220208_Subtype_Macrophages_Heatmap_5000.pdf", width = 12, height = 3)
p
dev.off()

# UMAP
codex_F_com_F_macrophages@meta.data$macrophages_phenograph = PhenoGraph_macrophages_result_k25
p = DimPlot(subset(codex_F_com_F_macrophages,downsample=500), 
            reduction = "umap", 
            label = FALSE, 
            repel = T,
            split.by = "macrophages_phenograph",
            ncol = 10,
            pt.size = 0.4) + NoLegend()

pdf("20220208_Subtype_Macrophages_Umap_Split.pdf", width = 12, height = 6)
p
dev.off()

# Cluster and Define subtype via pvclust and pheatmap
test_php_bind = data.frame()
for (i in 1:25) {
  temp = subset(codex_F_com_F_macrophages,idents = i,features=c("c-Myc","CD107A","CD4","CD44", 
                                                                "HIF1a","HLA-DR","KI67","p-AMPK","p-mTOR", 
                                                                "p-S6","p53","PD-L1"))
  temp = as.data.frame(temp@assays$CODEX@data)
  test_php_bind = bind_rows(test_php_bind, rowMeans(temp))
}
test_php = as.matrix(test_php_bind)
rownames(test_php) = 1:25

p_dat = t(scale(test_php)) # scale column
#Then we conduct hierarchical cluster analysis with multiscale bootstrap with number of bootstrap 1000, using average method and correlation-based dissimilarity matrix as follows:
result <- pvclust(p_dat, method.dist="euclidean", method.hclust="average", nboot=10000, parallel=TRUE) # cluster col
plot(result)

pheatmap(t(test_php),
         scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         clustering_method = "average",
         angle_col = "0",
         name = "Expression",
         color = c(colorRampPalette(colors = c("blue","white","red"))(100)))

bk <- c(seq(-2,-0.1,by=0.01),seq(0,6,by=0.01))
pheatmap(t(test_php),
         scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         clustering_method = "average",
         angle_col = "0",
         name = "Expression",
         color = c(colorRampPalette(colors = c("blue","white","red"))(100)),
         legend_breaks = seq(-2,6,2),
         breaks = bk)

bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pheatmap(t(test_php),
         scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         clustering_method = "average",
         angle_col = "0",
         name = "Expression",
         color = c(colorRampPalette(colors = c("blue","white","red"))(100)),
         legend_breaks = seq(-2,2,1),
         breaks = bk)

# Definition
new.cluster.ids <- c("Macrophages_HIF1a_int",#1
                     "Macrophages_CD107a+",#2 
                     "Macrophages_CD107a_int",#3
                     "Macrophages_others",#4
                     "Macrophages_HLA-DR_int",#5
                     "Macrophages_p-mTOR+",#6
                     "Macrophages_p-AMPK+",#7
                     "Macrophages_others",#8 
                     "Macrophages_KI67+",#9 
                     "Macrophages_p-AMPK+",#10
                     "Macrophages_HLA-DR_int",#11
                     "Macrophages_CD4+",#12
                     "Macrophages_p-S6_int",#13 
                     "Macrophages_CD44+",#14
                     "Macrophages_p-AMPK_int",#15
                     "Macrophages_CD4_int",#16
                     "Macrophages_HIF1a+",#17
                     "Macrophages_CD44_int",#18
                     "Macrophages_HLA-DR+",#19
                     "Macrophages_PD-L1+",#20
                     "Macrophages_p-S6+",#21
                     "Macrophages_c-Myc+",#22
                     "Macrophages_p53+",#23
                     "Macrophages_p-mTOR_int",#24
                     "Macrophages_HIF1a+"#25
)

names(new.cluster.ids) <- levels(codex_F_com_F_macrophages)
codex_F_com_F_macrophages <- RenameIdents(codex_F_com_F_macrophages, new.cluster.ids)
codex_F_com_F_macrophages@meta.data$MacrophagesSubtype = Idents(codex_F_com_F_macrophages)

# bind macrophages subtypes with other cells to further check and confirm via MAV
codex_F_com_F_others <- subset(codex_F_com_F,
                              idents="Macrophages",
                              invert = TRUE) %>%
  Seurat::NormalizeData(verbose = FALSE,normalization.method = NULL) %>%
  ScaleData(verbose = FALSE)

codex_F_com_F_others@meta.data$macrophages_phenograph = 0
codex_F_com_F_others@meta.data$MacrophagesSubtype = "cell_others"

codex_F_com_F@meta.data$cellid = as.numeric(rownames(codex_F_com_F@meta.data))
codex_F_com_F_macrophages@meta.data$cellid = as.numeric(rownames(codex_F_com_F_macrophages@meta.data))
codex_F_com_F_others@meta.data$cellid = as.numeric(rownames(codex_F_com_F_others@meta.data))

temp = rbind(codex_F_com_F_macrophages@meta.data,codex_F_com_F_others@meta.data)
temp = temp[order(temp$cellid),]

identical(codex_F_com_F@meta.data$TMA_reg,temp$TMA_reg)
identical(codex_F_com_F@meta.data$clusters,temp$clusters)

codex_F_com_F@meta.data$macrophages_phenograph = temp$macrophages_phenograph
codex_F_com_F@meta.data$MacrophagesSubtype = temp$MacrophagesSubtype

#### generate csv to MAV
#Read excel file from qupath into data frame
load("~/Desktop/Project/Codex/Analysis_07v/01_GenerateSeuratInput/Outputs/csd_raw.RData")
load("~/Desktop/Project/Codex/Analysis_07v/01_GenerateSeuratInput/Outputs/csd.RData")
load("~/Desktop/Project/Codex/Analysis_07v/01_GenerateSeuratInput/Outputs/FilterID.RData")

csd_F <- csd[FilterID,]
csd_raw_F <- csd_raw[FilterID,]

csd_F_F = csd_F[colnames(codex_F_com_F),]
csd_raw_F_F =  csd_raw_F[colnames(codex_F_com_F),]

Panel <- c("actin","aSMA","c_Myc","Caspase_3","CD107A","CD11C","CD163","CD20","CD21","CD3",
           "CD31","CD4","CD44","CD45","CD45RO","CD68","CD8","DAPI","E_Cadherin","FOXP3",
           "Glypican3","Hepar1","HIF1a","HistoneH3","HLA_DR","Keratin","KI67","p_AMPK","p_mTOR","p_S6",
           "p53","Pan_CK","PD_L1","PD1","Podoplanin","Twist1","Vimentin")

qupath_csd<-csd_raw_F_F
Panel_36 <- Panel[-18]
qupath_column_names<-paste0(Panel_36, " Nucleus Intensity")

#Isolate columns whose name we want to change I chose to use the entire cell. 
vxm2 = csd_F_F
#Replace data frame names with MAV compatible column names
names(vxm2)<- qupath_column_names
# Add columns to data frame for X Y coordinates following MAV required nomenclature. 
XMin=qupath_csd$`Centroid X px`#*2.65
XMax=qupath_csd$`Centroid X px`#*2.65
YMin=qupath_csd$`Centroid Y px`#*2.65
YMax=qupath_csd$`Centroid Y px`#*2.65
XMin=data.frame(XMin)
XMax=data.frame(XMax)
YMin=data.frame(YMin)
YMax=data.frame(YMax)
vxm2<- bind_cols(vxm2, XMin)
vxm2<- bind_cols(vxm2, XMax)
vxm2<- bind_cols(vxm2, YMin)
vxm2<- bind_cols(vxm2, YMax)

converttomav <- vxm2

#Combine dataframes
cellid = codex_F_com_F@meta.data$cellid
cellid <- data.frame(cellid)
names(cellid) = 'cellid'

MacrophagesSubtype<-codex_F_com_F@meta.data$MacrophagesSubtype
MacrophagesSubtype <- data.frame(MacrophagesSubtype)
names(MacrophagesSubtype) = 'MacrophagesSubtype'

macrophages_phenograph = codex_F_com_F@meta.data$macrophages_phenograph
macrophages_phenograph <- data.frame(macrophages_phenograph)
names(macrophages_phenograph) = 'macrophages_phenograph'

gc_csd_final4<-bind_cols(converttomav, cellid, macrophages_phenograph, MacrophagesSubtype)
gc_csd_final4<-gc_csd_final4 %>%  mutate_all(as.character)
gc_csd_final4$Class  <- csd_raw_F_F$Class

#name the cell seg file. You need one cell seg file per region. For region 1 the cell seg file should be name reg001_
write.table(gc_csd_final4, file = "testing_allReg_seurat_._analysis_20220210_MacrophagesFeatures_MacrophagesSubtype.csv", sep=",", col.names = NA, row.names = TRUE)

save.image("~/Desktop/Project/Codex/Analysis_07v/04_DefineMacrophageSubtype/Outputs_v1/07_DefineMacrophageSubtype_V1.RData")

#### V2: subtype: Macrophages: modification on V1 after checking back on MAV ####
load("~/Desktop/Project/Codex/Analysis_07v/04_DefineMacrophageSubtype/Outputs_v1/codex_F_com_F_macrophages_PhenoGraph_PCA_UMAP.RData")

codex_F_com_F_macrophages_F <- RunPCA(codex_F_com_F_macrophages, features=c("HLA-DR","CD4","CD44","CD107A","HIF1a","KI67","p-S6","c-Myc","CD11C","Vimentin"))
codex_F_com_F_macrophages_F <- RunUMAP(codex_F_com_F_macrophages_F, method = "umap-learn", metric = "correlation", dims = 1:length(codex_F_com_F_macrophages_F@reductions$pca))

codex_F_com_F_macrophages_F_phenograph <- subset(codex_F_com_F_macrophages_F,features=c("HLA-DR","CD4","CD44","CD107A","HIF1a","KI67","p-S6","c-Myc","CD11C","Vimentin"))
codex_F_com_F_macrophages_F_phenograph = ScaleData(codex_F_com_F_macrophages_F_phenograph)
codex_F_com_F_macrophages_F_phenograph <- codex_F_com_F_macrophages_F_phenograph@assays$CODEX@scale.data
# PG_elbow(t(codex_F_com_F_macrophages_F_phenograph), k_from = 10, k_to = 30, k_step=5) # cluster row
PhenoGraph_Ma_F_result_k25 <-as.numeric(membership(Rphenograph(t(codex_F_com_F_macrophages_F_phenograph),k=25)))

save(PhenoGraph_Ma_F_result_k25, file = "PhenoGraph_Ma_F_result_k25.RData")
save(codex_F_com_F_macrophages_F, file = "codex_F_com_F_macrophage_F_PCA_UMAP_beforeIdents.RData")

Idents(codex_F_com_F_macrophages_F) = PhenoGraph_Ma_F_result_k25
levels(codex_F_com_F_macrophages_F) = as.character(1:27)

# Heatmap
codex_F_com_F_macrophages_F_SelectedFeatures <- subset(codex_F_com_F_macrophages_F,features=c("HLA-DR","CD4","CD44","CD107A","HIF1a","KI67","p-S6","c-Myc","CD11C","Vimentin"))
codex_F_com_F_macrophages_F_SelectedFeatures = ScaleData(codex_F_com_F_macrophages_F_SelectedFeatures)

p = DoHeatmap(subset(codex_F_com_F_macrophages_F_SelectedFeatures, downsample=100), slot = "scale.data", assay = "CODEX",
              features = rownames(codex_F_com_F_macrophages_F_SelectedFeatures), angle = 90, size = 4,disp.min = -2, disp.max = 2) + NoLegend()
pdf("20220308_SFig1_Subtype_Ma_Heatmap_100_R.pdf", width = 12, height = 3)
p
dev.off()

p = DoHeatmap(subset(codex_F_com_F_macrophages_F_SelectedFeatures, downsample=5000), slot = "scale.data", assay = "CODEX",
              features = rownames(codex_F_com_F_macrophages_F_SelectedFeatures), angle = 90, size = 4,disp.min = -2, disp.max = 2) + NoLegend()
pdf("20220308_SFig1_Subtype_Ma_Heatmap_5000_R.pdf", width = 12, height = 3)
p
dev.off()

# UMAP
codex_F_com_F_macrophages_F@meta.data$macrophages_phenograph = PhenoGraph_Ma_F_result_k25
p = DimPlot(subset(codex_F_com_F_macrophages_F,downsample=500), 
            reduction = "umap", 
            label = FALSE, 
            repel = T,
            split.by = "macrophages_phenograph",
            ncol = 9,
            pt.size = 0.4) + NoLegend()

pdf("20220308_SFig1_Subtype_Ma_Umap_Split_R.pdf", width = 12, height = 6)
p
dev.off()

# Cluster and Define subtype via pvclust and pheatmap
test_php_bind = data.frame()
for (i in 1:27) {
  temp = subset(codex_F_com_F_macrophages_F,idents = i,features=c("HLA-DR","CD4","CD44","CD107A","HIF1a","KI67","p-S6","c-Myc","CD11C","Vimentin"))
  temp = as.data.frame(temp@assays$CODEX@data)
  test_php_bind = bind_rows(test_php_bind, rowMeans(temp))
}
test_php = as.matrix(test_php_bind)
rownames(test_php) = 1:27

p_dat = t(scale(test_php)) # scale column
#Then we conduct hierarchical cluster analysis with multiscale bootstrap with number of bootstrap 1000, using average method and correlation-based dissimilarity matrix as follows:
result <- pvclust(p_dat, method.dist="euclidean", method.hclust="average", nboot=10000, parallel=TRUE) # cluster col
plot(result)

saveRDS(test_php,file = "TrainingV3_macrosub_test_php.Rds")

pheatmap(t(test_php),
         scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         clustering_method = "average",
         angle_col = "0",
         name = "Expression",
         color = c(colorRampPalette(colors = c("blue","white","red"))(100)))

bk <- c(seq(-2,-0.1,by=0.01),seq(0,6,by=0.01))
pheatmap(t(test_php),
         scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         clustering_method = "average",
         angle_col = "0",
         name = "Expression",
         color = c(colorRampPalette(colors = c("blue","white","red"))(100)),
         legend_breaks = seq(-2,6,2),
         breaks = bk)

bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pheatmap(t(test_php),
         scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         clustering_method = "average",
         angle_col = "0",
         name = "Expression",
         color = c(colorRampPalette(colors = c("blue","white","red"))(100)),
         legend_breaks = seq(-2,2,1),
         breaks = bk)


bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pheatmap(t(test_php),
         scale = "none",
         cluster_rows = T,
         cluster_cols = T,
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         clustering_method = "average",
         angle_col = "0",
         name = "Expression",
         color = c(colorRampPalette(colors = c("blue","white","red"))(100)),
         legend_breaks = seq(-2,2,1),
         breaks = bk)

# Definition
new.cluster.ids <- c("Macrophages_Vimentin+",#1
                     "Macrophages_HLA-DR+",#2 
                     "Macrophages_CD11c+_HLA-DR+",#3
                     "Macrophages_CD11c+_HLA-DR+",#4
                     "Macrophages_p-S6+",#5
                     "Macrophages_KI67+",#6
                     "Macrophages_CD11c+_CD44+_KI67+",#7
                     "Macrophages_others",#8 
                     "Macrophages_c-Myc+",#9 
                     "Macrophages_others",#10
                     "Macrophages_Vimentin+",#11
                     "Macrophages_CD4+",#12
                     "Macrophages_CD44+",#13 
                     "Macrophages_HLA-DR+",#14
                     "Macrophages_CD44+",#15
                     "Macrophages_CD4+",#16
                     "Macrophages_others",#17
                     "Macrophages_others",#18
                     "Macrophages_HIF1a+",#19
                     "Macrophages_CD107a+",#20
                     "Macrophages_HIF1a+_c-Myc+",#21
                     "Macrophages_HIF1a+",#22
                     "Macrophages_CD107a+",#23
                     "Macrophages_others",#24
                     "Macrophages_HIF1a+_Vimentin+",#25
                     "Macrophages_others",#26
                     "Macrophages_CD11c+_CD44+"#27
)

names(new.cluster.ids) <- levels(codex_F_com_F_macrophages_F)
codex_F_com_F_macrophages_F <- RenameIdents(codex_F_com_F_macrophages_F, new.cluster.ids)
codex_F_com_F_macrophages_F@meta.data$MacrophagesSubtype = Idents(codex_F_com_F_macrophages_F)

# bind tumor subtypes with all stromal cells to further check and confirm via MAV
load("~/Desktop/Project/Codex/Analysis_07v/02_DefineCelltype/Outputs/codex_F_com_F.RData") # Only scale the object after subset; Not redo PCA and UMAP!
codex_F_com_F_others <- subset(codex_F_com_F,
                                idents="Macrophages",
                                invert = TRUE) %>%
  Seurat::NormalizeData(verbose = FALSE,normalization.method = NULL) %>%
  ScaleData(verbose = FALSE)

codex_F_com_F_others@meta.data$macrophages_phenograph = 0
codex_F_com_F_others@meta.data$MacrophagesSubtype = "cell_others"

codex_F_com_F@meta.data$cellid = as.numeric(rownames(codex_F_com_F@meta.data))
codex_F_com_F_macrophages_F@meta.data$cellid = as.numeric(rownames(codex_F_com_F_macrophages_F@meta.data))
codex_F_com_F_others@meta.data$cellid = as.numeric(rownames(codex_F_com_F_others@meta.data))

temp = rbind(codex_F_com_F_macrophages_F@meta.data,codex_F_com_F_others@meta.data)
temp = temp[order(temp$cellid),]

identical(codex_F_com_F@meta.data$TMA_reg,temp$TMA_reg)
identical(codex_F_com_F@meta.data$clusters,temp$clusters)

codex_F_com_F@meta.data$macrophages_phenograph = temp$macrophages_phenograph
codex_F_com_F@meta.data$MacrophagesSubtype = temp$MacrophagesSubtype


#Read excel file from qupath into data frame
Panel <- c("actin","aSMA","c_Myc","Caspase_3","CD107A","CD11C","CD163","CD20","CD21","CD3",
           "CD31","CD4","CD44","CD45","CD45RO","CD68","CD8","DAPI","E_Cadherin","FOXP3",
           "Glypican3","Hepar1","HIF1a","HistoneH3","HLA_DR","Keratin","KI67","p_AMPK","p_mTOR","p_S6",
           "p53","Pan_CK","PD_L1","PD1","Podoplanin","Twist1","Vimentin")
load("~/Desktop/Project/Codex/Analysis_07v/01_GenerateSeuratInput/Outputs/csd_raw.RData")
load("~/Desktop/Project/Codex/Analysis_07v/01_GenerateSeuratInput/Outputs/csd.RData")
load("~/Desktop/Project/Codex/Analysis_07v/01_GenerateSeuratInput/Outputs/FilterID.RData")

csd_F <- csd[FilterID,]
csd_raw_F <- csd_raw[FilterID,]

csd_F_F = csd_F[colnames(codex_F_com_F),]
csd_raw_F_F =  csd_raw_F[colnames(codex_F_com_F),]

qupath_csd<-csd_raw_F_F
Panel_36 <- Panel[-18]
qupath_column_names<-paste0(Panel_36, " Nucleus Intensity")

#Isolate columns whose name we want to change I chose to use the entire cell. 
vxm2 = csd_F_F
#Replace data frame names with MAV compatible column names
names(vxm2)<- qupath_column_names
# Add columns to data frame for X Y coordinates following MAV required nomenclature. 
XMin=qupath_csd$`Centroid X px`#*2.65
XMax=qupath_csd$`Centroid X px`#*2.65
YMin=qupath_csd$`Centroid Y px`#*2.65
YMax=qupath_csd$`Centroid Y px`#*2.65
XMin=data.frame(XMin)
XMax=data.frame(XMax)
YMin=data.frame(YMin)
YMax=data.frame(YMax)
vxm2<- bind_cols(vxm2, XMin)
vxm2<- bind_cols(vxm2, XMax)
vxm2<- bind_cols(vxm2, YMin)
vxm2<- bind_cols(vxm2, YMax)

converttomav <- vxm2

#Combine dataframes
cellid = codex_F_com_F@meta.data$cellid
cellid <- data.frame(cellid)
names(cellid) = 'cellid'

MacrophagesSubtype<-codex_F_com_F@meta.data$MacrophagesSubtype
MacrophagesSubtype <- data.frame(MacrophagesSubtype)
names(MacrophagesSubtype) = 'MacrophagesSubtype'

macrophages_phenograph = codex_F_com_F@meta.data$macrophages_phenograph
macrophages_phenograph <- data.frame(macrophages_phenograph)
names(macrophages_phenograph) = 'macrophages_phenograph'

gc_csd_final4<-bind_cols(converttomav, cellid, macrophages_phenograph, MacrophagesSubtype)
gc_csd_final4<-gc_csd_final4 %>%  mutate_all(as.character)
gc_csd_final4$Class  <- csd_raw_F_F$Class

#name the cell seg file. You need one cell seg file per region. For region 1 the cell seg file should be name reg001_
write.table(gc_csd_final4, file = "testing_allReg_seurat_._analysis_20220309_MacroFeatures_MacroSubtype_v2.csv", sep=",", col.names = NA, row.names = TRUE)

#### V3: subtype: Macrophages: modification on V2 after checking back on MAV ####
load("~/Desktop/Project/Codex/Analysis_07v/04_DefineMacrophageSubtype/Outputs_v2/07_DefineMacrophageSubtype_V2.RData")
Idents(codex_F_com_F_macrophages_F) = codex_F_com_F_macrophages_F@meta.data$macrophages_phenograph

new.cluster.ids <- c("Macrophages_Vimentin+",#1
                     "Macrophages_HLA-DR+",#2 
                     "Macrophages_CD11c+_HLA-DR+",#3
                     "Macrophages_CD11c+_HLA-DR+",#4
                     "Macrophages_p-S6+",#5
                     "Macrophages_KI67+",#6
                     "Macrophages_artifact",#7
                     "Macrophages_others",#8 
                     "Macrophages_c-Myc+",#9 
                     "Macrophages_others",#10
                     "Macrophages_Vimentin+",#11
                     "Macrophages_CD4+",#12
                     "Macrophages_CD44+",#13 
                     "Macrophages_HLA-DR+",#14
                     "Macrophages_CD44+",#15
                     "Macrophages_CD4+",#16
                     "Macrophages_p-S6+",#17
                     "Macrophages_others",#18
                     "Macrophages_HIF1a+",#19
                     "Macrophages_CD107a+",#20
                     "Macrophages_HIF1a+_c-Myc+",#21
                     "Macrophages_HIF1a+",#22
                     "Macrophages_CD107a+",#23
                     "Macrophages_HIF1a+",#24
                     "Macrophages_HIF1a+_Vimentin+",#25
                     "Macrophages_others",#26
                     "Macrophages_artifact"#27
)

names(new.cluster.ids) <- levels(codex_F_com_F_macrophages_F)
codex_F_com_F_macrophages_F <- RenameIdents(codex_F_com_F_macrophages_F, new.cluster.ids)
codex_F_com_F_macrophages_F@meta.data$MacrophagesSubtype = Idents(codex_F_com_F_macrophages_F)

# bind tumor subtypes with all stromal cells to further check and confirm via MAV
codex_F_com_F_macrophages_F@meta.data$cellid = as.numeric(rownames(codex_F_com_F_macrophages_F@meta.data))

temp = rbind(codex_F_com_F_macrophages_F@meta.data,codex_F_com_F_others@meta.data)
temp = temp[order(temp$cellid),]

identical(codex_F_com_F@meta.data$TMA_reg,temp$TMA_reg)
identical(codex_F_com_F@meta.data$clusters,temp$clusters)

codex_F_com_F@meta.data$macrophages_phenograph = temp$macrophages_phenograph
codex_F_com_F@meta.data$MacrophagesSubtype = temp$MacrophagesSubtype

#### generate csv to MAV
#Read excel file from qupath into data frame
Panel <- c("actin","aSMA","c_Myc","Caspase_3","CD107A","CD11C","CD163","CD20","CD21","CD3",
           "CD31","CD4","CD44","CD45","CD45RO","CD68","CD8","DAPI","E_Cadherin","FOXP3",
           "Glypican3","Hepar1","HIF1a","HistoneH3","HLA_DR","Keratin","KI67","p_AMPK","p_mTOR","p_S6",
           "p53","Pan_CK","PD_L1","PD1","Podoplanin","Twist1","Vimentin")
load("~/Desktop/Project/Codex/Analysis_07v/01_GenerateSeuratInput/Outputs/csd_raw.RData")
load("~/Desktop/Project/Codex/Analysis_07v/01_GenerateSeuratInput/Outputs/csd.RData")
load("~/Desktop/Project/Codex/Analysis_07v/01_GenerateSeuratInput/Outputs/FilterID.RData")

csd_F <- csd[FilterID,]
csd_raw_F <- csd_raw[FilterID,]

csd_F_F = csd_F[colnames(codex_F_com_F),]
csd_raw_F_F =  csd_raw_F[colnames(codex_F_com_F),]

qupath_csd<-csd_raw_F_F
Panel_36 <- Panel[-18]
qupath_column_names<-paste0(Panel_36, " Nucleus Intensity")

#Isolate columns whose name we want to change I chose to use the entire cell. 
vxm2 = csd_F_F
#Replace data frame names with MAV compatible column names
names(vxm2)<- qupath_column_names
# Add columns to data frame for X Y coordinates following MAV required nomenclature. 
XMin=qupath_csd$`Centroid X px`#*2.65
XMax=qupath_csd$`Centroid X px`#*2.65
YMin=qupath_csd$`Centroid Y px`#*2.65
YMax=qupath_csd$`Centroid Y px`#*2.65
XMin=data.frame(XMin)
XMax=data.frame(XMax)
YMin=data.frame(YMin)
YMax=data.frame(YMax)
vxm2<- bind_cols(vxm2, XMin)
vxm2<- bind_cols(vxm2, XMax)
vxm2<- bind_cols(vxm2, YMin)
vxm2<- bind_cols(vxm2, YMax)

converttomav <- vxm2

#Combine dataframes
cellid = codex_F_com_F@meta.data$cellid
cellid <- data.frame(cellid)
names(cellid) = 'cellid'

MacrophagesSubtype<-codex_F_com_F@meta.data$MacrophagesSubtype
MacrophagesSubtype <- data.frame(MacrophagesSubtype)
names(MacrophagesSubtype) = 'MacrophagesSubtype'

macrophages_phenograph = codex_F_com_F@meta.data$macrophages_phenograph
macrophages_phenograph <- data.frame(macrophages_phenograph)
names(macrophages_phenograph) = 'macrophages_phenograph'

gc_csd_final4<-bind_cols(converttomav, cellid, macrophages_phenograph, MacrophagesSubtype)
gc_csd_final4<-gc_csd_final4 %>%  mutate_all(as.character)
gc_csd_final4$Class  <- csd_raw_F_F$Class

#name the cell seg file. You need one cell seg file per region. For region 1 the cell seg file should be name reg001_
write.table(gc_csd_final4, file = "testing_allReg_seurat_._analysis_20220314_MacroFeatures_MacroSubtype_v3.csv", sep=",", col.names = NA, row.names = TRUE)

MacrophageSubtype = gc_csd_final4
save(MacrophageSubtype, file = "20220315_MacrophageSubtype.RData")


#### Subtype: UMAP and Heatmap ####
## UMAP
codex_F_com_F_macrophages_F <- subset(codex_F_com_F_macrophages_F,idents=c("Macrophages_artifact"),invert = TRUE) %>%
  Seurat::NormalizeData(verbose = FALSE,normalization.method = NULL) %>%
  ScaleData(verbose = FALSE)

codex_F_com_F_macrophages_F <- RunPCA(codex_F_com_F_macrophages_F, features=c("HLA-DR","CD4","CD44","CD107A","HIF1a","KI67","p-S6","c-Myc","CD11C","Vimentin"))
codex_F_com_F_macrophages_F <- RunUMAP(codex_F_com_F_macrophages_F, method = "umap-learn", metric = "correlation", dims = 1:length(codex_F_com_F_macrophages_F@reductions$pca))

subset_list = list()
downsample_cycle = round(table(Idents(codex_F_com_F_macrophages_F))/20,0)
for (i in 1:13) {
  temp = subset(codex_F_com_F_macrophages_F, idents = names(downsample_cycle)[i])
  subset_list[[i]] = subset(temp, downsample = downsample_cycle[i])
}
codex_F_com_F_subset_sample_prob = merge(x = subset_list[[1]], y = unlist(subset_list)[2:13], merge.data = T) # no scale data
codex_F_com_F_subset_sample_prob = subset(codex_F_com_F_macrophages_F, cells = colnames(codex_F_com_F_subset_sample_prob))

cols = c("#FA339A", # 1
         "#F6D026", # 2
         "#74F882", # 3
         "#7E8EE3", # 4
         "#8358CF", # 5
         "#27A7FA", # 6
         "#EE6600", # 7
         "#D083FF", # 8
         "#62FFFF", # 9
         "#AE2519", # 10
         "#9C9B99", # 11
         "#162C9B", # 12
         "#FBF88B"# 13
         #"#B58794" # 14
         #"#149684", # 15
         #"#D92B00", # 16
         #"#73904D", # 17
         #"#406E30", # 18
         #"#E3D78C", # 19
         #"#0F0F0F", # 20
         # "#336699", # 21
         # "#669999", # 22
         # "#6600ff", # 23
         # "#330066", # 24
         #"#cc99cc", # 25
         # "#666633", # 26
         #"#ffcccc", # 27
         # "#cc6600", # 28 
         # "#5b7192", # 29
         #"#99cc99" # 30
)

orders = names(table(Idents(codex_F_com_F_macrophages_F)))

pal<-colorRampPalette(cols)
image(x=1:13,y=1,z=as.matrix(1:13),col=pal(13))

p = DimPlot(codex_F_com_F_subset_sample_prob, 
            reduction = "umap", 
            group.by = "MacrophagesSubtype",
            label = FALSE, 
            repel = T,
            cols = cols,
            order = rev(orders),
            pt.size = 0.4) + NoLegend()
# theme(legend.position = "right") + 
# guides(color = guide_legend(ncol = 2, byrow = T))

pdf("20220417_MacrophagesSubtype_Umap_Merged_1in20.pdf", width = 5, height = 5)
p
dev.off()

## Freq
Freq = round((table(Idents(codex_F_com_F_macrophages_F))/nrow(codex_F_com_F_macrophages_F@meta.data))*100,2)
Freq = as.data.frame(Freq)
colnames(Freq)[1] = "Subtype"

p = ggplot(data = Freq, aes(x = Subtype, y = Freq, fill = Subtype)) + 
  geom_col(width = 0.8, position = position_dodge(0.65)) +
  scale_fill_manual(values = c("#FA339A", # 1
                               "#F6D026", # 2
                               "#74F882", # 3
                               "#7E8EE3", # 4
                               "#8358CF", # 5
                               "#27A7FA", # 6
                               "#EE6600", # 7
                               "#D083FF", # 8
                               "#62FFFF", # 9
                               "#AE2519", # 10
                               "#9C9B99", # 11
                               "#162C9B", # 12
                               "#FBF88B" # 13
                               #"#B58794" # 14
                               #"#149684", # 15
                               #"#D92B00", # 16
                               #"#73904D", # 17
                               #"#406E30", # 18
                               #"#E3D78C", # 19
                               #"#0F0F0F", # 20
                               # "#336699", # 21
                               # "#669999", # 22
                               # "#6600ff", # 23
                               # "#330066", # 24
                               #"#cc99cc", # 25
                               # "#666633", # 26
                               #"#ffcccc", # 27
                               # "#cc6600", # 28 
                               # "#5b7192", # 29
                               #"#99cc99" # 30
  )) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + NoLegend()+
  geom_text(aes(label = paste0(Freq,"%")), vjust = -.5) + 
  ylim(0,15)

pdf("20220417_heatmap_barplot_NoLegend.pdf", width = 10, height = 5)
p
dev.off()


## Heatmap
temp_subset = subset(codex_F_com_F_macrophages_F, features=c("HLA-DR","CD4","CD44","CD107A","HIF1a","KI67","p-S6","c-Myc","CD11C","Vimentin"))
temp_subset = subset(temp_subset, downsample=2000)
temp_subset = ScaleData(temp_subset)

toplot_heatmap = temp_subset@assays$CODEX@scale.data

annotation_col = as.data.frame(temp_subset@meta.data)
annotation_col = subset(annotation_col, select = "MacrophagesSubtype")
annotation_col = data.frame(cellid = rownames(annotation_col),
                            MacrophagesSubtype = annotation_col$MacrophagesSubtype)

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#6EC388", "white", "#C82033")) 


ha = HeatmapAnnotation(MacrophagesSubtype = annotation_col$MacrophagesSubtype,
                       col = list(MacrophagesSubtype = c("Macrophages_Vimentin+" = "#FA339A",
                                                   "Macrophages_HLA-DR+" = "#F6D026",
                                                   "Macrophages_CD11c+_HLA-DR+" = "#74F882",
                                                   "Macrophages_p-S6+" = "#7E8EE3",
                                                   "Macrophages_KI67+" = "#8358CF",
                                                   "Macrophages_others" = "#27A7FA",
                                                   "Macrophages_c-Myc+" = "#EE6600",
                                                   "Macrophages_CD4+" = "#D083FF",
                                                   "Macrophages_CD44+" = "#62FFFF",
                                                   "Macrophages_HIF1a+" = "#AE2519",
                                                   "Macrophages_CD107a+" = "#9C9B99",
                                                   "Macrophages_HIF1a+_c-Myc+" = "#162C9B",
                                                   "Macrophages_HIF1a+_Vimentin+" = "#FBF88B")),
                       show_legend = F)

ht = Heatmap(toplot_heatmap,
             name = "Expression",
             col = col_fun,
             width = nrow(toplot_heatmap)*unit(12, "mm"), 
             height = nrow(toplot_heatmap)*unit(7, "mm"),
             column_split = annotation_col$MacrophagesSubtype,
             column_title = NULL,
             cluster_rows = F,
             cluster_columns = F,
             show_column_names = F,
             top_annotation = ha,
             column_gap = unit(2,"mm"))

pdf("20220417_heatmap_2000.pdf", width = 13, height = 5)
ht
dev.off()





