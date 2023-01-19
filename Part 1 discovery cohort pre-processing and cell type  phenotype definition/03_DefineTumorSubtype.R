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


#### V1: subtype: tumor cells ####
load("~/Desktop/Project/Codex/Analysis_07v/02_DefineCelltype/Outputs/codex_F_com_F.RData")

codex_F_com_F_tumor <- subset(codex_F_com_F,idents="Tumor cells") %>%
  Seurat::NormalizeData(verbose = FALSE,normalization.method = NULL) %>%
  ScaleData(verbose = FALSE)

codex_F_com_F_tumor <- RunPCA(codex_F_com_F_tumor, features=c("c-Myc","CD107A","CD44","E-Cadherin",
                                                      "Glypican3","Hepar1","HIF1a","HLA-DR","KI67",
                                                      "p-AMPK","p-mTOR","p-S6","p53","PD-L1","Twist1","Vimentin"))
codex_F_com_F_tumor <- RunUMAP(codex_F_com_F_tumor, method = "umap-learn", metric = "correlation", dims = 1:length(codex_F_com_F_tumor@reductions$pca))

codex_F_com_F_tumor_phenograph <- subset(codex_F_com_F_tumor,features=c("c-Myc","CD107A","CD44","E-Cadherin",
                                                     "Glypican3","Hepar1","HIF1a","HLA-DR","KI67",
                                                     "p-AMPK","p-mTOR","p-S6","p53","PD-L1","Twist1","Vimentin"))
codex_F_com_F_tumor_phenograph = ScaleData(codex_F_com_F_tumor_phenograph)
codex_F_com_F_tumor_phenograph <- codex_F_com_F_tumor_phenograph@assays$CODEX@scale.data
# PG_elbow(t(codex_F_com_F_tumor_phenograph), k_from = 10, k_to = 30, k_step=5) # cluster row
PhenoGraph_tumor_result_k25 <-as.numeric(membership(Rphenograph(t(codex_F_com_F_tumor_phenograph),k=25)))

save(PhenoGraph_tumor_result_k25, file = "PhenoGraph_tumor_result_k25.RData")

Idents(codex_F_com_F_tumor) = PhenoGraph_tumor_result_k25
levels(codex_F_com_F_tumor) = as.character(1:39)

save(codex_F_com_F_tumor,file = "codex_F_com_F_tumor_PhenoGraph_PCA_UMAP.RData")

# Heatmap
codex_F_com_F_tumor_SelectedFeatures <- subset(codex_F_com_F_tumor,features=c("c-Myc","CD107A","CD44","E-Cadherin",
                                                                        "Glypican3","Hepar1","HIF1a","HLA-DR","KI67",
                                                                        "p-AMPK","p-mTOR","p-S6","p53","PD-L1","Twist1","Vimentin"))
codex_F_com_F_tumor_SelectedFeatures = ScaleData(codex_F_com_F_tumor_SelectedFeatures)

p = DoHeatmap(subset(codex_F_com_F_tumor_SelectedFeatures, downsample=100), slot = "scale.data", assay = "CODEX",
              features = rownames(codex_F_com_F_tumor_SelectedFeatures), angle = 90, size = 4,disp.min = -2, disp.max = 2) + NoLegend()
pdf("20220124_SFig1_Subtype_Tumor_Heatmap_100.pdf", width = 12, height = 3)
p
dev.off()

p = DoHeatmap(subset(codex_F_com_F_tumor_SelectedFeatures, downsample=5000), slot = "scale.data", assay = "CODEX",
              features = rownames(codex_F_com_F_tumor_SelectedFeatures), angle = 90, size = 4,disp.min = -2, disp.max = 2) + NoLegend()
pdf("20220124_SFig1_Subtype_Tumor_Heatmap_5000.pdf", width = 12, height = 3)
p
dev.off()

# UMAP
codex_F_com_F_tumor@meta.data$tumor_phenograph = PhenoGraph_tumor_result_k25
p = DimPlot(subset(codex_F_com_F_tumor,downsample=500), 
            reduction = "umap", 
            #group.by = "tumor_phenograph",
            label = FALSE, 
            repel = T,
            split.by = "tumor_phenograph",
            ncol = 10,
            pt.size = 0.4) + NoLegend()

pdf("20220124_SFig1_Subtype_Tumor_Umap_Split.pdf", width = 12, height = 6)
p
dev.off()

# Cluster and Define subtype via pvclust and pheatmap
test_php_bind = data.frame()
for (i in 1:39) {
  temp = subset(codex_F_com_F_tumor,idents = i,features=c("c-Myc","CD107A","CD44","E-Cadherin",
                                                          "Glypican3","Hepar1","HIF1a","HLA-DR","KI67",
                                                          "p-AMPK","p-mTOR","p-S6","p53","PD-L1","Twist1","Vimentin"))
  temp = as.data.frame(temp@assays$CODEX@data)
  test_php_bind = bind_rows(test_php_bind, rowMeans(temp))
}
test_php = as.matrix(test_php_bind)
rownames(test_php) = 1:39

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
new.cluster.ids <- c("Tumor_others",#1
                     "Tumor_c-Myc+",#2 
                     "Tumor_others",#3
                     "Tumor_others",#4
                     "Tumor_Hepar1+",#5
                     "Tumor_CD44+_E-Cadherin+",#6
                     "Tumor_Hepar1_int_Twist1_int",#7
                     "Tumor_Twist1+",#8 
                     "Tumor_p53+_p-mTOR+_PD-L1+",#9 
                     "Tumor_PD-L1_int",#10
                     "Tumor_others",#11
                     "Tumor_p-AMPK+_c-Myc+_Glypican3+_Hepar1+",#12
                     "Tumor_Hepar1_int_p-AMPK_int",#13 
                     "Tumor_KI67+",#14
                     "Tumor_others",#15
                     "Tumor_Vimentin+",#16
                     "Tumor_p53+",#17
                     "Tumor_Glypican3+_Hepar1+",#18
                     "Tumor_HLA-DR_int_Vimentin_int",#19
                     "Tumor_others",#20
                     "Tumor_p-AMPK_int",#21
                     "Tumor_p-AMPK+",#22
                     "Tumor_p-mTOR_int",#23
                     "Tumor_E-Cadherin+",#24
                     "Tumor_HLA-DR",#25
                     "Tumor_CD107A+",#26
                     "Tumor_CD44+",#27
                     "Tumor_E-Cadherin+",#28
                     "Tumor_HIF1a+",#29
                     "Tumor_p-S6+",#30
                     "Tumor_others",#31
                     "Tumor_others",#32
                     "Tumor_p-S6_int",#33
                     "Tumor_CD107A_int",#34
                     "Tumor_Glypican3+",#35
                     "Tumor_others",#36
                     "Tumor_c-Myc+",#37
                     "Tumor_HIF1a_int",#38    
                     "Tumor_others"#39
)

names(new.cluster.ids) <- levels(codex_F_com_F_tumor)
codex_F_com_F_tumor <- RenameIdents(codex_F_com_F_tumor, new.cluster.ids)
codex_F_com_F_tumor@meta.data$TumorSubtype = Idents(codex_F_com_F_tumor)

# bind tumor subtypes with all stromal cells to further check and confirm via MAV
codex_F_com_F_stromal <- subset(codex_F_com_F,
                              idents="Tumor cells",
                              invert = TRUE) %>%
  Seurat::NormalizeData(verbose = FALSE,normalization.method = NULL) %>%
  ScaleData(verbose = FALSE)

codex_F_com_F_stromal@meta.data$tumor_phenograph = 0
codex_F_com_F_stromal@meta.data$TumorSubtype = "Stromal"

codex_F_com_F@meta.data$cellid = as.numeric(rownames(codex_F_com_F@meta.data))
codex_F_com_F_tumor@meta.data$cellid = as.numeric(rownames(codex_F_com_F_tumor@meta.data))
codex_F_com_F_stromal@meta.data$cellid = as.numeric(rownames(codex_F_com_F_stromal@meta.data))

temp = rbind(codex_F_com_F_tumor@meta.data,codex_F_com_F_stromal@meta.data)
temp = temp[order(temp$cellid),]

identical(codex_F_com_F@meta.data$TMA_reg,temp$TMA_reg)
identical(codex_F_com_F@meta.data$clusters,temp$clusters)

codex_F_com_F@meta.data$tumor_phenograph = temp$tumor_phenograph
codex_F_com_F@meta.data$TumorSubtype = temp$TumorSubtype

#### generate csv to MAV
#Read excel file from qupath into data frame
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

TumorSubtype<-codex_F_com_F@meta.data$TumorSubtype
TumorSubtype <- data.frame(TumorSubtype)
names(TumorSubtype) = 'TumorSubtype'

tumor_phenograph = codex_F_com_F@meta.data$tumor_phenograph
tumor_phenograph <- data.frame(tumor_phenograph)
names(tumor_phenograph) = 'tumor_phenograph'

gc_csd_final4<-bind_cols(converttomav, cellid, tumor_phenograph, TumorSubtype)
gc_csd_final4<-gc_csd_final4 %>%  mutate_all(as.character)
gc_csd_final4$Class  <- csd_raw_F_F$Class

#name the cell seg file. You need one cell seg file per region. For region 1 the cell seg file should be name reg001_
write.table(gc_csd_final4, file = "testing_allReg_seurat_._analysis_20220124_TumorFeatures_TumorSubtype.csv", sep=",", col.names = NA, row.names = TRUE)

save.image("~/Desktop/Project/Codex/Analysis_07v/03_DefineTumorSubtype/Outputs_v1/07_DefineTumorSubtype_V1.RData")



#### V2: subtype of tumor cells : modification on V1 after checking back on MAV ####
load("~/Desktop/Project/Codex/Analysis_07v/02_DefineCelltype/Outputs/codex_F_com_F.RData")
load("~/Desktop/Project/Codex/Analysis_07v/03_DefineTumorSubtype/Outputs_v1/codex_F_com_F_tumor_PhenoGraph_PCA_UMAP.RData")
load("~/Desktop/Project/Codex/Analysis_07v/03_DefineTumorSubtype/Outputs_v1/PhenoGraph_tumor_result_k25.RData")

codex_F_com_F_tumor_F <- subset(codex_F_com_F_tumor,idents=c(9,12),invert = TRUE) %>%
  Seurat::NormalizeData(verbose = FALSE,normalization.method = NULL) %>%
  ScaleData(verbose = FALSE)

codex_F_com_F_tumor_F <- RunPCA(codex_F_com_F_tumor_F, features=c("c-Myc","CD107A","CD44","E-Cadherin",
                                                              "Glypican3","Hepar1","HIF1a","HLA-DR","KI67",
                                                              "p-mTOR","p-S6","p53","Twist1","Vimentin"))
codex_F_com_F_tumor_F <- RunUMAP(codex_F_com_F_tumor_F, method = "umap-learn", metric = "correlation", dims = 1:length(codex_F_com_F_tumor_F@reductions$pca))

codex_F_com_F_tumor_F_phenograph <- subset(codex_F_com_F_tumor_F,features=c("c-Myc","CD107A","CD44","E-Cadherin",
                                                                        "Glypican3","Hepar1","HIF1a","HLA-DR","KI67",
                                                                        "p-mTOR","p-S6","p53","Twist1","Vimentin"))
codex_F_com_F_tumor_F_phenograph = ScaleData(codex_F_com_F_tumor_F_phenograph)
codex_F_com_F_tumor_F_phenograph <- codex_F_com_F_tumor_F_phenograph@assays$CODEX@scale.data
# PG_elbow(t(codex_F_com_F_tumor_F_phenograph), k_from = 10, k_to = 30, k_step=5) # cluster row
PhenoGraph_tumor_F_result_k25 <-as.numeric(membership(Rphenograph(t(codex_F_com_F_tumor_F_phenograph),k=25)))

save(PhenoGraph_tumor_F_result_k25, file = "PhenoGraph_tumor_F_result_k25.RData")
save(codex_F_com_F_tumor_F, file = "codex_F_com_F_tumor_F_PCA_UMAP_beforeIdents.RData")

Idents(codex_F_com_F_tumor_F) = PhenoGraph_tumor_F_result_k25
levels(codex_F_com_F_tumor_F) = as.character(1:43)

# Heatmap
codex_F_com_F_tumor_F_SelectedFeatures <- subset(codex_F_com_F_tumor_F,features=c("c-Myc","CD107A","CD44","E-Cadherin",
                                                                              "Glypican3","Hepar1","HIF1a","HLA-DR","KI67",
                                                                              "p-mTOR","p-S6","p53","Twist1","Vimentin"))
codex_F_com_F_tumor_F_SelectedFeatures = ScaleData(codex_F_com_F_tumor_F_SelectedFeatures)

p = DoHeatmap(subset(codex_F_com_F_tumor_F_SelectedFeatures, downsample=100), slot = "scale.data", assay = "CODEX",
              features = rownames(codex_F_com_F_tumor_F_SelectedFeatures), angle = 90, size = 4,disp.min = -2, disp.max = 2) + NoLegend()
pdf("20220130_SFig1_Subtype_Tumor_Heatmap_100_R.pdf", width = 12, height = 3)
p
dev.off()

p = DoHeatmap(subset(codex_F_com_F_tumor_F_SelectedFeatures, downsample=5000), slot = "scale.data", assay = "CODEX",
              features = rownames(codex_F_com_F_tumor_F_SelectedFeatures), angle = 90, size = 4,disp.min = -2, disp.max = 2) + NoLegend()
pdf("20220130_SFig1_Subtype_Tumor_Heatmap_5000_R.pdf", width = 12, height = 3)
p
dev.off()

# UMAP
codex_F_com_F_tumor_F@meta.data$tumor_phenograph = PhenoGraph_tumor_F_result_k25
p = DimPlot(subset(codex_F_com_F_tumor_F,downsample=500), 
            reduction = "umap", 
            #group.by = "tumor_phenograph",
            label = FALSE, 
            repel = T,
            split.by = "tumor_phenograph",
            ncol = 10,
            pt.size = 0.4) + NoLegend()

pdf("20220130_SFig1_Subtype_Tumor_Umap_Split_R.pdf", width = 12, height = 6)
p
dev.off()

# Cluster and Define subtype via pvclust and pheatmap
test_php_bind = data.frame()
for (i in 1:43) {
  temp = subset(codex_F_com_F_tumor_F,idents = i,features=c("c-Myc","CD107A","CD44","E-Cadherin",
                                                          "Glypican3","Hepar1","HIF1a","HLA-DR","KI67",
                                                          "p-mTOR","p-S6","p53","Twist1","Vimentin"))
  temp = as.data.frame(temp@assays$CODEX@data)
  test_php_bind = bind_rows(test_php_bind, rowMeans(temp))
}
test_php = as.matrix(test_php_bind)
rownames(test_php) = 1:43

p_dat = t(scale(test_php)) # scale column
#Then we conduct hierarchical cluster analysis with multiscale bootstrap with number of bootstrap 1000, using average method and correlation-based dissimilarity matrix as follows:
result <- pvclust(p_dat, method.dist="euclidean", method.hclust="average", nboot=10000, parallel=TRUE) # cluster col
plot(result)

saveRDS(test_php, file = "TrainingV2_tumorsub_test_php.Rds")

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
new.cluster.ids <- c("Tumor_others",#1                     
                     "Tumor_others",#2
                     "Tumor_others",#3
                     "Tumor_others",#4
                     "Tumor_Hepar1+_Glypican3+",#5   
                     "Tumor_CD44+_E-Cadherin+",#6
                     "Tumor_Hepar1_int_c-Myc_int",#7
                     "Tumor_c-Myc+",#8
                     "Tumor_others",#9
                     "Tumor_Glypican3+",#10
                     "Tumor_p53+",#11
                     "Tumor_Hepar1_int",#12
                     "Tumor_E-Cadherin+",#13
                     "Tumor_CD107a+",#14
                     "Tumor_others",#15
                     "Tumor_others",#16
                     "Tumor_p-S6+_c_Myc+",#17
                     "Tumor_others",#18
                     "Tumor_others",#19
                     "Tumor_Vimentin_int_HLA-DR_int",#20
                     "Tumor_p-mTOR+",#21
                     "Tumor_Glypican3+_Ki67+",#22
                     "Tumor_Glypican3+",#23
                     "Tumor_E-Cadherin+",#24
                     "Tumor_Hepar1_int_c-Myc_int",#25
                     "Tumor_CD107a+",#26
                     "Tumor_CD107a+",#27
                     "Tumor_HLA-DR+",#28
                     "Tumor_Twist1+",#29
                     "Tumor_HIF1a+_Vimentin+",#30
                     "Tumor_Ki67+",#31
                     "Tumor_p-S6_int",#32
                     "Tumor_p-S6+_Ki67_int",#33
                     "Tumor_others",#34
                     "Tumor_others",#35
                     "Tumor_Hepar1+",#36
                     "Tumor_CD44+",#37
                     "Tumor_others",#38
                     "Tumor_c-Myc+",#39
                     "Tumor_others",#40
                     "Tumor_Vimentin+",#41
                     "Tumor_c-Myc+",#42   
                     "Tumor_others"#43
)

names(new.cluster.ids) <- levels(codex_F_com_F_tumor_F)
codex_F_com_F_tumor_F <- RenameIdents(codex_F_com_F_tumor_F, new.cluster.ids)
codex_F_com_F_tumor_F@meta.data$TumorSubtype = Idents(codex_F_com_F_tumor_F)


# bind tumor subtypes with all stromal cells to further check and confirm via MAV
#load("~/OneDrive/01_Project/01_Codex/Bioinfo_Analysis/RealRun/Analysis_07v/02_DefineCelltype/Outputs/codex_F_com_F.RData") # Only scale the object after subset; Not redo PCA and UMAP!
codex_F_com_F_stromal <- subset(codex_F_com_F,
                                idents="Tumor cells",
                                invert = TRUE) %>%
  Seurat::NormalizeData(verbose = FALSE,normalization.method = NULL) %>%
  ScaleData(verbose = FALSE)

codex_F_com_F_stromal@meta.data$tumor_phenograph = 0
codex_F_com_F_stromal@meta.data$TumorSubtype = "Stromal"

codex_F_com_F@meta.data$cellid = as.numeric(rownames(codex_F_com_F@meta.data))
codex_F_com_F_tumor_F@meta.data$cellid = as.numeric(rownames(codex_F_com_F_tumor_F@meta.data))
codex_F_com_F_stromal@meta.data$cellid = as.numeric(rownames(codex_F_com_F_stromal@meta.data))

TumorSub_DeletedCellID = WhichCells(codex_F_com_F_tumor, 
                                idents = c(9,12))

codex_F_com_F_F = subset(codex_F_com_F,
                       cells = TumorSub_DeletedCellID,
                       invert = TRUE)  %>%
  ScaleData(verbose = FALSE)

temp = rbind(codex_F_com_F_tumor_F@meta.data,codex_F_com_F_stromal@meta.data)
temp = temp[order(temp$cellid),]

identical(codex_F_com_F_F@meta.data$TMA_reg,temp$TMA_reg)
identical(codex_F_com_F_F@meta.data$clusters,temp$clusters)

codex_F_com_F_F@meta.data$tumor_phenograph = temp$tumor_phenograph
codex_F_com_F_F@meta.data$TumorSubtype = temp$TumorSubtype


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

csd_F_F = csd_F[colnames(codex_F_com_F_F),]
csd_raw_F_F =  csd_raw_F[colnames(codex_F_com_F_F),]

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
cellid = codex_F_com_F_F@meta.data$cellid
cellid <- data.frame(cellid)
names(cellid) = 'cellid'

TumorSubtype<-codex_F_com_F_F@meta.data$TumorSubtype
TumorSubtype <- data.frame(TumorSubtype)
names(TumorSubtype) = 'TumorSubtype'

tumor_phenograph = codex_F_com_F_F@meta.data$tumor_phenograph
tumor_phenograph <- data.frame(tumor_phenograph)
names(tumor_phenograph) = 'tumor_phenograph'

gc_csd_final4<-bind_cols(converttomav, cellid, tumor_phenograph, TumorSubtype)
gc_csd_final4<-gc_csd_final4 %>%  mutate_all(as.character)
gc_csd_final4$Class  <- csd_raw_F_F$Class

#name the cell seg file. You need one cell seg file per region. For region 1 the cell seg file should be name reg001_
write.table(gc_csd_final4, file = "testing_allReg_seurat_._analysis_20220201_TumorFeatures_TumorSubtype_v2_20220308.csv", sep=",", col.names = NA, row.names = TRUE)

# TumorSubtype = gc_csd_final4
# save(TumorSubtype,file="20220315_TumorSubtype.RData")

#### Subtype: UMAP and Heatmap ####
## UMAP
subset_list = list()
downsample_cycle = round(table(Idents(codex_F_com_F_tumor_F))/500,0)
for (i in 1:23) {
  temp = subset(codex_F_com_F_tumor_F, idents = names(downsample_cycle)[i])
  subset_list[[i]] = subset(temp, downsample = downsample_cycle[i])
}
codex_F_com_F_subset_sample_prob = merge(x = subset_list[[1]], y = unlist(subset_list)[2:23], merge.data = T) # no scale data
codex_F_com_F_subset_sample_prob = subset(codex_F_com_F_tumor_F, cells = colnames(codex_F_com_F_subset_sample_prob))

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
        "#FBF88B", # 13
        "#B58794", # 14
        "#149684", # 15
        "#D92B00", # 16
        "#73904D", # 17
        "#406E30", # 18
        "#E3D78C", # 19
        "#0F0F0F", # 20
        # "#336699", # 21
        # "#669999", # 22
        # "#6600ff", # 23
        # "#330066", # 24
        "#cc99cc", # 25
        # "#666633", # 26
        "#ffcccc", # 27
        # "#cc6600", # 28 
        # "#5b7192", # 29
        "#99cc99" # 30
        )

orders = names(table(Idents(codex_F_com_F_tumor_F)))

pal<-colorRampPalette(cols)
image(x=1:23,y=1,z=as.matrix(1:23),col=pal(23))

p = DimPlot(codex_F_com_F_subset_sample_prob, 
            reduction = "umap", 
            group.by = "TumorSubtype",
            label = FALSE, 
            repel = T,
            cols = cols,
            order = rev(orders),
            pt.size = 0.4) + NoLegend()
  # theme(legend.position = "right") + 
  # guides(color = guide_legend(ncol = 2, byrow = T))

pdf("20220417_TumorSubtype_Umap_Merged_1in500.pdf", width = 5, height = 5)
p
dev.off()

## Freq
Freq = round((table(codex_F_com_F_tumor_F@meta.data$TumorSubtype)/nrow(codex_F_com_F_tumor_F@meta.data))*100,2)
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
                                      "#FBF88B", # 13
                                      "#B58794", # 14
                                      "#149684", # 15
                                      "#D92B00", # 16
                                      "#73904D", # 17
                                      "#406E30", # 18
                                      "#E3D78C", # 19
                                      "#0F0F0F", # 20
                                      # "#336699", # 21
                                      # "#669999", # 22
                                      # "#6600ff", # 23
                                      # "#330066", # 24
                                      "#cc99cc", # 25
                                      # "#666633", # 26
                                      "#ffcccc", # 27
                                      # "#cc6600", # 28 
                                      # "#5b7192", # 29
                                      "#99cc99" # 30
                                      )) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + NoLegend()+
  scale_y_break(c(15,25)) + 
  geom_text(aes(label = paste0(Freq,"%")), vjust = -.5) + 
  ylim(0,35)

pdf("20220417_heatmap_barplot_NoLegend.pdf", width = 15, height = 5)
p
dev.off()


## Heatmap
temp_subset = subset(codex_F_com_F_tumor_F, features=c("c-Myc","CD107A","CD44","E-Cadherin",
                                                          "Glypican3","Hepar1","HIF1a","HLA-DR","KI67",
                                                          "p-mTOR","p-S6","p53","Twist1","Vimentin"))
temp_subset = subset(temp_subset, downsample=2000)
temp_subset = ScaleData(temp_subset)

toplot_heatmap = temp_subset@assays$CODEX@scale.data

annotation_col = as.data.frame(temp_subset@meta.data)
annotation_col = subset(annotation_col, select = "TumorSubtype")
annotation_col = data.frame(cellid = rownames(annotation_col),
                            TumorSubtype = annotation_col$TumorSubtype)

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#6EC388", "white", "#C82033")) 


ha = HeatmapAnnotation(TumorSubtype = annotation_col$TumorSubtype,
                       col = list(TumorSubtype = c("Tumor_others" = "#FA339A",
                                                   "Tumor_Hepar1+_Glypican3+" = "#F6D026",
                                                   "Tumor_CD44+_E-Cadherin+" = "#74F882",
                                                   "Tumor_Hepar1_int_c-Myc_int" = "#7E8EE3",
                                                   "Tumor_c-Myc+" = "#8358CF",
                                                   "Tumor_Glypican3+" = "#27A7FA",
                                                   "Tumor_p53+" = "#EE6600",
                                                   "Tumor_Hepar1_int" = "#D083FF",
                                                   "Tumor_E-Cadherin+" = "#62FFFF",
                                                   "Tumor_CD107a+" = "#AE2519",
                                                   "Tumor_p-S6+_c_Myc+" = "#9C9B99",
                                                   "Tumor_Vimentin_int_HLA-DR_int" = "#162C9B",
                                                   "Tumor_p-mTOR+" = "#FBF88B",
                                                   "Tumor_Glypican3+_Ki67+" = "#B58794",
                                                   "Tumor_HLA-DR+" = "#149684",
                                                   "Tumor_Twist1+" = "#D92B00",
                                                   "Tumor_HIF1a+_Vimentin+" = "#73904D",
                                                   "Tumor_Ki67+" = "#406E30",
                                                   "Tumor_p-S6_int" = "#E3D78C",
                                                   "Tumor_p-S6+_Ki67_int" = "#0F0F0F",
                                                   "Tumor_Hepar1+" = "#cc99cc",
                                                   "Tumor_CD44+" = "#ffcccc",
                                                   "Tumor_Vimentin+" = "#99cc99")),
                       show_legend = F)

ht = Heatmap(toplot_heatmap,
             name = "Expression",
             col = col_fun,
             width = nrow(toplot_heatmap)*unit(12, "mm"), 
             height = nrow(toplot_heatmap)*unit(7, "mm"),
             column_split = annotation_col$TumorSubtype,
             column_title = NULL,
             cluster_rows = F,
             cluster_columns = F,
             show_column_names = F,
             top_annotation = ha,
             column_gap = unit(2,"mm"))

pdf("20220417_heatmap_2000.pdf", width = 13, height = 5)
ht
dev.off()
