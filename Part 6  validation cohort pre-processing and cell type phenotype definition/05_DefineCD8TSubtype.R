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

#### V1: subtype: CD8T ####
load("~/Desktop/Project/Codex/Analysis_07v_Validation/02_DefineCelltype/Outputs/codex_F_com_F.RData")
codex_F_com_F_CD8T <- subset(codex_F_com_F,idents="CD8+ T cells") %>%
  Seurat::NormalizeData(verbose = FALSE,normalization.method = NULL) %>%
  ScaleData(verbose = FALSE)

codex_F_com_F_CD8T <- RunPCA(codex_F_com_F_CD8T, features=c("CD107A","HLA-DR","CD44","CD45RO","KI67","p-S6","HIF1a","c-Myc","PD1"))
codex_F_com_F_CD8T <- RunUMAP(codex_F_com_F_CD8T, method = "umap-learn", metric = "correlation", dims = 1:length(codex_F_com_F_CD8T@reductions$pca))

codex_F_com_F_CD8T_phenograph <- subset(codex_F_com_F_CD8T,features=c("CD107A","HLA-DR","CD44","CD45RO","KI67","p-S6","HIF1a","c-Myc","PD1"))
codex_F_com_F_CD8T_phenograph = ScaleData(codex_F_com_F_CD8T_phenograph)
codex_F_com_F_CD8T_phenograph <- codex_F_com_F_CD8T_phenograph@assays$CODEX@scale.data
# PG_elbow(t(codex_F_com_F_CD8T_phenograph), k_from = 10, k_to = 30, k_step=5) # cluster row
PhenoGraph_CD8T_result_k25 <-as.numeric(membership(Rphenograph(t(codex_F_com_F_CD8T_phenograph),k=25)))

save(PhenoGraph_CD8T_result_k25, file = "PhenoGraph_CD8T_result_k25.RData")

Idents(codex_F_com_F_CD8T) = PhenoGraph_CD8T_result_k25
levels(codex_F_com_F_CD8T) = as.character(1:25)

save(codex_F_com_F_CD8T,file = "codex_F_com_F_CD8T_PhenoGraph_PCA_UMAP.RData")

# Heatmap
codex_F_com_F_CD8T_SelectedFeatures <- subset(codex_F_com_F_CD8T,
                                                     features=c("CD107A","HLA-DR","CD44","CD45RO","KI67","p-S6","HIF1a","c-Myc","PD1"))
codex_F_com_F_CD8T_SelectedFeatures = ScaleData(codex_F_com_F_CD8T_SelectedFeatures)

p = DoHeatmap(subset(codex_F_com_F_CD8T_SelectedFeatures, downsample=100), slot = "scale.data", assay = "CODEX",
              features = rownames(codex_F_com_F_CD8T_SelectedFeatures), angle = 90, size = 4,disp.min = -2, disp.max = 2) + NoLegend()
pdf("20220411_Subtype_CD8T_Heatmap_100.pdf", width = 12, height = 3)
p
dev.off()

p = DoHeatmap(subset(codex_F_com_F_CD8T_SelectedFeatures, downsample=5000), slot = "scale.data", assay = "CODEX",
              features = rownames(codex_F_com_F_CD8T_SelectedFeatures), angle = 90, size = 4,disp.min = -2, disp.max = 2) + NoLegend()
pdf("20220411_Subtype_CD8T_Heatmap_5000.pdf", width = 12, height = 3)
p
dev.off()

# UMAP
codex_F_com_F_CD8T@meta.data$CD8T_phenograph = PhenoGraph_CD8T_result_k25
p = DimPlot(subset(codex_F_com_F_CD8T,downsample=500), 
            reduction = "umap", 
            label = FALSE, 
            repel = T,
            split.by = "CD8T_phenograph",
            ncol = 9,
            pt.size = 0.4) + NoLegend()

pdf("20220411_Subtype_CD8T_Umap_Split.pdf", width = 12, height = 6)
p
dev.off()

# Cluster and Define subtype via pvclust and pheatmap
test_php_bind = data.frame()
for (i in 1:25) {
  temp = subset(codex_F_com_F_CD8T,idents = i,features=c("CD107A","HLA-DR","CD44","CD45RO","KI67","p-S6","HIF1a","c-Myc","PD1"))
  temp = as.data.frame(temp@assays$CODEX@data)
  test_php_bind = bind_rows(test_php_bind, rowMeans(temp))
}
test_php = as.matrix(test_php_bind)
rownames(test_php) = 1:25

test_php_validation = test_php
test_php_training = readRDS("~/Desktop/Project/Codex/Analysis_07v/06_DefineCD8TSubtype/Outputs_v2/TrainingV2_CD8Tsub_test_php.Rds")
identical(colnames(test_php_training),colnames(test_php_validation))
rownames(test_php_training) = paste0(rownames(test_php_training),"_T")
test_php = rbind(test_php_training,test_php_validation)

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
         angle_col = "90",
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
         angle_col = "90",
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
         angle_col = "90",
         name = "Expression",
         color = c(colorRampPalette(colors = c("blue","white","red"))(100)),
         legend_breaks = seq(-2,2,1),
         breaks = bk)

# Definition
new.cluster.ids <- c("CD8T_others",#1
                     "CD8T_HLA-DR+",#2 
                     "CD8T_p-S6+",#3
                     "CD8T_KI67+",#4
                     "CD8T_CD45RO+",#5
                     "CD8T_HIF1a+",#6
                     "CD8T_HLA-DR+",#7
                     "CD8T_p-S6+",#8 
                     "CD8T_p-S6+",#9 
                     "CD8T_CD45RO+",#10
                     "CD8T_HIF1a+",#11
                     "CD8T_p-S6+_CD45RO+",#12
                     "CD8T_CD45RO+",#13 
                     "CD8T_others",#14
                     "CD8T_c-Myc+",#15
                     "CD8T_p-S6+_CD44+",#16
                     "CD8T_CD44+",#17
                     "CD8T_others",#18
                     "CD8T_c-Myc+",#19
                     "CD8T_CD107a+",#20
                     "CD8T_HIF1a+",#21 
                     "CD8T_CD44+",#22 #
                     "CD8T_CD107a+",#23
                     "CD8T_CD107a+",#24
                     "CD8T_others"#25
)
names(new.cluster.ids) <- levels(codex_F_com_F_CD8T)
codex_F_com_F_CD8T <- RenameIdents(codex_F_com_F_CD8T, new.cluster.ids)
codex_F_com_F_CD8T@meta.data$CD8TSubtype = Idents(codex_F_com_F_CD8T)

# bind macrophages subtypes with other cells to further check and confirm via MAV
codex_F_com_F_others <- subset(codex_F_com_F,
                              idents="CD8+ T cells",
                              invert = TRUE) %>%
  Seurat::NormalizeData(verbose = FALSE,normalization.method = NULL) %>%
  ScaleData(verbose = FALSE)

codex_F_com_F_others@meta.data$CD8T_phenograph = 0
codex_F_com_F_others@meta.data$CD8TSubtype = "cell_others"

codex_F_com_F@meta.data$cellid = as.numeric(rownames(codex_F_com_F@meta.data))
codex_F_com_F_CD8T@meta.data$cellid = as.numeric(rownames(codex_F_com_F_CD8T@meta.data))
codex_F_com_F_others@meta.data$cellid = as.numeric(rownames(codex_F_com_F_others@meta.data))

temp = rbind(codex_F_com_F_CD8T@meta.data,codex_F_com_F_others@meta.data)
temp = temp[order(temp$cellid),]

identical(codex_F_com_F@meta.data$TMA_reg,temp$TMA_reg)
identical(codex_F_com_F@meta.data$clusters,temp$clusters)

codex_F_com_F@meta.data$CD8T_phenograph = temp$CD8T_phenograph
codex_F_com_F@meta.data$CD8TSubtype = temp$CD8TSubtype

#### generate csv to MAV
#Read excel file from qupath into data frame
csd_raw <- readRDS("~/Desktop/Project/Codex/Analysis_07v_Validation/01_GenerateSeuratInput/Outputs/csd_raw.Rds")
csd <- readRDS("~/Desktop/Project/Codex/Analysis_07v_Validation/01_GenerateSeuratInput/Outputs/csd.Rds")
FilterID <- readRDS("~/Desktop/Project/Codex/Analysis_07v_Validation/01_GenerateSeuratInput/Outputs/FilterID.Rds")

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

CD8TSubtype<-codex_F_com_F@meta.data$CD8TSubtype
CD8TSubtype <- data.frame(CD8TSubtype)
names(CD8TSubtype) = 'CD8TSubtype'

CD8T_phenograph = codex_F_com_F@meta.data$CD8T_phenograph
CD8T_phenograph <- data.frame(CD8T_phenograph)
names(CD8T_phenograph) = 'CD8T_phenograph'

gc_csd_final4<-bind_cols(converttomav, cellid, CD8T_phenograph, CD8TSubtype)
gc_csd_final4<-gc_csd_final4 %>%  mutate_all(as.character)
gc_csd_final4$Class  <- csd_raw_F_F$Class

# #name the cell seg file. You need one cell seg file per region. For region 1 the cell seg file should be name reg001_
# write.table(gc_csd_final4, file = "testing_allReg_seurat_._analysis_20220310_CD8TFeatures_CD8TSubtype_V1.csv", sep=",", col.names = NA, row.names = TRUE)

CD8TSubtype = gc_csd_final4
save(CD8TSubtype,file = "20220411_CD8TSubtype.RData")

# #### Subtype: UMAP and Heatmap ####
## UMAP
subset_list = list()
downsample_cycle = round(table(Idents(codex_F_com_F_CD8T))/20,0)
for (i in 1:11) {
  temp = subset(codex_F_com_F_CD8T, idents = names(downsample_cycle)[i])
  subset_list[[i]] = subset(temp, downsample = downsample_cycle[i])
}
codex_F_com_F_subset_sample_prob = merge(x = subset_list[[1]], y = unlist(subset_list)[2:11], merge.data = T) # no scale data
codex_F_com_F_subset_sample_prob = subset(codex_F_com_F_CD8T, cells = colnames(codex_F_com_F_subset_sample_prob))

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
         "#9C9B99" # 11
         #"#162C9B", # 12
         #"#FBF88B", # 13
         #"#B58794", # 14
         #"#149684" # 15
         #"#D92B00" # 16
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

orders = names(table(Idents(codex_F_com_F_CD8T)))

pal<-colorRampPalette(cols)
image(x=1:11,y=1,z=as.matrix(1:11),col=pal(11))

p = DimPlot(codex_F_com_F_subset_sample_prob, 
            reduction = "umap", 
            group.by = "CD8TSubtype",
            label = FALSE, 
            repel = T,
            cols = cols,
            order = rev(orders),
            pt.size = 0.4) + NoLegend()
# theme(legend.position = "right") + 
# guides(color = guide_legend(ncol = 2, byrow = T))

pdf("20220422_CD8TSubtype_Umap_Merged_1in20.pdf", width = 5, height = 5)
p
dev.off()

## Freq
Freq = round((table(Idents(codex_F_com_F_CD8T))/nrow(codex_F_com_F_CD8T@meta.data))*100,2)
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
                               "#9C9B99" # 11
                               #"#162C9B", # 12
                               #"#FBF88B", # 13
                               #"#B58794", # 14
                               #"#149684", # 15
                               #"#D92B00" # 16
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
        axis.ticks.x=element_blank()) + NoLegend() +
  #scale_y_break(c(15,20)) + 
  geom_text(aes(label = paste0(Freq,"%")), vjust = -.5) + 
  ylim(0,20)

pdf("20220422_heatmap_barplot_NoLegend.pdf", width = 8, height = 5)
p
dev.off()

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
                               "#9C9B99" # 11
                               #"#162C9B", # 12
                               #"#FBF88B", # 13
                               #"#B58794", # 14
                               #"#149684", # 15
                               #"#D92B00" # 16
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
        axis.ticks.x=element_blank()) + 
  #scale_y_break(c(15,20)) + 
  geom_text(aes(label = paste0(Freq,"%")), vjust = -.5) + 
  ylim(0,20)

pdf("20220422_heatmap_barplot.pdf", width = 10, height = 5)
p
dev.off()

## Heatmap
temp_subset = subset(codex_F_com_F_CD8T, features=c("CD107A","HLA-DR","CD44","CD45RO","KI67","p-S6","HIF1a","c-Myc","PD1"))
temp_subset = subset(temp_subset, downsample=2000)
temp_subset = ScaleData(temp_subset)

toplot_heatmap = temp_subset@assays$CODEX@scale.data

annotation_col = as.data.frame(temp_subset@meta.data)
annotation_col = subset(annotation_col, select = "CD8TSubtype")
annotation_col = data.frame(cellid = rownames(annotation_col),
                            CD8TSubtype = annotation_col$CD8TSubtype)

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#6EC388", "white", "#C82033")) 


ha = HeatmapAnnotation(CD8TSubtype = annotation_col$CD8TSubtype,
                       col = list(CD8TSubtype = c("CD8T_others" = "#FA339A",
                                                  "CD8T_HLA-DR+" = "#F6D026",
                                                  "CD8T_p-S6+" = "#74F882",
                                                  "CD8T_KI67+" = "#7E8EE3",
                                                  "CD8T_CD45RO+" = "#8358CF",
                                                  "CD8T_HIF1a+" = "#27A7FA",
                                                  "CD8T_p-S6+_CD45RO+" = "#EE6600",
                                                  "CD8T_c-Myc+" = "#D083FF",
                                                  "CD8T_p-S6+_CD44+" = "#62FFFF",
                                                  "CD8T_CD44+" = "#AE2519",
                                                  "CD8T_CD107a+" = "#9C9B99")),
                       show_legend = F)

ht = Heatmap(toplot_heatmap,
             name = "Expression",
             col = col_fun,
             width = nrow(toplot_heatmap)*unit(12, "mm"), 
             height = nrow(toplot_heatmap)*unit(7, "mm"),
             column_split = annotation_col$CD8TSubtype,
             column_title = NULL,
             cluster_rows = F,
             cluster_columns = F,
             show_column_names = F,
             top_annotation = ha,
             column_gap = unit(2,"mm"))

pdf("20220422_heatmap_2000.pdf", width = 10, height = 5)
ht
dev.off()








