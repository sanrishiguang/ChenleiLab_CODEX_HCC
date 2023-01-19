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

load("/Users/taozhou/Desktop/Project/Codex/Analysis_07v/01_GenerateSeuratInput/Outputs/csd.RData")
load("/Users/taozhou/Desktop/Project/Codex/Analysis_07v/01_GenerateSeuratInput/Outputs/csd_raw.RData")
load("/Users/taozhou/Desktop/Project/Codex/Analysis_07v/01_GenerateSeuratInput/Outputs/FilterID.RData")
load("/Users/taozhou/Desktop/Project/Codex/Analysis_07v/01_GenerateSeuratInput/Outputs/TrainingData.RData")

Panel <- c("actin","aSMA","c_Myc","Caspase_3","CD107A","CD11C","CD163","CD20","CD21","CD3",
           "CD31","CD4","CD44","CD45","CD45RO","CD68","CD8","DAPI","E_Cadherin","FOXP3",
           "Glypican3","Hepar1","HIF1a","HistoneH3","HLA_DR","Keratin","KI67","p_AMPK","p_mTOR","p_S6",
           "p53","Pan_CK","PD_L1","PD1","Podoplanin","Twist1","Vimentin")

CellMarkers_13 <- c("aSMA","CD11C","CD20","CD3",
                    "CD31","CD4","CD45","CD68","CD8",
                    "Glypican3","Hepar1","Podoplanin","Pan-CK")

#### Create the filtered original Seurat object ####
csd_F <- csd[FilterID,]
csd_raw_F <- csd_raw[FilterID,]

t_csd_F <- t(as.matrix(csd_F))
colnames(t_csd_F) <- rownames(csd_F)
rownames(t_csd_F) <- colnames(csd_F)

codex_F <- CreateSeuratObject(t_csd_F, project = "new", assay="CODEX") %>%
  Seurat::NormalizeData(verbose = FALSE,normalization.method = NULL) %>%
  ScaleData(verbose = FALSE)

codex_F <- AddMetaData(codex_F,csd_raw_F$TMA,col.name="TMA")

save(codex_F, file = "codex_F.RData")

#### Create the Seurat object removed batch effect with ComBat ####
codex_F_com = codex_F

dat = as.matrix(codex_F@assays$CODEX@counts)
dat = as.data.frame(t(dat))
dat$TMA = codex_F@meta.data$TMA

edata <- t(as.matrix(dat[,1:36]))
pheno <- subset(dat, select = "TMA")
combat_edata <- ComBat(dat = edata, batch = pheno$TMA)

dat_rbe = as.data.frame(t(combat_edata))
dat_rbe = as.data.frame(cbind(pheno$TMA,dat_rbe),stringsAsFactors = FALSE)
colnames(dat_rbe)[1]="TMA"

codex_F_com@assays$CODEX@counts = t(as.matrix(dat_rbe[,2:37]))
codex_F_com@assays$CODEX@data = t(as.matrix(dat_rbe[,2:37]))
codex_F_com = ScaleData(codex_F_com)

save(codex_F_com, file = "codex_F_com.RData")

#### RunPCA and RunUMAP to detect if the batch effect was removed ####
codex_F <- RunPCA(codex_F, features=CellMarkers_13)
codex_F <- RunUMAP(codex_F, method = "umap-learn", metric = "correlation", dims = 1:length(codex_F@reductions$pca))

p = DimPlot(codex_F, reduction = "umap", group.by = "TMA", order = rev(unique(codex_F@meta.data$TMA))) + NoLegend()
pdf("20220124_SFig1_ComBat_codex_F_merged.pdf", width = 4, height = 4)
p
dev.off()

p = DimPlot(codex_F, reduction = "umap", group.by = "TMA", split.by = "TMA", ncol = 5,
            order = rev(unique(codex_F@meta.data$TMA))) + NoLegend()
pdf("20220124_SFig1_ComBat_codex_F_split.pdf", width = 8, height = 4)
p
dev.off()

codex_F_com <- RunPCA(codex_F_com, features=CellMarkers_13)
codex_F_com <- RunUMAP(codex_F_com, method = "umap-learn", metric = "correlation", dims = 1:length(codex_F_com@reductions$pca))

p = DimPlot(codex_F_com, reduction = "umap", group.by = "TMA", order = rev(unique(codex_F_com@meta.data$TMA))) + NoLegend()
pdf("20220124_SFig1_ComBat_codex_F_com_merged.pdf", width = 4, height = 4)
p
dev.off()

p = DimPlot(codex_F_com, reduction = "umap", group.by = "TMA", split.by = "TMA", ncol = 5,
            order = rev(unique(codex_F_com@meta.data$TMA))) + NoLegend()
pdf("20220124_SFig1_ComBat_codex_F_com_split.pdf", width = 8, height = 4)
p
dev.off()

#### Cluster cells via PhenoGraph using the batch effect removed scale.data ####
codex_F_com_phenograph <- as.matrix(codex_F_com@assays$CODEX@scale.data)
# PG_elbow(t(codex_F_com_phenograph), k_from = 10, k_to = 30, k_step=5) # cluster row
PhenoGraph_codexFcom_scale_result_k25 <- as.numeric(membership(Rphenograph(t(codex_F_com_phenograph),k=25)))
save(PhenoGraph_codexFcom_scale_result_k25, file = "PhenoGraph_codexFcom_scale_result_k25.RData")

codex_F_com_phenograph.seurat = codex_F_com
Idents(codex_F_com_phenograph.seurat) = PhenoGraph_codexFcom_scale_result_k25
levels(codex_F_com_phenograph.seurat) = as.character(1:54) # to order the top annotation order of the following heatmap

p = DoHeatmap(subset(codex_F_com_phenograph.seurat, downsample=100), slot = "scale.data", assay = "CODEX",
              features = rownames(codex_F_com_phenograph.seurat), angle = 90, size = 4, disp.min = -2, disp.max = 2)# + NoLegend()
pdf("20220124_SFig1_Celltype_AllMarkers_heatmap_100.pdf", width = 12, height = 7)
p
dev.off()

p = DoHeatmap(subset(codex_F_com_phenograph.seurat, downsample=5000), slot = "scale.data", assay = "CODEX",
              features = rownames(codex_F_com_phenograph.seurat), angle = 90, size = 4, disp.min = -2, disp.max = 2) + NoLegend()
pdf("20220124_SFig1_Celltype_AllMarkers_heatmap_5000.pdf", width = 12, height = 5)
p
dev.off()

# Define cluster 44 as Unidentified
# Use other cells and 13 cell markers to further define cell types
codex_F_com_rmNeg_13cm.seurat = subset(codex_F_com_phenograph.seurat, 
                                     idents = c(44), 
                                     invert = TRUE,
                                     features = CellMarkers_13)
codex_F_com_rmNeg_13cm.seurat = ScaleData(codex_F_com_rmNeg_13cm.seurat)

codex_F_com_rmNeg_13cm.phenograph = codex_F_com_rmNeg_13cm.seurat@assays$CODEX@scale.data
# PG_elbow(t(codex_F_com_rmNeg_13cm.phenograph), k_from = 10, k_to = 30, k_step=5) # cluster row
PhenoGraph_13cm_result_k25 <-as.numeric(membership(Rphenograph(t(codex_F_com_rmNeg_13cm.phenograph),k=25)))
save(PhenoGraph_13cm_result_k25, file = "PhenoGraph_13cm_result_k25.RData")

Idents(codex_F_com_rmNeg_13cm.seurat) = PhenoGraph_13cm_result_k25
levels(codex_F_com_rmNeg_13cm.seurat) = as.character(1:55)

p = DoHeatmap(subset(codex_F_com_rmNeg_13cm.seurat, downsample=100), slot = "scale.data", assay = "CODEX",
              features = rownames(codex_F_com_rmNeg_13cm.seurat), angle = 90, size = 4,disp.min = -2, disp.max = 2)# + NoLegend()
pdf("20220124_SFig1_Celltype_CellMarkers_heatmap_100.pdf", width = 12, height = 7)
p
dev.off()

p = DoHeatmap(subset(codex_F_com_rmNeg_13cm.seurat, downsample=5000), slot = "scale.data", assay = "CODEX",
              features = rownames(codex_F_com_rmNeg_13cm.seurat), angle = 90, size = 4,disp.min = -2, disp.max = 2) + NoLegend()
pdf("20220124_SFig1_Celltype_CellMarkers_heatmap_5000.pdf", width = 12, height = 3)
p
dev.off()

new.cluster.ids <- c("Tumor cells",#1
                     "Tumor cells",#2 
                     "Tumor cells",#3
                     "Tumor cells",#4
                     "Tumor cells",#5
                     "Tumor cells",#6
                     "Tumor cells",#7
                     "Tumor cells",#8 
                     "Biliary tract cells",#9 
                     "Tumor cells",#10
                     "CD4+ T cells",#11
                     "Tumor cells",#12
                     "Tumor cells",#13
                     "Tumor cells",#14
                     "Tumor cells",#15
                     "Tumor cells",#16
                     "CD8+ T cells",#17
                     "Tumor cells",#18
                     "Tumor cells",#19
                     "Tumor cells",#20
                     "Unidentified",#21
                     "Tumor cells",#22
                     "Tumor cells",#23
                     "Tumor cells",#24
                     "Tumor cells",#25
                     "Tumor cells",#26
                     "Tumor cells",#27
                     "Tumor cells",#28
                     "B cells",#29
                     "Tumor cells",#30
                     "Fibroblasts",#31
                     "Lymphatic endothelial cells",#32
                     "Macrophages",#33
                     "Tumor cells",#34
                     "CD4+ T cells",#35
                     "Tumor cells",#36
                     "Tumor cells",#37
                     "Tumor cells",#38
                     "Tumor cells",#39
                     "Tumor cells",#40
                     "Tumor cells",#41
                     "Fibroblasts",#42
                     "Tumor cells",#43
                     "Macrophages",#44
                     "Endothelial cells",#45
                     "Macrophages",#46
                     "Endothelial cells",#47 
                     "Fibroblasts",#48 
                     "Tumor cells",#49 
                     "Tumor cells",#50
                     "Tumor cells",#51
                     "Tumor cells",#52
                     "Tumor cells",#53
                     "Tumor cells",#54
                     "Tumor cells"#55
)

names(new.cluster.ids) <- levels(codex_F_com_rmNeg_13cm.seurat)
codex_F_com_rmNeg_13cm.seurat <- RenameIdents(codex_F_com_rmNeg_13cm.seurat, new.cluster.ids)
codex_F_com_rmNeg_13cm.seurat@meta.data$cluster =  PhenoGraph_13cm_result_k25

# Add the results of cell clustering into metadata of codex_F_2
codex_F_com_keepNeg_13cm.seurat = subset(codex_F_com_phenograph.seurat, 
                                       idents = c(44),
                                       features = CellMarkers_13)

codex_F_com_keepNeg_13cm.seurat@meta.data$cluster = 0
codex_F_com_keepNeg_13cm.seurat@meta.data$cellid = as.numeric(rownames(codex_F_com_keepNeg_13cm.seurat@meta.data))
codex_F_com_rmNeg_13cm.seurat@meta.data$cellid = as.numeric(rownames(codex_F_com_rmNeg_13cm.seurat@meta.data))
codex_F_com_keepNeg_13cm.seurat@meta.data$celltype = "Unidentified"
codex_F_com_rmNeg_13cm.seurat@meta.data$celltype = Idents(codex_F_com_rmNeg_13cm.seurat)
identical(colnames(codex_F_com_keepNeg_13cm.seurat@meta.data),colnames(codex_F_com_rmNeg_13cm.seurat@meta.data))

codex_F_com_phenograph_celltype = rbind(codex_F_com_keepNeg_13cm.seurat@meta.data,
                                      codex_F_com_rmNeg_13cm.seurat@meta.data)
codex_F_com_phenograph_celltype = codex_F_com_phenograph_celltype[order(codex_F_com_phenograph_celltype$cellid),]
identical(as.numeric(rownames(codex_F_com@meta.data)),codex_F_com_phenograph_celltype$cellid)

codex_F_com@meta.data$celltype = codex_F_com_phenograph_celltype$celltype
codex_F_com@meta.data$clusters = codex_F_com_phenograph_celltype$cluster

identical(codex_F_com@meta.data$TMA,csd_raw_F$TMA)
codex_F_com@meta.data$reg = csd_raw_F$reg
Idents(codex_F_com) = codex_F_com@meta.data$celltype
codex_F_com@meta.data$TMA_reg = paste0(codex_F_com@meta.data$TMA,"_",codex_F_com@meta.data$reg)


#### Delete Unidentified cells and cells in  regions with poor quality ####
# Second round QC, pick passed cells to proceed following analysis
UnidentifiedCellID = WhichCells(codex_F_com, 
                       idents = "Unidentified")
BadRegCellID = rownames(codex_F_com@meta.data[codex_F_com@meta.data$TMA_reg %in% c("TMA9_3_reg011","TMA6_11_reg018","TMA7_10_reg006"),])
DeletedCellID = unique(c(UnidentifiedCellID,BadRegCellID))

codex_F_com_F = subset(codex_F_com,
                       cells = DeletedCellID,
                       invert = TRUE)

codex_F_com_F = ScaleData(codex_F_com_F)

save(codex_F_com_F,file = "codex_F_com_F.RData")

codex_F_com_F <- RunPCA(codex_F_com_F, features=CellMarkers_13)
codex_F_com_F <- RunUMAP(codex_F_com_F, method = "umap-learn", metric = "correlation", dims = 1:length(codex_F_com_F@reductions$pca))

#### Plot UMAP to show the defined cell types and Plot heatmap to show 13 markers expression of these cell types ####
cols = c("#8358CF", # B cells
         "#D083FF", # CD4+ T cells
         "#27A7FA", # CD8+ T cells
         "#62FFFF", # Macrophages
         "#74F882", # Fibroblasts
         "#F6D026", # Endothelial cells
         "#EE6600", # Lymphatic endothelial cells
         "#AE2519", # Biliary tract cells
         "#FA339A" # Tumor cells
)

orders = c("B cells",
           "CD4+ T cells",
           "CD8+ T cells",
           "Macrophages",
           "Fibroblasts",
           "Endothelial cells",
           "Lymphatic endothelial cells",
           "Biliary tract cells",
           "Tumor cells"
)

pal<-colorRampPalette(cols)
image(x=1:9,y=1,z=as.matrix(1:9),col=pal(9))

# Modify the UMAP of codex_F_com_F with defined cell types. Sample 1/500 for each cell type.
subset_list = list()
downsample_cycle = round(table(Idents(codex_F_com_F))/500,0)
for (i in 1:9) {
  temp = subset(codex_F_com_F, idents = names(downsample_cycle)[i])
  subset_list[[i]] = subset(temp, downsample = downsample_cycle[i])
}
codex_F_com_F_subset_sample_prob = merge(x = subset_list[[1]], y = unlist(subset_list)[2:9], merge.data = T) # no scale data
codex_F_com_F_subset_sample_prob = subset(codex_F_com_F, cells = colnames(codex_F_com_F_subset_sample_prob))

p = DimPlot(codex_F_com_F_subset_sample_prob, 
        reduction = "umap", 
        group.by = "celltype",
        label = FALSE, 
        repel = T,
        cols = cols,
        order = rev(orders),
        pt.size = 0.4) +
  theme(legend.position = "top") + 
  guides(color = guide_legend(ncol = 2, byrow = T))

pdf("20220124_Fig1_Celltype_Umap_Merged.pdf", width = 6, height = 8)
p
dev.off()

p = DimPlot(codex_F_com_F_subset_sample_prob, 
            reduction = "umap", 
            group.by = "celltype",
            label = FALSE, 
            repel = T,
            cols = cols,
            order = rev(orders),
            split.by = "celltype",
            ncol = 5,
            pt.size = 0.4) + NoLegend()

pdf("20220124_Fig1_Celltype_Umap_Split.pdf", width = 13, height = 8)
p
dev.off()

# Modify the HEATMAP of codex_F_com_F with defined cell types. Sample 5000 for each cell type.
temp_subset = subset(codex_F_com_F, features = CellMarkers_13)
temp_subset = subset(temp_subset, downsample=7000)
temp_subset = ScaleData(temp_subset)

toplot_heatmap = temp_subset@assays$CODEX@scale.data

annotation_col = as.data.frame(temp_subset@meta.data)
annotation_col = subset(annotation_col, select = "celltype")
annotation_col = data.frame(cellid = rownames(annotation_col),
                            celltype = annotation_col$celltype)

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#6EC388", "white", "#C82033")) 


ha = HeatmapAnnotation(celltype = annotation_col$celltype,
                       col = list(celltype = c("B cells" = "#8358CF",
                                               "CD4+ T cells" = "#D083FF",
                                               "CD8+ T cells" = "#27A7FA",
                                               "Macrophages" = "#62FFFF",
                                               "Fibroblasts" = "#74F882",
                                               "Endothelial cells" = "#F6D026",
                                               "Lymphatic endothelial cells" = "#EE6600",
                                               "Biliary tract cells" = "#AE2519",
                                               "Tumor cells" = "#FA339A")),
                       show_legend = F)

ht = Heatmap(toplot_heatmap,
             name = "Expression",
             col = col_fun,
             width = nrow(toplot_heatmap)*unit(12, "mm"), 
             height = nrow(toplot_heatmap)*unit(7, "mm"),
             column_split = annotation_col$celltype,
             column_title = NULL,
             cluster_rows = F,
             cluster_columns = F,
             show_column_names = F,
             top_annotation = ha,
             column_gap = unit(2,"mm"))

pdf("20220124_Fig1_heatmap_7000.pdf", width = 10, height = 5)
ht
dev.off()

# FeaturePlot to detect the expression of all features
codex_F_com_F_subset_sample_prob_scale = ScaleData(codex_F_com_F_subset_sample_prob)
p = FeaturePlot(codex_F_com_F_subset_sample_prob_scale,
                slot = "scale.data",
                features = rownames(codex_F_com),
                ncol = 6,
                min.cutoff = 0,
                max.cutoff = 2,
                cols = c("#AFDFD5", "#D96B0C"))

pdf("20220124_Celltype_Umap_FeaturesExpres.pdf", width = 20, height = 17)
p
dev.off()

## Percentage of each celltype as annotation of HEATMAP
Freq = round((table(codex_F_com_F@meta.data$celltype)/nrow(codex_F_com_F@meta.data))*100,2)
Freq = as.data.frame(Freq)
colnames(Freq)[1] = "Celltype"

p = ggplot(data = Freq, aes(x = Celltype, y = Freq, fill = Celltype)) + 
  geom_col(width = 0.8, position = position_dodge(0.65)) +
  scale_fill_manual(values = c("#8358CF", "#AE2519", "#D083FF", "#27A7FA", "#F6D026", "#74F882", "#EE6600", "#62FFFF", "#FA339A")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_break(c(10,70)) + 
  geom_text(aes(label = paste0(Freq,"%")), vjust = -.5) + 
  ylim(0,75)

pdf("20220124_Fig1_heatmap_barplot.pdf", width = 10, height = 3)
p
dev.off()

















































