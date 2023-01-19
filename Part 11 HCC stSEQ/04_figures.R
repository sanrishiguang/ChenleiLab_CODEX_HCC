##R4.1.0

# Packages and signatures -------------------------------------------------


library(Seurat)  #‘4.1.1’
library(SeuratData) #‘0.2.2’
library(ggplot2)   #‘3.3.6’
library(patchwork) #‘1.1.1’
library(dplyr)   #‘1.0.9’
library(SPATA2)   #‘0.1.0’
library(ggpubr)  #‘0.4.0’


# Gene Signature 
Mac_VIM_high <- readRDS(file='./MyGeneSignature/Mac_VIM_hi_genesets_R_FCtop50.RDS') 

T.reg = c("FOXP3", "TNFRSF18", "TNFRSF9", "TIGIT", "IKZF2", "CTLA4",
          "TNFRSF4", "CCR8", "IL2RA", "BATF", "IL2RB", "ACP5", "CD27",
          "ICOS", "CTSC", "PTTG1", "VDR", "IL1R2", "CD7", "ERI1", "CSF1",
          "TFRC", "CSF2RB", "GLRX", "NCF4", "DPYSL2", "MAGEH1", "CALM3",
          "FCRL3", "PARK7", "HTATIP2", "CRADD", "IL32", "RGS1", "STAT3", "PTPN22")  ##zhang zm



# 225847-T (HCC-1-N):  cluster4----------------------------------------------------------------
load('./20220823 整理Figure/225847_T.Rdata')

p1 <- plotSurfaceAverage(object = spata_obj, 
                         color_by = Mac_VIM_high, 
                         smooth = TRUE, 
                         pt_size = 2,
                         smooth_span = 0.2, 
                         pt_clrsp = "inferno"
) +  labs(title='Mac_VIM_high')

p2 <- plotSurfaceAverage(object = spata_obj, 
                         color_by = T.reg, 
                         smooth = TRUE, 
                         pt_size = 2,
                         smooth_span = 0.2, 
                         pt_clrsp = "inferno")+  labs(title='Treg')
p3 <- plotSurface(object = spata_obj, 
                  color_by = "seurat_clusters", 
                  display_image = F, 
                  pt_alpha = 1) # nolint

p4 <- plotDistributionAcross2(df = joined_df,
                              variables = mySig_gs,
                              across = "seurat_clusters",
                              plot_type = "violin",
                              clrp = 'jco',
                              nrow=2)

temp <- joined_df %>% dplyr::filter(seurat_clusters == 4) %>% select(mySig_Mac_VIM_high, mySig_T.reg)
p5 <- ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
                color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = T, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
                title =  'seurat-cluster 4')

gseaDF <- getGseaDf(
  object = spata_obj, 
  across = "seurat_clusters",
  across_subset = '4',
  method_de = "wilcox", 
  signif_threshold = 0.05 # extract top 20 most significant gene sets
) 
gseaDF_HM <- gseaDF %>% dplyr::filter(grepl("^HM|BP|RCTM", label))
gseaDF_HM <- gseaDF_HM[order(gseaDF_HM$fdr, decreasing = T),]
gseaDF_HM$label <- factor(gseaDF_HM$label, levels = gseaDF_HM$label)
gseaDF_HM <- tail(gseaDF_HM, 7)
p6 <- ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')

pdf(file = './20220823 整理Figure/225847_T.pdf', height = 5, width = 30)
p1|p2 |p3 |p4 |p5|p6
dev.off()

plotViolinplot(spata_obj, 'IL1B',
               across = 'seurat_clusters',
               test_groupwise = 'kruskal.test'
)

####IL1B
df <- getExpressionMatrix(spata_obj) %>% t() %>% as.data.frame()
df_int <- df %>% select('IL1B')
df_int$barcodes <- rownames(df_int)
joined_df <- merge(joined_df, df_int, by='barcodes')


temp <- joined_df %>% dplyr::filter(seurat_clusters == 4) %>% select(mySig_Mac_VIM_high, IL1B)
ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'IL1B',
          fill ='#8491B4FF',color='darkgrey',size = 2,shape = 21,  # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#7E6148FF", fill = "#B09C85FF"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = T, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
          title =  'seurat-cluster 4')


# 227312_T: cluster3----------------------------------------------------------------
rm(list=ls())
gc()

load('./20220823 整理Figure/227312_T.Rdata')

p1 <- plotSurfaceAverage(object = spata_obj, 
                         color_by = Mac_VIM_high, 
                         smooth = TRUE, 
                         pt_size = 2,
                         smooth_span = 0.2, 
                         pt_clrsp = "inferno"
) +  labs(title='Mac_VIM_high')

p2 <- plotSurfaceAverage(object = spata_obj, 
                         color_by = T.reg, 
                         smooth = TRUE, 
                         pt_size = 2,
                         smooth_span = 0.2, 
                         pt_clrsp = "inferno")+  labs(title='Treg')
p3 <- plotSurface(object = spata_obj, 
                  color_by = "seurat_clusters", 
                  display_image = F, 
                  pt_alpha = 1) # nolint

p4 <- plotDistributionAcross2(df = joined_df,
                              variables = mySig_gs,
                              across = "seurat_clusters",
                              plot_type = "violin",
                              clrp = 'jco',
                              nrow=2)

temp <- joined_df %>% dplyr::filter(seurat_clusters == 3) %>% select(mySig_Mac_VIM_high, mySig_T.reg)
p5 <- ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
                color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = T, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
                title =  'seurat-cluster 3')

gseaDF <- getGseaDf(
  object = spata_obj, 
  across = "seurat_clusters",
  across_subset = '3',
  method_de = "wilcox", 
  signif_threshold = 0.05 # extract top 20 most significant gene sets
) 
gseaDF_HM <- gseaDF %>% dplyr::filter(grepl("^HM|BP", gseaDF$label))
gseaDF_HM <- gseaDF_HM[order(gseaDF_HM$fdr, decreasing = T),]
gseaDF_HM$label <- factor(gseaDF_HM$label, levels = gseaDF_HM$label)
gseaDF_HM <- tail(gseaDF_HM, 7)
p6 <- ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')

pdf(file = './20220823 整理Figure/227312_T.pdf', height = 5, width = 30)
p1|p2 |p3 |p4 |p5|p6
dev.off()

plotViolinplot(spata_obj, 'IL1B',
               across = 'seurat_clusters',
               test_groupwise = 'kruskal.test') 

plotViolinplot(spata_obj, 'IL1B',
               across = 'seurat_clusters',
               #test_groupwise = 'kruskal.test',
               test_pairwise = 'wilcox.test',
               ref_group = '3',
               step_increase =0.10) 

####IL1B
df <- getExpressionMatrix(spata_obj) %>% t() %>% as.data.frame()
df_int <- df %>% select('IL1B')
df_int$barcodes <- rownames(df_int)
joined_df <- merge(joined_df, df_int, by='barcodes')


temp <- joined_df %>% dplyr::filter(seurat_clusters == 3) %>% select(mySig_Mac_VIM_high, IL1B)
ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'IL1B',
          fill ='#8491B4FF',color='darkgrey',size = 2,shape = 21, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#7E6148FF", fill = "#B09C85FF"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = T, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
          title =  'seurat-cluster 4')





# 084717_A_orig (HCC-5-A) cluster 0&2&6 ----------------------------------------------------------------



rm(list = ls())
gc()

load(file='./20220823 整理Figure/084717_A.Rdata')




p1 <- plotSurfaceAverage(object = spata_obj, 
                         color_by = Mac_VIM_high, 
                         smooth = TRUE, 
                         pt_size = 2,
                         smooth_span = 0.2, 
                         pt_clrsp = "inferno")
p2 <- plotSurfaceAverage(object = spata_obj, 
                         color_by = T.reg, 
                         smooth = TRUE, 
                         pt_size = 2,
                         smooth_span = 0.2, 
                         pt_clrsp = "inferno")

p3 <- plotSurface(object = spata_obj, color_by = "seurat_clusters", display_image = F, pt_alpha = 1) # nolint


p4 <- plotDistributionAcross2(df = joined_df,
                              variables = mySig_gs,
                              across = "seurat_clusters",
                              plot_type = "violin",
                              clrp = 'jco',
                              nrow=2)


########correlation between Mac-vim-hi and Treg score


temp <- joined_df %>% dplyr::filter(seurat_clusters == 0) %>% select(mySig_Mac_VIM_high, mySig_T.reg)

p5.0 <- ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
                  color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = T, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
                  title =  'seurat-cluster 0')




temp <- joined_df %>% dplyr::filter(seurat_clusters == 2) %>% select(mySig_Mac_VIM_high, mySig_T.reg)

p5.2 <- ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
                  color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = T, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
                  title =  'seurat-cluster 2')

temp <- joined_df %>% dplyr::filter(seurat_clusters == 6) %>% select(mySig_Mac_VIM_high, mySig_T.reg)

p5.6 <- ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
                  color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = T, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
                  title =  'seurat-cluster 6')



###cluster0
gseaDF <- getGseaDf(
  object = spata_obj, 
  across = "seurat_clusters",
  across_subset = '0',
  method_de = "wilcox", 
  signif_threshold = 0.05 # extract top 20 most significant gene sets
) 
gseaDF_HM <- gseaDF %>% dplyr::filter(grepl("^HM|BP", gseaDF$label))
gseaDF_HM <- gseaDF_HM[order(gseaDF_HM$fdr, decreasing = T),]
gseaDF_HM$label <- factor(gseaDF_HM$label, levels = gseaDF_HM$label)
gseaDF_HM <- tail(gseaDF_HM, 7)

p6.0 <- ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')



###cluster2
gseaDF <- getGseaDf(
  object = spata_obj, 
  across = "seurat_clusters",
  across_subset = '2',
  method_de = "wilcox", 
  signif_threshold = 0.05 # extract top 20 most significant gene sets
) 
gseaDF_HM <- gseaDF %>% dplyr::filter(grepl("^HM|BP", gseaDF$label))
gseaDF_HM <- gseaDF_HM[order(gseaDF_HM$fdr, decreasing = T),]
gseaDF_HM$label <- factor(gseaDF_HM$label, levels = gseaDF_HM$label)
gseaDF_HM <- tail(gseaDF_HM, 7)

p6.2 <- ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')


###cluster6
gseaDF <- getGseaDf(
  object = spata_obj, 
  across = "seurat_clusters",
  across_subset = '6',
  method_de = "wilcox", 
  signif_threshold = 0.05 # extract top 20 most significant gene sets
) 
gseaDF_HM <- dplyr::filter(gseaDF, grepl("^HM|BP", gseaDF$label))
gseaDF_HM <- gseaDF_HM[-1,]
gseaDF_HM <- gseaDF_HM[order(gseaDF_HM$fdr, decreasing = T),]
gseaDF_HM$label <- factor(gseaDF_HM$label, levels = gseaDF_HM$label)
gseaDF_HM <- tail(gseaDF_HM, 7)

p6.6 <- ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')


pdf(file = './20220823 整理Figure/084717_A.pdf', height = 5, width = 20)
p1|p2 |p3 |p4
dev.off()

pdf(file = './20220823 整理Figure/084717_B.pdf', height = 5, width = 15)
p5.0|p5.2|p5.6
dev.off()

pdf(file = './20220823 整理Figure/084717_C.pdf', height = 4, width = 26)
p6.0|p6.2|p6.6
dev.off()


plotViolinplot(spata_obj, 'IL1B',
               across = 'seurat_clusters',
               test_groupwise = 'kruskal.test') 


####IL1B
df <- getExpressionMatrix(spata_obj) %>% t() %>% as.data.frame()
df_int <- df %>% select('IL1B')
df_int$barcodes <- rownames(df_int)
joined_df <- merge(joined_df, df_int, by='barcodes')


temp <- joined_df %>% dplyr::filter(seurat_clusters == 0) %>% select(mySig_Mac_VIM_high, IL1B)
ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'IL1B',
          fill ='#8491B4FF',color='darkgrey',size = 2,shape = 21, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#7E6148FF", fill = "#B09C85FF"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = T, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
          title =  'seurat-cluster 0')


# YangFu HCC2 (HCC-R1) cluster 5 & 7 ----------------------------------------------------------------


rm(list = ls())
gc()

load(file='./20220823 整理Figure/YangFU-HCC2.Rdata')



p1 <- plotSurfaceAverage(object = spata_obj, 
                         color_by = Mac_VIM_high, 
                         smooth = TRUE, 
                         pt_size = 2,
                         smooth_span = 0.2, 
                         pt_clrsp = "inferno")
p2 <- plotSurfaceAverage(object = spata_obj, 
                         color_by = T.reg, 
                         smooth = TRUE, 
                         pt_size = 2,
                         smooth_span = 0.2, 
                         pt_clrsp = "inferno")

p3 <- plotSurface(object = spata_obj, color_by = "seurat_clusters", display_image = F, pt_alpha = 1) # nolint


p4 <- plotDistributionAcross2(df = joined_df,
                              variables = mySig_gs,
                              across = "seurat_clusters",
                              plot_type = "violin",
                              clrp = 'jco',
                              nrow=2)


########correlation between Mac-vim-hi and Treg score


temp <- joined_df %>% dplyr::filter(seurat_clusters == 5) %>% select(mySig_Mac_VIM_high, mySig_T.reg)

p5.5 <- ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
                  color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = T, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
                  title =  'seurat-cluster 5')




temp <- joined_df %>% dplyr::filter(seurat_clusters == 7) %>% select(mySig_Mac_VIM_high, mySig_T.reg)

p5.7 <- ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
                  color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = T, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
                  title =  'seurat-cluster 7')





###cluster5
gseaDF <- getGseaDf(
  object = spata_obj, 
  across = "seurat_clusters",
  across_subset = '5',
  method_de = "wilcox", 
  signif_threshold = 0.05 # extract top 20 most significant gene sets
) 
gseaDF_HM <- gseaDF %>% dplyr::filter(grepl("^HM|BP", gseaDF$label))
gseaDF_HM <- gseaDF_HM[order(gseaDF_HM$fdr, decreasing = T),]
gseaDF_HM$label <- factor(gseaDF_HM$label, levels = gseaDF_HM$label)
gseaDF_HM <- tail(gseaDF_HM, 7)

p6.5 <- ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')



###cluster7
gseaDF <- getGseaDf(
  object = spata_obj, 
  across = "seurat_clusters",
  across_subset = '7',
  method_de = "wilcox", 
  signif_threshold = 0.05 # extract top 20 most significant gene sets
) 
gseaDF_HM <- gseaDF %>% dplyr::filter(grepl("^HM|BP", gseaDF$label))
gseaDF_HM <- gseaDF_HM[order(gseaDF_HM$fdr, decreasing = T),]
gseaDF_HM$label <- factor(gseaDF_HM$label, levels = gseaDF_HM$label)
gseaDF_HM <- tail(gseaDF_HM, 7)

p6.7 <- ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')


pdf(file = './20220823 整理Figure/YangFu HCC2 A.pdf', height = 5, width = 20)
p1|p2 |p3 |p4 
dev.off()

pdf(file = './20220823 整理Figure/YangFu HCC2 B.pdf', height = 5, width = 30)
p5.5|p5.7|p6.5|p6.7
dev.off()


plotViolinplot(spata_obj, 'IL1B',
               across = 'seurat_clusters',
               test_groupwise = 'kruskal.test') 


####IL1B
df <- getExpressionMatrix(spata_obj) %>% t() %>% as.data.frame()
df_int <- df %>% select('IL1B')
df_int$barcodes <- rownames(df_int)
joined_df <- merge(joined_df, df_int, by='barcodes')


# temp <- joined_df %>% dplyr::filter(seurat_clusters == 5) %>% select(mySig_Mac_VIM_high, IL1B)
# ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'IL1B',
#           color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
#           add = "reg.line",  # Add regressin line
#           add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
#           conf.int = TRUE, # Add confidence interval
#           cor.coef = T, # Add correlation coefficient. see ?stat_cor
#           cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
#           title =  'seurat-cluster 5')

temp <- joined_df %>% dplyr::filter(seurat_clusters == 7) %>% select(mySig_Mac_VIM_high, IL1B)
ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'IL1B',
          fill ='#8491B4FF',color='darkgrey',size = 2,shape = 21, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#7E6148FF", fill = "#B09C85FF"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = T, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
          title =  'seurat-cluster 7')


# YangFu HCC3 (HCC-R2) cluster 6 ----------------------------------------------------------------

##cluster 6

rm(list = ls())
gc()



load(file='./20220823 整理Figure/YangFU-HCC3.Rdata')




p1 <- plotSurfaceAverage(object = spata_obj, 
                         color_by = Mac_VIM_high, 
                         smooth = TRUE, 
                         pt_size = 2,
                         smooth_span = 0.2, 
                         pt_clrsp = "inferno")
p2 <- plotSurfaceAverage(object = spata_obj, 
                         color_by = T.reg, 
                         smooth = TRUE, 
                         pt_size = 2,
                         smooth_span = 0.2, 
                         pt_clrsp = "inferno")


p3 <- plotSurface(object = spata_obj, color_by = "seurat_clusters", display_image = F, pt_alpha = 1) # nolint

p4 <- plotDistributionAcross2(df = joined_df,
                              variables = mySig_gs,
                              across = "seurat_clusters",
                              plot_type = "violin",
                              clrp = 'jco',
                              nrow=2)


########correlation between Mac-vim-hi and Treg score


temp <- joined_df %>% dplyr::filter(seurat_clusters == 6) %>% select(mySig_Mac_VIM_high, mySig_T.reg)

p5 <- ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
                color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = T, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
                title =  'seurat-cluster 6')


# 
# 
# temp <- joined_df %>% dplyr::filter(seurat_clusters == 7) %>% select(mySig_Mac_VIM_high, mySig_T.reg)
# 
# p5.7 <- ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
#                   color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
#                   add = "reg.line",  # Add regressin line
#                   add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
#                   conf.int = TRUE, # Add confidence interval
#                   cor.coef = T, # Add correlation coefficient. see ?stat_cor
#                   cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
#                   title =  'seurat-cluster 7')





###cluster6
gseaDF <- getGseaDf(
  object = spata_obj, 
  across = "seurat_clusters",
  across_subset = '6',
  method_de = "wilcox", 
  signif_threshold = 0.05 # extract top 20 most significant gene sets
) 
gseaDF_HM <- gseaDF %>% dplyr::filter(grepl("^HM|BP", gseaDF$label))
gseaDF_HM <- gseaDF_HM[order(gseaDF_HM$fdr, decreasing = T),]
gseaDF_HM$label <- factor(gseaDF_HM$label, levels = gseaDF_HM$label)
gseaDF_HM <- tail(gseaDF_HM, 7)

p6 <- ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')


# 
# ###cluster7
# gseaDF <- getGseaDf(
#   object = spata_obj, 
#   across = "seurat_clusters",
#   across_subset = '7',
#   method_de = "wilcox", 
#   signif_threshold = 0.05 # extract top 20 most significant gene sets
# ) 
# gseaDF_HM <- gseaDF %>% dplyr::filter(grepl("^HM|BP", gseaDF$label))
# gseaDF_HM <- gseaDF_HM[order(gseaDF_HM$fdr, decreasing = T),]
# gseaDF_HM$label <- factor(gseaDF_HM$label, levels = gseaDF_HM$label)
# gseaDF_HM <- tail(gseaDF_HM, 7)
# 
# p6.7 <- ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')


pdf(file = './20220823 整理Figure/YangFu HCC3.pdf', height = 5, width = 30)
p1|p2 |p3 |p4 |p5|p6
dev.off()

plotViolinplot(spata_obj, 'IL1B',
               across = 'seurat_clusters',
               test_groupwise = 'kruskal.test'
)


####IL1B 
df <- getExpressionMatrix(spata_obj) %>% t() %>% as.data.frame()
df_int <- df %>% select('IL1B')
df_int$barcodes <- rownames(df_int)
joined_df <- merge(joined_df, df_int, by='barcodes')


temp <- joined_df %>% dplyr::filter(seurat_clusters == 6) %>% select(mySig_Mac_VIM_high, IL1B)
ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'IL1B',
          fill ='#8491B4FF',color='darkgrey',size = 2,shape = 21, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#7E6148FF", fill = "#B09C85FF"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = T, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
          title =  'seurat-cluster 6')


# 226766_C_L_orig  (HCC-4-L) cluster 4 & 6 ----------------------------------------------------------------

##cluster 4 EMT
##cluster 6 immune

rm(list = ls())
gc()


load(file='./20220823 整理Figure/226766_C_L.Rdata')



p1 <- plotSurfaceAverage(object = spata_obj, 
                         color_by = Mac_VIM_high, 
                         smooth = TRUE, 
                         pt_size = 2,
                         smooth_span = 0.2, 
                         pt_clrsp = "inferno")
p2 <- plotSurfaceAverage(object = spata_obj, 
                         color_by = T.reg, 
                         smooth = TRUE, 
                         pt_size = 2,
                         smooth_span = 0.2, 
                         pt_clrsp = "inferno")

p3 <- plotSurface(object = spata_obj, color_by = "seurat_clusters", display_image = F, pt_alpha = 1) # nolint


p4 <- plotDistributionAcross2(df = joined_df,
                              variables = mySig_gs,
                              across = "seurat_clusters",
                              plot_type = "violin",
                              clrp = 'jco',
                              nrow=2)


########correlation between Mac-vim-hi and Treg score


temp <- joined_df %>% dplyr::filter(seurat_clusters == 4) %>% select(mySig_Mac_VIM_high, mySig_T.reg)

p5.4 <- ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
                  color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = T, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
                  title =  'seurat-cluster 4')




temp <- joined_df %>% dplyr::filter(seurat_clusters == 6) %>% select(mySig_Mac_VIM_high, mySig_T.reg)

p5.6 <- ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
                  color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = T, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
                  title =  'seurat-cluster 6')





###cluster4
gseaDF <- getGseaDf(
  object = spata_obj, 
  across = "seurat_clusters",
  across_subset = '4',
  method_de = "wilcox", 
  signif_threshold = 0.05 # extract top 20 most significant gene sets
) 
gseaDF_HM <- gseaDF %>% dplyr::filter(grepl("^HM|BP", gseaDF$label))
gseaDF_HM <- gseaDF_HM[order(gseaDF_HM$fdr, decreasing = T),]
gseaDF_HM$label <- factor(gseaDF_HM$label, levels = gseaDF_HM$label)
gseaDF_HM <- tail(gseaDF_HM, 7)

p6.4 <- ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')



###cluster6
gseaDF <- getGseaDf(
  object = spata_obj, 
  across = "seurat_clusters",
  across_subset = '6',
  method_de = "wilcox", 
  signif_threshold = 0.05 # extract top 20 most significant gene sets
) 
gseaDF_HM <- gseaDF %>% dplyr::filter(grepl("^HM|BP", gseaDF$label))
gseaDF_HM <- gseaDF_HM[order(gseaDF_HM$fdr, decreasing = T),]
gseaDF_HM$label <- factor(gseaDF_HM$label, levels = gseaDF_HM$label)
gseaDF_HM <- tail(gseaDF_HM, 7)

p6.6 <- ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')


pdf(file = './20220823 整理Figure/226766L A.pdf', height = 5, width = 20)
p1|p2 |p3 |p4
dev.off()

pdf(file = './20220823 整理Figure/226766L B.pdf', height = 5, width = 30)
p5.4|p5.6|p6.4|p6.6
dev.off()



plotViolinplot(spata_obj, 'IL1B',
               across = 'seurat_clusters',
               test_groupwise = 'kruskal.test'
)

####IL1B correlation
df <- getExpressionMatrix(spata_obj) %>% t() %>% as.data.frame()
df_int <- df %>% select('IL1B')
df_int$barcodes <- rownames(df_int)
joined_df <- merge(joined_df, df_int, by='barcodes')


temp <- joined_df %>% dplyr::filter(seurat_clusters == 4) %>% select(mySig_Mac_VIM_high, IL1B)
ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'IL1B',
          fill ='#8491B4FF',color='darkgrey',size = 2,shape = 21, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#7E6148FF", fill = "#B09C85FF"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = T, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
          title =  'seurat-cluster 4')

# temp <- joined_df %>% dplyr::filter(seurat_clusters == 6) %>% select(mySig_Mac_VIM_high, IL1B)
# ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'IL1B',
#           color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
#           add = "reg.line",  # Add regressin line
#           add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
#           conf.int = TRUE, # Add confidence interval
#           cor.coef = T, # Add correlation coefficient. see ?stat_cor
#           cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
#           title =  'seurat-cluster 6')

