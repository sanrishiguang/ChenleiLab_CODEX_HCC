

###R 4.1.0
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






# 227312-T : cluster3----------------------------------------------------------------

spata_obj <- readRDS(file = './227312_T solo/spata-obj-227312-T-DSEA.RDS' )

##gene signature
spata_obj <- addGeneSet(spata_obj, class_name = 'mySig', gs_name = 'Mac_VIM_high', genes=Mac_VIM_high, overwrite = T)
spata_obj <- addGeneSet(spata_obj, class_name = 'mySig', gs_name = 'T.reg', genes=T.reg, overwrite = T)

getGeneSetOverview(spata_obj)

spata_obj <- adjustDefaultInstructions(object = spata_obj,pt_clrp = "jco") #jco, npg, milo, nejm

plotSurface(object = spata_obj, color_by = "seurat_clusters", display_image = F, pt_alpha = 1) # nolint

plotSurfaceAverage(object = spata_obj, 
                   color_by = Mac_VIM_high, 
                   smooth = TRUE, 
                   pt_size = 2,
                   smooth_span = 0.2, 
                   pt_clrsp = "inferno")
plotSurfaceAverage(object = spata_obj, 
                   color_by = T.reg, 
                   smooth = TRUE, 
                   pt_size = 2,
                   smooth_span = 0.2, 
                   pt_clrsp = "inferno")

# get only names of discrete features
discrete_features <-
  getFeatureNames(object = spata_obj, of_class = c("character", "factor"))

# output
discrete_features

# create a spata-data.frame
mySig_gs <- getGeneSets(object = spata_obj, of_class = "mySig")
joined_df <- joinWith(object = spata_obj,
                      spata_df = getCoordsDf(spata_obj),
                      gene_sets = mySig_gs,
                      features = discrete_features)

#output
joined_df

#violin
plotDistributionAcross2(df = joined_df,
                        variables = mySig_gs,
                        across = "seurat_clusters",
                        plot_type = "violin",
                        clrp = 'jco',
                        nrow=2)


########correlation between Mac-vim-hi and Treg score

##seurat-cluster 3 is the co-occurrence cluster

temp <- joined_df %>% dplyr::filter(seurat_clusters == 3) %>% select(mySig_Mac_VIM_high, mySig_T.reg)

ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
          color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = T, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
          title =  'seurat-cluster 4')



########enriched terms in seurat-cluster 4

# spata_obj <- 
#   runGsea(
#     object = spata_obj, 
#     across = "seurat_clusters",
#     methods_de = "wilcox" 
#   )

gseaDF <- getGseaDf(
  object = spata_obj, 
  across = "seurat_clusters",
  across_subset = '3',
  method_de = "wilcox", 
  signif_threshold = 0.05 # extract top 20 most significant gene sets
) 

gseaDF_HM <- gseaDF %>% dplyr::filter(grepl("^HM|BP|RCTM", label))
gseaDF_HM <- gseaDF_HM[order(gseaDF_HM$fdr, decreasing = T),]
gseaDF_HM$label <- factor(gseaDF_HM$label, levels = gseaDF_HM$label)
gseaDF_HM <- tail(gseaDF_HM, 7)
library(ggplot2)
ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')


save(list = ls(), file='./20220823 整理Figure/227312_T.Rdata')











# 084717_A_orig (HCC-5-A): cluster 0&2&6 ----------------------------------------------------------------

##cluster 0 & 2 & 6 


pathname <- '084717_A_orig'
spata_obj <- readRDS(file = paste0('./',pathname ,'/',pathname, '_DEA.RDS' ))

##gene signature
spata_obj <- addGeneSet(spata_obj, class_name = 'mySig', gs_name = 'Mac_VIM_high', genes=Mac_VIM_high, overwrite = T)
spata_obj <- addGeneSet(spata_obj, class_name = 'mySig', gs_name = 'T.reg', genes=T.reg, overwrite = T)

getGeneSetOverview(spata_obj)

spata_obj <- adjustDefaultInstructions(object = spata_obj,pt_clrp = "jco") #jco, npg, milo, nejm

plotSurface(object = spata_obj, color_by = "seurat_clusters", display_image = F, pt_alpha = 1) # nolint

plotSurfaceAverage(object = spata_obj, 
                   color_by = Mac_VIM_high, 
                   smooth = TRUE, 
                   pt_size = 2,
                   smooth_span = 0.2, 
                   pt_clrsp = "inferno")
plotSurfaceAverage(object = spata_obj, 
                   color_by = T.reg, 
                   smooth = TRUE, 
                   pt_size = 2,
                   smooth_span = 0.2, 
                   pt_clrsp = "inferno")

# get only names of discrete features
discrete_features <-
  getFeatureNames(object = spata_obj, of_class = c("character", "factor"))

# output
discrete_features

# create a spata-data.frame
mySig_gs <- getGeneSets(object = spata_obj, of_class = "mySig")
joined_df <- joinWith(object = spata_obj,
                      spata_df = getCoordsDf(spata_obj),
                      gene_sets = mySig_gs,
                      features = discrete_features)

#output
joined_df

#violin
plotDistributionAcross2(df = joined_df,
                        variables = mySig_gs,
                        across = "seurat_clusters",
                        plot_type = "violin",
                        clrp = 'jco',
                        nrow=2)


########correlation between Mac-vim-hi and Treg score

##seurat-cluster 2 is the co-occurrence cluster

temp <- joined_df %>% dplyr::filter(seurat_clusters == 2) %>% select(mySig_Mac_VIM_high, mySig_T.reg)
library(ggpubr)

ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
          color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = T, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
          title =  'seurat-cluster 2')



########enriched terms in seurat-cluster 2

spata_obj <- 
  runGsea(
    object = spata_obj, 
    across = "seurat_clusters",
    methods_de = "wilcox" 
  )

gseaDF <- getGseaDf(
  object = spata_obj, 
  across = "seurat_clusters",
  across_subset = '2',
  method_de = "wilcox", 
  signif_threshold = 0.05 # extract top 20 most significant gene sets
) 

gseaDF_HM <- gseaDF %>% dplyr::filter(grepl("^HM|BP|RCTM", gseaDF$label))
gseaDF_HM <- gseaDF_HM[order(gseaDF_HM$fdr, decreasing = T),]
gseaDF_HM$label <- factor(gseaDF_HM$label, levels = gseaDF_HM$label)
gseaDF_HM <- tail(gseaDF_HM, 7)
library(ggplot2)
ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')


save(list = ls(), file='./20220823 整理Figure/084717_A.Rdata')





# YangFu HCC2 (HCC-R1) : cluster 5 & 7-------------------------------------------------------------


spata_obj <- readRDS(file = './YangFu HCC2/spata-obj-YangFu_HCC2-DSEA.RDS' )

##gene signature
spata_obj <- addGeneSet(spata_obj, class_name = 'mySig', gs_name = 'Mac_VIM_high', genes=Mac_VIM_high, overwrite = T)
spata_obj <- addGeneSet(spata_obj, class_name = 'mySig', gs_name = 'T.reg', genes=T.reg, overwrite = T)

getGeneSetOverview(spata_obj)

spata_obj <- adjustDefaultInstructions(object = spata_obj,pt_clrp = "jco") #jco, npg, milo, nejm

plotSurface(object = spata_obj, color_by = "seurat_clusters", display_image = F, pt_alpha = 1) # nolint

plotSurfaceAverage(object = spata_obj, 
                   color_by = Mac_VIM_high, 
                   smooth = TRUE, 
                   pt_size = 2,
                   smooth_span = 0.2, 
                   pt_clrsp = "inferno")
plotSurfaceAverage(object = spata_obj, 
                   color_by = T.reg, 
                   smooth = TRUE, 
                   pt_size = 2,
                   smooth_span = 0.2, 
                   pt_clrsp = "inferno")

# get only names of discrete features
discrete_features <-
  getFeatureNames(object = spata_obj, of_class = c("character", "factor"))

# output
discrete_features

# create a spata-data.frame
mySig_gs <- getGeneSets(object = spata_obj, of_class = "mySig")
joined_df <- joinWith(object = spata_obj,
                      spata_df = getCoordsDf(spata_obj),
                      gene_sets = mySig_gs,
                      features = discrete_features)

#output
joined_df

#violin
plotDistributionAcross2(df = joined_df,
                        variables = mySig_gs,
                        across = "seurat_clusters",
                        plot_type = "violin",
                        clrp = 'jco',
                        nrow=2)


########correlation between Mac-vim-hi and Treg score

##seurat-cluster 7 is the co-occurrence cluster

temp <- joined_df %>% dplyr::filter(seurat_clusters == 5) %>% select(mySig_Mac_VIM_high, mySig_T.reg)

ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
          color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = T, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
          title =  'seurat-cluster 7')



########enriched terms in seurat-cluster 7

# spata_obj <- 
#   runGsea(
#     object = spata_obj, 
#     across = "seurat_clusters",
#     methods_de = "wilcox" 
#   )

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
library(ggplot2)
ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')


save(list = ls(), file='./20220823 整理Figure/YangFU-HCC2.Rdata')



# YangFu HCC3 (HCC-R2) cluster 5 & 7-------------------------------------------------------------


spata_obj <- readRDS(file = './YangFu HCC3/spata-obj-YangFu_HCC3-DSEA.RDS' )

##gene signature
spata_obj <- addGeneSet(spata_obj, class_name = 'mySig', gs_name = 'Mac_VIM_high', genes=Mac_VIM_high, overwrite = T)
spata_obj <- addGeneSet(spata_obj, class_name = 'mySig', gs_name = 'T.reg', genes=T.reg, overwrite = T)

getGeneSetOverview(spata_obj)

spata_obj <- adjustDefaultInstructions(object = spata_obj,pt_clrp = "jco") #jco, npg, milo, nejm

plotSurface(object = spata_obj, color_by = "seurat_clusters", display_image = F, pt_alpha = 1) # nolint

plotSurfaceAverage(object = spata_obj, 
                   color_by = Mac_VIM_high, 
                   smooth = TRUE, 
                   pt_size = 2,
                   smooth_span = 0.2, 
                   pt_clrsp = "inferno")
plotSurfaceAverage(object = spata_obj, 
                   color_by = T.reg, 
                   smooth = TRUE, 
                   pt_size = 2,
                   smooth_span = 0.2, 
                   pt_clrsp = "inferno")

# get only names of discrete features
discrete_features <-
  getFeatureNames(object = spata_obj, of_class = c("character", "factor"))

# output
discrete_features

# create a spata-data.frame
mySig_gs <- getGeneSets(object = spata_obj, of_class = "mySig")
joined_df <- joinWith(object = spata_obj,
                      spata_df = getCoordsDf(spata_obj),
                      gene_sets = mySig_gs,
                      features = discrete_features)

#output
joined_df

#violin
plotDistributionAcross2(df = joined_df,
                        variables = mySig_gs,
                        across = "seurat_clusters",
                        plot_type = "violin",
                        clrp = 'jco',
                        nrow=2)


########correlation between Mac-vim-hi and Treg score

##seurat-cluster 6 is the co-occurrence cluster

temp <- joined_df %>% dplyr::filter(seurat_clusters == 6) %>% select(mySig_Mac_VIM_high, mySig_T.reg)

ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
          color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = T, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
          title =  'seurat-cluster 6')



########enriched terms in seurat-cluster 6

# spata_obj <- 
#   runGsea(
#     object = spata_obj, 
#     across = "seurat_clusters",
#     methods_de = "wilcox" 
#   )

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
library(ggplot2)
ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')


save(list = ls(), file='./20220823 整理Figure/YangFU-HCC3.Rdata')




# 226766_C_L_orig (HCC-4-L) -----------------------------------------------------------------------

##cluster 4 EMT
##cluster 6 immune


pathname <- '226766_C_L_orig'
spata_obj <- readRDS(file = paste0('./',pathname ,'/',pathname, '_DEA.RDS' ))

##gene signature
spata_obj <- addGeneSet(spata_obj, class_name = 'mySig', gs_name = 'Mac_VIM_high', genes=Mac_VIM_high, overwrite = T)
spata_obj <- addGeneSet(spata_obj, class_name = 'mySig', gs_name = 'T.reg', genes=T.reg, overwrite = T)

getGeneSetOverview(spata_obj)

spata_obj <- adjustDefaultInstructions(object = spata_obj,pt_clrp = "jco") #jco, npg, milo, nejm

plotSurface(object = spata_obj, color_by = "seurat_clusters", display_image = F, pt_alpha = 1) # nolint

plotSurfaceAverage(object = spata_obj, 
                   color_by = Mac_VIM_high, 
                   smooth = TRUE, 
                   pt_size = 2,
                   smooth_span = 0.2, 
                   pt_clrsp = "inferno")
plotSurfaceAverage(object = spata_obj, 
                   color_by = T.reg, 
                   smooth = TRUE, 
                   pt_size = 2,
                   smooth_span = 0.2, 
                   pt_clrsp = "inferno")

# get only names of discrete features
discrete_features <-
  getFeatureNames(object = spata_obj, of_class = c("character", "factor"))

# output
discrete_features

# create a spata-data.frame
mySig_gs <- getGeneSets(object = spata_obj, of_class = "mySig")
joined_df <- joinWith(object = spata_obj,
                      spata_df = getCoordsDf(spata_obj),
                      gene_sets = mySig_gs,
                      features = discrete_features)

#output
joined_df

#violin
plotDistributionAcross2(df = joined_df,
                        variables = mySig_gs,
                        across = "seurat_clusters",
                        plot_type = "violin",
                        clrp = 'jco',
                        nrow=2)


########correlation between Mac-vim-hi and Treg score

##seurat-cluster 4 is the co-occurrence cluster

temp <- joined_df %>% dplyr::filter(seurat_clusters == 4) %>% select(mySig_Mac_VIM_high, mySig_T.reg)
library(ggpubr)

ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
          color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = T, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
          title =  'seurat-cluster 4')

##seurat-cluster 6 is the co-occurrence cluster

temp <- joined_df %>% dplyr::filter(seurat_clusters == 6) %>% select(mySig_Mac_VIM_high, mySig_T.reg)
library(ggpubr)

ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
          color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = T, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
          title =  'seurat-cluster 6')





########enriched terms in seurat-cluster 4

spata_obj <- 
  runGsea(
    object = spata_obj, 
    across = "seurat_clusters",
    methods_de = "wilcox" 
  )

#cluster 4
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
library(ggplot2)
ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')



#cluster 6
gseaDF <- getGseaDf(
  object = spata_obj, 
  across = "seurat_clusters",
  across_subset = '6',
  method_de = "wilcox", 
  signif_threshold = 0.05 # extract top 20 most significant gene sets
) 

gseaDF_HM <- gseaDF %>% dplyr::filter(grepl("^HM|BP|RCTM", label))
gseaDF_HM <- gseaDF_HM[order(gseaDF_HM$fdr, decreasing = T),]
gseaDF_HM$label <- factor(gseaDF_HM$label, levels = gseaDF_HM$label)
gseaDF_HM <- tail(gseaDF_HM, 7)
library(ggplot2)
ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')

save(list = ls(), file='./20220823 整理Figure/226766_C_L.Rdata')


# 225847-T (HCC-1-N) cluster4----------------------------------------------------------------

pathname <- '225847_A_T_orig'
spata_obj <- readRDS(file = paste0('./',pathname ,'/',pathname, '_DEA.RDS' ))

##gene signature
spata_obj <- addGeneSet(spata_obj, class_name = 'mySig', gs_name = 'Mac_VIM_high', genes=Mac_VIM_high, overwrite = T)
spata_obj <- addGeneSet(spata_obj, class_name = 'mySig', gs_name = 'T.reg', genes=T.reg, overwrite = T)

getGeneSetOverview(spata_obj)

spata_obj <- adjustDefaultInstructions(object = spata_obj,pt_clrp = "jco") #jco, npg, milo, nejm

plotSurface(object = spata_obj, color_by = "seurat_clusters", display_image = F, pt_alpha = 1) # nolint

plotSurfaceAverage(object = spata_obj, 
                   color_by = Mac_VIM_high, 
                   smooth = TRUE, 
                   pt_size = 2,
                   smooth_span = 0.2, 
                   pt_clrsp = "inferno")
plotSurfaceAverage(object = spata_obj, 
                   color_by = T.reg, 
                   smooth = TRUE, 
                   pt_size = 2,
                   smooth_span = 0.2, 
                   pt_clrsp = "inferno")

# get only names of discrete features
discrete_features <-
  getFeatureNames(object = spata_obj, of_class = c("character", "factor"))

# output
discrete_features

# create a spata-data.frame
mySig_gs <- getGeneSets(object = spata_obj, of_class = "mySig")
joined_df <- joinWith(object = spata_obj,
                      spata_df = getCoordsDf(spata_obj),
                      gene_sets = mySig_gs,
                      features = discrete_features)

#output
joined_df

#violin
plotDistributionAcross2(df = joined_df,
                        variables = mySig_gs,
                        across = "seurat_clusters",
                        plot_type = "violin",
                        clrp = 'jco',
                        nrow=2)


########correlation between Mac-vim-hi and Treg score

##seurat-cluster 4 is the co-occurrence cluster

temp <- joined_df %>% dplyr::filter(seurat_clusters == 4) %>% select(mySig_Mac_VIM_high, mySig_T.reg)
library(ggpubr)

ggscatter(temp, x = 'mySig_Mac_VIM_high', y = 'mySig_T.reg',
          color = "darkgrey", shape = 21, fill='grey',size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = T, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"),
          title =  'seurat-cluster 4')



########enriched terms in seurat-cluster 4

spata_obj <- 
  runGsea(
    object = spata_obj, 
    across = "seurat_clusters",
    methods_de = "wilcox" 
  )

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
library(ggplot2)
ggplot(gseaDF_HM, mapping = aes(x=-log(fdr), y=label, fill=pval)) + geom_bar(stat = 'identity')


save(list = ls(), file='./20220823 整理Figure/225847_T.Rdata')





