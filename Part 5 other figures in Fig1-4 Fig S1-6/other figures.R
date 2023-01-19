


###   R 4.1.0


# library(data.table)
library(tidyverse)  #‘1.3.2’
library(tidyr) #‘1.2.0’
library(ggpubr) # ‘0.4.0’
library(ggthemes) #‘4.2.4’
library(Rtsne)  # ‘0.16’
library(Rcpp)  #‘1.0.8.3’
library(igraph)  #‘1.3.2’
library(ggthemes)  # ‘4.2.4’
library(survival)   #‘3.3.1’
library(survminer)   # ‘0.4.9’
library(pheatmap)  # ‘1.0.12’
library(ggsci)  # ‘2.9’
library(reshape2) #‘1.4.4’


mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                 legend.key = element_rect(fill = "white", colour = "white"), #图标
                 legend.background = (element_rect(colour= "white", fill = "white")))



# Fig. S6C-D characteristics of SP ----------------------------------------

load('import.Rdata')

###CN: Fig.S6A
temp <- mydata_pheno_anno[c(var_allcn,'SP')] 
temp_agg <- aggregate(.~ SP,data=temp,  FUN = mean, na.rm=T)
library(reshape2)
temp_agg_melt <- melt(temp_agg, id='SP')
# ggplot(temp_agg_melt, aes(x=variable, y=SP)) +
#   geom_point(aes(size=value, fill=value), shape=21, colour='black') +
#   scale_size_area(max_size = 10) +
#   scale_fill_gradient2(low = "blue", high = 'red', mid = "white", midpoint = 5)+
#   theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1)) +
#   coord_flip()

ggplot(temp_agg_melt, aes(x=variable, y=SP)) +
  geom_point(aes(size=value, fill=value), shape=21, colour='black') +
  scale_size_area(max_size = 8) +
  scale_fill_gradient2(low = "blue", high = 'red', mid = "white", midpoint = 5)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1)) 

##subtype
temp <- mydata_pheno_anno[c(var_subtypes,'SP')] 
temp_agg <- aggregate(.~ SP,data=temp,  FUN = mean, na.rm=T)
library(reshape2)
temp_agg_melt <- melt(temp_agg, id='SP')

ggplot(temp_agg_melt, aes(x=variable, y=SP)) +
  geom_point(aes(size=value, fill=value), shape=21, colour='black') +
  scale_size_area(max_size = 8) +
  scale_fill_gradient2(low = "blue", high = 'red', mid = "white", midpoint = 5)+
  theme(axis.text.x=element_text(angle=90,  hjust=1, vjust=0)) 
# coord_flip()
# ggsave(p, file='./random forest/5SP_ct_w.png',width = 12, height = 4)


#CT: Fig.S6B
temp <- mydata_pheno_anno[c(var_ct,'SP')] 
temp_agg <- aggregate(.~ SP,data=temp,  FUN = mean, na.rm=T)

temp_agg[,2:12] = scale(temp_agg[,2:12])

library(reshape2)
temp_agg_melt <- melt(temp_agg, id='SP')

ggplot(temp_agg_melt, aes(x=variable, y=SP)) +
  geom_point(aes(size=value, fill=value), shape=21, colour='black') +
  scale_size(range=c(-2,10)) +
  scale_fill_gradient2(low = "blue", high = 'red', mid = "white", midpoint = 0.2)+
  theme(axis.text.x=element_text(angle=90,  hjust=1, vjust=0)) 





# Fig. S4D boxplot of CP-based-CN-----------------------------------------------------------------

mydata_pheno_anno_TP <- readRDS("Class_Pheno_Anno.Rds")

mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), 
                 legend.key = element_rect(fill = "white", colour = "white"), 
                 legend.background = (element_rect(colour= "white", fill = "white")))

colnames(mydata_pheno_anno_TP)

df <- mydata_pheno_anno_TP %>% select(Class, SP,CN0_Hepar1_LII : CN9_CD44)
df_melt <- melt(df,id=c('Class',"SP"))

p <- ggplot(df_melt, aes(x=variable, y=value,color=variable)) + 
  geom_boxplot(outlier.alpha = 0.1)+
  #coord_flip()+
  ggtitle("The percentage of CN in samples")+ 
  xlab('')+
  ylab('Percentage')+ 
  mytheme+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=.5))+
  guides(color=F)+scale_color_manual(values=c(pal_npg('nrc',alpha = 1.0)(10),pal_npg('nrc',alpha = 0.9)(10),pal_npg('nrc',alpha = 0.7)(10),pal_npg('nrc',alpha = 0.8)(10)))+
  ylim(0, 20)
p






# #Fig. S4B boxplot of CT-based-CN----------------------------------------------------------------

mydata_pheno_anno_TP <- read_csv(file = 'E:\\myPythonProject\\TrainingData\\20220423 celltype CN\\cells_r=50_CN=10.csv')

mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), 
                 legend.key = element_rect(fill = "white", colour = "white"), 
                 legend.background = (element_rect(colour= "white", fill = "white")))

colnames(mydata_pheno_anno_TP)

df <- mydata_pheno_anno_TP %>% select(Class, neighborhood10)
df$neighborhood10 = paste0("CN_",df$neighborhood10)

cycle = unique(df$Class)
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(df, Class == cycle[i])
  Freq = bind_rows(Freq, round((table(sub$neighborhood10)/nrow(sub))*100,4))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0
Freq$Class = rownames(Freq)
Freq = melt(Freq,id="Class")
Freq$value = as.numeric(Freq$value)
Freq$variable = as.character(Freq$variable)

Freq$variable[Freq$variable=='CN_0'] <- 'CN0_ImmMix_B' 
Freq$variable[Freq$variable=='CN_1'] <- 'CN1_Tumor'
Freq$variable[Freq$variable=='CN_2'] <- 'CN2_ImmMix_CD8T' 
Freq$variable[Freq$variable=='CN_3'] <- 'CN3_Tumor_ImmMix_Endo'
Freq$variable[Freq$variable=='CN_4'] <- 'CN4_Tumor_ImmMix_Fibr'
Freq$variable[Freq$variable=='CN_5'] <- 'CN5_ImmMix_Fibr'
Freq$variable[Freq$variable=='CN_6'] <- 'CN6_ImmMix_Lym' 
Freq$variable[Freq$variable=='CN_7'] <- 'CN7_ImmMix_Bil' 
Freq$variable[Freq$variable=='CN_8'] <- 'CN8_Tumor_ImmMix_CD8T'  
Freq$variable[Freq$variable=='CN_9'] <- 'CN9_Tumor_ImmMix_Mac'  



ggplot(Freq,aes(x=variable, y=value,color=variable)) + 
  geom_boxplot(outlier.alpha = 0.1)+
  #coord_flip()+
  ggtitle("The percentage of CN in samples")+ 
  xlab('')+
  ylab('Percentage')+ 
  mytheme+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=.5))+
  guides(color=F)+scale_color_manual(values=c(pal_npg('nrc',alpha = 1.0)(10)))




# Barplor in Fig. S4C CP-based-CN  cell type constitution------------------------------------------------------


load("All_CN_TP_0424.Rdata")
Freq = subset(All_CN_TP, select = c("celltype","All_CN"))
Freq_toplot = c()
cycle = unique(Freq$All_CN)
for (i in 1:length(cycle)) {
  temp = subset(Freq, All_CN == cycle[i])
  Freq_toplot = bind_rows(Freq_toplot, table(temp$celltype))
}
Freq_toplot = as.data.frame(Freq_toplot)
rownames(Freq_toplot) = cycle
Freq_toplot[is.na(Freq_toplot)]=0
Freq_toplot$CN = rownames(Freq_toplot)
Freq_toplot = melt(Freq_toplot, id.vars = "CN")
Freq_toplot = Freq_toplot[order(Freq_toplot$CN),]
CN = str_split_fixed(Freq_toplot$CN, "_",3)
CN = as.data.frame(CN)
Freq_toplot$CN_r = CN$V1
Freq_toplot$CN_r = substr(Freq_toplot$CN_r, 3,4)
Freq_toplot$CN_r = as.numeric(Freq_toplot$CN_r)
Freq_toplot = Freq_toplot[order(Freq_toplot$CN_r),]
# Freq_toplot$CN_r = factor(Freq_toplot$CN_r, levels = rev(c(29,20,28,37,19,
#                                                           15,39,21,14,27,
#                                                           7,18,26,9,12,
#                                                           3,33,4,36,30,
#                                                           35,8,31,22,24,
#                                                           17,32,11,23,13,
#                                                           6,25,38,0,16,
#                                                           2,34,10,1,5)))

Freq_toplot$CN_r = factor(Freq_toplot$CN_r, levels = rev(c(20,28,15,39,29,
                                                           14,38,27,1,21,
                                                           7,18,26,12,30,
                                                           35,4,36,3,33,
                                                           31,8,11,23,13,
                                                           24,17,32,10,5,
                                                           34,2,16,0,25,
                                                           6,22,19,9,37)))

ggplot(Freq_toplot, aes(x = CN_r, y = value)) + 
  geom_bar(aes(fill = factor(variable)), stat = "identity", position = "fill") + 
  scale_fill_manual(values = c(pal_npg('nrc',alpha = 0.6)(11))) +
  scale_y_continuous(labels = scales::percent) + 
  coord_flip()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())





# Barplot in Fig. S4A CT-based-CN  cell type constitution------------------------------------------------------


mydata_pheno_anno_TP <- read_csv(file = 'E:\\myPythonProject\\TrainingData\\20220423 celltype CN\\cells_r=50_CN=10.csv')
Freq = subset(mydata_pheno_anno_TP, select = c("celltype","neighborhood10"))
Freq_toplot = c()
cycle = unique(Freq$neighborhood10)
for (i in 1:length(cycle)) {
  temp = subset(Freq, neighborhood10 == cycle[i])
  Freq_toplot = bind_rows(Freq_toplot, table(temp$celltype))
}
Freq_toplot = as.data.frame(Freq_toplot)
rownames(Freq_toplot) = cycle
Freq_toplot[is.na(Freq_toplot)]=0
Freq_toplot$CN = rownames(Freq_toplot)
Freq_toplot = melt(Freq_toplot, id.vars = "CN")
Freq_toplot = Freq_toplot[order(Freq_toplot$CN),]

# Freq_toplot$CN_r = factor(Freq_toplot$CN_r, levels = rev(c(29,20,28,37,19,
#                                                           15,39,21,14,27,
#                                                           7,18,26,9,12,
#                                                           3,33,4,36,30,
#                                                           35,8,31,22,24,
#                                                           17,32,11,23,13,
#                                                           6,25,38,0,16,
#                                                           2,34,10,1,5)))

Freq_toplot$CN = factor(Freq_toplot$CN, levels = rev(c(7,2,0,9,4,5,6,8,3,1)))

ggplot(Freq_toplot, aes(x = CN, y = value)) + 
  geom_bar(aes(fill = factor(variable)), stat = "identity", position = "fill") + 
  scale_fill_manual(values = c(pal_npg('nrc',alpha = 0.6)(11))) +
  scale_y_continuous(labels = scales::percent) + 
  coord_flip()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



# CP-based-CN CP constitution heatmap in Fig. S4C ------------------------------------------------------


library(tidyverse)
df <- read_csv(file = 'E:\\myPythonProject\\TrainingData\\20220317 AllsubtypeCN\\cells_r=50_CN=40.csv')
gc_csd_CN <- df %>% select(127, 55:126)

gc_csd_CN_anno <- gc_csd_CN %>% mutate(All_CN=paste0('CN',neighborhood10))

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN1'] <- 'CN1_HLA_DR'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN2'] <- 'CN2_others_LII'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN3'] <- 'CN3_c_Myc'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN4'] <- 'CN4_p53'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN5'] <- 'CN5_others'

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN6'] <- 'CN6_E_Cad_LII'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN7'] <- 'CN7_Fibro_HII'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN8'] <- 'CN8_CD107a'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN9'] <- 'CN9_CD44'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN10'] <- 'CN10_others'

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN11'] <- 'CN11_p_S6'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN12'] <- 'CN12_Glypican3'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN13'] <- 'CN13_Hepar1_LII'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN14'] <- 'CN14_ImmMix_Tumor'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN15'] <- 'CN15_ImmMix_B'

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN16'] <- 'CN16_Hepar1'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN17'] <- 'CN17_p_mTOR'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN18'] <- 'CN18_Fibro'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN19'] <- 'CN19_Biliary'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN20'] <- 'CN20_ImmMix_Mac'

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN21'] <- 'CN21_ImmMix_Lym'
# gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN22'] <- ''
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN23'] <- 'CN23_p_S6'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN24'] <- 'CN24_Twist1'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN25'] <- 'CN25_E_Cad'

# gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN26'] <- ''
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN27'] <- 'CN27_Endo'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN28'] <- 'CN28_CD44_ImmMix'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN29'] <- 'CN29_Hif1a'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN30'] <- 'CN30_Ki67'

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN31'] <- 'CN31_CD107a'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN32'] <- 'CN32_p_mTOR'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN33'] <- 'CN33_c_Myc_LII'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN34'] <- 'CN34_others_LII'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN35'] <- 'CN35_Ki67_LII'

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN36'] <- 'CN36_p53'
# gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN37'] <- ''
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN38'] <- 'CN38_Hepar1_HII' #HII, High Immune Infiltrate
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN39'] <- 'CN39_ImmuneMix_B'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN0'] <- 'CN0_Hepar1_LII'



head(gc_csd_CN_anno)
gc_csd_CN_anno <- gc_csd_CN_anno %>%  select(All_CN, everything()) 
gc_csd_CN_anno <- gc_csd_CN_anno[-2]

df_cn <- aggregate(. ~ All_CN, gc_csd_CN_anno, mean)
df_cn_scale <- df_cn[-1] %>% scale()
rownames(df_cn_scale) <- df_cn$All_CN

library(ComplexHeatmap)
library(circlize)
Heatmap(df_cn_scale, col = colorRamp2(c(-2, 2, 6), c("#35978f", "white", "#bf812d")))




# CT-based-CN CT constitution heatmap in Fig. S4A ------------------------------------------------------


df <- read_csv(file = 'E:\\myPythonProject\\TrainingData\\20220423 celltype CN\\cells_r=50_CN=10.csv')
gc_csd_CN <- df %>% select(67, 56:66)

gc_csd_CN_anno <- gc_csd_CN %>% mutate(All_CN=paste0('CN',neighborhood10))

head(gc_csd_CN_anno)
gc_csd_CN_anno <- gc_csd_CN_anno %>%  select(All_CN, everything()) 
gc_csd_CN_anno <- gc_csd_CN_anno[-2]

df_cn <- aggregate(. ~ All_CN, gc_csd_CN_anno, mean)
df_cn_scale <- df_cn[-1] %>% scale()
rownames(df_cn_scale) <- df_cn$All_CN

library(ComplexHeatmap)
library(circlize)
rownames(df_cn_scale) = c('CN0_ImmMix_B',
                          'CN1_Tumor',
                          'CN2_ImmMix_CD8T',
                          'CN3_Tumor_ImmMix_Endo',
                          'CN4_Tumor_ImmMix_Fibr',
                          'CN5_ImmMix_Fibr',
                          'CN6_ImmMix_Lym',
                          'CN7_ImmMix_Bil',
                          'CN8_Tumor_ImmMix_CD8T',
                          'CN9_Tumor_ImmMix_Mac')
Heatmap(df_cn_scale, col = colorRamp2(c(-4, 0, 4), c("#35978f", "white", "#bf812d")))






# Fig.S4E-F overall CD4T vs CN2 specific CD4T in prognosis -----------------------------------------------------------------

####CN2 specific CD4T
load('ct-based-cn-specific CT.Rdata')

sur.cut <- surv_cutpoint(df_percent_anno_cn2, time= 'RFSday',event = 'RFS01' , variables = 'CN2_ImmMix_CD8T_CD4T' , minprop = 0.1)
sur.cat <- surv_categorize(sur.cut)
names(sur.cat) <- c("RFSday",'RFS01','group')
fit <- survfit(Surv(RFSday, RFS01) ~ group, data = sur.cat)
ggsurvplot(fit,palette = "npg",
           risk.table = TRUE, pval = TRUE,
           conf.int = FALSE, xlab="Time in Days",
           ggtheme = theme_bw(),
           title = 'RFS: CN2_ImmMix_CD8T-specific CD4T enrichment',
           surv.median.line = "hv",  # add the median survival pointer.
           # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
           tables.y.text = T,
           risk.table.pos = "in",
           risk.table.col = "strata",
           fontsize = 3.5,
           pval.size = 4,
           surv.plot.height = 0.8,
           tables.height = 0.2,
           pval.coord = c(2000, 0.9),
           legend = "right"
)

sur.cut <- surv_cutpoint(df_percent_anno_cn2, time= 'Osday',event = 'OS01' , variables = 'CN2_ImmMix_CD8T_CD4T' , minprop = 0.1)
sur.cat <- surv_categorize(sur.cut)
names(sur.cat) <- c("RFSday",'RFS01','group')
fit <- survfit(Surv(RFSday, RFS01) ~ group, data = sur.cat)
ggsurvplot(fit,palette = "npg",
           risk.table = TRUE, pval = TRUE,
           conf.int = FALSE, xlab="Time in Days",
           ggtheme = theme_bw(),
           title = 'OS: CN2_ImmMix_CD8T-specific CD4T enrichment',
           surv.median.line = "hv",  # add the median survival pointer.
           # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
           tables.y.text = T,
           risk.table.pos = "in",
           risk.table.col = "strata",
           fontsize = 3.5,
           pval.size = 4,
           surv.plot.height = 0.8,
           tables.height = 0.2,
           pval.coord = c(2000, 0.9),
           legend = "right"
)


#### overall CD4T 
mydata_pheno_anno_TP <- readRDS(file='Class_Pheno_Anno.Rds')

sur.cut <- surv_cutpoint(mydata_pheno_anno_TP, time= 'RFSday',event = 'RFS01' , variables = 'celltype_CD4T' , minprop = 0.1)
sur.cat <- surv_categorize(sur.cut)
names(sur.cat) <- c("RFSday",'RFS01','group')
fit <- survfit(Surv(RFSday, RFS01) ~ group, data = sur.cat)
ggsurvplot(fit,palette = "npg",
           risk.table = TRUE, pval = TRUE,
           conf.int = FALSE, xlab="Time in Days",
           ggtheme = theme_bw(),
           title = 'RFS: Overall CD4T',
           surv.median.line = "hv",  # add the median survival pointer.
           # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
           tables.y.text = T,
           risk.table.pos = "in",
           risk.table.col = "strata",
           fontsize = 3.5,
           pval.size = 4,
           surv.plot.height = 0.8,
           tables.height = 0.2,
           pval.coord = c(2000, 0.9),
           legend = "right"
)

sur.cut <- surv_cutpoint(mydata_pheno_anno_TP, time= 'Osday',event = 'OS01' , variables = 'celltype_CD4T' , minprop = 0.1)
sur.cat <- surv_categorize(sur.cut)
names(sur.cat) <- c("RFSday",'RFS01','group')
fit <- survfit(Surv(RFSday, RFS01) ~ group, data = sur.cat)
ggsurvplot(fit,palette = "npg",
           risk.table = TRUE, pval = TRUE,
           conf.int = FALSE, xlab="Time in Days",
           ggtheme = theme_bw(),
           title = 'OS: Overall CD4T',
           surv.median.line = "hv",  # add the median survival pointer.
           # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
           tables.y.text = T,
           risk.table.pos = "in",
           risk.table.col = "strata",
           fontsize = 3.5,
           pval.size = 4,
           surv.plot.height = 0.8,
           tables.height = 0.2,
           pval.coord = c(2000, 0.9),
           legend = "right"
)













# Fig. S2 cell type constitution groups -----------------------------------------------------------------

library(Seurat) #‘4.1.1’
library(dplyr)  #‘1.0.9’
library(patchwork) # ‘1.1.1’
library(plotly) #‘4.10.0’
library(ggplot2)  #‘3.3.6’
library(clusterProfiler)  #‘4.5.2’
library(org.Hs.eg.db)  #‘3.14.0’
library(enrichplot)  #‘1.14.2’
library(ggpubr)  # ‘0.4.0’
library(tidyverse)  #‘1.3.2’
library(GSVA)  #‘1.42.0’
library(fgsea)   #‘1.20.0’
library(msigdbr)  #‘7.5.1’
library(survminer)  # ‘0.4.9’
library(survival)   #‘3.3.1’

library(ggrepel) #‘0.9.1’
library(corrplot)  # ‘0.92’
library(gmodels)  #‘2.18.1.1’
library(dendextend)  #‘1.16.0’

Class_Pheno_Anno_ctCN <- readRDS("Class_Pheno_Anno_ctCN.Rds")
#celltype_percent = read.csv("~/Desktop/Project/Codex/Analysis_07v/07_GenerateAllSubtypesDF/testing_allReg_seurat_._analysis_20220324_AllSubtypes_withoutartifact_addct.csv")
celltype_percent = read.csv(file="E:\\11. CODEX\\TrainingData\\2022\\testing_allReg_seurat_._analysis_20220324_AllSubtypes_withoutartifact_addct.csv",header = T,row.names = 1,stringsAsFactors = F)

#### celltype: with tumor cells: percentage ####
celltype_percent_WithoutTumor = celltype_percent
cycle = unique(celltype_percent_WithoutTumor$Class)

Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(celltype_percent_WithoutTumor, Class == cycle[i])
  Freq = bind_rows(Freq, round(table(sub$celltype)/nrow(sub),4))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0

## order
hclust = hclust(dist(Freq))

patient_order = rownames(Freq)[hclust$order]

Freq = Freq[patient_order,]

## heatmap for annotation
rownames(Freq) = factor(rownames(Freq), levels = rownames(Freq))

annotation_col = subset(Class_Pheno_Anno_ctCN, select = c("Class", "Age", "Gender", "Pathology", "Differentiation", "TumorNumber", "TumorSizeLagestTumor", "ExtrahepaticMetastasis", "LymphaticMetastasis", "LiverCirrhosis", 
                                                          "Caspsule", "HepatitisB", "HepatitisC", "BileDuctThrombi", "VascularTumorEmboli"))
rownames(annotation_col) = annotation_col$Class
annotation_col = annotation_col[,-1]
annotation_col = annotation_col[patient_order,]

colnames(annotation_col) = c("Age", "Gender", "Pathology", "Differentiation", "Tumor Number", "Tumor Size", "Extrahepatic Metastasis", "Lymphatic Metastasis", "Liver Cirrhosis", 
                             "Caspsule", "Hepatitis B", "Hepatitis C", "Bile Duct Thrombi", "Vascular Tumor Emboli")

## re-classify patient groups
annotation_col$Age = ifelse(annotation_col$Age < 35, "<35", ifelse(annotation_col$Age > 65, ">65", "35~65"))
annotation_col$Gender = ifelse(annotation_col$Gender == 1, "M", "F")
annotation_col$Differentiation = ifelse(annotation_col$Differentiation == 1, "I", 
                                        ifelse(annotation_col$Differentiation == 2, "II", 
                                               ifelse(annotation_col$Differentiation == "2~3", "III", 
                                                      ifelse(annotation_col$Differentiation == 3, "III", 
                                                             ifelse(annotation_col$Differentiation == "3~4", "IV",
                                                                    ifelse(annotation_col$Differentiation == 4, "IV", ""))))))
annotation_col$`Tumor Number` = ifelse(annotation_col$`Tumor Number` == 1, "=1", ">1")
annotation_col$`Tumor Size` = ifelse(annotation_col$`Tumor Size` < 3, "[0,3)", 
                                     ifelse(annotation_col$`Tumor Size` < 5, "[3,5)",
                                            ifelse(annotation_col$`Tumor Size` < 10, "[5,10)", "[10,20]")))
annotation_col$`Extrahepatic Metastasis` = ifelse(annotation_col$`Extrahepatic Metastasis` == 0, "No", "Yes")
annotation_col$`Lymphatic Metastasis` = ifelse(annotation_col$`Lymphatic Metastasis` == 0, "No", "Yes")
annotation_col$`Lymphatic Metastasis` = ifelse(is.na(annotation_col$`Lymphatic Metastasis`), "No", annotation_col$`Lymphatic Metastasis`)
annotation_col$`Liver Cirrhosis` = ifelse(annotation_col$`Liver Cirrhosis` == 0, "No", "Yes")
annotation_col$Caspsule = ifelse(annotation_col$Caspsule == 1, "Yes", "No")
annotation_col$Caspsule = ifelse(is.na(annotation_col$Caspsule), "No", annotation_col$Caspsule)
annotation_col$`Hepatitis B` = ifelse(annotation_col$`Hepatitis B` == 0, "No", "Yes")
annotation_col$`Hepatitis C` = ifelse(annotation_col$`Hepatitis C` == 0, "No", "Yes")
annotation_col$`Bile Duct Thrombi` = ifelse(annotation_col$`Bile Duct Thrombi` == 0, "No", "Yes")
annotation_col$`Vascular Tumor Emboli` = ifelse(annotation_col$`Vascular Tumor Emboli` == 0, "No", "Yes")

annotation_col$Age = factor(annotation_col$Age, levels = c("<35", "35~65", ">65"))
annotation_col$Gender = factor(annotation_col$Gender, levels = c("M","F"))
annotation_col$`Tumor Size` = factor(annotation_col$`Tumor Size`, levels = c("[0,3)", "[3,5)", "[5,10)", "[10,20]"))

#scales::show_col(ggsci::pal_aaas("default", alpha = 0.5)(9))
mycolor = RColorBrewer::brewer.pal(12, "Paired")
scales::show_col(mycolor)

annotation_col_ha = ComplexHeatmap::HeatmapAnnotation(Age = annotation_col$Age,
                                                      Gender = annotation_col$Gender,
                                                      Pathology = annotation_col$Pathology,
                                                      Differentiation = annotation_col$Differentiation,
                                                      `Tumor Number` = annotation_col$`Tumor Number`,
                                                      `Tumor Size` = annotation_col$`Tumor Size`,
                                                      `Extrahepatic Metastasis` = annotation_col$`Extrahepatic Metastasis`,
                                                      `Lymphatic Metastasis` = annotation_col$`Lymphatic Metastasis`,
                                                      `Liver Cirrhosis` = annotation_col$`Liver Cirrhosis`,
                                                      Caspsule = annotation_col$Caspsule,
                                                      `Hepatitis B` = annotation_col$`Hepatitis B`,
                                                      `Hepatitis C` = annotation_col$`Hepatitis C`,
                                                      `Bile Duct Thrombi` = annotation_col$`Bile Duct Thrombi`,
                                                      `Vascular Tumor Emboli` = annotation_col$`Vascular Tumor Emboli`,
                                                      col = list(Age = c("<35" = "#e0e0e0",
                                                                         "35~65" = "#878787",
                                                                         ">65" = "#1a1a1a"),
                                                                 Gender = c("M" = "#878787", "F" = "#1a1a1a"),
                                                                 Pathology = c("1" = "#878787",
                                                                               "1,7" = mycolor[3],
                                                                               "2" = mycolor[2],
                                                                               "3" = mycolor[5],
                                                                               "4" = mycolor[6],
                                                                               "5" = mycolor[10],
                                                                               "6" = mycolor[11],
                                                                               "7" = mycolor[4]),
                                                                 Differentiation = c("I" = mycolor[5],
                                                                                     "II" = mycolor[4],
                                                                                     "III" = "#878787",
                                                                                     "IV" = mycolor[10]),
                                                                 `Tumor Number` = c("=1" = "#878787",
                                                                                    ">1" = "#1a1a1a"),
                                                                 `Tumor Size` = c("[0,3)" = mycolor[3],
                                                                                  "[3,5)" = mycolor[9],
                                                                                  "[5,10)" = "#878787",
                                                                                  "[10,20]" = mycolor[10]),
                                                                 `Extrahepatic Metastasis` = c("No" = "#878787", "Yes" = "#1a1a1a"),
                                                                 `Lymphatic Metastasis` = c("No" = "#878787", "Yes" = "#1a1a1a"),
                                                                 `Liver Cirrhosis` = c("No" = "#878787", "Yes" = "#1a1a1a"),
                                                                 Caspsule = c("No" = "#878787", "Yes" = "#1a1a1a"),
                                                                 `Hepatitis B` = c("No" = "#878787", "Yes" = "#1a1a1a"),
                                                                 `Hepatitis C` = c("No" = "#878787", "Yes" = "#1a1a1a"),
                                                                 `Bile Duct Thrombi` = c("No" = "#878787", "Yes" = "#1a1a1a"),
                                                                 `Vascular Tumor Emboli` = c("No" = "#878787", "Yes" = "#1a1a1a")
                                                      ))
ComplexHeatmap::Heatmap(t(as.matrix(Freq)), bottom_annotation = annotation_col_ha,
                        cluster_columns = F,
                        cluster_rows = F,
                        show_column_names = F,
                        width = ncol(t(as.matrix(Freq)))*unit(0.8, "mm")) # 16x12

## barplot
Freq_melt = Freq
Freq_melt$Patient = rownames(Freq_melt)
# Freq_melt = Freq_melt[order(Freq_melt$`Tumor cells`, decreasing = T),]
Freq_melt$Patient = factor(Freq_melt$Patient, levels = Freq_melt$Patient)
Freq_melt <- reshape2::melt(Freq_melt, id.vars=c("Patient"), variable.name="Celltype", value.name="Count")
Freq_melt = Freq_melt[order(Freq_melt$Patient),]

ggplot(data = Freq_melt, aes(x = Patient, y = Count, fill = Celltype)) + 
  geom_bar(stat = "identity", position = "stack", width = 1) + 
  theme(panel.grid.major.y = element_line(colour = "black"),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_x_discrete(labels = NULL) +
  #theme(axis.text.x = element_text(angle=40, hjust=1, vjust=1)) + 
  scale_fill_manual(values=c("B cells" = "#8358CF",
                             "Biliary tract cells" = "#AE2519",
                             "CD4+ T cells" = "#D083FF",
                             "CD8+ T cells" = "#27A7FA",
                             "Endothelial cells" = "#F6D026",
                             "Fibroblasts" = "#74F882",
                             "Lymphatic endothelial cells" = "#EE6600",
                             "Macrophages" = "#62FFFF",
                             "Tumor cells" = "#FA339A"))

#### relevance analysis ####
hclust = color_branches(hclust,h=0.25,groupLabels = F,col = ggsci::pal_npg("nrc")(9))
plot(hclust) # 14x2.5 inch

temp = as.data.frame(cutree(hclust, h = 0.25, order_clusters_as_data = F))
colnames(temp) = "hclust"

annotation_col$patientid = rownames(annotation_col)
temp$patientid = rownames(temp)

annotation_col = left_join(annotation_col, temp, by = "patientid")

annotation_col$hclust[annotation_col$hclust %in% 1:8] = '1'
annotation_col$hclust[annotation_col$hclust %in% 9] = '2'
annotation_col$hclust[annotation_col$hclust %in% 10] = '3'
annotation_col$hclust[annotation_col$hclust %in% 11:19] = '4'

CrossTable(annotation_col$hclust, annotation_col$`Vascular Tumor Emboli`, fisher=T, chisq =T, format="SPSS",
           prop.c = T,prop.t = T,prop.chisq = T,expected=T)

#### survival ####
plot(hclust)
abline(h = 0.25, col = "red")
group_surv = data.frame(cutree(hclust, h = 0.25,order_clusters_as_data=F))
group_surv$Class = rownames(group_surv)
group_surv_os = subset(Class_Pheno_Anno_ctCN, select = c("Class","RFSday","RFS01","Osday","OS01"))
group_surv = left_join(group_surv, group_surv_os, by = "Class")
colnames(group_surv)[1] = "group"

group_surv$group[group_surv$group %in% 1:8] = '1'
group_surv$group[group_surv$group == 9] = '2'
group_surv$group[group_surv$group == 10] = '3'
group_surv$group[group_surv$group %in% 11:19] = '4'

# res <- pairwise_survdiff(Surv(RFSday, RFS01) ~ group,
#                          data = group_surv, p.adjust.method = "none")
# View(res$p.value)
# 
# res <- pairwise_survdiff(Surv(Osday, OS01) ~ group,
#                          data = group_surv, p.adjust.method = "none")
# View(res$p.value)

fit<-survfit(Surv(RFSday, RFS01)~group, data=group_surv)
# filename <- "./Figures/F3_SP_survival.png"
# png(file=filename,width=2000,height=1750,res=300)
# p <- 
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = group_surv,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #conf.int = TRUE,         # show confidence intervals for 
  palette = "npg",
  xlab = "Time in days",   # customize X axis label.
  ggtheme = theme_classic(), 
  #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
  surv.median.line = "hv"  # add the median survival pointer.
  # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
)
# print(p)
# dev.off()

fit<-survfit(Surv(Osday, OS01)~group, data=group_surv)

ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = group_surv,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #conf.int = TRUE,         # show confidence intervals for 
  palette = "npg",
  xlab = "Time in days",   # customize X axis label.
  ggtheme = theme_classic(), 
  #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
  surv.median.line = "hv"  # add the median survival pointer.
  # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
)

res <- pairwise_survdiff(Surv(Osday, OS01) ~ group,
                         data = group_surv, p.adjust.method = "none")
View(res$p.value)

color = c('#fcad03','#8f876f', '#a4e0b0','#ab66b0')

temp = group_surv
temp$group[!(temp$group %in% 4)] = 'others'
fit<-survfit(Surv(Osday, OS01)~group, data=temp)

p = ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = temp,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #conf.int = TRUE,         # show confidence intervals for 
  palette = c(color[4],"grey"),
  xlab = "Time in days",   # customize X axis label.
  ggtheme = theme_classic(), 
  #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
  surv.median.line = "hv"  # add the median survival pointer.
  # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
)

pdf(file = "Surv4.pdf", width =5, height = 5)
print(p)
dev.off()

save(list = ls(),file = "cell_type_constitution_groups.Rdata")
