###R 4.1.0
rm(list=ls())
gc()

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


### Fig.S7A: CP-based-CN in validation dataset#####
library(tidyverse)
df <- read_csv(file = 'E:\\myPythonProject\\ValidationData\\20220412 AllsubtypeCN\\cells_r=50_CN=40.csv')
gc_csd_CN <- df %>% select(117, 56:116)

gc_csd_CN_anno <- gc_csd_CN %>% mutate(All_CN=paste0('CN',neighborhood10))

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN1'] <- 'CN1_Endo_Tumor'  ###?  Tumor Immmix
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN2'] <- 'CN2_others_LII'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN3'] <- 'CN3_Endo'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN4'] <- 'CN4_p_S6'
# gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN5'] <- ''

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN6'] <- 'CN6_others'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN7'] <- 'CN7_c_Myc_LII'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN8'] <- 'CN8_Biliary_ImmMix'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN9'] <- 'CN9_ImmMix_CD4T'   ###?
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN10'] <- 'CN10_CD107a'

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN11'] <- 'CN11_E_Cad_LII'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN12'] <- 'CN12_Hepar1_LII'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN13'] <- 'CN13_Fibro_HII'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN14'] <- 'CN14_Glypican3'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN15'] <- 'CN15_ImmMix_Lym'

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN16'] <- 'CN16_Twist1'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN17'] <- 'CN17_ImmMix_B'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN18'] <- 'CN18_CD44'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN19'] <- 'CN19_Hepar1_LII'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN20'] <- 'CN20_E_Cad'

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN21'] <- 'CN21_p_mTOR'
# gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN22'] <- ''
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN23'] <- 'CN23_Fibro_HII'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN24'] <- 'CN24_Hepar1' #
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN25'] <- 'CN25_p53'

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN26'] <- 'CN26_Biliary'
# gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN27'] <- ''
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN28'] <- 'CN28_Hepar1_HII'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN29'] <- 'CN29_p_S6'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN30'] <- 'CN30_ImmMix_Mac'

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN31'] <- 'CN31_CD107a'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN32'] <- 'CN32_Ki67'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN33'] <- 'CN33_Ki67_LII'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN34'] <- 'CN34_others'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN35'] <- 'CN35_c_Myc'

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN36'] <- 'CN36_c_Myc_LII'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN37'] <- 'CN37_Fibro'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN38'] <- 'CN38_Hif1a'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN39'] <- 'CN39_others_LII'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN0'] <- 'CN0_HLA_DR'

head(gc_csd_CN_anno)
gc_csd_CN_anno <- gc_csd_CN_anno %>%  select(All_CN, everything()) 
gc_csd_CN_anno <- gc_csd_CN_anno[-2]

df_cn <- aggregate(. ~ All_CN, gc_csd_CN_anno, mean)
df_cn_scale <- df_cn[-1] %>% scale()
rownames(df_cn_scale) <- df_cn$All_CN



library(ComplexHeatmap)
library(circlize)
Heatmap(df_cn_scale, col = colorRamp2(c(-2, 2, 6), c("#35978f", "white", "#bf812d")))


####barplot
All_CN_TP = read_csv(file = 'E:\\myPythonProject\\ValidationData\\20220412 AllsubtypeCN\\cells_r=50_CN=40.csv')
Freq = subset(All_CN_TP, select = c("celltype","neighborhood10"))
colnames(Freq)[2] = "All_CN"
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
library(reshape)
Freq_toplot = melt(Freq_toplot, id.vars = "CN")
Freq_toplot = Freq_toplot[order(Freq_toplot$CN),]
Freq_toplot$CN = as.numeric(Freq_toplot$CN)
Freq_toplot = Freq_toplot[order(Freq_toplot$CN),]
# Freq_toplot$CN_r = factor(Freq_toplot$CN_r, levels = rev(c(29,20,28,37,19,
#                                                           15,39,21,14,27,
#                                                           7,18,26,9,12,
#                                                           3,33,4,36,30,
#                                                           35,8,31,22,24,
#                                                           17,32,11,23,13,
#                                                           6,25,38,0,16,
#                                                           2,34,10,1,5)))

Freq_toplot$CN = factor(Freq_toplot$CN, levels = rev(c(30,38,9,17,23,
                                                       3,6,0,22,18,
                                                       15,13,37,4,29,
                                                       25,32,33,31,10,
                                                       28,8,34,5,1,
                                                       35,16,27,14,21,
                                                       20,11,24,12,39,
                                                       2,36,7,26,19)))

Freq_toplot <- Freq_toplot[-3]
names(Freq_toplot)[3] <- 'value'

ggplot(Freq_toplot, aes(x = CN, y = value)) + 
  geom_bar(aes(fill = factor(variable)), stat = "identity", position = "fill") + 
  scale_fill_manual(values = c(pal_npg('nrc',alpha = 0.6)(12))) +
  scale_y_continuous(labels = scales::percent) + #纵坐标变为百分比 
  coord_flip()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())





# Fig.S7B CN comparison in discovery and validation cohort ------------------------


gc_csd_CN_val <- read_csv(file = 'E:\\myPythonProject\\ValidationData\\20220412 AllsubtypeCN\\cells_r=50_CN=40.csv')
gc_csd_CN_anno <- gc_csd_CN_val %>% mutate(All_CN=paste0('CN',neighborhood10))

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN1'] <- 'CN1_Endo_Tumor'  ###?  Tumor Immmix
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN2'] <- 'CN2_others_LII' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN3'] <- 'CN3_Endo' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN4'] <- 'CN4_p_S6' 
# gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN5'] <- ''

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN6'] <- 'CN6_others' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN7'] <- 'CN7_c_Myc_LII' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN8'] <- 'CN8_Biliary_ImmMix' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN9'] <- 'CN9_ImmMix_CD4T'   ###?
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN10'] <- 'CN10_CD107a' 

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN11'] <- 'CN11_E_Cad_LII' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN12'] <- 'CN12_Hepar1_LII' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN13'] <- 'CN13_Fibro_HII' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN14'] <- 'CN14_Glypican3' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN15'] <- 'CN15_ImmMix_Lym' 

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN16'] <- 'CN16_Twist1' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN17'] <- 'CN17_ImmMix_B' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN18'] <- 'CN18_CD44' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN19'] <- 'CN19_Hepar1_LII' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN20'] <- 'CN20_E_Cad' 

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN21'] <- 'CN21_p_mTOR' 
# gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN22'] <- ''
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN23'] <- 'CN23_Fibro_HII' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN24'] <- 'CN24_Hepar1' #
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN25'] <- 'CN25_p53' 

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN26'] <- 'CN26_Biliary' 
# gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN27'] <- ''
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN28'] <- 'CN28_Hepar1_HII' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN29'] <- 'CN29_p_S6' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN30'] <- 'CN30_ImmMix_Mac' 

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN31'] <- 'CN31_CD107a' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN32'] <- 'CN32_Ki67' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN33'] <- 'CN33_Ki67_LII' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN34'] <- 'CN34_others' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN35'] <- 'CN35_c_Myc' 

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN36'] <- 'CN36_c_Myc_LII' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN37'] <- 'CN37_Fibro' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN38'] <- 'CN38_Hif1a' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN39'] <- 'CN39_others_LII' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN0'] <- 'CN0_HLA_DR' 

#saveRDS(gc_csd_CN_anno, file='cells_All_CN_annotaion.RDS')


##Training data with CN annotaion
gc_csd_CN_Tra <- read_csv('./TrainingData AllCN_r50_CN40_0331/All_CN_TP_0424.csv')
#gc_csd_CN_anno <- readRDS('cells_All_CN_annotaion.RDS')


pre_CN_subtype <- function(x,nx='CN',cn=neighborhood10,ct=Allsubtypes){
  x1 <- x %>% dplyr::select(CN={{cn}},{{ct}}) %>% dplyr::group_by(CN,{{ct}}) %>% 
    dplyr::summarise(count=n()) %>% tidyr::spread({{ct}},count)
  x2 <- as.data.frame(x1)
  row.names(x2) <- paste0(nx,x2$CN)
  df_CN_subtype <- x2[,-1]
  return(df_CN_subtype)
}

df_CN_subtype_val <- pre_CN_subtype(gc_csd_CN_anno,nx="Val_",cn=All_CN,ct=Allsubtypes)
df_CN_subtype_Tra <- pre_CN_subtype (gc_csd_CN_Tra,nx="Tra_",cn=All_CN,ct=Allsubtypes)

###subset common subtype
common_va <- intersect(colnames(df_CN_subtype_val),colnames(df_CN_subtype_Tra))

df_CN_subtype_val[is.na(df_CN_subtype_val)] <- 0
df_CN_subtype_Tra[is.na(df_CN_subtype_Tra)] <- 0


val_scale <- scale(df_CN_subtype_val)
tra_scale <- scale(df_CN_subtype_Tra)

val <- val_scale[,common_va]
tra <- tra_scale[,common_va]


df_CN_subtype <- rbind(val,tra)

bk <- c(seq(-5,-0.1,by=0.01),seq(0,5,by=0.01))
dev.off()
pdf(file='TrainingCN vs ValidationCN with annotaion.pdf',width = 20,height = 25)
pheatmap(df_CN_subtype,
         cellwidth = 15,
         cellheight = 15,
         border=FALSE,
         color=c(colorRampPalette(colors=c("#35978f", "white"))(length(bk)/2),colorRampPalette(color=c("white", "#bf812d"))(length(bk)/2)),
         scale = "row",
         cluster_rows = T,
         show_rownames = T,
         cluster_cols = T,
         legend_breaks=seq(-5,5,1),
         breaks=bk,
         cutree_rows=25)
dev.off()


# Fig. S8A TP in validation  -----------------------------------------------------------
###TMA2_8_reg007, TMA7_10_reg013 were removed from downstream analysis since the majority  was pan-ck+ cells.

gc_csd_CN_anno_filter <-gc_csd_CN_anno %>% 
  dplyr::filter(!(Class=='TMA2_8_reg007'|Class=='TMA7_10_reg013'))

cells_CN <- gc_csd_CN_anno_filter %>% dplyr::select(Class, All_CN)


dat <- as.data.frame(with(cells_CN, table(Class, All_CN)))
cells_CN_Freq <- spread(dat, All_CN, Freq)

dat_CN_freq <- cells_CN_Freq[,-1]
row.names(dat_CN_freq) <- cells_CN_Freq$Class


dat_CN_percent <- dat_CN_freq / rowSums(dat_CN_freq) * 100
row.names(dat_CN_percent)<- cells_CN_Freq$Class

save(gc_csd_CN_anno_filter, dat_CN_freq, dat_CN_percent, file='data for TPSP definition.Rdata')


load('data for TPSP definition.Rdata')

colnames(dat_CN_percent)
dat_CN_percent_TuCN <- dplyr::select(dat_CN_percent,
                                     CN10_CD107a,
                                     CN11_E_Cad_LII,
                                     CN12_Hepar1_LII,
                                     CN14_Glypican3,
                                     CN16_Twist1,
                                     # CN18_CD44+
                                     CN2_others_LII,
                                     CN20_E_Cad,
                                     CN21_p_mTOR,
                                     CN24_Hepar1,
                                     CN25_p53,
                                     # CN26_Biliary+
                                     CN29_p_S6,
                                     CN31_CD107a,
                                     # CN32_Ki67+
                                     # CN33_Ki67_LII+
                                     CN35_c_Myc,
                                     CN36_c_Myc_LII,
                                     # CN37_Fibro+
                                     CN4_p_S6,
                                     CN7_c_Myc_LII,
                                     CN34_others,
                                     CN39_others_LII,
                                     CN28_Hepar1_HII,
                                     # CN38_Hif1a,
                                     # CN18_CD44,
)


bk <- c(seq(0,1.9,by=0.01),seq(2,10,by=0.01))
p <- pheatmap(dat_CN_percent_TuCN,
              cellwidth = 10,
              cellheight = 2,
              # annotation_row = df_TP,
              border=FALSE,
              color=c(colorRampPalette(colors=c("blue","white"))(length(bk)/5),colorRampPalette(color=c("white","red"))(length(bk)*4/5)),
              scale = 'row',
              cluster_rows = T,
              show_rownames = F,
              show_colnames = T,
              legend_breaks=seq(0,10,2),
              breaks=bk,
              cutree_rows = 10)



row_cluster <- as.data.frame(cutree(p$tree_row,k=10))
table(row_cluster)
names(row_cluster) <- 'PhenoGraph'
row_cluster$Class <-  row.names(row_cluster)

dat_CN_percent_class <- dat_CN_percent_TuCN
dat_CN_percent_class$Class <- row.names(dat_CN_percent_TuCN)

df_merge <- merge(row_cluster,dat_CN_percent_class,by='Class')


df_merge <- df_merge[order(df_merge$PhenoGraph),]
temp <- df_merge[,c(-1,-2)]
row.names(temp) <- df_merge$Class 

anno <- df_merge[2]
row.names(anno) <- df_merge$Class
anno$PhenoGraph <- as.factor(anno$PhenoGraph)

p1 <- pheatmap(temp,
               cellwidth = 10,
               cellheight = 2,
               annotation_row = anno,
               border=FALSE,
               color=c(colorRampPalette(colors=c("blue","white"))(length(bk)/5),colorRampPalette(color=c("white","red"))(length(bk)*4/5)),
               scale = 'row',
               cluster_rows = T,
               show_rownames = F,
               show_colnames = T,
               legend_breaks=seq(0,10,2),
               breaks=bk,
               cutree_rows = 10)

# pdf(file = "7_pheatmap_AllCN.pdf", width =10, height = 15)
# print(p1)
# dev.off()

row_cluster$Cluster <- row_cluster$PhenoGraph

row_cluster$Cluster[row_cluster$Cluster== 1] <- 'TP_p_mTOR'
row_cluster$Cluster[row_cluster$Cluster== 2] <- 'TP_others'
row_cluster$Cluster[row_cluster$Cluster== 3] <- 'TP_Hepar1'
row_cluster$Cluster[row_cluster$Cluster== 4] <- 'TP_c_Myc'
row_cluster$Cluster[row_cluster$Cluster== 5] <- 'TP_CD107a'
row_cluster$Cluster[row_cluster$Cluster== 6] <- 'TP_E_cadherin'
row_cluster$Cluster[row_cluster$Cluster== 7] <- 'TP_Glypican3'
row_cluster$Cluster[row_cluster$Cluster== 8] <- 'TP_p53'
row_cluster$Cluster[row_cluster$Cluster== 9] <- 'TP_p_S6'
row_cluster$Cluster[row_cluster$Cluster== 10] <- 'TP_Twist1'

# row_cluster$Cluster[row_cluster$Cluster== 13] <- 'TP13'
# row_cluster$Cluster[row_cluster$Cluster== 14] <- 'TP_ImmuneMix'
# row_cluster$Cluster[row_cluster$Cluster== 15] <- 'TP_p53'
# row_cluster$Cluster[row_cluster$Cluster== 16] <- 'TP_Twist1'
# row_cluster$Cluster[row_cluster$Cluster== 17] <- 'TP_Hepar1'



df_merge <- merge(row_cluster,dat_CN_percent_class,by='Class')

write.csv(df_merge, file='class_TP.csv')


##
temp <- df_merge[,c(-1,-2,-3)]
row.names(temp) <- df_merge$Class 

anno <- df_merge[3]
row.names(anno) <- df_merge$Class


bk <- c(seq(0,1.9,by=0.01),seq(2,10,by=0.01))
ann_colors = list(
  Cluster = c("TP_c_Myc" = pal_npg('nrc',alpha = 0.6)(10)[1],
              "TP_CD107a" = pal_npg('nrc',alpha = 0.6)(10)[2],
              "TP_E_cadherin" = pal_npg('nrc',alpha = 0.6)(10)[3],
              "TP_Glypican3" = pal_npg('nrc',alpha = 0.6)(10)[4],
              "TP_Hepar1" = pal_npg('nrc',alpha = 0.6)(10)[5],
              "TP_others" = pal_npg('nrc',alpha = 0.6)(10)[6],
              "TP_p_mTOR" = pal_npg('nrc',alpha = 0.6)(10)[7],
              "TP_p_S6" = pal_npg('nrc',alpha = 0.6)(10)[8],
              "TP_p53" = pal_npg('nrc',alpha = 0.6)(10)[9],
              "TP_Twist1" = pal_npg('nrc',alpha = 0.6)(10)[10])
)
pheatmap(temp,
         cellwidth = 15,
         cellheight = 3,
         annotation_row = anno,
         border=FALSE,
         color=c(colorRampPalette(colors=c("#35978f", "white"))(length(bk)/5),colorRampPalette(color=c("white", "#bf812d"))(length(bk)*4/5)),
         scale = 'row',
         cluster_rows = T,
         show_rownames = F,
         show_colnames = T,
         legend_breaks=seq(0,10,2),
         breaks=bk,
         cutree_rows = 10,
         annotation_colors = ann_colors)




# Fig. S8B TP-p53 survival in validation ----------------------------------

mydata_pheno <- df_merge
row.names(mydata_pheno) <- df_merge$Class
mydata_pheno$Class0 <- rownames(mydata_pheno)
mydata_pheno$Class0 <- as.character(mydata_pheno$Class0)

load("HCC401_new.Rdata")
HCCpathology <- df_validation_new
HCCpathology$Class <- paste0(HCCpathology$TMA, "_reg",sprintf("%03d",HCCpathology$Reg))

for (i in 1:length(mydata_pheno$Class)) {
  temp<-strsplit(mydata_pheno$Class,'_')[[i]][1]
  temp1 <-strsplit(temp,"A")[[1]][2]
  temp2 <- strsplit(mydata_pheno$Class,'_')[[i]][3]
  mydata_pheno$Class[i]<- paste0(temp1,"_",temp2)
}

mydata_pheno_anno <- merge(mydata_pheno, HCCpathology, by="Class")


# for (i in seq_along(cy)) {
a <- mydata_pheno_anno
# a$Cluster[a$Cluster != cy[i]] <- 'other clusters'
a$Cluster[a$Cluster != 'TP_p53'] <- 'other clusters'
fit<-survfit(Surv(RFSday, RFS01)~Cluster, data=a)
# filename <- paste0('./TP survival/TP vs others/',cy[i],"_vs_other clusters.png")
# png(file=filename,width=1500,height=1750,res=300)
p <- ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = a,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #conf.int = TRUE,         # show confidence intervals for 
  palette = "npg",
  xlab = "Time in days",   # customize X axis label.
  ggtheme = theme_classic(), 
  #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = c('TP_p53','other clusters'),    # change legend labels.
  title = 'RFS',
  tables.y.text = T,
  risk.table.pos = "in",
  risk.table.col = "strata",
  fontsize = 3.5,
  pval.size = 4,
  surv.plot.height = 0.7,
  tables.height = 0.3,
  pval.coord = c(2000, 0.9),
  legend = "right"
)
print(p)
#   dev.off()
#   print(paste0(cy[i]))
# }














