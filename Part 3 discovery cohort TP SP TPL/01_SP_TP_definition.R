

##### R 3.6.3


rm(list=ls())
gc()


library(tidyverse) #v1.3.0
library(tidyr)  #v1.1.4
library(Rphenograph)  #v0.99.1
library(ggpubr)   #v0.4.0
library(ggthemes) #v4.2.0
library(Rtsne)  #v0.15
library(flowCore)   #v1.52.1
library(Rcpp)  #v1.0.5
library(cytofkit)   #v0.99.0
library(igraph) #v1.2.6
library(cytofexplorer)  #v0.1.0
 library(survival)   #v3.2.7
library(survminer)  #v0.4.8
library(pheatmap)  #v1.0.12
 library(ggsci)  #v2.9

# setwd("E:/11. CODEX/TrainingData/20220308/AllCN_r50_CN40_0331")

# load CP-based-CN40  ---------------------------------------------------------------

gc_csd_CN <- read_csv(file = 'E:\\myPythonProject\\TrainingData\\20220317 AllsubtypeCN\\cells_r=50_CN=40.csv')


# CN annotation -----------------------------------------------------------

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



# CN percentages in each region/class -----------------------------------------------------------

cells_CN <- gc_csd_CN_anno %>% dplyr::select(Class, All_CN)


dat <- as.data.frame(with(cells_CN, table(Class, All_CN)))
cells_CN_Freq <- spread(dat, All_CN, Freq)

dat_CN_freq <- cells_CN_Freq[,-1]
row.names(dat_CN_freq) <- cells_CN_Freq$Class


dat_CN_percent <- dat_CN_freq / rowSums(dat_CN_freq) * 100
row.names(dat_CN_percent)<- cells_CN_Freq$Class





# TP definition -----------------------------------------------------

####
bk <- c(seq(0,1.9,by=0.01),seq(2,10,by=0.01))
p <- pheatmap(dat_CN_percent,
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
         cutree_rows = 17)

##

row_cluster <- as.data.frame(cutree(p$tree_row,k=17))
table(row_cluster)
names(row_cluster) <- 'PhenoGraph'
row_cluster$Class <-  row.names(row_cluster)

dat_CN_percent_class <- dat_CN_percent
dat_CN_percent_class$Class <- row.names(dat_CN_percent)

df_merge <- merge(row_cluster,dat_CN_percent_class,by='Class')

##
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
               cluster_rows = F,
               show_rownames = F,
               show_colnames = T,
               legend_breaks=seq(0,10,2),
               breaks=bk,
               cutree_rows = 17)

pdf(file = "7_pheatmap_AllCN.pdf", width =10, height = 15)
print(p1)
dev.off()



####

row_cluster$Cluster <- row_cluster$PhenoGraph

row_cluster$Cluster[row_cluster$Cluster== 1] <- 'TP_Others'
row_cluster$Cluster[row_cluster$Cluster== 2] <- 'TP_p_mTOR'
row_cluster$Cluster[row_cluster$Cluster== 3] <- 'TP_Hepar1'
row_cluster$Cluster[row_cluster$Cluster== 4] <- 'TP_Ki67'
row_cluster$Cluster[row_cluster$Cluster== 5] <- 'TP_CD107a'
row_cluster$Cluster[row_cluster$Cluster== 6] <- 'TP_p_S6'
row_cluster$Cluster[row_cluster$Cluster== 7] <- 'TP_E_cadherin'
row_cluster$Cluster[row_cluster$Cluster== 8] <- 'TP_Glypican3'
row_cluster$Cluster[row_cluster$Cluster== 9] <- 'TP_c_Myc'
row_cluster$Cluster[row_cluster$Cluster== 10] <- 'TP_CN22'
row_cluster$Cluster[row_cluster$Cluster== 11] <- 'TP_ImmuneMix'
row_cluster$Cluster[row_cluster$Cluster== 12] <- 'TP_ImmuneMix'
row_cluster$Cluster[row_cluster$Cluster== 13] <- 'TP13'
row_cluster$Cluster[row_cluster$Cluster== 14] <- 'TP_ImmuneMix'
row_cluster$Cluster[row_cluster$Cluster== 15] <- 'TP_p53'
row_cluster$Cluster[row_cluster$Cluster== 16] <- 'TP_Twist1'
row_cluster$Cluster[row_cluster$Cluster== 17] <- 'TP_Hepar1'

dat_CN_percent_class <- dat_CN_percent
dat_CN_percent_class$Class <- row.names(dat_CN_percent)

df_merge <- merge(row_cluster,dat_CN_percent_class,by='Class')

##
temp <- df_merge[,c(-1,-2,-3)]
row.names(temp) <- df_merge$Class 

anno <- df_merge[3]
row.names(anno) <- df_merge$Class

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
               cutree_rows = 17)

pdf(file = "8_anno_pheatmap_AllCN.pdf", width =10, height = 15)
print(p1)
dev.off()


# TP combined with clinical information

mydata_pheno <- merge(row_cluster,dat_CN_percent_class,by='Class')
row.names(mydata_pheno) <- dat_CN_percent_class$Class
mydata_pheno$Class0 <- rownames(mydata_pheno)
mydata_pheno$Class0 <- as.character(mydata_pheno$Class0)

load("../HCC401_new.Rdata")
HCCpathology <- df_training_new
HCCpathology$Class <- paste0(HCCpathology$TMA, "_reg",sprintf("%03d",HCCpathology$Reg))

for (i in 1:length(mydata_pheno$Class)) {
  temp<-strsplit(mydata_pheno$Class,'_')[[i]][1]
  temp1 <-strsplit(temp,"A")[[1]][2]
  temp2 <- strsplit(mydata_pheno$Class,'_')[[i]][3]
  mydata_pheno$Class[i]<- paste0(temp1,"_",temp2)
}

mydata_pheno_anno <- merge(mydata_pheno, HCCpathology, by="Class")

write.csv(mydata_pheno_anno, file = 'mydata_pheno_anno.csv')



############SP definition ######
####

colnames(dat_CN_freq)
dat_selCN_freq <- dat_CN_freq %>% mutate(CN_tumor=(CN0_Hepar1_LII +
                                                   CN10_others +
                                                   CN11_p_S6 +
                                                   CN12_Glypican3 +
                                                   CN13_Hepar1_LII +
                                                   CN16_Hepar1 +
                                                   CN17_p_mTOR +
                                                   CN25_E_Cad +
                                                   CN3_c_Myc +
                                                   CN31_CD107a +
                                                   CN32_p_mTOR +
                                                   CN33_c_Myc_LII +
                                                   CN34_others_LII +
                                                   CN36_p53 +
                                                   CN38_Hepar1_HII +
                                                   CN4_p53 +
                                                   CN5_others +
                                                   CN6_E_Cad_LII +
                                                   CN8_CD107a +
                                                   CN2_others_LII +
                                                   CN23_p_S6 +
                                                   CN24_Twist1 +
                                                   CN22 +
                                                   CN26 +
                                                   CN37)) %>%  select(CN_tumor,
                                                                                  starts_with('CN18_'),
                                                                                   starts_with('CN7_'),
                                                                                   starts_with('CN15_'),
                                                                                   starts_with('CN20_'),
                                                                                   starts_with('CN21_'),
                                                                                   starts_with('CN9_'),
                                                                                   # starts_with('CN26_'),
                                                                                   starts_with('CN29_'),
                                                                                   starts_with('CN27_'),
                                                                                   starts_with('CN39_'),
                                                                                   # starts_with('CN37_'),
                                                                                  starts_with('CN19_'),
                                                                                  starts_with('CN28_'),
                                                                                  # starts_with('CN22_'),
                                                                                  starts_with('CN30_'),
                                                                                  starts_with('CN35'),
                                                                                  starts_with('CN1_'),
                                                                                  starts_with('CN14_'))

dat_selCN_percent <- dat_selCN_freq / rowSums(dat_selCN_freq) * 100
dat_selCN_percent_scale <- scale(dat_selCN_percent[,-1])


# 
# #######PhenoGraph clustering

# 
#PG_elbow(dat_selCN_percent_scale, k_from = 20, k_to = 50, k_step=5)
# 

#PhenoGraph parameters：
seed=123 
k <-30 
PhenoGraph_result <-as.numeric(membership(Rphenograph(dat_selCN_percent_scale,k=k)))

#Checkpoint
hist(PhenoGraph_result,unique(PhenoGraph_result))

#Input: tSNE parameters

max_iter=1000   
perplexity=10  
seed=123      
theta=0.5     
dims = 2      


if (exists('seed')) set.seed(seed)
tsne_result <- Rtsne(dat_selCN_percent_scale,
                     initial_dims = ncol(dat_selCN_percent_scale),
                     pca = FALSE,
                     dims = dims,
                     check_duplicates = FALSE,
                     perplexity=perplexity,
                     max_iter=max_iter,
                     theta=theta)$Y
row.names(tsne_result)<-row.names(dat_selCN_percent_scale)
colnames(tsne_result)<-c("tsne_1","tsne_2")
plot(tsne_result)

combined_data_analysed <- cbind(dat_selCN_percent_scale,
                                tsne_result,
                                PhenoGraph = PhenoGraph_result)

#Checkpoint
head(combined_data_analysed) 
combined_data_analysed <- as.data.frame(combined_data_analysed)

combined_data_analysed$PhenoGraph  <-as.factor(combined_data_analysed$PhenoGraph)
centers<-combined_data_analysed %>%
  group_by(PhenoGraph)  %>%
  summarise(tsne_1=median(tsne_1),tsne_2=median(tsne_2))

mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                 legend.key = element_rect(fill = "white", colour = "white"), #图标
                 legend.background = (element_rect(colour= "white", fill = "white")))


ggplot(combined_data_analysed)+
  geom_point(aes(x=tsne_1,y=tsne_2,colour=PhenoGraph))+
  mytheme+
  scale_color_manual(values=c(pal_npg('nrc',alpha = 1.0)(6)))

##tsne
ggplot(combined_data_analysed)+
  geom_point(aes(x=tsne_1,y=tsne_2,colour=PhenoGraph))+
  mytheme+
  scale_color_manual(values=c(pal_npg('nrc',alpha = 1.0)(6)))



##################phenograph heatmap, selected CN
sortdata <- combined_data_analysed[order(combined_data_analysed$PhenoGraph),]
unique(sortdata$PhenoGraph)
# sortdata$PhenoGraph <- as.character(sortdata$PhenoGraph)
annotation_row<-select(sortdata, PhenoGraph)

temp1 <- sortdata[,1:dim(dat_selCN_percent_scale)[2]]

bk <- c(seq(-2,-0.1,by=0.01),seq(0,5,by=0.01))

pheatmap(temp1,
         cellwidth = 10,
         cellheight = 1,
         annotation_row = annotation_row,
         border=FALSE,
         color=c(colorRampPalette(colors=c("blue","white"))(length(bk)/2),colorRampPalette(color=c("white","red"))(length(bk)/2)),
         scale = "none",
         cluster_rows = F,
         show_rownames = F,
         breaks = bk)

##annotation
combined_data_analysed$PhenoGraph <- as.character(combined_data_analysed$PhenoGraph)
combined_data_analysed$PhenoGraph[combined_data_analysed$PhenoGraph=='1'] <- 'SP_HII'
combined_data_analysed$PhenoGraph[combined_data_analysed$PhenoGraph=='2'] <- 'SP_Pf'
combined_data_analysed$PhenoGraph[combined_data_analysed$PhenoGraph=='3'] <- 'SP_MII_Tumor'
combined_data_analysed$PhenoGraph[combined_data_analysed$PhenoGraph=='4'] <- 'SP_MII_Fibro'
combined_data_analysed$PhenoGraph[combined_data_analysed$PhenoGraph=='5'] <- 'SP_LII'
combined_data_analysed$PhenoGraph[combined_data_analysed$PhenoGraph=='6'] <- 'SP_LII'

##tsne
ggplot(combined_data_analysed)+
  geom_point(aes(x=tsne_1,y=tsne_2,colour=PhenoGraph))+
  mytheme+
  scale_color_manual(values=c(pal_npg('nrc',alpha = 1.0)(6)))


##################phenograph heatmap, selected CN
sortdata <- combined_data_analysed[order(combined_data_analysed$PhenoGraph),]
unique(sortdata$PhenoGraph)
# sortdata$PhenoGraph <- as.character(sortdata$PhenoGraph)
annotation_row<-select(sortdata, PhenoGraph)

temp1 <- sortdata[,1:dim(dat_selCN_percent_scale)[2]]

bk <- c(seq(0,1.9,by=0.01),seq(2,4,by=0.01))

pdf(file = "9_pheatmape_phenograph_K=30_5SP_anno.pdf", width =10, height = 15)
pheatmap(temp1,
         cellwidth = 10,
         cellheight = 1,
         annotation_row = annotation_row,
         border=FALSE,
         color=c(colorRampPalette(colors=c("blue","white"))(length(bk)/2),colorRampPalette(color=c("white","red"))(length(bk)/2)),
         scale = "none",
         cluster_rows = F,
         show_rownames = F,
         breaks = bk)
dev.off()


pdf(file = "10_pheatmape_phenograph_K=30_5SP_class.pdf", width =10, height = 45)
pheatmap(temp1,
         cellwidth = 10,
         cellheight = 10,
         annotation_row = annotation_row,
         border=FALSE,
         color=c(colorRampPalette(colors=c("blue","white"))(length(bk)/2),colorRampPalette(color=c("white","red"))(length(bk)/2)),
         scale = "none",
         cluster_rows = F,
         show_rownames = T,
         breaks = bk)
dev.off()


############# SP combined with clinical information 
mydata_pheno <- combined_data_analysed
mydata_pheno$Class <- rownames(combined_data_analysed)
mydata_pheno$Class <- as.character(mydata_pheno$Class)
mydata_pheno$Class0 <- rownames(combined_data_analysed)
mydata_pheno$Class0 <- as.character(mydata_pheno$Class0)

load("../HCC401_new.Rdata")
HCCpathology <- df_training_new
HCCpathology$Class <- paste0(HCCpathology$TMA, "_reg",sprintf("%03d",HCCpathology$Reg))

for (i in 1:length(mydata_pheno$Class)) {
  temp<-strsplit(mydata_pheno$Class,'_')[[i]][1]
  temp1 <-strsplit(temp,"A")[[1]][2]
  temp2 <- strsplit(mydata_pheno$Class,'_')[[i]][3]
  mydata_pheno$Class[i]<- paste0(temp1,"_",temp2)
}

mydata_pheno_anno <- merge(mydata_pheno, HCCpathology, by="Class")
dat_CN_percent$Class <- row.names(dat_CN_percent)
mydata_pheno_anno <-  merge(dat_CN_percent,mydata_pheno_anno[-c(2:10)], by.x = 'Class', by.y = 'Class0')
mydata_pheno_anno <- mydata_pheno_anno %>% dplyr::rename(Class0 = Class.y)
write.csv(mydata_pheno_anno, file = 'selCN_k=30_5SP_anno.csv')


##########export information for each cell#####
##all cells tumor subtype CN = 15
##all cells immune subtype CN = 15

library(tidyverse)
library(patchwork) 
library(reshape) #melt，cast
library(tidyr) #spread，gather
library(stringr)

#TP pattern
tissue_pattern <- read.csv('mydata_pheno_anno.csv',header = T,row.names = 1,stringsAsFactors = F)
#SP pattern
ImmuneInf <- read.csv('selCN_k=30_5SP_anno.csv',header = T,row.names = 1,stringsAsFactors = F)



tissue_pattern <- dplyr::select(tissue_pattern, Class0,Cluster)
colnames(tissue_pattern) <- c('Class','TP')
unique(tissue_pattern$TP)


ImmuneInf <- ImmuneInf %>% select(Class, SP=PhenoGraph)


TP_SP <- merge(tissue_pattern, ImmuneInf,by='Class')


colnames(gc_csd_CN_anno)
allCN <- select(gc_csd_CN_anno,All_CN,3:54,neighborhood10)
colnames(allCN)

allCN <- select(allCN,cellid,Class, All_CN,Allsubtypes,TumorSubtype,StromalSubtype,everything())



All_CN_TP <- merge(TP_SP,allCN, by='Class')
head(All_CN_TP)
colnames(All_CN_TP)
unique(All_CN_TP$TP)


unique(All_CN_TP$StromalSubtype)
All_CN_TP$celltype <- All_CN_TP$StromalSubtype
All_CN_TP$celltype[with(All_CN_TP,grepl("CD4T", celltype))] <-'CD4T'
All_CN_TP$celltype[with(All_CN_TP,grepl("CD8T", celltype))] <-'CD8T'
All_CN_TP$celltype[with(All_CN_TP,grepl("Macrophages", celltype))] <-'Macrophages'

All_CN_TP$celltype[All_CN_TP$StromalSubtype=="CD4T_FOXP3+_KI67+"|All_CN_TP$StromalSubtype=="CD4T_FOXP3+"] <- "Treg"
All_CN_TP$celltype[All_CN_TP$StromalSubtype=="Macrophages_CD11c+_HLA-DR+"|All_CN_TP$StromalSubtype=="Macrophages_HLA-DR+"] <- "APC"

unique(All_CN_TP$celltype)

save(All_CN_TP,file = 'All_CN_TP_0424.Rdata')
write.csv(All_CN_TP, file = 'All_CN_TP_0424.csv')



######################## export information of each region/class###################################

##########All CN
cells_CN <-dplyr::select(All_CN_TP, Class, All_CN)
dat <- as.data.frame(with(cells_CN, table(Class, All_CN)))
dat1 <- spread(dat, All_CN, Freq)
dat2 <- dat1[,-1]
dat_percent <- dat2 / rowSums(dat2) * 100
# names(dat_percent) <- paste0("AllCN",c(0:39))
dat_percent$Class <- dat1$Class
write.csv(dat_percent,file = 'class_allCN.csv')

####All subtype
cells_subtype <-dplyr::select(All_CN_TP, Class, Allsubtypes)
dat <- as.data.frame(with(cells_subtype, table(Class, Allsubtypes)))
dat1 <- spread(dat, Allsubtypes, Freq)
dat2 <- dat1[,-1]
dat_percent <- dat2 / rowSums(dat2) * 100
dat_percent$Class <- dat1$Class
write.csv(dat_percent,file = 'class_allsubtypes.csv')

###All celltype
cells_celltype <-dplyr::select(All_CN_TP, Class, celltype)
dat <- as.data.frame(with(cells_celltype, table(Class, celltype)))
dat1 <- spread(dat, celltype, Freq)
dat2 <- dat1[,-1]
dat_percent <- dat2 / rowSums(dat2) * 100
dat_percent$Class <- dat1$Class
write.csv(dat_percent,file = 'class_celltype.csv')


####merge files
# rm(list=ls())
# gc()
allCN <- read.csv('class_allCN.csv',header = T,row.names = 1,stringsAsFactors = FALSE)
allsubtypes <- read.csv('class_allsubtypes.csv', header = T,row.names = 1,stringsAsFactors = FALSE)

allcelltype <- read.csv('class_celltype.csv', header = T,row.names = 1,stringsAsFactors = FALSE)
names(allcelltype)[1:11] <- paste0('celltype_',names(allcelltype)[1:11])



temp <- merge(allCN,allsubtypes,by='Class')
df_merge <- merge(temp, allcelltype,by='Class')


SP_anno <- read.csv(file='selCN_k=30_5SP_anno.csv',header = T,row.names = 1,stringsAsFactors = F)
SP_sel <- dplyr::select(SP_anno,Class ,SP=PhenoGraph)

TP_anno <- read.csv(file='mydata_pheno_anno.csv',header = T,row.names = 1,stringsAsFactors = F)
TP_sel <- dplyr::select(TP_anno,Class =  Class0, PhenoGraph, TP=Cluster, 45:106)

TP_SP_merge <- merge(SP_sel,TP_sel,by='Class')

mydata_pheno_anno_TP <- merge(df_merge,TP_SP_merge,by='Class')
colnames(mydata_pheno_anno_TP)
dim(mydata_pheno_anno_TP)
write.csv(mydata_pheno_anno_TP,file = 'Class_Pheno_Anno.csv')










