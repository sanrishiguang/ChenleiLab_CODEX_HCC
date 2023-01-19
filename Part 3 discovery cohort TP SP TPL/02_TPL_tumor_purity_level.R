
##### R 3.6.3

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



#load output data from Python
gc_csd_CN <- read_csv(file = 'E:\\myPythonProject\\TrainingData\\20220423 celltype CN\\cells_r=50_CN=10.csv')
colnames(gc_csd_CN)
# CN annotation -----------------------------------------------------------

gc_csd_CN_anno <- gc_csd_CN %>% mutate(All_CN=paste0('CN',neighborhood10))

gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN0'] <- 'CN0_ImmMix_B' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN1'] <- 'CN1_Tumor' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN2'] <- 'CN2_ImmMix_CD8T' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN3'] <- 'CN3_Tumor_ImmMix_Endo' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN4'] <- 'CN4_Tumor_ImmMix_Fibr' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN5'] <- 'CN5_ImmMix_Fibr'
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN6'] <- 'CN6_ImmMix_Lym' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN7'] <- 'CN7_ImmMix_Bil' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN8'] <- 'CN8_Tumor_ImmMix_CD8T' 
gc_csd_CN_anno$All_CN[gc_csd_CN_anno$All_CN=='CN9'] <- 'CN9_Tumor_ImmMix_Mac' 



# ct-based-cn in each region/class -----------------------------------------------------------

cells_CN <- gc_csd_CN_anno %>% dplyr::select(Class, All_CN)


dat <- as.data.frame(with(cells_CN, table(Class, All_CN)))
cells_CN_Freq <- spread(dat, All_CN, Freq)

dat_CN_freq <- cells_CN_Freq[,-1]
row.names(dat_CN_freq) <- cells_CN_Freq$Class


dat_CN_percent <- dat_CN_freq / rowSums(dat_CN_freq) * 100
row.names(dat_CN_percent)<- cells_CN_Freq$Class




# clustering in heatmap ---------------------------------------------------------


# pdf(file = "1_pheatmap_allCN.pdf", width =10, height = 15)
p <- pheatmap(dat_CN_percent,
              cellwidth = 10,
              cellheight = 1,
              show_rownames = F,
              cutree_rows = 3)
# dev.off() 

# pdf(file = "2_pheatmap_allCN.pdf", width =10, height = 15)
# bk <- c(seq(-4,-0.1,by=0.01),seq(0,4,by=0.01))
# p <- pheatmap(dat_CN_percent,
#               cellwidth = 10,
#               cellheight = 1,
#               # annotation_row = df_TP,
#               border=FALSE,
#               color=c(colorRampPalette(colors=c("blue","white"))(length(bk)/2),colorRampPalette(color=c("white","red"))(length(bk)/2)),
#               scale = 'row',
#               cluster_rows = T,
#               show_rownames = F,
#               show_colnames = T,
#               legend_breaks=seq(-4,4,2),
#               breaks=bk,
#               cutree_rows = 14)
# print(p)
# dev.off() 




row_cluster <- as.data.frame(cutree(p$tree_row,k=3))
table(row_cluster)
names(row_cluster) <- 'PhenoGraph'
row_cluster$Class <-  row.names(row_cluster)

dat_CN_percent_class <- dat_CN_percent
dat_CN_percent_class$Class <- row.names(dat_CN_percent)

df_merge <- merge(row_cluster,dat_CN_percent_class,by='Class')

df_merge <- df_merge[order(df_merge$PhenoGraph),]
temp <- df_merge[,c(-1,-2)]
row.names(temp) <- df_merge$Class 

anno <- df_merge[2]
row.names(anno) <- df_merge$Class
anno$PhenoGraph <- as.factor(anno$PhenoGraph)

p1 <- pheatmap(temp,
               cellwidth = 10,
               cellheight = 1,
               show_rownames = F,
               cutree_rows = 3,
               annotation_row = anno)

df_merge$PhenoGraph[df_merge$PhenoGraph==1] <- 'H'
df_merge$PhenoGraph[df_merge$PhenoGraph==2] <- 'L'
df_merge$PhenoGraph[df_merge$PhenoGraph==3] <- 'M'

df_merge <- df_merge %>% dplyr::rename(TumorPurityLevel=PhenoGraph)

anno <- df_merge[2]
row.names(anno) <- df_merge$Class
temp <- df_merge[,c(-1,-2)]
row.names(temp) <- df_merge$Class 
p1 <- pheatmap(temp,
               cellwidth = 10,
               cellheight = 1,
               show_rownames = F,
               cutree_rows = 3,
               annotation_row = anno)



# combination of TPL with clinical information -----------------------------------------------------

mydata_pheno <- df_merge
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

write.csv(mydata_pheno_anno, file = 'mydata_pheno_anno_ctCN.csv')



# TPL, TP, SP and clinical information for each region/class -----------------------------------------------------------------------


CN_basedon_celltype <- read.csv('mydata_pheno_anno_ctCN.csv',header = T,row.names = 1,stringsAsFactors = F)
colnames(CN_basedon_celltype)[3:12] <- paste0('ct_based_',colnames(CN_basedon_celltype)[3:12])
CN_basedon_celltype <- CN_basedon_celltype[2:13]
CN_basedon_celltype <- dplyr::rename(CN_basedon_celltype, Class=Class0)

Class_Pheno_Anno <- read.csv('Class_Pheno_Anno.csv',header = T,row.names = 1,stringsAsFactors = F)

Class_Pheno_Anno_ctCN <- merge(CN_basedon_celltype,Class_Pheno_Anno,by='Class')


write.csv(Class_Pheno_Anno_ctCN,file = 'Class_Pheno_Anno_ctCN.csv')
