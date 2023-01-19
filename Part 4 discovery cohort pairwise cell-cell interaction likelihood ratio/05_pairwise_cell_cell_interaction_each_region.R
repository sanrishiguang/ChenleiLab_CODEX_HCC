

########R 4.1.0

library(circlize)  #‘0.4.15’
library(dplyr)  #‘1.0.9’
library(ggplot2) #‘3.3.6’
library(reshape) # ‘0.8.9’
library(RColorBrewer) #‘1.1.3’
library(pheatmap)  #‘1.0.12’
library (plyr)  #‘1.8.7’
library(ggpubr) # ‘0.4.0’
library(tidyr)  #‘1.2.0’
library(ggsignif)  # ‘0.6.3’
library(ggsci)   # ‘2.9’
library(survival)   #‘3.3.1’
library(survminer)   # ‘0.4.9’
library(randomForest)  #‘4.7.1.1’



#######computation of pairwise cell-cell contacts in each class/region######


# All subtypes/phenotypes -------------------------------------------------------------


#定义likehood ratio
fun_LikelihoodRatio <- function(x){
  x1 <- x[,-c(1,2)]
  row.names(x1) <- x$Allsubtypes
  #colsums=0
  x2 <- x1[,colSums(x1)!=0]
  x2 <- x2[order(rownames(x2)),order(colnames(x2))]
  dat_filter <- as.matrix(x2)
  #rowSums(dat_filter)
  sum(dat_filter==Inf)
  dat_filter[dat_filter==Inf] <- 100000
  
  dat_new<-matrix(nrow = length(rownames(dat_filter)), ncol=length(colnames(dat_filter)), dimnames = list(rownames(dat_filter),colnames(dat_filter)))
  
  dat_sum <- sum(dat_filter)
  for (i in 1:length(rownames(dat_filter))) {
    for (j in 1:length(colnames(dat_filter))) 
    {
      a0 <- dat_sum/sum(dat_filter[i,])/sum(dat_filter[,j])
      a<- dat_filter[i,j]*a0
      b<- log(a+0.0001)
      dat_new[i,j]<-b
    }
  }
  return(dat_new)
}

#输入文件

group2 <- read.csv("./spatial analysis/spatialanalysis_output_subtype_class.csv",header = T,stringsAsFactors = F)

dim(group2)

names(group2) <- gsub('[.]','_',colnames(group2))
group2$Allsubtypes <- gsub(' ','_', group2$Allsubtypes)
group2$Allsubtypes <- gsub('[+]','_', group2$Allsubtypes)
group2$Allsubtypes <- gsub('[-]','_', group2$Allsubtypes)

group_by_class <- split(group2,f=group2$Class)

group_by_class_order <- lapply(group_by_class, fun_LikelihoodRatio)

group_by_class_df <- lapply(group_by_class_order, as.data.frame)

group_by_class_df <- lapply(group_by_class_df, function(x){
  x$Allsubytpes <- row.names(x)
  return(x)})

group_by_class_df_melt <- lapply(group_by_class_df, function(x){
  x3 <- reshape::melt(x,id='Allsubytpes')
  return(x3)})

df_in_one <- ldply (group_by_class_df_melt, data.frame)
df_in_one <- dplyr::rename(df_in_one,Class=.id)
df_in_one <- dplyr::mutate(df_in_one,index_surrounding = paste0(Allsubytpes,"_",variable))
df_in_one_new <- df_in_one


df_in_one_subtype_spread <- df_in_one_new %>% select(Class,index_surrounding,value) %>% spread(index_surrounding, value)





# celltype ----------------------------------------------------------------

fun_x <- function(x){
  x1 <- x[,-c(1,2)]
  row.names(x1) <- x$celltype
  #colsums=0的列删掉
  x2 <- x1[,colSums(x1)!=0]
  x2 <- x2[order(rownames(x2)),order(colnames(x2))]
  dat_filter <- as.matrix(x2)
  #rowSums(dat_filter)
  sum(dat_filter==Inf)
  dat_filter[dat_filter==Inf] <- 100000
  
  dat_new<-matrix(nrow = length(rownames(dat_filter)), ncol=length(colnames(dat_filter)), dimnames = list(rownames(dat_filter),colnames(dat_filter)))
  
  dat_sum <- sum(dat_filter)
  for (i in 1:length(rownames(dat_filter))) {
    for (j in 1:length(colnames(dat_filter))) 
    {
      a0 <- dat_sum/sum(dat_filter[i,])/sum(dat_filter[,j])
      a<- dat_filter[i,j]*a0
      b<- log(a+0.0001)
      dat_new[i,j]<-b
    }
  }
  return(dat_new)
}


df_in_one_new <- data.frame()

input_path <- './spatial analysis/spatialanalysis_output_celltype_class.csv'

group2 <- read.csv(input_path,header = T,stringsAsFactors = F)
dim(group2)
names(group2) <- gsub('[.]','_',colnames(group2))
# rownames(group2) <- gsub(' ','_', rownames(group2))
# rownames(group2) <- gsub('[+]','_', rownames(group2))
# rownames(group2) <- gsub('[/]','_', rownames(group2))
# 
group_by_class <- split(group2,f=group2$Class)

group_by_class_order <- lapply(group_by_class, fun_x)

group_by_class_df <- lapply(group_by_class_order, as.data.frame)

group_by_class_df <- lapply(group_by_class_df, function(x){
  x$celltype <- row.names(x)
  return(x)})

group_by_class_df_melt <- lapply(group_by_class_df, function(x){
  x3 <- reshape::melt(x,id='celltype')
  return(x3)})

df_in_one <- ldply (group_by_class_df_melt, data.frame)
df_in_one <- dplyr::rename(df_in_one,Class=.id)
df_in_one <- dplyr::mutate(df_in_one,index_surrounding = paste0(celltype,"_",variable))

df_in_one_new <- df_in_one


df_in_one_celltype_spread <- df_in_one_new %>% select(Class,index_surrounding,value) %>% spread(index_surrounding, value)



##combine clinical information

TP_anno <- read.csv(file='Class_Pheno_Anno.csv',header = T,row.names = 1,stringsAsFactors = F)

mydata_pheno_anno_subtype <- merge(TP_anno,df_in_one_subtype_spread,by='Class')
mydata_pheno_anno_subtype_celltype <- merge(mydata_pheno_anno_subtype,df_in_one_celltype_spread,by='Class')


write.csv(mydata_pheno_anno_subtype_celltype,file = 'Class_Pheno_Anno_spatialInteraction.csv')

View(t(mydata_pheno_anno_subtype_celltype))




#detailed information for each clas/region -------------------------------
rm(list=ls())
gc()

mydata_pheno_anno1 <- read.csv(file='Class_Pheno_Anno_ctCN.csv',header = T,row.names = 1,stringsAsFactors = F)
colnames(mydata_pheno_anno1) <- gsub('[.]','_',colnames(mydata_pheno_anno1))
colnames(mydata_pheno_anno1)


mydata_pheno_anno2 <- read.csv(file='Class_Pheno_Anno_spatialInteraction.csv', header = T, row.names = 1, stringsAsFactors = F)
colnames(mydata_pheno_anno2) <- gsub('[.]','_',colnames(mydata_pheno_anno2))
colnames(mydata_pheno_anno2)

mydata_pheno_anno <- merge(select(mydata_pheno_anno1, 1:12),mydata_pheno_anno2, by='Class')
colnames(mydata_pheno_anno)
dim(mydata_pheno_anno)

row.names(mydata_pheno_anno) <- mydata_pheno_anno$Class

#整理变量
mydata_pheno_anno$TNM_LH[mydata_pheno_anno$TNM=='I A'|mydata_pheno_anno$TNM=='I B'|mydata_pheno_anno$TNM=='II'] <- "low"
mydata_pheno_anno$TNM_LH[mydata_pheno_anno$TNM=='III A'|mydata_pheno_anno$TNM=='IV A'|mydata_pheno_anno$TNM=='IV B'] <- "high"

mydata_pheno_anno$Differentiation_LH [mydata_pheno_anno$Differentiation=='1'|mydata_pheno_anno$Differentiation=='2'|mydata_pheno_anno$Differentiation=='2~3'] <- 'high'
mydata_pheno_anno$Differentiation_LH [mydata_pheno_anno$Differentiation=='3'|mydata_pheno_anno$Differentiation=='4'] <- 'low'

mydata_pheno_anno$Age_60[mydata_pheno_anno$Age <60] <- "0"
mydata_pheno_anno$Age_60[mydata_pheno_anno$Age >=60] <- "1"


#
unique(mydata_pheno_anno$TP)
unique(mydata_pheno_anno$SP)

temp <- select(mydata_pheno_anno, Class, SP)
temp$Freq <- 1
temp <- spread(temp,SP,Freq)

temp1 <- select(mydata_pheno_anno, Class, TP)
temp1$Freq <- 1
temp1 <- spread(temp1,TP,Freq)

temp0 <- merge(temp,temp1, by='Class')
temp0[is.na(temp0)] <- 0

mydata_pheno_anno <- merge(mydata_pheno_anno,temp0, by='Class')


#spatial interaction variables
var_spaInt <- colnames(mydata_pheno_anno)[201 : 5505]
temp <- mydata_pheno_anno[var_spaInt]
temp_t <- t(temp)
nalist <- !(rowSums(is.na(temp_t)) > 200)
var_spaInt <- var_spaInt[nalist]

mydata_pheno_anno[var_spaInt] <- na.roughfix(mydata_pheno_anno[var_spaInt])


#ct_based_CN
var_ctBcn <- colnames(mydata_pheno_anno)[3: 12]

#allCN40
var_allcn <- colnames(mydata_pheno_anno)[13: 52]

#subtypes
var_subtypes <- colnames(mydata_pheno_anno)[53:124]

#celltypes
var_ct <- colnames(mydata_pheno_anno)[125:135]

#
var_anno <- colnames(mydata_pheno_anno)[c(157:198,5506:5508)]
#TP, SP
var_TP <- unique(mydata_pheno_anno$TP)
var_SP <- unique(mydata_pheno_anno$SP)

save(list=ls(), file='import.Rdata')



# Fig.2C-D harzard_ratio --------------------------------------------------------

rm(list=ls())
gc()

load('import.Rdata')

# covariates <- c(var_ct,var_subtypes)
# covariates <- c(var_spaInt)
covariates <- c(var_spaInt, var_subtypes,var_allcn,var_ct, var_ctBcn)

#
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(RFSday, RFS01)~', x)))
univ_formulas

#
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = mydata_pheno_anno)})
univ_models


univ_results <- lapply(univ_models,
                       function(x){
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         #获取HR
                         HR <-signif(x$coef[2], digits=3);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                         res<-c(p.value,HR,HR.confint.lower,HR.confint.upper)
                         names(res)<-c("p.value","HR", 'HR.confint.lower','HR.confint.upper')
                         return(res)
                       })

res <- t(as.data.frame(univ_results, check.names = FALSE))
View(as.data.frame(res))



dysfx <- as.data.frame(res)
dysfx$index <- row.names(dysfx)
dysfx$Category <- row.names(dysfx)
dysfx$Category[dysfx$Category %in% var_ct] <- 'Cell type'
dysfx$Category[dysfx$Category %in% var_subtypes] <- 'Cell state'
dysfx$Category[dysfx$Category %in% var_spaInt] <- 'cell-cell spatial interaction'
dysfx$Category[dysfx$Category %in% var_allcn] <- 'CN based on cell state'
dysfx$Category[dysfx$Category %in% var_ctBcn] <- 'CN based on cell type'


write.csv(dysfx, file='./Figures/F8_harzard_ratio/Tables S2.csv')
# # dysfx_p0.01 <- dplyr::filter(dysfx, p.value < 0.05 & (HR > 1.1 | HR < 0.9))
# dysfx <- read.csv('./Figures/F8_harzard_ratio/uni_cox_all_index.csv', header = T, row.names = 1)
dysfx $Category <- factor(dysfx $Category, levels = c('Cell type',
                                                      'Cell state',
                                                      'CN based on cell type',
                                                      'CN based on cell state',
                                                      'cell-cell spatial interaction'
))
dysfx_p0.01 <- dplyr::filter(dysfx, p.value < 0.01 )



dysfx_plot <- as.data.frame(table(dysfx$Category))
dysfx_p0.01_plot <-  as.data.frame(table(dysfx_p0.01$Category))

dysfx_p0.01_plot$Percentage <- round ( (dysfx_p0.01_plot$Freq / dysfx_plot$Freq) *100, 2)
dysfx_p0.01_plot

mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                 legend.key = element_rect(fill = "white", colour = "white"), #图标
                 legend.background = (element_rect(colour= "white", fill = "white")))

dysfx_p0.01_plot = dysfx_p0.01_plot[order(dysfx_p0.01_plot$Percentage),]
dysfx_p0.01_plot$Var1 = factor(dysfx_p0.01_plot$Var1, levels = c('Cell type',
                                                                 'CN based on cell type',
                                                                 'cell-cell spatial interaction',
                                                                 'Cell state',
                                                                 'CN based on cell state'
))
p <- ggplot(data=dysfx_p0.01_plot,mapping = aes(x=Var1, y=Percentage, fill=Var1)) + 
  geom_bar(stat = 'identity') +
  scale_fill_manual(values=c("Cell type" = "#edfcc2",
                             "Cell state" = "#455655",
                             "CN based on cell type" = "#60bdaf",
                             "CN based on cell state" = "#f88aaf",
                             "cell-cell spatial interaction" = "#a1d8b1"))+
  mytheme+
  coord_flip() +
  guides(fill=FALSE) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))
p
ggsave(p, filename = './Figures/F8_harzard_ratio/Perc_cox_0.01_0730.pdf', width = 4, height = 2)  



##Fig. 2C
names(dysfx_p0.01)
library(RColorBrewer)
cols=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)

dysfx_p0.01_sort <- dysfx_p0.01[order(dysfx_p0.01$HR,decreasing = T),]
dysfx_p0.01_sort$index <- factor(dysfx_p0.01_sort$index, levels = dysfx_p0.01$index[order(dysfx_p0.01$HR,decreasing = T)])

dysfx_p0.01_sort = dysfx_p0.01_sort[c("Macrophages_Vimentin_","CN35_Ki67_LII","Tumor_Ki67_","CN7_Fibro_HII","Tumor_CD8T"),]

p <- ggplot(dysfx_p0.01_sort, aes(x = index, y = HR, color = Category)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymax = HR.confint.upper, ymin = HR.confint.lower))+
  scale_color_manual(values=c("Cell type" = "#edfcc2",
                              "Cell state" = "#455655",
                              "CN based on cell type" = "#60bdaf",
                              "CN based on cell state" = "#f88aaf",
                              "cell-cell spatial interaction" = "#a1d8b1"))+
  coord_flip()+
  # ylim(0,20)+
  geom_hline(yintercept = 1,color = 'black')+
  guides(color=F)+
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))
p
ggsave(p, file='./Figures/F8_harzard_ratio/cox_0.01_HR_0730.pdf', width = 5, height = 2)



