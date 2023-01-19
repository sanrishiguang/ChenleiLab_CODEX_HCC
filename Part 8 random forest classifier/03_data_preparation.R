

###R 4.1.0


# 1. cell type, cell phenotype, and CN in each class/region ---------------------------------------------------------------------


rm(list=ls())
gc()
library(dplyr)  #‘1.0.9’
library(tidyverse)  #‘1.3.2’
library(reshape) #‘0.8.9’

load('data for TPSP definition.Rdata')


##########All CN
cells_CN <-dplyr::select(gc_csd_CN_anno_filter, Class, All_CN)  
dat <- as.data.frame(with(cells_CN, table(Class, All_CN)))
dat1 <- spread(dat, All_CN, Freq)
dat2 <- dat1[,-1]
dat_percent <- dat2 / rowSums(dat2) * 100
# names(dat_percent) <- paste0("AllCN",c(0:39))
dat_percent$Class <- dat1$Class
write.csv(dat_percent,file = 'class_allCN.csv')

####All subtype
cells_subtype <-dplyr::select(gc_csd_CN_anno_filter, Class, Allsubtypes)
dat <- as.data.frame(with(cells_subtype, table(Class, Allsubtypes)))
dat1 <- spread(dat, Allsubtypes, Freq)
dat2 <- dat1[,-1]
dat_percent <- dat2 / rowSums(dat2) * 100
dat_percent$Class <- dat1$Class
write.csv(dat_percent,file = 'class_allsubtypes.csv')

###All celltype
cells_celltype <-dplyr::select(gc_csd_CN_anno_filter, Class, celltype)
dat <- as.data.frame(with(cells_celltype, table(Class, celltype)))
dat1 <- spread(dat, celltype, Freq)
dat2 <- dat1[,-1]
dat_percent <- dat2 / rowSums(dat2) * 100
dat_percent$Class <- dat1$Class
write.csv(dat_percent,file = 'class_celltype.csv')


##merge tables
# rm(list=ls())
# gc()
allCN <- read.csv('class_allCN.csv',header = T,row.names = 1,stringsAsFactors = FALSE)
allsubtypes <- read.csv('class_allsubtypes.csv', header = T,row.names = 1,stringsAsFactors = FALSE)
allcelltype <- read.csv('class_celltype.csv', header = T,row.names = 1,stringsAsFactors = FALSE)


names(allcelltype)[1:12] <- paste0('celltype_',names(allcelltype)[1:12])



temp <- merge(allCN,allsubtypes,by='Class')
df_merge <- merge(temp, allcelltype,by='Class')



write.csv(df_merge,file = 'Class_CN_CT_ST.csv')




# 2. pairwise cell-cell interation likelihood ratio in validation dataset -------------------------------------------------------------


#######computation of pairwise cell-cell contacts


rm(list=ls())
gc()
# setwd("E:\\11. CODEX\\TrainingData\\20220308\\AllCN_r50_CN40_0331\\spatial analysis")
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

#
fun_LikelihoodRatio <- function(x){
  x1 <- x[,-c(1,2)]
  row.names(x1) <- x$Allsubtypes
  
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

#subtype
input_path <- paste0("./spatial analysis/spatialanalysis_output_subtype_5TPs_class.csv")

group2 <- read.csv(input_path,header = T,stringsAsFactors = F)
dim(group2)

names(group2) <- gsub('[.]','_',colnames(group2))

group2$Allsubtypes <- gsub(' ','_', group2$Allsubtypes)
group2$Allsubtypes <- gsub('[+]','_', group2$Allsubtypes)
group2$Allsubtypes <- gsub('[-]','_', group2$Allsubtypes)
group2$Allsubtypes <- gsub('[/]','_', group2$Allsubtypes)

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

write.csv(df_in_one_subtype_spread, file="Class_Spatial.csv")




# 3. combine data: class ct st spatial ---------------------------------------

df_spatial <- read.csv('Class_Spatial.csv', header = T, row.names = 1)
df_ctst <- read.csv("Class_CN_CT_ST.csv", header = T, row.names = 1)

df_merge <- merge(df_spatial, df_ctst, by='Class')
df_merge$Class0 <- df_merge$Class
mydata_pheno <- df_merge
mydata_pheno$Class <- as.character(mydata_pheno$Class)


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
View(as.data.frame(mydata_pheno_anno[1:2,]))

write.csv(mydata_pheno_anno, file='mydata_pheno_anno.csv')
save(mydata_pheno_anno,df_merge, file = 'mydata_pheno_anno.Rdata')



# 4. common variables between discovery and validation dataset ------------



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
library(randomForest)  #‘4.7.1.1’


###Validation dataset
load(file='mydata_pheno_anno.Rdata')
Val_mydata <- df_merge
#View(as.data.frame(names(mydata_pheno_anno)))
names(Val_mydata) <- gsub('[.]','_',colnames(Val_mydata))



#spatial interaction variables
var_spaInt <- colnames(Val_mydata)[2:3722]
temp <- Val_mydata[var_spaInt]
temp_t <- t(temp)
nalist <- !(rowSums(is.na(temp_t)) > 60)
var_spaInt <- var_spaInt[nalist]

Val_mydata[var_spaInt] <- na.roughfix(Val_mydata[var_spaInt])


Val_anno <- names(mydata_pheno_anno)[3836:3898]
Val_ct <- names(Val_mydata)[3824:3835]
Val_spaInt <- names(Val_mydata)[2:3722]
Val_subtypes <- names(Val_mydata)[3763:3823]

###training dataset
load(file='./random forest classifier in TrainingData/import.Rdata')
Tra_mydata <- mydata_pheno_anno

##common terms
com_ct <- intersect(Val_ct, var_ct)
com_subtypes <- intersect(Val_subtypes, var_subtypes)
com_spaInt <- intersect(Val_spaInt, var_spaInt)

com_ct <- com_ct[c(-1,-2, -3, -6, -7,-8,-10)]


####discovery/training spit into 8:2
train <- sample(nrow(Tra_mydata), 0.8*nrow(Tra_mydata))
Tra_mydata_df.train <- Tra_mydata[train,]
Tra_mydata_df.validate <- Tra_mydata[-train,]
table(Tra_mydata_df.train$Class)
table(Tra_mydata_df.validate$Class)



drop_var <- paste0(com_subtypes,"_", com_subtypes)
drop_var <- intersect(drop_var, com_spaInt)

Tra_mydata_df.train <- Tra_mydata_df.train %>% dplyr::select(!drop_var)
Tra_mydata_df.validate <- Tra_mydata_df.validate %>% dplyr::select(!drop_var)

`%notin%` <- Negate(`%in%`)
com_spaInt_filter <- com_spaInt [com_spaInt %notin% drop_var]

save(list=ls(), file='importData for random forest.Rdata')


