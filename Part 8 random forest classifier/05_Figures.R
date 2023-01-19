
### R 4.1.0


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
library(pROC)  #‘1.18.0’
set.seed(1234)




# Fig. 5A ROC -------------------------------------------------------------

######SP_Pf
load(file = "./Figures/SP_Pf_all.Rdata")

###Training set: ROC
roc1 = roc(df.train$class, as.numeric(res$N),
           main = "SP-Pf: Confidence intervals (training)", 
           percent=TRUE,
           ci = TRUE,                  # compute AUC (of AUC by default)
           print.auc = TRUE,
           auc.polygon = T,
           grid=c(1, 1),
           grid.col=c("grey"),
           max.auc.polygon=TRUE,
           auc.polygon.col="skyblue")  

plot(roc1)
sp.obj1 <- ci.sp(roc1, sensitivities=seq(0, 100, 1), boot.n=100) 
plot(sp.obj1, type="shape", col = "#EBB98A")
dev.off()


###testing set: ROC
roc1 = roc(df_roc$observed, as.numeric(df_roc$prop),
           main = "SP-Pf: Confidence intervals (testing)", 
           percent=TRUE,
           ci = TRUE,                  # compute AUC (of AUC by default)
           print.auc = TRUE,
           auc.polygon = T,
           grid=c(1, 1),
           grid.col=c("grey"),
           max.auc.polygon=TRUE,
           auc.polygon.col="skyblue")           # print the AUC (will contain the CI)
plot(roc1)
sp.obj1 <- ci.sp(roc1, sensitivities=seq(0, 100, 1), boot.n=100) 
plot(sp.obj1, type="shape", col = "#C7DFDA")
dev.off()



####### SP-HII 

load(file = "./Figures/SP_HII_all.Rdata")

###Training set: ROC
roc1 = roc(df.train$class, as.numeric(res$N),
           main = "SP-HII: Confidence intervals (training)", 
           percent=TRUE,
           ci = TRUE,                  # compute AUC (of AUC by default)
           print.auc = TRUE,
           auc.polygon = T,
           grid=c(1, 1),
           grid.col=c("grey"),
           max.auc.polygon=TRUE,
           auc.polygon.col="skyblue")

plot(roc1)
sp.obj1 <- ci.sp(roc1, sensitivities=seq(0, 100, 1), boot.n=100) 
plot(sp.obj1, type="shape", col = "#EBB98A")
dev.off()

###testing set: ROC
roc1 = roc(df_roc$observed, as.numeric(df_roc$prop),
           main = "SP-HII: Confidence intervals (testing)", 
           percent=TRUE,
           ci = TRUE,                  # compute AUC (of AUC by default)
           print.auc = TRUE,
           auc.polygon = T,
           grid=c(1, 1),
           grid.col=c("grey"),
           max.auc.polygon=TRUE,
           auc.polygon.col="skyblue")           # print the AUC (will contain the CI)

plot(roc1)
sp.obj1 <- ci.sp(roc1, sensitivities=seq(0, 100, 1), boot.n=100) 
plot(sp.obj1, type="shape", col = "#C7DFDA")
dev.off()




###### SP_LII 

load(file = "./Figures/SP_LII_all.Rdata")

###Training set: ROC
roc1 = roc(df.train$class, as.numeric(res$N),
           main = "SP-LII: Confidence intervals (training)", 
           percent=TRUE,
           ci = TRUE,                  # compute AUC (of AUC by default)
           print.auc = TRUE,
           auc.polygon = T,
           grid=c(1, 1),
           grid.col=c("grey"),
           max.auc.polygon=TRUE,
           auc.polygon.col="skyblue")  
plot(roc1)
sp.obj1 <- ci.sp(roc1, sensitivities=seq(0, 100, 1), boot.n=100) 
plot(sp.obj1, type="shape", col = "#EBB98A")
dev.off()


###Testing set : ROC
roc1 = roc(df_roc$observed, as.numeric(df_roc$prop),
           main = "SP-LII: Confidence intervals (testing)", 
           percent=TRUE,
           ci = TRUE,                  # compute AUC (of AUC by default)
           print.auc = TRUE,
           auc.polygon = T,
           grid=c(1, 1),
           grid.col=c("grey"),
           max.auc.polygon=TRUE,
           auc.polygon.col="skyblue")           # print the AUC (will contain the CI)

plot(roc1)
sp.obj1 <- ci.sp(roc1, sensitivities=seq(0, 100, 1), boot.n=100) 
plot(sp.obj1, type="shape", col = "#C7DFDA")
dev.off()




# Fig. 5D Gini importance -------------------------------------------------

rm(list=ls())
gc()
library(tidyverse)
library(randomForest)  
load(file='./Figures/F4/import.Rdata')

temp <- mydata_pheno_anno
temp$SP_LII <- as.character(temp$SP_LII)
temp$SP_HII <- as.character(temp$SP_HII)
temp$SP_Pf <- as.character(temp$SP_Pf)

library(ggpubr)


var_to_plot <-  c('SP_HII', 'SP_LII','SP_Pf')


for (sp in var_to_plot){
  inputpath <- paste0('./',sp,'_GINI.csv')  ##改一下
  df_imp <- read.csv(inputpath,header = T, row.names = 1, stringsAsFactors = F)
  df_imp_top <- df_imp[order(df_imp$MeanDecreaseGini,decreasing = T)[1:30],]
  cyc <- df_imp_top$index
  
  FeatureEnrichment <- c()
  for (i in seq_along(cyc)) {
    temp[,cyc[i]] <- scale(temp[[cyc[i]]])
    a <- temp[[cyc[i]]][temp[[sp]] == '0']
    b <- temp[[cyc[i]]][temp[[sp]] == '1']
    a0 <- mean(a)
    b0 <- mean(b)
    fc <- b0-a0
    FeatureEnrichment <- append(FeatureEnrichment,fc)
  }
  
  df_imp_top$FeatureEnrichment <-  FeatureEnrichment
  
  p <- ggplot(data=df_imp_top, mapping = aes(x=reorder(index,MeanDecreaseGini),y=MeanDecreaseGini, fill=FeatureEnrichment)) +
    geom_bar(stat = 'identity')+
    coord_flip() + 
    scale_fill_gradient2(low = "blue", high = 'red', mid = "white", midpoint = 0)
  
  ggsave(p, file= paste0('./Figures/F4/',sp,'_GINI.pdf'), height = 5, width = 7)
  #ggsave(p, file= paste0('./Figures/F4/',sp,'_GINI.png'), height = 5, width = 7)
  
}



sp='SP_Pf'
inputpath <- paste0('./',sp,'_GINI.csv')  ##改一下
df_imp <- read.csv(inputpath,header = T, row.names = 1, stringsAsFactors = F)
df_imp_top <- df_imp[order(df_imp$MeanDecreaseGini,decreasing = T)[1:30],]
cyc <- df_imp_top$index

FeatureEnrichment <- c()
for (i in seq_along(cyc)) {
  temp[,cyc[i]] <- scale(temp[[cyc[i]]])
  a <- temp[[cyc[i]]][temp[[sp]] == '0']
  b <- temp[[cyc[i]]][temp[[sp]] == '1']
  a0 <- mean(a)
  b0 <- mean(b)
  fc <- b0-a0
  FeatureEnrichment <- append(FeatureEnrichment,fc)
}

df_imp_top$FeatureEnrichment <-  FeatureEnrichment

library(ggbreak)

p <- ggplot(data=df_imp_top, mapping = aes(x=reorder(index,MeanDecreaseGini),y=MeanDecreaseGini, fill=FeatureEnrichment)) +
  geom_bar(stat = 'identity')+
  coord_flip() + 
  #scale_fill_gradient2(low = "cyan", high = 'red', mid = "black", midpoint = 0)
  scale_fill_gradientn(colours = c("#4DBBD5FF",  "#E64B35FF"),  values = scales::rescale(c(-0.5, -0.05, 0, 0.05, 0.5, 1.0,1.8))) 
p1 <- p+scale_y_break(c(3,25),scales = 0.5)



print(p1)
ggsave(p1, file= paste0('./Figures/F4/',sp,'_GINI.pdf'), height = 5, width = 7)
#ggsave(p, file= paste0('./Figures/F4/',sp,'_GINI.png'), height = 5, width = 7)


library("scales")
library(ggsci)
show_col(pal_npg("nrc")(10))
show_col(pal_npg("nrc", alpha = 0.8)(10))

sp='SP_HII'
inputpath <- paste0('./',sp,'_GINI.csv')  ##改一下
df_imp <- read.csv(inputpath,header = T, row.names = 1, stringsAsFactors = F)
df_imp_top <- df_imp[order(df_imp$MeanDecreaseGini,decreasing = T)[1:30],]
cyc <- df_imp_top$index

FeatureEnrichment <- c()
for (i in seq_along(cyc)) {
  temp[,cyc[i]] <- scale(temp[[cyc[i]]])
  a <- temp[[cyc[i]]][temp[[sp]] == '0']
  b <- temp[[cyc[i]]][temp[[sp]] == '1']
  a0 <- mean(a)
  b0 <- mean(b)
  fc <- b0-a0
  FeatureEnrichment <- append(FeatureEnrichment,fc)
}

df_imp_top$FeatureEnrichment <-  FeatureEnrichment

p <- ggplot(data=df_imp_top, mapping = aes(x=reorder(index,MeanDecreaseGini),y=MeanDecreaseGini, fill=FeatureEnrichment)) +
  geom_bar(stat = 'identity')+
  coord_flip() + 
  #scale_fill_gradient2(low = "cyan", high = 'red', mid = "black", midpoint = 0)
  scale_fill_gradientn(colours = c("#4DBBD5FF",  "#E64B35FF"),  values = scales::rescale(c(-1.0,-0.5, -0.05, 0, 0.05,0.5, 0.9)))
print(p)
ggsave(p, file= paste0('./Figures/F4/',sp,'_GINI.pdf'), height = 5, width = 7)
#ggsave(p, file= paste0('./Figures/F4/',sp,'_GINI.png'), height = 5, width = 7)


sp='SP_LII'
inputpath <- paste0('./',sp,'_GINI.csv')  ##改一下
df_imp <- read.csv(inputpath,header = T, row.names = 1, stringsAsFactors = F)
df_imp_top <- df_imp[order(df_imp$MeanDecreaseGini,decreasing = T)[1:30],]
cyc <- df_imp_top$index

FeatureEnrichment <- c()
for (i in seq_along(cyc)) {
  temp[,cyc[i]] <- scale(temp[[cyc[i]]])
  a <- temp[[cyc[i]]][temp[[sp]] == '0']
  b <- temp[[cyc[i]]][temp[[sp]] == '1']
  a0 <- mean(a)
  b0 <- mean(b)
  fc <- b0-a0
  FeatureEnrichment <- append(FeatureEnrichment,fc)
}

df_imp_top$FeatureEnrichment <-  FeatureEnrichment

p <- ggplot(data=df_imp_top, mapping = aes(x=reorder(index,MeanDecreaseGini),y=MeanDecreaseGini, fill=FeatureEnrichment)) +
  geom_bar(stat = 'identity')+
  coord_flip() + 
  #scale_fill_gradient2(low = "cyan", high = 'red', mid = "black", midpoint = 0)
  scale_fill_gradientn(colours = c("#4DBBD5FF",  "#E64B35FF"),  values = scales::rescale(c(-1.0,-0.5, -0.05, 0, 0.05,0.5, 1.2)))
print(p)
ggsave(p, file= paste0('./Figures/F4/',sp,'_GINI.pdf'), height = 5, width = 7)
#ggsave(p, file= paste0('./Figures/F4/',sp,'_GINI.png'), height = 5, width = 7)



# Fig. 5B evaluate classifiers on validation dataset ----------------------

rm(list=ls())
gc()

load(file ='mydata_pheno_anno.Rdata')  ##validation data
Val_mydata_anno <- mydata_pheno_anno

Val_SP_LII_predict <- readRDS(file = 'Validation _predict_SP_LII.Rds')
names(Val_SP_LII_predict) <- c('SP_LII', 'SP_LII_N', 'SP_LII_Y')

Val_SP_HII_predict <- readRDS(file = 'Validation _predict_SP_HII.Rds')
names(Val_SP_HII_predict) <- c('SP_HII', 'SP_HII_N', 'SP_HII_Y')

Val_SP_Pf_predict <- readRDS(file = 'Validation _predict_SP_Pf.Rds')
names(Val_SP_Pf_predict) <- c('SP_Pf', 'SP_Pf_N', 'SP_Pf_Y')


Val_SP_predict <- cbind(Val_SP_LII_predict, Val_SP_HII_predict, Val_SP_Pf_predict)

rownames(Val_SP_predict) <- Val_mydata_anno$Class0



Val_SP_predict_filter <- Val_SP_predict[!(Val_SP_predict$SP_LII == 'N' & Val_SP_predict$SP_HII == 'N' & Val_SP_predict$SP_Pf == 'N'), ]

Val_SP_predict_filter$SP <- 'SP_LII'
Val_SP_predict_filter$SP[Val_SP_predict_filter$SP_HII == 'Y'] <- 'SP_HII'
Val_SP_predict_filter$SP[Val_SP_predict_filter$SP_Pf == 'Y'] <- 'SP_Pf'

Val_SP_predict_filter$Class <- rownames(Val_SP_predict_filter)


###
df_anno <- mydata_pheno_anno %>% select(Class0, RFS01, RFSday, OS01, Osday)

df_anno_SP <-  merge(Val_SP_predict_filter, df_anno, by.x = 'Class', by.y = 'Class0')




######heatmap

df_anno_SP_cn <- df_anno_SP %>% select(Class, SP)

class_allCN <- read.csv(file='class_allCN.csv', header = T, row.names = 1)

SP_class_allCN <- merge(df_anno_SP_cn, class_allCN, by='Class')


SP_class_allCN <- SP_class_allCN[order(SP_class_allCN$SP),]

# save(df_anno_SP,SP_class_allCN,file = './Figures/F5/classifier on validationset.Rdata')

library(pheatmap)
annotation <- SP_class_allCN %>% select(SP)
rownames(annotation) <- SP_class_allCN$Class
rownames(SP_class_allCN) <- SP_class_allCN$Class

temp <- SP_class_allCN %>% select( "CN30_ImmMix_Mac",
                                   "CN17_ImmMix_B",
                                   "CN15_ImmMix_Lym",
                                   "CN8_Biliary_ImmMix",
                                   "CN9_ImmMix_CD4T" ,
                                   
                                   "CN13_Fibro_HII"  ,
                                   "CN23_Fibro_HII"  ,
                                   "CN37_Fibro"  ,
                                   
                                   "CN3_Endo",
                                   # "CN1_Endo_Tumor",
                                   
                                   # "CN0_HLA_DR",
                                   #"CN6_others",
                                   
                                   # "CN26_Biliary" ,
                                   
                                   "CN32_Ki67"   ,
                                   "CN33_Ki67_LII"  ,
                                   #"CN5"
)

pheatmap(t(scale(temp)),
         annotation_col = annotation,
         cluster_cols = F,
         color = colorRampPalette(c("#35978f", "white", "#bf812d"))(100))




# Fig. 5C classifiers on all samples 400 ----------------------------------


rm(list=ls())
gc()

# library(data.table)
library(tidyverse)
# library(patchwork)
library(tidyr)
#library(ComplexHeatmap)
library(ggthemes)
library(Rcpp)
library(igraph)
library(ggthemes)
library(pheatmap)
library(ggsci)
library(survival)
library(survminer)
library(randomForest)

######choose common features between discovery and validation set

###Validation dataset
load(file='mydata_pheno_anno.Rdata')
Val_mydata <- df_merge
View(as.data.frame(names(mydata_pheno_anno)))
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

Val_com <- Val_mydata[c('Class',com_ct,com_subtypes, com_spaInt )]
Tra_com <- Tra_mydata[c('Class',com_ct,com_subtypes, com_spaInt )]

All_class <- rbind(Tra_com, Val_com)


#####SP-HII classifier
load(file='SP_HII.Rdata')

forest.pred <- predict (fit.forest,All_class)
pre_df <- as.data.frame(predict(fit.forest, All_class, type = 'prob'))
table(forest.pred)

Val_SP_HII_predict <- cbind(forest.pred, pre_df)

#####SP-LII classifier
load(file='SP_LII.Rdata')

forest.pred <- predict (fit.forest,All_class)
pre_df <- as.data.frame(predict(fit.forest, All_class, type = 'prob'))
table(forest.pred)

Val_SP_LII_predict <- cbind(forest.pred, pre_df)


#####SP-Pf classifier
load(file='SP_Pf.Rdata')

forest.pred <- predict (fit.forest,All_class)
pre_df <- as.data.frame(predict(fit.forest, All_class, type = 'prob'))
table(forest.pred)

Val_SP_Pf_predict <- cbind(forest.pred, pre_df)


###combine

names(Val_SP_LII_predict) <- c('SP_LII', 'SP_LII_N', 'SP_LII_Y')

names(Val_SP_HII_predict) <- c('SP_HII', 'SP_HII_N', 'SP_HII_Y')

names(Val_SP_Pf_predict) <- c('SP_Pf', 'SP_Pf_N', 'SP_Pf_Y')


Val_SP_predict <- cbind(Val_SP_LII_predict, Val_SP_HII_predict, Val_SP_Pf_predict)

rownames(Val_SP_predict) <- All_class$Class

All_class_predict <- cbind(Val_SP_predict,All_class )

#
All_class_predict_filter <- All_class_predict
All_class_predict_filter$SP [All_class_predict_filter$SP_LII == 'N' & All_class_predict_filter$SP_HII == 'N' & All_class_predict_filter$SP_Pf == 'N' ] <- 'SP_others'

All_class_predict_filter$SP [All_class_predict_filter$SP_LII == 'Y']<- 'SP_LII'
All_class_predict_filter$SP[All_class_predict_filter$SP_HII == 'Y'] <- 'SP_HII'
All_class_predict_filter$SP[All_class_predict_filter$SP_Pf == 'Y'] <- 'SP_Pf'


####clinical information
mydata_pheno <- All_class_predict_filter
mydata_pheno$Class0 <- rownames(mydata_pheno)
mydata_pheno$Class0 <- as.character(mydata_pheno$Class0)
mydata_pheno$Class <- as.character(mydata_pheno$Class0)

load("HCC401_new.Rdata")
HCCpathology <- rbind(df_validation_new, df_training_new)
HCCpathology$Class <- paste0(HCCpathology$TMA, "_reg",sprintf("%03d",HCCpathology$Reg))

for (i in 1:length(mydata_pheno$Class)) {
  temp<-strsplit(mydata_pheno$Class,'_')[[i]][1]
  temp1 <-strsplit(temp,"A")[[1]][2]
  temp2 <- strsplit(mydata_pheno$Class,'_')[[i]][3]
  mydata_pheno$Class[i]<- paste0(temp1,"_",temp2)
}

mydata_pheno_anno <- merge(mydata_pheno, HCCpathology, by="Class")


fit <- survfit(Surv(Osday, OS01) ~ SP, data = mydata_pheno_anno)
ggsurvplot(fit,
           #palette = "npg",
           palette= c('#E64B35FF', '#4DBBD5FF','lightgrey','#F39B7FFF'),
           risk.table = F, 
           pval = TRUE,
           #conf.int = FALSE, 
           xlab="Time in Days",
           ggtheme = theme_bw(),
           surv.median.line = "hv", # add the median survival pointer.
           tables.y.text = T,
           #risk.table.pos = "in",
           #risk.table.col = "strata",
           fontsize = 3.5,
           pval.size = 4,
           surv.plot.height = 0.7,
           tables.height = 0.3,
           pval.coord = c(2000, 0.9)
)

fit <- survfit(Surv(RFSday, RFS01) ~ SP, data = mydata_pheno_anno)
ggsurvplot(fit,
           #palette = "npg",
           palette= c('#E64B35FF', '#4DBBD5FF','lightgrey','#F39B7FFF'),
           risk.table = F, 
           pval = TRUE,
           #conf.int = FALSE, 
           xlab="Time in Days",
           ggtheme = theme_bw(),
           surv.median.line = "hv", # add the median survival pointer.
           tables.y.text = T,
           #risk.table.pos = "in",
           #risk.table.col = "strata",
           fontsize = 3.5,
           pval.size = 4,
           surv.plot.height = 0.7,
           tables.height = 0.3,
           pval.coord = c(2000, 0.9)
)

# Fig. 5E TP-specific HR --------------------------------------------------

rm(list=ls())
gc()

# library(data.table)
library(tidyverse)
# library(patchwork)
library(tidyr)
#library(ComplexHeatmap)
library(ggthemes)
library(Rcpp)
library(igraph)
library(ggthemes)
library(pheatmap)
library(ggsci)
library(survival)
library(survminer)
###training dataset
load(file='./random forest classifier in TrainingData/import.Rdata')
Tra_mydata <- mydata_pheno_anno

library(RColorBrewer)
cols=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)


mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                 legend.key = element_rect(fill = "white", colour = "white"), #图标
                 legend.background = (element_rect(colour= "white", fill = "white")))

#####validation,需要用到common指标
load(file='importData for random forest.Rdata')


#####TP-specific
sp='SP_LII'


for (sp in unique(Tra_mydata$SP)) {
  mydata_pheno_anno_sp <- dplyr::filter(Tra_mydata,SP==sp)
  

  df <- mydata_pheno_anno_sp[ ,c(com_ct, com_subtypes)]
  mads <- apply(df,2,mad)
  useful_v <- colnames(df)[mads>0.0000001]
  
 
  covariates <- useful_v
  
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(RFSday, RFS01)~', x)))
  univ_formulas
  
 
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = mydata_pheno_anno_sp)})
  univ_models
  
  
  
  univ_results <- lapply(univ_models,
                         function(x){
                           x <- summary(x)
                          
                           p.value<-signif(x$wald["pvalue"], digits=3)
                          
                           HR <-signif(x$coef[2], digits=3);
                           
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                           res<-c(p.value,HR,HR.confint.lower,HR.confint.upper)
                           names(res)<-c("p.value","HR", 'HR.confint.lower','HR.confint.upper')
                           return(res)
                         })
  
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  #View(as.data.frame(res))
  
  dysfx <- as.data.frame(res)
  dysfx$index <- row.names(dysfx)
  dysfx$Category <- row.names(dysfx)
  dysfx$Category[dysfx$Category %in% var_ct] <- 'Cell type'
  dysfx$Category[dysfx$Category %in% var_subtypes] <- 'Cell subgroup'
  # dysfx$Category[dysfx$Category %in% var_spaInt] <- 'cell-cell spatial interaction'
  #dysfx$Category[dysfx$Category %in% var_allcn] <- 'CN based on cell subgroup'
  #dysfx$Category[dysfx$Category %in% var_ctBcn] <- 'CN based on cell type'
  
  dysfx_p0.05 <- subset(dysfx, `p.value` <0.05)
  
  #write.csv(dysfx_p0.05,file=paste0('./Figures/F16_SP-sepcific prognosis/cox_0.05_', sp, '.csv') )
  
  names(dysfx_p0.05)
  
  dysfx_p0.05_sort <- dysfx_p0.05[order(dysfx_p0.05$HR,decreasing = T),]
  dysfx_p0.05_sort$index <- factor(dysfx_p0.05_sort$index, levels = dysfx_p0.05$index[order(dysfx_p0.05$HR,decreasing = T)])
  p <- ggplot(dysfx_p0.05_sort, aes(x = index, y = HR, color= Category)) +
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymax = HR.confint.upper, ymin = HR.confint.lower))+
    scale_color_manual(values=c(pal_npg('nrc',alpha = 1.0)(6)))+
    # coord_flip()+
    # ylim(0,20)+
    geom_hline(yintercept = 1,color = 'black')+
    # guides(color=F)+
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))+
    coord_flip()
  
  ggsave(p, file=paste0('./Figures/F8/pvalue0.05_HR_', sp, '.pdf'), width = 5, height = 2)
  
}




#####TP-specific
sp='SP_LII'


for (sp in unique(Tra_mydata$SP)) {
  mydata_pheno_anno_sp <- dplyr::filter(Tra_mydata,SP==sp)
  
  
  df <- mydata_pheno_anno_sp[ ,c(com_ct, com_subtypes,com_spaInt_filter)]
  mads <- apply(df,2,mad)
  useful_v <- colnames(df)[mads>0.0000001]
  
  
  covariates <- useful_v
  
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(RFSday, RFS01)~', x)))
  univ_formulas
  
  
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = mydata_pheno_anno_sp)})
  univ_models
  
  
  
  univ_results <- lapply(univ_models,
                         function(x){
                           x <- summary(x)
                           
                           p.value<-signif(x$wald["pvalue"], digits=3)
                           
                           HR <-signif(x$coef[2], digits=3);
                           
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                           res<-c(p.value,HR,HR.confint.lower,HR.confint.upper)
                           names(res)<-c("p.value","HR", 'HR.confint.lower','HR.confint.upper')
                           return(res)
                         })
  
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  #View(as.data.frame(res))
  
  dysfx <- as.data.frame(res)
  dysfx$index <- row.names(dysfx)
  dysfx$Category <- row.names(dysfx)
  dysfx$Category[dysfx$Category %in% var_ct] <- 'Cell type'
  dysfx$Category[dysfx$Category %in% var_subtypes] <- 'Cell subgroup'
  dysfx$Category[dysfx$Category %in% var_spaInt] <- 'PCI-LR'
  #dysfx$Category[dysfx$Category %in% var_allcn] <- 'CN based on cell subgroup'
  #dysfx$Category[dysfx$Category %in% var_ctBcn] <- 'CN based on cell type'
  
  dysfx_p0.05 <- subset(dysfx, `p.value` <0.05)
  
  write.csv(dysfx_p0.05,file=paste0('./Figures/1121_cox_0.05_', sp, '.csv') )
  
  names(dysfx_p0.05)
  
  dysfx_p0.05_sort <- dysfx_p0.05[order(dysfx_p0.05$HR,decreasing = T),]
  dysfx_p0.05_sort$index <- factor(dysfx_p0.05_sort$index, levels = dysfx_p0.05$index[order(dysfx_p0.05$HR,decreasing = T)])
  p <- ggplot(dysfx_p0.05_sort, aes(x = index, y = HR, color= Category)) +
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymax = HR.confint.upper, ymin = HR.confint.lower))+
    scale_color_manual(values=c(pal_npg('nrc',alpha = 1.0)(6)))+
    # coord_flip()+
    # ylim(0,20)+
    geom_hline(yintercept = 1,color = 'black')+
    # guides(color=F)+
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))+
    coord_flip()
  
  ggsave(p, file=paste0('./Figures/1121_pvalue0.05_HR_', sp, '.pdf'), width = 5, height = 2)
  
}








