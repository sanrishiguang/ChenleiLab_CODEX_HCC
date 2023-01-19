


### R 4.1.0

#####12.6  SP_LII -------------------------------------------------------------------

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


load(file='importData for random forest.Rdata')

####Train
df.train <- Tra_mydata_df.train[c(com_ct,com_subtypes,com_spaInt_filter, 'SP_LII')]

df.train <- df.train %>% rename(class = SP_LII)
df.train$class <- factor(df.train$class, levels=c(0,1),
                         labels=c("N", "Y"))

sum(is.na(df.train))



###Val

df.validate <- Tra_mydata_df.validate[c(com_ct,com_subtypes,com_spaInt_filter, 'SP_LII')]

df.validate <- df.validate %>% rename(class = SP_LII)
df.validate$class <- factor(df.validate$class, levels=c(0,1),
                            labels=c("N", "Y"))

sum(is.na(df.validate))

dim(df.validate)
dim(df.train)


library (randomForest)
set.seed(1234)




######Kolmogorov-Smirnov test, filtering features
colnames(df.train)[dim(df.train)[2]]
cyc <- colnames(df.train)[-dim(df.train)[2]]

library("dgof")
sig_predictors <- c()
for (i in seq_along(cyc)) {
  a <- df.train[[cyc[i]]][df.train$class== 'Y']
  b <- df.train[[cyc[i]]][df.train$class== 'N']
  kstest <- ks.test(a,b)
  if (kstest$p.value < 0.05) {
    sig_predictors <- append(sig_predictors,cyc[i])
  }
  print(sig_predictors)
}

df.train.sig <- df.train %>% select(class, sig_predictors)


##spearman correlation between var_spaInt
df.cor <- df.train.sig [,colnames(df.train.sig) %in% var_spaInt]


cyc1 <- colnames(df.cor)

hypercorindex <- data.frame(index1=character(), index2=character())
for (i in seq_along(cyc1)) {
  for (j in seq_along(cyc1)) {
    if (j>i) {
      mytest <- cor.test(df.cor[[i]],df.cor[[j]], method =  "spearman")
      coreff <- abs(mytest$estimate)
      pval <- mytest$p.value
      if (coreff > 0.85 & pval < 0.05) {
        print(paste0(cyc1[i], ' vs ', cyc1[j]))
        add_var <- data.frame(index1 =  cyc1[i], index2 = cyc1[j] )
        hypercorindex <- rbind(hypercorindex, add_var)
      }
    }
  }
}


##keep one 
keep_var <- !(colnames(df.train.sig) %in% hypercorindex[[1]])
df.train.sig_fil <- df.train.sig[keep_var]  
colnames(df.train.sig_fil)



df.train.sig_fil2 <- df.train.sig_fil

df.cor <- df.train.sig_fil2
df.cor <- df.cor[,-1]

cyc1 <- colnames(df.cor)

hypercorindex <- data.frame(index1=character(), index2=character())
for (i in seq_along(cyc1)) {
  for (j in seq_along(cyc1)) {
    if (j>i) {
      mytest <- cor.test(df.cor[[i]],df.cor[[j]], method =  "spearman")
      coreff <- abs(mytest$estimate)
      pval <- mytest$p.value
      if (coreff > 0.85 & pval < 0.05) {
        print(paste0(cyc1[i], ' vs ', cyc1[j]))
        add_var <- data.frame(index1 =  cyc1[i], index2 = cyc1[j] )
        hypercorindex <- rbind(hypercorindex, add_var)
      }
    }
  }
}


keep_var <- !(colnames(df.train.sig_fil2) %in% hypercorindex[[2]])
df.train.sig_fil3 <- df.train.sig_fil2[keep_var]  
colnames(df.train.sig_fil3)



df.train.sig_fil2 <- df.train.sig_fi3

df.cor <- df.train.sig_fil2
df.cor <- df.cor[,-1]

cyc1 <- colnames(df.cor)

hypercorindex <- data.frame(index1=character(), index2=character())
for (i in seq_along(cyc1)) {
  for (j in seq_along(cyc1)) {
    if (j>i) {
      mytest <- cor.test(df.cor[[i]],df.cor[[j]], method =  "spearman")
      coreff <- abs(mytest$estimate)
      pval <- mytest$p.value
      if (coreff > 0.85 & pval < 0.05) {
        print(paste0(cyc1[i], ' vs ', cyc1[j]))
        add_var <- data.frame(index1 =  cyc1[i], index2 = cyc1[j] )
        hypercorindex <- rbind(hypercorindex, add_var)
      }
    }
  }
}


keep_var <- !(colnames(df.train.sig_fil2) %in% hypercorindex[[2]])
df.train.sig_fil3 <- df.train.sig_fil2[keep_var]  
colnames(df.train.sig_fil3)





##############establish random forest classifier

#choose optimal mtry value
n <- ncol(df.train.sig_fil3) -1
errRate <- c(1)
for (i in 1:100){
  m <- randomForest(class~.,data=df.train.sig_fil3,mtry=i,proximity=TRUE,na.action=na.roughfix)
  err<-mean(m$err.rate)
  errRate[i] <- err
}
print(errRate)


m= which.min(errRate)
print(m)


#choose optimal ntree value
set.seed(4321)
rf_ntree <- randomForest(class~.,data=df.train.sig_fil3,na.action=na.roughfix , ntree=5000)
plot(rf_ntree)



set.seed(4321)
fit.forest <- randomForest (class~., data=df.train.sig_fil3,    
                            na.action=na.roughfix,
                            importance=TRUE,
                            ntree=1500,
                            mtry=94,
                            proximity=T,
                            nodesize=10,   
                            maxnodes=200,
                            nPerm=1)
fit.forest

max(treesize(fit.forest, terminal=TRUE))  
MDSplot(fit.forest, df.train$class)


##GINI
varImpPlot(fit.forest,n.var=20)

df_imp <- as.data.frame(importance (fit.forest, type=2))  ##GINI
df_imp$index <- row.names(df_imp)
df_imp_top <- df_imp[order(df_imp$MeanDecreaseGini,decreasing = T)[1:30],]
ggplot(data=df_imp_top, mapping = aes(x=reorder(index,MeanDecreaseGini),y=MeanDecreaseGini, fill="#E64B35FF")) +
  geom_bar(stat = 'identity')+
  coord_flip()+
  guides (fill=F)

write.csv(df_imp,file='./SP_LII_GINI.csv')

##predict on tesing set

forest.pred <- predict (fit.forest,df.validate)
forest.perf <- table (df.validate$class, forest.pred,
                      dnn=c("Actual","Predicted"))
forest.perf


###roc based on training
library(pROC)
res <- as.data.frame(fit.forest$votes)

plot.roc(df.train$class, as.numeric(res$N),
         main = "Confidence intervals", 
         percent=TRUE,
         ci = TRUE,                  # compute AUC (of AUC by default)
         print.auc = TRUE)  



####classifier on testing
pre_df <- as.data.frame(predict(fit.forest, df.validate, type = 'prob'))
predictions <- pre_df$N

forest.pred <- predict (fit.forest,df.validate)

df_roc <- as.data.frame(cbind(predictions,forest.pred,df.validate$class))
names(df_roc) <- c('prop','predicted','observed')

plot.roc(df_roc$observed, as.numeric(df_roc$prop),
         main = "Confidence intervals", 
         percent=TRUE,
         ci = TRUE,                  # compute AUC (of AUC by default)
         print.auc = TRUE)           # print the AUC (will contain the CI)




####classifier on validation dataset 


ValSet98 <- Val_mydata[colnames(df.train.sig_fil3)[-1]]
dim(ValSet98)
sum(is.na(ValSet98))
forest.pred <- predict (fit.forest,ValSet98)
table(forest.pred)

forest.pred <- predict (fit.forest,ValSet98)
pre_df <- as.data.frame(predict(fit.forest, ValSet98, type = 'prob'))

load(file='mydata_pheno_anno.Rdata')
Val_mydata_anno <- mydata_pheno_anno
Val_mydata_anno$SP_LII <- forest.pred
table(Val_mydata_anno$SP_LII)

Val_SP_LII_predict <- cbind(forest.pred,pre_df )


library(survival)
library(survminer)
fit <- survfit(Surv(RFSday, RFS01) ~ SP_LII, data = Val_mydata_anno)
ggsurvplot(fit,palette = "npg",
           risk.table = TRUE, pval = TRUE,
           conf.int = FALSE, xlab="Time in Days",
           ggtheme = theme_classic())

save(Val_SP_LII_predict,Val_mydata, fit.forest, Tra_mydata_df.train,df.train, df.validate, Tra_mydata_df.validate, file = 'SP_LII.Rdata')
saveRDS(Val_SP_LII_predict, file='Validation _predict_SP_LII.Rds')
save(list=ls(), file = "SP_LII_all.Rdata")
