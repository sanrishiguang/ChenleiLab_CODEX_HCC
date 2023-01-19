

### R 4.1.0

library(tidyverse)   #‘1.3.2’
library(tidyr)  #‘1.2.0’
library(ggpubr)  # ‘0.4.0’
library(ggthemes)  #‘4.2.4’
library(Rtsne)  # ‘0.16’
library(Rcpp)  #‘1.0.8.3’
library(igraph)  #‘1.3.2’
library(ggthemes)  # ‘4.2.4’
library(survival)   #‘3.3.1’
library(survminer)   # ‘0.4.9’
library(pheatmap)  # ‘1.0.12’
library(ggsci)  # ‘2.9’


mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                 legend.key = element_rect(fill = "white", colour = "white"), #图标
                 legend.background = (element_rect(colour= "white", fill = "white")))



load('import.Rdata')  

mydata <- mydata_pheno_anno

mydata$Metastasis <- mydata$ExtrahepaticMetastasis
mydata$Metastasis[mydata$LymphaticMetastasis == 1] <- 'P_LN'
mydata$Metastasis[mydata$ExtrahepaticMetastasis == 1] <- 'P_LM'
mydata$Metastasis[mydata$Metastasis != 'P_LN' & mydata$Metastasis != 'P_LM'] <- 'P'

mydata <- mydata %>% filter(Metastasis != 'P_LN')

mydata$Metastasis <- factor(mydata$Metastasis,c('P','P_LM'))
table(mydata$Metastasis)

(per_anno <- 55/(55+228))


#####Fig. 6B  Macrophages_vimentin+ percentage variation#######

temp <-  dplyr::select (mydata, Metastasis, Percentage = Macrophages_Vimentin_)
# p<- ggplot(temp,aes(x=Metastasis,y=Percentage))+
#   geom_boxplot(aes(color = Metastasis),outlier.colour=NA)+
#   geom_point(aes(color = Metastasis),position = 'jitter',size=0.5,alpha=0.5)+
#   #stat_compare_means(test="kruskal.test",label.x=1.5,label.y.npc='bottom',size=2,vjust=2)+
#   stat_compare_means(#comparisons = list(c('P','P_LM')), 
#                      aes(label="p.signif"),
#                      method = "wilcox.test",
#                      label = "p.signif")+
#   # geom_signif(comparisons = list(c('P','P_LM')), 
#   #             test = "wilcox.test",
#   #             map_signif_level=T,
#   #             step_increase = 0.1,
#   #             size = 0.5,
#   #             textsize = 2)+
#   labs(title='Macrophages_Vimentin+')+
#   mytheme+
#   scale_color_manual(values=rev(c(pal_npg('nrc',alpha = 1.0)(2))))+
#   theme(plot.title = element_text(size=10))+ 
#   theme(axis.text.x = element_blank())+ylim(0,2)
# p

ggplot(temp,aes(x=Metastasis,y=Percentage))+
  geom_boxplot(aes(color = Metastasis),outlier.colour=NA)+
  geom_point(aes(color = Metastasis),position = 'jitter',alpha=0.5)+
  stat_compare_means(#comparisons = list(c('P','P_LM')), 
    aes(label="p.signif"),
    method = "wilcox.test",
    label = "p.signif")+
  # geom_signif(comparisons = list(c('P','P_LM')), 
  #             test = "wilcox.test",
  #             map_signif_level=T,
  #             step_increase = 0.1) +
  labs(title='Mac_Vim+')+
  mytheme+
  scale_color_manual(values=rev(c(pal_npg('nrc',alpha = 1.0)(2))))+
  #theme(plot.title = element_text(size=15)) + 
  guides(color=F)+
  #coord_flip()+
  ylim(0, 2)


temp <-  dplyr::select (mydata, TNM_LH, Percentage = Macrophages_Vimentin_)
temp$TNM_LH = factor(temp$TNM_LH,levels = c("low","high"))
# ggplot(temp,aes(x=TNM_LH,y=Percentage))+
#   geom_boxplot(aes(color = TNM_LH),outlier.colour=NA)+
#   geom_point(aes(color = TNM_LH),position = 'jitter',size=0.5,alpha=0.5)+
#   # stat_compare_means(test="wilcox.test",label.x=1.5,label.y.npc='bottom',size=2,vjust=2)+
#   stat_compare_means(#comparisons = list(c('P','P_LM')), 
#     aes(label="p.signif"),
#     method = "wilcox.test",
#     label = "p.signif")+
#   # geom_signif(comparisons = list(c('low','high')), 
#   #             test = "wilcox.test",
#   #             map_signif_level=T,
#   #             step_increase = 0.1,
#   #             size = 0.5,
#   #             textsize = 2)+
#   mytheme+
#   scale_color_manual(values=rev(c(pal_npg('nrc',alpha = 1.0)(2))))+
#   theme(plot.title = element_text(size=10))+ 
#   theme(axis.text.x = element_blank())+ylim(0,2)




###TNM stage
ggplot(temp,aes(x=TNM_LH,y=Percentage))+
  geom_boxplot(aes(color = TNM_LH),outlier.colour=NA)+
  geom_point(aes(color = TNM_LH),position = 'jitter',alpha=0.5)+
  stat_compare_means(#comparisons = list(c('P','P_LM')), 
    aes(label="p.signif"),
    method = "wilcox.test",
    label = "p.signif")+
  # geom_signif(comparisons = list(c('P','P_LM')),
  #             test = "wilcox.test",
  #             map_signif_level=T,
  #             step_increase = 0.1) +
  labs(title='Mac_Vim+')+
  mytheme+
  scale_color_manual(values=rev(c(pal_npg('nrc',alpha = 1.0)(2))))+
  #theme(plot.title = element_text(size=15)) + 
  guides(color=F)+
  #coord_flip()+
  ylim(0, 2)



############Fig. 6A: survival###################

##In all samples
sur.cut <- surv_cutpoint(mydata_pheno_anno, time= 'Osday',event = 'OS01' , variables = 'Macrophages_Vimentin_', minprop = 0.5)
#summary(sur.cut)
sur.cat <- surv_categorize(sur.cut)
#head(sur.cat)
#table(sur.cat$AllCN32)
names(sur.cat) <- c("Osday",'OS01','Metastasis') 
fit <- survfit(Surv(Osday, OS01) ~ Metastasis, data = sur.cat)
ggsurvplot(fit,palette = "npg",
           risk.table = TRUE, pval = TRUE,
           conf.int = FALSE, xlab="Time in Days",
           title = 'Macrophages_Vimentin+ in all samples (OS)',
           ggtheme = theme_bw(),
           #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
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
           legend = "right")


sur.cut <- surv_cutpoint(mydata_pheno_anno, time= 'RFSday',event = 'RFS01' , variables = 'Macrophages_Vimentin_', minprop = 0.5)
#summary(sur.cut)
sur.cat <- surv_categorize(sur.cut)
#head(sur.cat)
#table(sur.cat$AllCN32)
names(sur.cat) <- c("RFSday",'RFS01','Metastasis') 
fit <- survfit(Surv(RFSday, RFS01) ~ Metastasis, data = sur.cat)
ggsurvplot(fit,palette = "npg",
           risk.table = TRUE, pval = TRUE,
           conf.int = FALSE, xlab="Time in Days",
           title = 'Macrophages_Vimentin+ in all samples (OS)',
           ggtheme = theme_bw(),
           #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
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
           legend = "right")


rm(list=ls())
gc()

# library(tidyverse)
# spat_class <- read.csv('./spatial analysis/spatialanalysis_output_subtype_class.csv')
# spat_class <- spat_class[-1]
# 
# group2 <- aggregate(. ~ Allsubtypes, spat_class, sum)
# row.names(group2) <- group2$Allsubtypes
# group2 <- group2[,-1]
# dim(group2)
# names(group2) <- gsub('[.]','_',colnames(group2))
# rownames(group2) <- gsub(' ','_', rownames(group2))
# rownames(group2) <- gsub('[+]','_', rownames(group2))
# rownames(group2) <- gsub('-','_', rownames(group2))
# 
# group2 <- group2[order(rownames(group2)),order(colnames(group2))]
# 
# dat <- as.matrix(group2)
# dat_filter <- dat
# #dat_filter <- dat_filter[rowSums(dat_filter)>100,colSums(dat_filter)>100]
# rowSums(dat_filter)
# sum(dat_filter==Inf)
# # dat_filter[dat_filter==Inf] <- 50000000     
# dat_filter
# #按likelihood ratios 
# dat_new<-matrix(nrow = length(rownames(dat_filter)), ncol=length(colnames(dat_filter)), dimnames = list(rownames(dat_filter),colnames(dat_filter)))
# 
# dat_sum <- sum(dat_filter)
# for (i in 1:length(rownames(dat_filter))) {
#   for (j in 1:length(colnames(dat_filter)))
#   {
#     a0 <- dat_sum/sum(dat_filter[i,])/sum(dat_filter[,j])
#     a<- dat_filter[i,j]*a0
#     b<- log(a+0.0001)
#     dat_new[i,j]<-b
#   }
# }
# 
# dat_LR <-as.data.frame (dat_new)
# dat_LR$Allsubtypes <- row.names(group2)
# save(dat_LR,file='./Figures/F19-1.Rdata')
# 







# Fig. 6C macrophage-vimentin interacting cells ---------------------------

rm(list=ls())
gc()

library(tidyverse)
spat_class <- read.csv('./spatial analysis/spatialanalysis_output_subtype_class.csv')
spat_class <- spat_class[-1]

group2 <- aggregate(. ~ Allsubtypes, spat_class, sum)
row.names(group2) <- group2$Allsubtypes
group2 <- group2[,-1]
dim(group2)
names(group2) <- gsub('[.]','_',colnames(group2))
rownames(group2) <- gsub(' ','_', rownames(group2))
rownames(group2) <- gsub('[+]','_', rownames(group2))
rownames(group2) <- gsub('-','_', rownames(group2))

group2 <- group2[order(rownames(group2)),order(colnames(group2))]

dat <- as.matrix(group2)
dat_filter <- dat
#dat_filter <- dat_filter[rowSums(dat_filter)>100,colSums(dat_filter)>100]
rowSums(dat_filter)
sum(dat_filter==Inf)
# dat_filter[dat_filter==Inf] <- 50000000    
dat_filter
#按likelihood ratios 公式进行计算
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

dat_LR <-as.data.frame (dat_new)
dat_LR$Allsubtypes <- row.names(group2)
#save(dat_LR,file='./Figures/F19-1.Rdata')


library(tidyverse)
library(ggsci)
#load('./Figures/F19-1.Rdata')

dat1 <- dat_LR %>% filter(Allsubtypes=="Macrophages_Vimentin_") %>% select(-Allsubtypes)
dat1_t <- as.data.frame(t(dat1))
dat1_t$Allsubtypes <- rownames(dat1_t)
dat1_t_top10 <- dat1_t[order(dat1_t$Macrophages_Vimentin_,decreasing = T),][2:50, ]

write.csv(dat1_t_top10, file='Table S6.csv')

dat1_t_top10$Color = ifelse(dat1_t_top10$Macrophages_Vimentin_>0, "H", "L")

ggplot(data=dat1_t_top10, mapping = aes(x=reorder(Allsubtypes,Macrophages_Vimentin_,decreasing =T),y=Macrophages_Vimentin_,fill=Color)) +
  geom_bar(stat = 'identity',colour = "grey")+
  #coord_flip()+
  guides (fill=F)+
  xlab('Cell Subtypes') +
  ylab('Likelihood Ratio') +
  scale_color_manual(values=c(pal_npg('nrc',alpha = 1.0)(2)))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5))

ggplot(data=dat1_t_top10, mapping = aes(x=reorder(Allsubtypes,Macrophages_Vimentin_),y=Macrophages_Vimentin_, fill=Macrophages_Vimentin_)) +
  geom_bar(stat = 'identity')+
  coord_flip()+
  guides (fill=F)+
  xlab('Cell Subtypes') +
  ylab('Likelihood Ratio')


##### Fig. 6D-E: likelihood ratio >0  for prognosis analysis ########

load('mydata_pheno_anno.Rdata')  


covariates <- paste0('Macrophages_Vimentin__',dat1_t_top10$Allsubtypes[dat1_t_top10$Macrophages_Vimentin_ > 0])

library(survival)
library(survminer)

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(RFSday, RFS01)~', x)))
univ_formulas


univ_models <- lapply( univ_formulas, function(x){coxph(x, data = mydata_pheno_anno)})
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
View(as.data.frame(res))

dysfx <- as.data.frame(res)
dysfx$index <- row.names(dysfx)
dysfx_sort <- dysfx[order(dysfx$HR,decreasing = T),]
dysfx_sort$index <- factor(dysfx_sort$index, levels = dysfx$index[order(dysfx$HR)])
dysfx_sort$Significance <- 'ns'
dysfx_sort$Significance [dysfx_sort$p.value <0.05] <- 'p < 0.05'
save(dysfx_sort, file='./Figures/F19-2.Rdata')


load(file='./Figures/F19-2.Rdata')
#load('./random forest/import.Rdata')
load('import.Rdata')


dysfx_sort$index = factor(dysfx_sort$index, levels = dysfx_sort$index)
##Harzard ratio
p <- ggplot(dysfx_sort, aes(x = index, y = HR,color=Significance)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymax = HR.confint.upper, ymin = HR.confint.lower))+
  scale_color_manual(values=rev(c(pal_npg('nrc',alpha = 1.0)(2))))+
  # coord_flip()+
  # ylim(0,20)+
  geom_hline(yintercept = 1,color = 'black')+
  # guides(color=F)+
  theme_bw() +
  theme(axis.text.x = element_text(angle=-90, hjust=0, vjust=0.5))
#coord_flip()
p


##survival curve for cox p < 0.05
cycle_variable <- as.character(dysfx_sort$index[dysfx_sort$p.value < 0.05])
cycle_variable= cycle_variable[1]

# for (i in 1:length(cycle_variable)) {
#   sur.cut <- surv_cutpoint(mydata_pheno_anno, time= 'RFSday',event = 'RFS01' , variables = cycle_variable[i], minprop = 0.3)
#   #summary(sur.cut)
#   sur.cat <- surv_categorize(sur.cut)
#   #head(sur.cat)
#   #table(sur.cat$AllCN32)
#   names(sur.cat) <- c("RFSday",'RFS01','group') 
#   fit <- survfit(Surv(RFSday, RFS01) ~ group, data = sur.cat)
#   filename <- paste0('./Figures/F19_Macrphage_Vimentin_spatial_interacting_cells/opti_RFS',"_",cycle_variable[i],".png")
#   png(file=filename,width=1500,height=1500,res=300)
#   p <- ggsurvplot(fit,palette = "npg",
#                   risk.table = TRUE, pval = TRUE,
#                   conf.int = FALSE, xlab="Time in Days",
#                   ggtheme = theme_classic(),
#                   title = cycle_variable[i])
#   print(p)
#   dev.off()
#   print(cycle_variable[i])
# }

sur.cut <- surv_cutpoint(mydata_pheno_anno, time= 'RFSday',event = 'RFS01' , variables = cycle_variable, minprop = 0.3)
#summary(sur.cut)
sur.cat <- surv_categorize(sur.cut)
#head(sur.cat)
#table(sur.cat$AllCN32)
names(sur.cat) <- c("RFSday",'RFS01','group') 
fit <- survfit(Surv(RFSday, RFS01) ~ group, data = sur.cat)
ggsurvplot(fit,
           palette = ggsci::pal_npg("nrc")(10)[c(5,9)],
           risk.table = TRUE, 
           pval = TRUE,
           #conf.int = FALSE, 
           xlab="Time in Days",
           ggtheme = theme_bw(),
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
           legend = "right",
           title = paste0(cycle_variable,"RFS"))



sur.cut <- surv_cutpoint(mydata_pheno_anno, time= 'Osday',event = 'OS01' , variables = cycle_variable, minprop = 0.3)
#summary(sur.cut)
sur.cat <- surv_categorize(sur.cut)
#head(sur.cat)
#table(sur.cat$AllCN32)
names(sur.cat) <- c("Osday",'OS01','group') 
fit <- survfit(Surv(Osday, OS01) ~ group, data = sur.cat)
ggsurvplot(fit,
           palette = ggsci::pal_npg("nrc")(10)[c(5,9)],
           risk.table = TRUE, 
           pval = TRUE,
           #conf.int = FALSE, 
           xlab="Time in Days",
           ggtheme = theme_bw(),
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
           legend = "right",
           title = paste0(cycle_variable,"OS"))



# Fig. 6F density plot ----------------------------------------------------

load(file='./Figures/F19-2.Rdata')
load('./random forest/import.Rdata') 


library(ggsignif)
mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), 
                 legend.key = element_rect(fill = "white", colour = "white"), #图标
                 legend.background = (element_rect(colour= "white", fill = "white")))


mydata <- mydata_pheno_anno
table(mydata$TNM_LH)
mydata$TNM_LH <- factor(mydata$TNM_LH,c('high','low'))

mydata$Metastasis <- mydata$ExtrahepaticMetastasis
mydata$Metastasis[mydata$LymphaticMetastasis == 1] <- 'P_LN'
mydata$Metastasis[mydata$ExtrahepaticMetastasis == 1] <- 'P_LM'
mydata$Metastasis[mydata$Metastasis != 'P_LN' & mydata$Metastasis != 'P_LM'] <- 'P'
mydata <- mydata %>% filter(Metastasis != 'P_LN')
mydata$Metastasis <- factor(mydata$Metastasis,c('P_LM','P'))
table(mydata$Metastasis)

temp <-  dplyr::select (mydata, Metastasis, LikelihoodRatio = 'CD4T_FOXP3__Macrophages_Vimentin_')
temp$LikelihoodRatio <- round(temp$LikelihoodRatio,digits = 2)   


ggplot(temp,aes(x=Metastasis,y=LikelihoodRatio))+
  geom_boxplot(aes(color = Metastasis),outlier.colour=NA)+
  geom_point(aes(color = Metastasis),position = 'jitter',size=0.5,alpha=0.5)+
  stat_compare_means(test="wilcox.test",label.x=1.5,label.y.npc='bottom',size=2,vjust=2)+
  geom_signif(comparisons = list(c('P','P_LM')), 
              test = "wilcox.test",
              map_signif_level=T,
              step_increase = 0.1,
              size = 0.5,
              textsize = 2)+
  labs(title='CD4T_FOXP3__Macrophages_Vimentin_')+
  mytheme+
  scale_color_manual(values=c(pal_npg('nrc',alpha = 1.0)(5)))+
  theme(plot.title = element_text(size=10))+ 
  theme(axis.text.x = element_blank())

temp$Metastasis = factor(temp$Metastasis, levels = c("P","P_LM"))
ggdensity(temp, x='LikelihoodRatio', color='Metastasis', fill = 'Metastasis', palette = ggsci::pal_npg("nrc")(10)[c(9,5)])  


temp <-  dplyr::select (mydata, TNM_LH, LikelihoodRatio = 'CD4T_FOXP3__Macrophages_Vimentin_')
temp$LikelihoodRatio <- round(temp$LikelihoodRatio,digits = 2)   


ggplot(temp,aes(x=TNM_LH,y=LikelihoodRatio))+
  geom_boxplot(aes(color = TNM_LH),outlier.colour=NA)+
  geom_point(aes(color = TNM_LH),position = 'jitter',size=0.5,alpha=0.5)+
  stat_compare_means(test="wilcox.test",label.x=1.5,label.y.npc='bottom',size=2,vjust=2)+
  geom_signif(comparisons = list(c('high','low')), 
              test = "wilcox.test",
              map_signif_level=T,
              step_increase = 0.1,
              size = 0.5,
              textsize = 2)+
  labs(title='CD4T_FOXP3__Macrophages_Vimentin_')+
  mytheme+
  scale_color_manual(values=c(pal_npg('nrc',alpha = 1.0)(5)))+
  theme(plot.title = element_text(size=10))+ 
  theme(axis.text.x = element_blank())

library(ggpubr)
temp$TNM_LH = factor(temp$TNM_LH, levels = c("low","high"))
ggdensity(temp, x='LikelihoodRatio', color='TNM_LH', fill = 'TNM_LH', palette = ggsci::pal_npg("nrc")(10)[c(9,5)]) 













