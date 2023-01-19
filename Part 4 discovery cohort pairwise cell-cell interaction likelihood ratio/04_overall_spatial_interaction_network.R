

##### R 4.1.0


rm(list=ls())
gc()


library(tidyverse) #‘1.3.2’
library(tidyr)  #‘1.2.0’
library(ggpubr)    # ‘0.4.0’
library(ggthemes)    #‘4.2.4’
library(pheatmap)    # ‘1.0.12’
library(corrplot)    #‘0.92’
library(psych)   #‘2.2.5’
library(reshape)  #‘0.8.9’
library(cowplot)  # ‘1.1.1’
library(gplots)   #‘3.1.3’
library(stringr)    # ‘1.4.0’
 


load(file = 'All_CN_TP_0424.Rdata')

cells_subtype <- All_CN_TP %>% dplyr::select(Class, Allsubtypes)
dat <- as.data.frame(with(cells_subtype, table(Class, Allsubtypes)))
cells_subtype_Freq <- spread(dat, Allsubtypes, Freq)
dat_subtype_freq <- cells_subtype_Freq[,-1]
row.names(dat_subtype_freq) <- cells_subtype_Freq$Class
dat_subtype_percent <- dat_subtype_freq / rowSums(dat_subtype_freq) * 100
row.names(dat_subtype_percent)<- cells_subtype_Freq$Class

MM<-corr.test(dat_subtype_percent, method="pearson")
temp <-  corrplot(MM$r,method = "square",order = 'hclust',
                p.mat = MM$p, sig.level = 0.05, insig = "label_sig")

#### Spacial Analysis ####
dat_SP_LII <- read.csv('./spatial analysis/spatialanalysis_output_subtype_SP_LII_class.csv', header = T,  stringsAsFactors = F, check.names = F)
dat_SP_HII <- read.csv('./spatial analysis/spatialanalysis_output_subtype_SP_HII_class.csv', header = T, stringsAsFactors = F, check.names = F)
dat_SP_MII_Tumor <- read.csv('./spatial analysis/spatialanalysis_output_subtype_SP_MII_Tumor_class.csv', header = T, stringsAsFactors = F, check.names = F)
dat_SP_MII_Fibro <- read.csv('./spatial analysis/spatialanalysis_output_subtype_SP_MII_Fibro_class.csv', header = T, stringsAsFactors = F, check.names = F)
dat_SP_Pf <- read.csv('./spatial analysis/spatialanalysis_output_subtype_SP_Pf_class.csv', header = T,stringsAsFactors = F, check.names = F)

functionX <- function(dat){
  dat <- dplyr::rename(dat, subtype=Allsubtypes)
  dat0 <- as.matrix(dat[c(-1,-2)])
  print(sum(dat0==Inf))
  dat0[dat0==Inf] <- 100000
  dat1 <- dat0/1000
  dat2 <- round(dat1,2)
  dat3 <- as.data.frame(dat2)
  dat4 <- cbind(dat[2],dat3)
  dat5 <- aggregate(. ~ subtype, dat4, sum)
  return(dat5)
}

dat_SP_LII_sum <- functionX(dat_SP_LII)
dat_SP_HII_sum <- functionX(dat_SP_HII)
dat_SP_MII_Tumor_sum <- functionX(dat_SP_MII_Tumor)
dat_SP_MII_Fibro_sum <- functionX(dat_SP_MII_Fibro)
dat_SP_Pf_sum <- functionX(dat_SP_Pf)

a <- dplyr::setdiff(dat_SP_HII_sum$subtype,dat_SP_LII_sum$subtype)
dat_SP_LII_sum$'Tumor_Glypican3+_Ki67+' <- 0

group <- rbind(dat_SP_LII_sum,dat_SP_HII_sum,dat_SP_MII_Tumor_sum,dat_SP_MII_Fibro_sum,dat_SP_Pf_sum)

group2 <- aggregate(. ~ subtype, group, sum)
# group2 %>% dplyr::group_by(subtype) %>% dplyr::summarise(across(everything(), sum))
row.names(group2) <- group2$subtype
group2 <- group2[,-1]
dim(group2)

group2 <- group2[order(rownames(group2)),order(colnames(group2))]

dat <- as.matrix(group2)
dat_filter <- dat
#dat_filter <- dat_filter[rowSums(dat_filter)>100,colSums(dat_filter)>100]
rowSums(dat_filter)
sum(dat_filter==Inf)
# dat_filter[dat_filter==Inf] <- 50000000     ##Tumor - Tumor 是Inf
dat_filter
#按likelihood ratios 
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
##画图
dat_LR <-dat_new


###########
mat <- t(round(dat_LR,2))

bk <- c(seq(-6,-0.1,by=0.01),seq(0,6,by=0.01))
pheatmap(mat,
         cellwidth = 15,
         cellheight = 15,
         border=FALSE,
         color=c(colorRampPalette(colors=c("blue","white"))(length(bk)/2),colorRampPalette(color=c("white","red"))(length(bk)/2)),
         scale = "none",
         cluster_rows = F,
         show_rownames = T,
         cluster_cols = F,
         legend_breaks=seq(-6,6,2),
         breaks=bk)

#### plot ####
temp1 = MM$r[rownames(temp$corr),colnames(temp$corr)]   
melt_corr = melt(temp1)
melt_corr$X1 <- factor(melt_corr$X1)
melt_corr$X2 <- factor(melt_corr$X2)

mat1 = mat[rownames(temp$corr),colnames(temp$corr)]
melted.dat_notTransposed <- melt(data.matrix(mat1))
melted.dat_notTransposed$X1 <- factor(melted.dat_notTransposed$X1)
melted.dat_notTransposed$X2 <- factor(melted.dat_notTransposed$X2)


p1 <-ggplot(melt_corr, aes(y = X1,
                           x = X2, fill=value)) +        
  geom_tile() +         
  scale_fill_gradient2(low = "#6aa84f", high = "#e69138", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation")+
  labs(x="",y="")+
  theme(legend.position = c(1.035,0.8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

melted.dat_notTransposed$color = ifelse(melted.dat_notTransposed$value >0, "#455655","#a1d8b1")
p2 <- ggplot(melted.dat_notTransposed,aes(y=X1,x=X2))+
  geom_point(aes(size = melted.dat_notTransposed$value), color = melted.dat_notTransposed$color)+   
  # scale_color_gradient("#455655" = "#455655",
  #                      "#a1d8b1" = "#a1d8b1",
  #                      name = "Likelihood\nRatio")+
  scale_size(range=c(-5,5),name = "Likelihood\nRatio")+
  labs(x="",y="")+
  theme(#legend.position = c(1.1,0.5),
    panel.background = element_rect(fill = "transparent"), 
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

p3 = ggdraw() +
  draw_plot(p1, 0, 0, 0.95, 1) +
  draw_plot(p2, 0, 0, 1, 1)+ggtitle('Overlay')



pdf('overlay_correlation_likelihood_0902.pdf',width = 20, height = 20)
p3
dev.off()

