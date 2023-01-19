
#R version 4.2.1

library('dplyr') #1.0.10
library(cowplot)  #1.1.1
library(readr)  #2.1.3
library(tidyverse) #1.3.2
library(ggplot2) #3.3.6
library(PupillometryR) #0.0.4
library(gghalves) #0.1.3
library(reshape2)  #1.4.4
library(Hmisc) #4.7-1

###Read your data
oneslice_final=read.csv('B:\\codex_density\\oneslice_final_Macrophages_Vimentin+.csv')
colnames(oneslice_final)='score'
oneslice_final$group='P'
oneslice_final$lable='1'


oneslice_LM_final=read.csv('B:\\codex_density\\oneslice_LM_final_Macrophages_Vimentin+.csv')
colnames(oneslice_LM_final)='score'
oneslice_LM_final$group='LM'
oneslice_LM_final$lable='0'
total=rbind(oneslice_final,oneslice_LM_final)


a=data.frame(score=hdquantile(oneslice_final$score,seq(0,1,0.1)),lable='1')
a=a[-c(1,11),]
hdquantile(oneslice_final$score,seq(0,1,0.1),se=TRUE)
a$se=c(  2.85 , 4.22  ,4.68 , 5.47 , 6.92 , 9.10 ,10.28, 16.08, 29.11)
a$max=a$score+1.96*a$se
a$min=a$score-1.96*a$se
a$group=c(1:9)
b=data.frame(score=hdquantile(oneslice_LM_final$score,seq(0,1,0.1)),lable='0')
b=b[-c(1,11),]
hdquantile(oneslice_LM_final$score,seq(0,1,0.1),se=TRUE)
b$se=c( 4.56 , 6.24 , 7.79 , 8.41,  9.15  ,9.13, 10.19, 15.51, 28.63 )
b$max=b$score+1.96*b$se
b$min=b$score-1.96*b$se
b$group=c(1:9)


for (i in c(1:9)){
  if (a$score[i]>b$score[i]){
    if (a$min[i]>b$max[i]){
      a$result[i]='p < 0.05'
      b$result[i]='p < 0.05'
    }
    else{
      a$result[i]='ns'
      b$result[i]='ns'
    }
  }
  else{
    if (b$max[i]>a$max[i]){
      a$result[i]='p < 0.05'
      b$result[i]='p < 0.05'
    }
    else{
      a$result[i]='ns'
      b$result[i]='ns'
    }
  }
}


info_line2=rbind(a,b)

library(RColorBrewer)
library(ggsci)
cols=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)
scale_color_manual(values=c(pal_npg('nrc',alpha = 1.0)(5)))
col=scale_color_manual(values=c(pal_npg('nrc',alpha = 1.0)(5)))

options(scipen = 0)
options(digits =3)
a=ks.test(oneslice_final$score,oneslice_LM_final$score)[["p.value"]]
a='<2e-16'
b=wilcox.test(oneslice_final$score,oneslice_LM_final$score)[["p.value"]]
b='1.41e-16'

out=ggplot(info_line2,color=result)+
  geom_point(aes(as.numeric(lable),score,color=result), size = 3)+
  geom_line(aes(as.numeric(lable),score,group=group,color=result), size = 1)+
  geom_half_violin(aes(x= as.numeric(lable)+0.2,y=score,fill=group),side = "r",data = oneslice_final)+
  geom_half_violin(aes(x= as.numeric(lable)-0.2,y=score,fill=group),side = "l",data = oneslice_LM_final)+
  geom_boxplot(aes(x =as.numeric(lable)+0.08, y = score),outlier.shape = NA,
               alpha = 0.3, width = .05, colour = "BLACK",data = oneslice_final) +
  geom_boxplot(aes(x =as.numeric(lable)-0.08, y = score),outlier.shape = NA,
               alpha = 0.3, width = .05, colour = "BLACK",data=oneslice_LM_final)+
  geom_jitter(aes(x= as.numeric(lable)+0.15,y=score), size = .25,data = oneslice_final,width = .03)+
  geom_jitter(aes(x= as.numeric(lable)-0.15,y=score), size = .25,data = oneslice_LM_final,width = .03)+
  scale_color_manual(values=c('black'))+
  ylim(0,1600)+
  ylab('Distance')+xlab('Group')+theme_cowplot() +
  #theme(axis.ticks = element_blank(), axis.text.y = element_blank())+
  coord_flip()+geom_text(x=0.5,y=1300,label=paste("KS test:P.value ", a, sep ="")  ,size = 4)+
  geom_text(x=0.3,y=1300,label=paste("Wilcox test:P.value =", b, sep ="")  ,size = 4)+
  scale_x_continuous(breaks = NULL)+
  scale_fill_manual(values = c('#E2A4A1','#AFB3CD'))

pdf(file='B:\\codex_density\\P&LM\\result.pdf',height=3,width = 7)
print(out)
dev.off()

