#Read in file
design <- read.table("design.txt",header = T)
qPCR_air <- read.table("qPCR_air.txt",header = T)
qPCR_design <- merge(design,qPCR_air,by="SampleID")
write.table(qPCR_design,"qPCR_design.txt",append = F,quote = F,sep = '\t',row.names = F,col.names = T)

qPCR_design$number <- log10((qPCR_design$copynumber_16S))
qPCR_design$number2 <- (qPCR_design$copynumber_16S)
qPCR_mean <- format(mean(qPCR_design$number2),scientific=TRUE)
G1_mean <- format(mean(qPCR_design[1:16,]$number2),scientific=TRUE)
G2_mean <- format(mean(qPCR_design[17:33,]$number2),scientific=TRUE)

library(ggplot2)
library(ggsci)
library(dplyr) 
library(ggpubr)
library(patchwork)

shapiro.test(qPCR_design$number)
bartlett.test(number~group ,data=qPCR_design)
bartlett.test(number~location ,data=qPCR_design)
bartlett.test(number~month ,data=qPCR_design)

my_comparisons1 <- list( c("G1", "G2"))

p1 <- ggplot(qPCR_design,aes(x=group,y=number))+
  geom_boxplot(aes(fill = group)) +
  geom_jitter(aes(shape = group), 
    position = position_jitter(0),
    size = 3)+
  scale_color_npg()+
  stat_compare_means(comparisons = my_comparisons1,method = "t.test",label.y = 13)+
  ylab(expression(paste("16S rRNA gene copy number/ ",m^3)))+
  xlab("Group")+
  theme_bw()
p1

my_comparisons2 <-list( c("CC", "PC"),c("PC", "SC"),c("CC", "SC")) 
p2 <- ggplot(qPCR_design,aes(x=location,y=number))+
  geom_boxplot(aes(fill = location)) +
  facet_wrap(~group)+
  geom_jitter(aes(shape = location), 
              position = position_jitter(0),
              size = 3)+
  scale_color_npg()+
  stat_compare_means(comparisons = my_comparisons2,method = "t.test")+
  ylab(expression(paste("16S rRNA gene copy number/ ",m^3)))+
  xlab("Location")+
  theme_bw()
p2
qPCR_design$month <- factor(qPCR_design$month,levels = c("3","4","5","10","11","12"))
p3 <- ggline(qPCR_design, x = "month", y = "number", point.size = 2,
       add = c("mean_sd", "jitter","boxplot"),color = "group")+
  ylab(expression(paste("16S rRNA gene copy number/ ",m^3)))+
  xlab("Month")+scale_color_npg()+
  theme_bw()
  
p3
p=p2 / (p1+p3)
p
ggsave("Fig_qPCR_air_result.pdf", p, 
       width = 8, height = 7)
