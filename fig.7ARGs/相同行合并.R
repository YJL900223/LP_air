#####读入数据
resfam <- read.delim("resfam.count2.txt",header = T)


####相同行名合并
resfam_sum <- aggregate(. ~ Resfam,data=resfam,FUN="sum")
write.table(resfam_sum,"resfam_sum.txt",append = F,quote = F,sep = '\t',row.names = F,col.names = T)

#读入物种数据
otu <- read.delim('resfam_sum.txt', row.names = 1, header = T)

otu <- t(otu)
#head(otu)
#加载包
library(vegan)

##物种丰富度 Richness 指数
richness <- rowSums(otu > 0)


##Shannon（以下 Shannon 公式的对数底数均设使用 e，在 R 中即表示为 exp(1)）
#Shannon 指数
shannon_index <- diversity(otu, index = 'shannon', base = exp(1))


alpha_data <- data.frame(shannon_index,richness)#写入自己需要的指数
alpha_data[is.na(alpha_data)]<-0  
write.csv(alpha_data, file=paste("alpha_merged","_sample.csv",sep=""))

#读取数据
library(reshape2)
library(ggplot2)
library(ggpubr)
library(patchwork)
alpha <- read.csv('alpha_merged_sample.csv', stringsAsFactors = FALSE)
colnames(alpha)[1] <- "SampleID"
mapping <- read.table("metadata.txt",header = T)
alpha <- merge(alpha,mapping,by="SampleID")
alpha <- alpha[,1:4]
alpha$group <- factor(alpha$group)

alpha1 <- melt(alpha, id = c('SampleID', 'group'))


#多变量情况，添加分面的箱线图
my_comparisons <- list(c("G1", "G2"))

p1 <- ggplot(alpha1[1:33,], aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.5, size = 0.5) +
  
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank()) +
  labs(x = 'Groups', y = 'Shannon index')+stat_compare_means(comparisons = my_comparisons,
                                                label = 'p.signif')+theme_light()
p1
mean(alpha1$value[1:16])
mean(alpha1$value[17:33])
p2 <- ggplot(alpha1[34:66,], aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.5, size = 0.5) +
  
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank()) +
  labs(x = 'Groups', y = 'Detected number ')+stat_compare_means(comparisons = my_comparisons,
                                                label = 'p.signif',method = "t.test")+theme_light()
p2

library(pheatmap)
library(ggplot2)
data <- read.delim("resfam_sum1.txt",header = T,row.names = 1)
data1 <- as.data.frame(t(data))
data1$SampleID <- rownames(data1)
data2 <- merge(data1,mapping,by="SampleID")
data2 <- data2[,-c(18,20:25)]

data3 <- melt(data2,id.vars = c("SampleID","group","site2"))
max(data3$value)
min(data3$value)
mean(data3$value)


p3<- ggplot(data=data3,aes(x=site2,y=variable))+
  facet_wrap(~group)
p3 <- p3+geom_tile(alpha=0.5,aes(fill=value))+       #geom_tile绘制热图，使用value列作为颜色填充
  theme(axis.text.x=element_text(angle=0,vjust=0.5,hjust = 0.5,size=14))+
  scale_fill_gradient(name=NULL,low = "white",high = "red",limits=c(0,300000))+ theme_bw()+
  xlab(NULL) + ylab(NULL)
   

p3
p= (p1+p2) /
  p3
p

ggsave("Fig_resfam_air_result.pdf", p, 
       width = 8, height = 7)





