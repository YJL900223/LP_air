# Fig.2a
#读入物种数据
otu <- read.delim('differ_env.txt', row.names = 1, header = T)

otu <- t(otu)
#head(otu)
#加载包
library(vegan)

##物种丰富度 Richness 指数
#richness <- rowSums(otu > 0)
#或
richness <- estimateR(otu)[1, ]

##Shannon（以下 Shannon 公式的对数底数均设使用 e，在 R 中即表示为 exp(1)）
#Shannon 指数
shannon_index <- diversity(otu, index = 'shannon', base = exp(1))

#Shannon 多样性
shannon_diversity <- exp(1)^shannon_index

#Shannon 均匀度（Pielou 均匀度）
pielou <- shannon_index / log(richness, exp(1))

##Simpson
#Gini-Simpson 指数（我们平时常用的 Simpson 指数即为 Gini-Simpson 指数）
gini_simpson_index <- diversity(otu, index = 'simpson')

#经典 Simpson 指数（使用频率比较低）
simpson_index <- 1 - gini_simpson_index

#Invsimpson 指数（Gini-Simpson 的倒数）
invsimpson_index <- 1 / gini_simpson_index
#或
invsimpson_index <- diversity(otu, index = 'invsimpson')

#Simpson 多样性
simpson_diversity <- 1 / (1 - gini_simpson_index)

#Simpson 均匀度（equitability 均匀度）
equitability <- 1 / (richness * (1 - gini_simpson_index))

##Chao1 & ACE
#Chao1 指数
chao1 <- estimateR(otu)[2, ]

#ACE 指数
ace <- estimateR(otu)[4, ]

##goods_coverage 指数
goods_coverage <- 1 - rowSums(otu == 1) / rowSums(otu)

alpha_data <- data.frame(pielou,shannon_index,richness,gini_simpson_index)#写入自己需要的指数
alpha_data[is.na(alpha_data)]<-0  
write.csv(alpha_data, file=paste("alpha_merged","_sample.csv",sep=""))

#读取数据
library(reshape2)
library(ggplot2)
library(ggpubr)
alpha <- read.csv('alpha_merged_sample.csv', stringsAsFactors = FALSE)
colnames(alpha)[1] <- "SampleID"
mapping <- read.table("differ_design.txt",header = T)
alpha <- merge(alpha,mapping,by="SampleID")
alpha$group <- factor(alpha$group)

alpha1 <- melt(alpha, id = c('SampleID', 'group'))


#多变量情况，添加分面的箱线图
my_comparisons <- list(c("ISS", "LP"),c("park","LP"),c("home_plant","LP"),c("room","LP"),c("outside","LP"))

p1 <- ggplot(alpha1, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.5, size = 0.5) +
  facet_wrap(~variable, 2, scales = 'free') +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank()) +
  labs(x = 'Groups', y = '')+stat_compare_means(comparisons = my_comparisons,
                                                label = 'p.signif')+theme_light()
p1
#ggsave(paste("fig2a_alpha_rdp", ".pdf", sep=""), p1, width =10, height = 8)


# Fig.2b
library(vegan)

##读取 OTU 丰度表，排序
otu <- read.delim('differ_env.txt', row.names = 1,header = T)
otu <- data.frame(t(otu))

#物种数据 Hellinger 预转化（处理包含很多 0 值的群落物种数据时，推荐使用）
otu<- decostand(otu, method = 'hellinger')

#根据物种组成计算样方距离，如 Bray-curtis "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" or "mahalanobis".距离，详情 ?vegdist
bray_dis <- vegdist(otu, method = 'bray')      #结果以 dist 数据类型存储

#输出距离矩阵
#write.table(as.matrix(bray_dis), 'bray_distance.txt', sep = '\t', col.names = NA, quote = FALSE)

#PCoA 排序，详情 ?cmdscale
pcoa <- cmdscale(bray_dis, k = (nrow(otu) - 1), eig = TRUE)




#各 PCoA 轴的特征值
pcoa_eig <- pcoa$eig
#先评估下负特征值（末尾几个轴）
barplot(pcoa_eig)



#各 PCoA 轴的解释量
pcoa_exp <- pcoa$eig/sum(pcoa$eig)

#样方排序坐标
site <- pcoa$point
#或
site <- scores(pcoa)

#断棍模型认为前 4 轴特征值具有代表性，这里就暂且先映射前 4 轴
species <- wascores(pcoa$points[,1:4], otu)


#断棍模型认为前 4 轴特征值具有代表性，这里就展示前 4 轴
#前 4 轴解释量
pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')
pcoa3 <- paste('PCoA axis3 :', round(100*pcoa_exp[3], 2), '%')
pcoa4 <- paste('PCoA axis4 :', round(100*pcoa_exp[4], 2), '%')



##################################ggplot2 作图，已知第3、4 轴无法区分样方差异，就只观测前 2 轴
#添加分组信息
site <- data.frame(pcoa$point)[1:2]
site$name <- rownames(site)
group <- read.table("differ_design.txt",header = T)
colnames(group)[1] <- "name"
#为样本点坐标添加分组信息
site <- merge(site, group, by = 'name', all.x = TRUE)

# species_top10 <- data.frame(species[abundance_top10,1:2])
# species_top10$name <- rownames(species_top10)

#ggplot2 作图
library(ggplot2)
library(ggsci)

dis <- as.matrix(bray_dis)

adonis_result <- adonis(dis~group, group, permutations = 999)
adonis_result
p2 = ggplot(site, aes(x=X1, y=X2)) + geom_point(alpha=.7, size=2,aes(color=group)) +scale_fill_npg()+
  labs(x=paste("PCoA 1 (", format(100 * pcoa_eig[1] / sum(pcoa_eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * pcoa_eig[2] / sum(pcoa_eig), digits=4), "%)", sep=""),
       title=paste("Bray_Curtis"," PCoA"," R^2=0.36 P=0.001",sep=""))  + 
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.6, show.legend = FALSE) +theme_light()
p2


#ggsave('Bray_Curtis3.pdf', p2, width = 5.5, height = 5.5)

# Fig.2c
library("reshape2")
library("ggplot2")
library("digest")
library("ggrepel")
library("ggpubr")
# 1. 读取输入文件

# # 读取样品分类学文件
tax_sample = read.table("tax_6Genus.txt", header=T, row.names= 1, sep="\t", comment.char="") 

# 读取实验设计
design = read.table("differ_design.txt", header=T, row.names= 1, sep="\t", comment.char="") 

# 提取样品组信息,默认为group可指定
sampFile = data.frame(group=design[,"group"],
                      sample=row.names(design), 
                      row.names = row.names(design))

# 数据筛选，筛选两文件中共有
idx = rownames(sampFile) %in% colnames(tax_sample) # match design with alpha
sampFile = sampFile[idx,]
tax_sample = tax_sample[,rownames(sampFile)] 

# 2. 统计与绘图

# 按丰度降序排序
mean_sort = tax_sample[(order(-rowSums(tax_sample))), ]
mean_sort = as.data.frame(mean_sort)

# 再次检验计算是否出错
# colSums(mean_sort)

# 保存变量备份，并输出至文件
merge_tax=mean_sort[c(1:10),]


# 按组合并求均值

# 转置样品名添加组名，并去除多余的两个样品列
mat_t = t(merge_tax)
mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,c(-1,-3)]

# 按组求均值，转置，再添加列名
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$group
colnames(mat_mean_final) = geno

# 保存变量备份，并输出至文件
mean_sort=as.data.frame(mat_mean_final)
# write.table("\t", file=paste(opts$output,"_group.txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
# write.table(mean_sort, file=paste(opts$output,"_group.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)

# 数据转换长表格并绘图
mean_sort$tax = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("tax")))
# 设置分类学顺序，默认字母，可选丰度或手动
# data_all$tax  = factor(data_all$tax, levels=rownames(mean_sort))   
data_all$value <- as.numeric(data_all$value)
p3 = ggplot(data_all, aes(x=variable, y = value, fill = tax )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + scale_fill_brewer(palette="Paired")+
  xlab("Groups")+ylab("Percentage (%)")+ theme_classic()+theme_light()
p3

# 保存pdf和png格式方便查看和编辑
#ggsave(paste("stackplot_genus", "_group.pdf", sep=""), p3, width = 10, height = 8)

# Fig.2d
data_vn <- read.csv('vennDiagram.csv',header = T)
A <- as.vector(as.matrix(data_vn[which(data_vn$sample=="LP"),]['taxonomy']))
B <- as.vector(as.matrix(data_vn[which(data_vn$sample=="ISS"),]['taxonomy']))
C <- as.vector(as.matrix(data_vn[which(data_vn$sample=="park"),]['taxonomy']))
D <- as.vector(as.matrix(data_vn[which(data_vn$sample=="room"),]['taxonomy']))
E <- as.vector(as.matrix(data_vn[which(data_vn$sample=="home_plant"),]['taxonomy']))
#install.packages('ggvenn') # 安装

library(ggvenn)#导入
x = list('LP' = A,'ISS' = B,'park' = C,'room' = D,'home_plant'= E)
p4 <- ggvenn(x,
             show_percentage = F,
             #fill_color = c('#F8766D', '#7CAE00', '#00BFC4','#C77CFF','green'),
             fill_alpha = 0.5,
             
             stroke_color = "black",
             stroke_alpha = 1,
             stroke_size = 1,
             stroke_linetype = "solid",
             #set_name_color = c('#F8766D', '#7CAE00', '#00BFC4','#C77CFF','green'),
             set_name_size = 6,
             text_color = "black",
             text_size = 5
)
p4

#组合图
#devtools::install_github("thomasp85/patchwork")

library(patchwork)
p=p1 /
  (p2+p3)+plot_annotation(tag_levels = "a")+plot_layout(widths = c(2,2))
p
ggsave(paste("Figure", "3.pdf", sep=""), p, width = 12, height = 10)




