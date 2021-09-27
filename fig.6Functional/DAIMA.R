
KEGG <- read.table("L2.txt",header = TRUE,sep = "\t")
#这里有一点需要说明，因为KEGG第一层分类中有几个分类的名字特别长，因此需要处理一下，使其能够实现自动换行。
library(ggplot2)
library(reshape2)
library(ggsci)
swr = function(string, nwrap = 12){
   paste(strwrap(string,width = nwrap),collapse = "\n")
}
swr = Vectorize(swr)
KEGG$L1 <- swr(KEGG$L1)
#画图
pathway2_plot <- ggplot(KEGG, aes(L2, value.sum.mean, fill = group)) +
   geom_col(position = 'dodge', width = 0.8, colour = 'black', size = 0.05) +
   geom_errorbar(aes(ymin = value.sum.mean - value.sum.sd, ymax = value.sum.mean + value.sum.sd), size = 0.05, width = .35, position = position_dodge(width = .8)) +  #?????????ߣ???ֵ?��?׼?
   #scale_fill_manual(values = c('red', 'blue')) +
   theme(legend.title = element_blank(), legend.position = c(0.9, 0.95)) +
   coord_flip() +  
   theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent',  color = 'black')) +
   geom_text(aes(label = '', y = value.sum.mean + value.sum.sd + 0.5), size = 4, position = position_dodge(0.8)) +
   labs(x = 'KEGG pathway2', y = 'Relative Abundance (%)')+
   facet_grid(L1~.,space = "free_y",scales = "free_y")+
   scale_fill_npg()
pathway2_plot
ggsave('L2_2.pdf', pathway2_plot, width = 8, height = 10)
