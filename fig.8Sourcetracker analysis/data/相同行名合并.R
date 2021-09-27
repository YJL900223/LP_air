#####读入数据
LP <- read.table("data//otutab.txt",header = T)

tax <- read.table("data//taxonomy.txt",header = T)
otu_tax <- merge(LP,tax,by="OTUID")
otu_tax1 <- otu_tax[,-c(1,35,36,37,38,39,41)]
####相同行名合并
LP_sum <- aggregate(. ~ Genus,data=otu_tax1,FUN="sum")
write.table(LP_sum,"sum_g.txt",append = F,quote = F,sep = '\t',row.names = F,col.names = T)

park_sum <- aggregate(. ~ taxonomy,data=park,FUN="sum")
write.table(park_sum,"park_sum.txt",append = F,quote = F,sep = '\t',row.names = F,col.names = T)

home_plant_sum <- aggregate(. ~ taxonomy,data=home_plant,FUN="sum")
write.table(home_plant_sum,"home_plant_sum.txt",append = F,quote = F,sep = '\t',row.names = F,col.names = T)

class_outside_sum <- aggregate(. ~ taxonomy,data=class_outside,FUN="sum")
write.table(class_outside_sum,"class_outside_sum.txt",append = F,quote = F,sep = '\t',row.names = F,col.names = T)

ISS_sum <- aggregate(. ~ taxonomy,data=ISS,FUN="sum")
write.table(ISS_sum,"ISS_sum.txt",append = F,quote = F,sep = '\t',row.names = F,col.names = T)

######合并数据集

####按照by.x列合并，只保留要溯源的sink数据里所包含的物种进行溯源，sink不包含的不管
a <- merge(LP_sum,park_sum,all.x = T,all.y = T)
b <- merge(a,home_plant_sum,all.x = T,all.y = T)
c <- merge(b,class_outside_sum,all.x = T,all.y = T)
all <- merge(c,ISS_sum,all.x = T,all.y = T)
####去除NA值
all[is.na(all)] <-0
###保存数据
write.table(all,"differ_env.txt",append = F,quote = F,sep = '\t',row.names = F,col.names = T)


