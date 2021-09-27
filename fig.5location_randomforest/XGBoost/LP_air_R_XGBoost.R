#读取数据
genus = as.data.frame(t(read.delim("sum_g.txt",row.names=1,header = T)))
genus$SampleID <- rownames(genus)
design = read.delim("design.txt",header = T)
genus_design <- merge(genus,design[,c(1,2)],by="SampleID")
genus_design1 <- genus_design[,-1]
#将数据集分为训练集和测试集,比例为7:3
train_sub = sample(nrow(genus_design1),7/10*nrow(genus_design1))
train_data = genus_design1[train_sub,]
test_data = genus_design1[-train_sub,]

library(xgboost)
library(Matrix)
####训练集的数据预处理
# 将自变量转化为矩阵
traindata1 <- data.matrix(train_data[,c(1:193)]) 
# 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
traindata2 <- Matrix(traindata1,sparse=T)
train_data[,194][which(train_data$group=="G1")]= 0
train_data[,194][which(train_data$group=="G2")]= 1
traindata3 <- train_data[,194]
# 将自变量和因变量拼接为list
traindata4 <- list(data=traindata2,label=traindata3) 
# 构造模型需要的xgb.DMatrix对象，处理对象为稀疏矩阵
dtrain <- xgb.DMatrix(data = traindata4$data, label = traindata4$label) 

####测试集的数据预处理
# 将自变量转化为矩阵
testset1 <- data.matrix(test_data[,c(1:193)]) 
# 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
testset2 <- Matrix(testset1,sparse=T) 
# 将因变量转化为numeric
test_data[,194][which(test_data$group=="G1")]= 0
test_data[,194][which(test_data$group=="G2")]= 1
testset3 <- test_data[,194]
# 将自变量和因变量拼接为list
testset4 <- list(data=testset2,label=testset3) 
# 构造模型需要的xgb.DMatrix对象，处理对象为稀疏矩阵
dtest <- xgb.DMatrix(data = testset4$data, label = testset4$label)


xgb <- xgboost(data = dtrain,max_depth=6, eta=0.5,  objective='binary:logistic', nround=25)
#重要重要性排序 
library(Rcpp)
library(Ckmeans.1d.dp)
importance <- xgb.importance(traindata2@Dimnames[[2]], model = xgb)  
head(importance)
xgb.ggplot.importance(importance)
library(ggplot2)
importance$Feature=factor(importance$Feature,levels = importance$Feature)
p=ggplot(data = importance, mapping = aes(x=Feature,y=Gain,fill=Feature)) + 
  geom_bar(stat="identity")+coord_flip()+theme_classic()
p
ggsave(paste("imp_feature",".pdf", sep=""), p, width = 4, height =4)
library(ggplot2)
library(dplyr)
library(forcats)

# 通过另一列的值来对因子重新排序（升序）
importance %>%
  mutate(name = fct_reorder(Feature, Gain)) %>%
  ggplot(aes(x=Feature, y=Gain)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()

# 通过另一列的值来对因子重新排序（降序）
importance %>%
  mutate(name = fct_reorder(Feature, Gain, .desc=T)) %>%
  ggplot(aes(x=Feature, y=Gain)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
#在测试集上预测
pre_xgb = round(predict(xgb,newdata = dtest))
#输出混淆矩阵
table(test_data$group,pre_xgb,dnn=c("真实值","预测值"))
#install.packages("pROC")
library(pROC)
xgboost_roc <- roc(test_data$group,as.numeric(pre_xgb))
#绘制ROC曲线和AUC值
plot(xgboost_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='xgboost模型ROC曲线')

