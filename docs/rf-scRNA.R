#library(randomForest)
set.seed(100)

load('t3.Seurat.RData')
data1 <- as.data.frame(as.matrix(t3@data))
cut_matrix <- readRDS('cut_matrix.rds')

cut_matrix <- as.data.frame(cut_matrix)
names(cut_matrix) <- c("Two", "Four", "Eight", "Sixteen")

data1 <- as.data.frame(t(data1))
row.names(data1) <- gsub("-", ".", row.names(data1))
data1.test <- subset(data1, row.names(data1) %in% row.names(cut_matrix))
cut_matrix <- subset(cut_matrix, row.names(cut_matrix) %in% row.names(data1.test))
data1.cbind <- cbind(data1.test, cut_matrix)

fwrite(data1.cbind, 't3.expression.cutree.csv',row.names = T, col.names = T, quote=F, sep="\t")

data1.cbind <- read.table('t3.expression.cutree.csv', row.names=1)
names(data1.cbind) <- gsub("-", "_", names(data1.cbind))
train <- sample(nrow(data1.cbind), 0.7*nrow(data1.cbind), replace = FALSE)

TrainSet <- data1.cbind[train,]
ValidSet <- data1.cbind[-train,]

TrainSet.backup <- data1.cbind[train,]
ValidSet.backup <- data1.cbind[-train,]

TrainSet.backup$Two <- as.factor(TrainSet.backup$Two)
ValidSet.backup$Two <- as.factor(ValidSet.backup$Two)

fwrite(TrainSet.backup, 't3.trainset.csv', row.names = T, col.names = T, quote=F, sep="\t")
fwrite(ValidSet.backup, 't3.validset.csv', row.names = T, col.names = T, quote=F, sep="\t")

TrainSet$Two <- as.factor(TrainSet$Two)
TrainSet$Four <- as.factor(TrainSet$Four)
TrainSet$Eight <- as.factor(TrainSet$Eight)
TrainSet$Sixteen <- as.factor(TrainSet$Sixteen)

TrainSet <- TrainSet[,c(1:1000, 33203)]
ValidSet <- ValidSet[,c(1:1000, 33203)]
model1 <- randomForest(TrainSet$Two ~ ., data= TrainSet, ntree = 500, mtry = 6, importance = TRUE)
model1

predTrain <- predict(model1, TrainSet, type = "class")
table(predTrain, TrainSet$Two)

predValid <- predict(model1, ValidSet, type = "class")
mean(predValid == ValidSet$Two) #0.9239654
table(predValid,ValidSet$Two)

#predValid   1   2
#1 500  39
#2  40 460

importance(model1)
varImpPlot(model1)

# install.packages('inTrees')
# library(inTrees)
# treeList <- RF2List(model1) 
# X <- TrainSet[,c(1:1000)]
# target <- TrainSet[,c(1001)]
# exec <- extractRules(treeList, X) 
# exec[1:2,]
# ruleMetric <- getRuleMetric(exec,X,target)
# ruleMetric[1:2,]
# as.data.table(ruleMetric)
# ruleMetric <- pruneRule(ruleMetric, X, target)
# ruleMetric[1:2,]
# (ruleMetric <- selectRuleRRF(ruleMetric, X, target))
# (learner <- buildLearner(ruleMetric, X, target))
# readableRules <- presentRules(ruleMetric, colnames(X))
# readableRules[1:2,]


## trying with top 1000 variable genes instead of first 1000-----
topgenes <- t3@var.genes[1:1000]
data2 <- as.data.frame(as.matrix(t3@data))
data2 <- subset(data2, row.names(data2) %in% topgenes)
data2 <- as.data.frame(t(data2))
row.names(data2) <- gsub("-", ".", row.names(data2))
data2.test <- subset(data2, row.names(data2) %in% row.names(cut_matrix))
cut_matrix2 <- subset(cut_matrix, row.names(cut_matrix) %in% row.names(data2.test))
data2.cbind <- cbind(data2.test, cut_matrix2)
names(data2.cbind) <- gsub("-", "_", names(data2.cbind))
train2 <- sample(nrow(data2.cbind), 0.7*nrow(data2.cbind), replace = FALSE)
TrainSet2 <- data2.cbind[train2,]
ValidSet2 <- data2.cbind[-train2,]
TrainSet2$Two <- as.factor(TrainSet2$Two)
TrainSet2$Four <- as.factor(TrainSet2$Four)
TrainSet2$Eight <- as.factor(TrainSet2$Eight)
TrainSet2$Sixteen <- as.factor(TrainSet2$Sixteen)

TrainSet2 <- TrainSet2[,c(1:1001)]
ValidSet2 <- ValidSet2[,c(1:1001)]
model2 <- randomForest(TrainSet2$Two ~ ., data= TrainSet2, ntree = 500, mtry = 6, importance = TRUE)
model2
save(model2, file='model2.rf.RData')
importance(model2)
varImpPlot(model2)
imp.model2 <- as.data.frame(importance(model2))
imp.model2$gene <- row.names(imp.model2)
imp.model2 <- as.data.table(imp.model2)
imp.model2 <- imp.model2[order(-MeanDecreaseAccuracy)]
cat(imp.model2[c(1:100),gene], sep="\n")

randomForest::MDSplot(model2, TrainSet2$Two)

#MDSplot(iris.rf, iris$Species)

# plotting specific trees----
library(reprtree)
tr <- reprtree:::plot.getTree(model2)

model

tree <- getTree(model2, k=1, labelVar=TRUE)
realtree <- reprtree:::as.tree(tree, model)

# evalutate full model----
load('rf.fullmodel.RData') # breaks my session

# run RF on only variable genes minus the cell cycle genes, ribosomal genes ----

# on cluster:
options(expressions = 5e5)
library(randomforest)
library(randomForest)
library(data.table)

genes.keep <- as.data.frame(t3@var.genes)
genes.keep <- as.data.table(genes.keep[!(genes.keep$`t3@var.genes` %in% c(g1.s.genes, g2.m.genes, m.genes, s.genes)),])
names(genes.keep)<- "gene"
genes.keep <-as.data.table(genes.keep[-(grep('^RP', genes.keep$gene)),])
genes.keep <- as.vector(unlist(genes.keep$gene))
#fwrite(genes.keep, 't3.vargenes.csv', sep='\t', col.names = F, row.names = F, quote=F)
#genes.keep <- read.csv('t3.vargenes.csv', header=F)
#genes.keep <- unlist(genes.keep$V1, use.names=F)
library(Seurat)
load('t3.Seurat.RData')

cut_matrix <- readRDS('cut_matrix.rds')
cut_matrix <- as.data.frame(cut_matrix)
names(cut_matrix) <- c("Two", "Four", "Eight", "Sixteen")
cut_matrix <- readRDS('cut_matrix.rds')
cut_matrix <- as.data.frame(cut_matrix)
names(cut_matrix) <- c("Two", "Four", "Eight", "Sixteen")
cut_matrix
head(cut_matrix)
cut_matrix <- subset(cut_matrix, select=c('Two'))
#row.names(cut_matrix) <- gsub("-", ".", row.names(cut_matrix))

data1 <- as.data.frame(as.matrix(t3@data))
data1 <- as.data.frame(t(data1))
row.names(data1) <- gsub("-", ".", row.names(data1))
data1.test <- subset(data1, row.names(data1) %in% row.names(cut_matrix))
cut_matrix <- subset(cut_matrix, row.names(cut_matrix) %in% row.names(data1.test))
data1.cbind <- cbind(data1.test, cut_matrix)

data1.cbind.sub <- data1.cbind[,c(genes.keep, "Two")]
names(data1.cbind.sub) <- gsub("-", "_", names(data1.cbind.sub))
train <- sample(nrow(data1.cbind.sub), 0.7*nrow(data1.cbind.sub), replace = FALSE)
TrainSet <- data1.cbind.sub[train,]
ValidSet <- data1.cbind.sub[-train,]
TrainSet$Two <- as.factor(TrainSet$Two)

#> dim(TrainSet)
#[1] 2424 9186

#> dim(ValidSet)
#[1] 1039 9186
#

model <- randomForest(TrainSet$Two ~ ., data= TrainSet, ntree = 500, mtry = 6, importance = TRUE) #started 10:10 am september 12
model
save(model, file='model.vargenes.rf.RData')
