# computes variance partition and plots it
library(variancePartition)
library(limma)
library(DESeq2)
## input
load(_deseq_object_here_) ## ddsmat ** requires input
load(_attribute_data_here_) ## clinData ** requires input
## removal of samples with expression estimates with counts in less than 20% of cases.
RNAseq <- ddsMat[apply(counts(ddsMat),1,function(x) sum(x==0))<ncol(ddsMat)*0.8,]
## voom transformation
RNAseq.voom <- voom(counts(RNAseq))$E
## prepareSeqData
expData <- RNAseq.voom
sampleOrder <- sapply(rownames(dataMat), FUN=function(x) which(colData(RNAseq)$biobankID%in%x))
expData <- expData[, sampleOrder]
## diagnosis of the samples
figdiag <- sapply(colnames(expData), FUN=function(x)colData(RNAseq)$figdiag[which(rownames(colData(RNAseq))%in%x)])
## 
rownames(dataMat) <- colnames(expData)
dataMat$figdiag <- figdiag
dataMat <- droplevels(dataMat)
rownames(dataMat) <- colnames(expData)
dataMat$Diag_PGA <- as.numeric(as.matrix(dataMat$Diag_PGA))

# formula
form <- ~ _your_formula_here # ** requires input

## variance estimation
varFitSdiag <- fitExtractVarPartModel(expData, form, dataMat)
varfitMat <- data.frame(do.call(cbind, varFitSdiag))
## plotting
## Fig1F
attrVar <- rowSums(varfitMat[, -which(names(varfitMat)%in%c("figdiag", "Residuals"))])
plotMat <- data.frame(cbind(attrVar, varfitMat$figdiag))
names(plotMat) <- c("attributes", "diagnosis")
plotMat <- reshape2::melt(plotMat)
plotMat$variable <- factor(plotMat$variable,
                           levels = c('diagnosis','attributes'),ordered = TRUE)
p <- ggplot(plotMat, aes(x=variable, y=value)) + 
  geom_boxplot()+
  theme_bw() + theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(p)

## Fig1G
plotMat <- varfitMat[, -which(names(varfitMat)%in%c("Residuals"))]

## top 10 attributes
test <- colSums(plotMat)
top10 <- names(sort(test, decreasing = T))[1:10]
plotMat <- reshape2::melt(plotMat)
testMat <- plotMat[which(plotMat$variable%in%top10), ]
testMat$variable <- factor(testMat$variable, levels=top10, ordered = T)

p <- ggplot(testMat, aes(x=variable, y=value)) + 
  geom_boxplot()+
  theme_bw() + theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(p)

## all
plotMat <- varfitMat[, -which(names(varfitMat)%in%c("Residuals"))]
plotMat <- reshape2::melt(plotMat)

p <- ggplot(plotMat, aes(x=variable, y=value)) + 
  geom_boxplot()+
  theme_bw() + theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(p)