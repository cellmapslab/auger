# plots transcriptomics heatmap
library(DESeq2)
library(RColorBrewer)
library(pheatmap)

## input
load(_data_here_) # ddsMat ** require input
load(_data_here) # genderGenes ** require input
## removal of samples with expression estimates with counts in less than 20% of cases.
RNAseq <- ddsMat[apply(counts(ddsMat),1,function(x) sum(x==0))<ncol(ddsMat)*0.8,]
## remove gender associated genes
geneNames <- data.frame(symbol=rownames(RNAseq))
RNAseq2 <- RNAseq[which(geneNames$symbol%in%genderGenes$symbol==F),]
## data tranformation
vsdD <-  varianceStabilizingTransformation(RNAseq2[, which(colData(RNAseq2)$sampleType=="D")], blind=F)
## 1000 most varrying genes for clustering
topVarGenes <- head(order(rowVars(assay(vsdD)), decreasing = TRUE), 1000)
mat  <- assay(vsdD)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
## annotations
annotMat <- data.frame(diagnosis=colData(vsdD)$newNames) 
rownames(annotMat) <- colnames(mat)
## colors of the annotation
colset <- c(brewer.pal(n = 12, name = "Paired"), "#999999")
names(colset) <- levels(annotMat$diagnosis)
colset.list <- list(diagnosis = colset)
## plot the clusters
pheatmap(mat, annotation_col = annotMat, clustering_distance_cols = "euclidean", annotation_colors = colset.list, labels_row = F, labels_col = F, clustering_method = "ward.D2")
