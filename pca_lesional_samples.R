## plots pca
library(DESeq2)
library(RColorBrewer)
library(ggplot2)
## input
load(_deseq_object_named_ddsMat_) #** requires input
## removal of samples with expression estimates with counts in less than 20% of cases.
RNAseq <- ddsMat[apply(counts(ddsMat),1,function(x) sum(x==0))<ncol(ddsMat)*0.8,]
## variance stabilization transformation
vsdD <-  varianceStabilizingTransformation(RNAseq[, which(colData(RNAseq)$sampleType=="D")], blind=F)
## colors for the diagnosis
n=12 #** requires input
colset <- c(brewer.pal(n = n, name = "Paired"), "#999999")

pcaData <- plotPCA(vsdD, intgroup = c("newNames"),ntop=1000) +
  scale_color_manual(values=colset) +
  labs(title="Lesional samples PCA")+
  theme_bw() + theme(legend.title=element_blank())+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

print(pcaData)
