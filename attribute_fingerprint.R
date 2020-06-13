## plots fingerprint of attributes per sample
library(ggplot2)
## input
load(_clinical_attributes_) # ** requires input
# plot per sample
for(i in c(1:nrow(dataMat))){
  plotMat <- data.frame(cbind(names(dataMat), t(dataMat[i,])))
  names(plotMat) <- c("attribute", "value")
  p <- ggplot(plotMat, aes(x=attribute, y=value, group=1)) + geom_line()+
    coord_polar()+ ggtitle(paste0("patientID: ", rownames(dataMat)[i]))+
    theme(axis.text.x = 
            element_text(
              size=5,
              vjust=20,
              hjust=20,
              angle=-90 -360 / length(unique(plotMat$attribute)) * seq_along(plotMat$attribute)
            ))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line())
  print(p)
}
