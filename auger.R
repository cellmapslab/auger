# computes associations between genes and clinical attributes
# --max-ppsize=500000
library(glmnet)
library(openxlsx)
## input: contains imputed clinical phenotypes, expression data of lesional and nonlesional data
load(_data_here) # ** require input
modSuperList <- list()
matSuperList <- list()
## loop over each imputed dataset
for(iset in c(1:length(datasets))){
  
  modList <- list()
  matList <- list()
  ## imputed dataset
  dataMat<- datasets[[iset]]
## loop over each clinical phenotype
for(i in c(1:ncol(dataMat))){
  names(dataMat)[i] <- "myY"
  myY <- dataMat$myY
  dataMat <- dataMat[, -which(names(dataMat)%in%"myY")]
  ## identify the datatype of the clinical phenotype numeric/ordinal/categorical
  iDataType <- dataTypes$dataType[which(dataTypes$attributeNames%in%names(dataMat)[i])]
  ## identify the family type based on the datatype of clinical phenotype
  if(iDataType%in%"ordinal" | iDataType%in%"numeric"){
    fam.type="gaussian"
    myY <- as.numeric(as.matrix(myY))
  } else if(iDataType%in%"binary"){
    fam.type="binomial"
  } else if(iDataType%in%"categorical"){
    fam.type="multinomial"
  }
  ## lesional data model
  DdataMat <- data.frame(cbind(dataMat, lesMat))
  mod.mat <- model.matrix(~0+.,DdataMat)
  dMod <- cv.glmnet(x=mod.mat, y=myY, alpha=0.5, family=fam.type)
  dMat <- data.matrix(coef(dMod, dMod$lambda.1se))
  ## non lesional data model
  HdataMat <- data.frame(cbind(dataMat, nlesMat))
  mod.mat <- model.matrix(~0+.,HdataMat)
  hMod <- cv.glmnet(x=mod.mat, y=myY, alpha=0.5, family=fam.type)
  hMat <- data.matrix(coef(hMod, hMod$lambda.1se))
  ## combine the results
  mat <- data.frame(cbind(rownames(dMat), dMat, hMat))
  names(mat) <- c("var", "D", "H")
  ## save the results
  ##write.xlsx(mat, row.names=F, col.names=T, 
  ##           file=paste0(res.dir, Sys.Date(), "_elasticNet(0.5)_allgenes_imputatinset-", 
  ##                       iset, "_", names(impClinData.matched)[i], ".xlsx"))
  matList[[i]] <- mat
  modList[[i]] <- list(dMod, hMod)
}
modSuperList[[iset]] <- modList
matSuperList[[iset]] <- matList
  
}

save(matSuperList, modSuperList, file=_outfile_) # ** require input