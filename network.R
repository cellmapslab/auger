# plots the network
library(igraph)
library(RCy3)
## input
load(_gene-attribute_association_) # ** require input
## sets which have a selected gene
geneSets <- isets[which(lapply(isets, FUN=function(x) nrow(x))>0)]
## convert to edge list
for(i in c(1:length(geneSets))){
  geneSets[[i]]$attribute <- names(geneSets)[i]
}
## edges
net <- do.call(rbind, geneSets)
names(net) <- c("target", "edgeweight", "source")
## node attribute
nodeattribute <- unique(reshape2::melt(t(net[, c("source", "target")]))[, -2])
names(nodeattribute) <- c("nodeType", "nodeID")
nodeAttr <- data.frame(nodeName=V(g)$name, nodeType=V(g)$type, row.names = V(g)$name, stringsAsFactors = F)
## visualization
g <- graph.data.frame(net[, c(1,3)], directed = F)
V(g)$type <- V(g)$name %in% net[, 3]

gc <- igraph.to.graphNEL(g)
createNetworkFromGraph(gc, title="phenotype-gene", collection="imcis")

loadTableData(nodeAttr)
## set node shapes
column <- 'nodeType'
values <- c("FALSE", "TRUE")
shapes <- c("ELLIPSE", "DIAMOND")
setNodeShapeMapping(column, values, shapes)
## set node sizes 
sizes <- rep(10, nrow(nodeAttr))
sizes [which(nodeAttr$nodeType=="TRUE")] <-100
lockNodeDimensions(FALSE)
setNodeSizeBypass(node.names=V(g)$name, new.sizes = sizes)
