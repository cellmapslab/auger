# plots various acanthosis related plots
library(RColorBrewer)
library(tidyverse) # arrange
library(pheatmap)
library(colorRamps) # blue2yellow
library(DESeq2)
colset <- c(brewer.pal(n = 12, name = "Paired"), "#999999")
colset.attr <- brewer.pal(n = 3, name = "Greys")
## input
load(_data_here) # ** require input

## heatmap
aka.mat <- assay(vsdD[which(rownames(vsdD)%in%geneSets$Hist_Aka_quant), ]) ## expression data
annotMat <- data.frame(diagnosis=colData(vsdD)$newNames) ## diagnosis annotation
annotMat$acanthosis <- dataMat$Hist_Aka_quant ## add clinical phenotype data
rownames(annotMat) <- colnames(aka.mat)
names(colset) <- levels(annotMat$diagnosis)
names(colset.attr) <- levels(annotMat$acanthosis)
colset.list <- list(diagnosis = colset, acanthosis=colset.attr)
## plotting
pheatmap(aka.mat, color = blue2yellow(100), scale = "row", annotation_col = annotMat, fontsize_row = 5, annotation_colors = colset.list, cutree_rows = 2, cutree_cols = 3)

## circular barplot 
plotMat <- data.frame(figdiag=annotMat$diagnosis, aka=dataMat$Hist_Aka_quant, id=1:nrow(annotMat))
plotMat$aka <- as.numeric(as.matrix(plotMat$aka))
# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 1
to_add <- data.frame( matrix(NA, empty_bar*nlevels(plotMat$figdiag), ncol(plotMat)) )
colnames(to_add) <- colnames(plotMat)
to_add$figdiag <- rep(levels(plotMat$figdiag), each=empty_bar)
plotMat <- rbind(plotMat, to_add)
plotMat <- plotMat %>% arrange(figdiag)
plotMat$id <- seq(1, nrow(plotMat))
# Get the name and the y position of each label
label_data <- plotMat
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- plotMat %>% 
  group_by(figdiag) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

p <- ggplot(plotMat, aes(x=as.factor(id), y=aka, fill=figdiag)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=aka, fill=figdiag), stat="identity", alpha=0.5) +
  # customize colors
  scale_fill_manual(values=colset)+
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 3, xend = start, yend = 3), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 2, xend = start, yend = 2), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 1, xend = start, yend = 1), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(plotMat$id),4), y = c(0, 1, 2, 3), label = c("0", "1", "2", "3") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=aka, fill=figdiag), stat="identity", alpha=0.5) +
  ylim(-10,3) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  #geom_text(data=label_data, aes(x=id, y=value+10, label="", hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -1, xend = end, yend = -1), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -2, label=figdiag), colour = "black", alpha=0.8, size=2, fontface="bold", inherit.aes = FALSE)

## plotting
print(p)

## barplot
aka.pathnet <- do.call(rbind, lapply(1:nrow(aka.path), FUN=function(i) {
  x<- aka.path$symbols[i]; 
  data.frame(genes=gsub("\\s", "", unlist(strsplit(unlist(x), ","))), 
  GO_Term =aka.path$Var1[i])}))

aka.path$Var1 <- factor(aka.path$Var1, levels = aka.path$Var1[order(aka.path$value)])

p <- ggplot(data=aka.path, aes(x=Var1, y=value)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=value), hjust=-0.5, size=5)+
  theme_classic()+ theme(text = element_text(size=20))+
  xlab("")+ ylab("#genes")+coord_flip()

print(p)

