#!/usr/local/bin/Rscript

#
# CTR_edw23_0003 ::: MeDiPS PCA
#
# 
# Copyright Russell S. Hamilton (rsh46@cam.ac.uk)
#


library('RColorBrewer')
library('ggplot2')
library("ggrepel")
library("cowplot")
library("matrixStats")
library("ggdendro")
library("ggforce")
library("ggalt")

wd <- "~/Documents/CTR-Groups/Erica_Watson/CTR_edw23_0003"
setwd(wd)

Project         <- "CTR_edw23_0003"
elementTextSize <- 8
WindowSize      <- 5
NumMostVariable <- 500

table        <- read.table(gzfile(paste0("CTR_edw23_0003_MeDIPS_Region_Coverage_Combo_", WindowSize,"Kb.table.gz")), header = TRUE, row.names=1)
table        <- table[!!rowSums(abs(table[])),]
names(table) <- gsub("^X", "", names(table))

sampleCondition <- c("+/gt",    "+/+",     "+/+",     "+/gt",    "+/gt",    "+/+",     "+/gt",    "+/gt", 
                     "+/gt",    "+/+",     "+/gt",    "+/gt",    "+/+",     "+/+",     "+/+",     "+/+",
                     "gt/gt",   "gt/gt",   "gt/gt",   "gt/gt",   "gt/gt",   "gt/gt",   "gt/gt",   "gt/gt",
                     "control", "control", "control", "control", "control", "control", "control", "control" ) 

message("+-------------------------------------------------------------------------------")
message("+ Create PCA Plots")
message("+-------------------------------------------------------------------------------")


comparison  <-  "All Samples"    
rv          <- rowVars(as.matrix(table))
select      <- order(rv, decreasing = TRUE)[seq_len(min(NumMostVariable, length(rv)))]
selection   <- paste0(NumMostVariable, "_MostVariable")
pca         <- prcomp(t(table), scale. = TRUE)
pc1var      <- round(summary(pca)$importance[2,1]*100, digits=1)
pc2var      <- round(summary(pca)$importance[2,2]*100, digits=1)
pc1lab      <- paste0("PC1 (",as.character(pc1var),"%)")
pc2lab      <- paste0("PC2 (",as.character(pc2var),"%)")
plot.tableA <- data.frame(names=colnames(table), pca$x, sampleCondition=sampleCondition)


#plot.tableA <- subset(plot.tableA, names!="C57_31")
#plot.tableA <- subset(plot.tableA, names!="173_1")

#C57Bl/6: green, ++: red ,  +gt: purple and gtgt: blue?


pltA.PCA    <- ggplot(plot.tableA, aes(x = PC1, y = PC2, col = (factor(sampleCondition))) ) +
              # geom_encircle(aes(group = factor(sampleCondition)), alpha=0.75, show.legend = FALSE) +
               geom_point(size = 4, alpha=.75 ) + geom_text_repel(aes(label=plot.tableA$names), size=4, nudge_x=1, nudge_y=1) +
               xlab(pc1lab) + ylab(pc2lab) + ggtitle(paste( "PCA ", comparison,  sep="")) +
               scale_colour_manual(name="Condition", values = c("gt/gt"="blue", "+/+"="red", "+/gt"="purple", "control"="#1B9E77")) +
               theme(text = element_text(size=elementTextSize)) 
pltA.PCA

pltA.PCA.nl    <- ggplot(plot.tableA, aes(x = PC1, y = PC2, col = (factor(sampleCondition))) ) +
 # geom_mark_ellipse(aes(fill = factor(sampleCondition)), alpha=0.1, show.legend = FALSE) +
  geom_encircle(aes(group = factor(sampleCondition)), alpha=0.75, show.legend = FALSE) +
  geom_point(size = 4, alpha=.75 ) + 
  xlab(pc1lab) + ylab(pc2lab) + ggtitle(paste( "PCA ", comparison,  sep="")) +
  scale_colour_manual(name="Condition", values = c("gt/gt"="blue", "+/+"="red", "+/gt"="purple", "control"="#1B9E77")) +
  scale_fill_manual(name="Condition", values = c("gt/gt"="blue", "+/+"="red", "+/gt"="purple", "control"="#1B9E77")) +
  theme(text = element_text(size=elementTextSize)) 
pltA.PCA.nl


pltA.hcl    <- ggdendrogram(hclust(dist(t(table[select,c(1,4,5,7,8,9,11,12,25:ncol(table))]))), rotate = FALSE, size = 2) + ggtitle(paste( "Hclust ", comparison,  sep=""))

loadings      <- as.data.frame(pca$rotation)
pca.1         <- loadings[ order(loadings$PC1,decreasing=TRUE), ]
pca.1.25      <- pca.1[c(1:25),]
pltA.pca.1.25 <- ggplot(data=pca.1.25, aes(x=factor(rownames(pca.1.25),levels=unique(rownames(pca.1.25))), y=PC1)) + 
  geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
pca.2         <- loadings[ order(loadings$PC2,decreasing=TRUE), ]
pca.2.25      <- pca.2[c(1:25),]
pltA.pca.2.25 <- ggplot(data=pca.2.25, aes(x=factor(rownames(pca.2.25),levels=unique(rownames(pca.2.25))), y=PC2)) + 
  geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 











comparison <-  "control vs +/gt"    # 1,4,5,7,8,9,11,12
rv          <- rowSums((table[,c(1,4,5,7,8,9,11,12,25:ncol(table))] - rowMeans(table[,c(1,4,5,7,8,9,11,12,25:ncol(table))]))^2)/(dim(table[,c(1,4,5,7,8,9,11,12,25:ncol(table))])[2] - 1)
select      <- order(rv, decreasing = TRUE)[seq_len(min(NumMostVariable, length(rv)))]
selection   <- paste0(NumMostVariable, "_MostVariable")
pca         <- prcomp(t(table[select,c(1,4,5,7,8,9,11,12,25:ncol(table))]), scale. = TRUE)
pc1var      <- round(summary(pca)$importance[2,1]*100, digits=1)
pc2var      <- round(summary(pca)$importance[2,2]*100, digits=1)
pc1lab      <- paste0("PC1 (",as.character(pc1var),"%)")
pc2lab      <- paste0("PC2 (",as.character(pc2var),"%)")
plot.table1 <- data.frame(names=colnames(table[,c(1,4,5,7,8,9,11,12,25:ncol(table))]), pca$x, sampleCondition=sampleCondition[c(1,4,5,7,8,9,11,12,25:ncol(table))])
plt1.PCA    <- ggplot(plot.table1, aes(x = PC1, y = PC2, col = (factor(sampleCondition))) ) +
               geom_point(size = 4, alpha=.75 ) + 
               geom_text_repel(aes(label=plot.table1$names), size=3, nudge_x=1, nudge_y=1) +
               xlab(pc1lab) + ylab(pc2lab) + ggtitle(paste( "PCA ", comparison,  sep="")) +
               scale_colour_manual(name="Condition", values = c("purple", "#1B9E77")) +
               theme(text = element_text(size=elementTextSize)) 
plt1.hcl    <- ggdendrogram(hclust(dist(t(table[select,c(1,4,5,7,8,9,11,12,25:ncol(table))]))), rotate = FALSE, size = 2) + ggtitle(paste( "Hclust ", comparison,  sep=""))

loadings      <- as.data.frame(pca$rotation)
pca.1         <- loadings[ order(loadings$PC1,decreasing=TRUE), ]
pca.1.25      <- pca.1[c(1:25),]
plt1.pca.1.25 <- ggplot(data=pca.1.25, aes(x=factor(rownames(pca.1.25),levels=unique(rownames(pca.1.25))), y=PC1)) + 
                 geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
pca.2         <- loadings[ order(loadings$PC2,decreasing=TRUE), ]
pca.2.25      <- pca.2[c(1:25),]
plt1.pca.2.25 <- ggplot(data=pca.2.25, aes(x=factor(rownames(pca.2.25),levels=unique(rownames(pca.2.25))), y=PC2)) + 
                 geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 



comparison  <- "control vs +/+"  # 2,3,6,10,13,14,15,16
rv          <- rowSums((table[,c(2,3,6,10,13,14,15,16,25:ncol(table))] - rowMeans(table[,c(2,3,6,10,13,14,15,16,25:ncol(table))]))^2)/(dim(table[,c(2,3,6,10,13,14,15,16,25:ncol(table))])[2] - 1)
select      <- order(rv, decreasing = TRUE)[seq_len(min(NumMostVariable, length(rv)))]
selection   <- paste0(NumMostVariable, "_MostVariable")
pca         <- prcomp(t(table[select,c(2,3,6,10,13,14,15,16,25:ncol(table))]), scale. = TRUE)
pc1var      <- round(summary(pca)$importance[2,1]*100, digits=1)
pc2var      <- round(summary(pca)$importance[2,2]*100, digits=1)
pc1lab      <- paste0("PC1 (",as.character(pc1var),"%)")
pc2lab      <- paste0("PC2 (",as.character(pc2var),"%)")
plot.table2 <- data.frame(names=colnames(table[,c(2,3,6,10,13,14,15,16,25:ncol(table))]), pca$x, sampleCondition=sampleCondition[c(2,3,6,10,13,14,15,16,25:ncol(table))])
plt2.PCA    <- ggplot(plot.table2, aes(x = PC1, y = PC2, col = (factor(sampleCondition))) ) +
               geom_point(size = 4, alpha=.75 ) + 
               geom_text_repel(aes(label=plot.table2$names), size=3, nudge_x=1, nudge_y=1) +
               xlab(pc1lab) + ylab(pc2lab) + ggtitle(paste( "PCA ", comparison,  sep="")) +
               scale_colour_manual(name="Condition", values = c("red", "#1B9E77")) +
               theme(text = element_text(size=elementTextSize)) 
plt2.hcl    <- ggdendrogram(hclust(dist(t(table[select,c(2,3,6,10,13,14,15,16,25:ncol(table))]))), rotate = FALSE, size = 2) + ggtitle(paste( "Hclust ", comparison,  sep=""))

loadings      <- as.data.frame(pca$rotation)
pca.1         <- loadings[ order(loadings$PC1,decreasing=TRUE), ]
pca.1.25      <- pca.1[c(1:25),]
plt2.pca.1.25 <- ggplot(data=pca.1.25, aes(x=factor(rownames(pca.1.25),levels=unique(rownames(pca.1.25))), y=PC1)) + 
                 geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
pca.2         <- loadings[ order(loadings$PC2,decreasing=TRUE), ]
pca.2.25      <- pca.2[c(1:25),]
plt2.pca.2.25 <- ggplot(data=pca.2.25, aes(x=factor(rownames(pca.2.25),levels=unique(rownames(pca.2.25))), y=PC2)) + 
                 geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


comparison  <-  "control vs gt/gt"      # 17-24
rv          <- rowSums((table[,c(17:24,25:ncol(table))] - rowMeans(table[,c(17:24,25:ncol(table))]))^2)/(dim(table[,c(17:24,25:ncol(table))])[2] - 1)
select      <- order(rv, decreasing = TRUE)[seq_len(min(NumMostVariable, length(rv)))]
selection   <- paste0(NumMostVariable, "_MostVariable")
pca         <- prcomp(t(table[select,c(17:24,25:ncol(table))]), scale. = TRUE)
pc1var      <- round(summary(pca)$importance[2,1]*100, digits=1)
pc2var      <- round(summary(pca)$importance[2,2]*100, digits=1)
pc1lab      <- paste0("PC1 (",as.character(pc1var),"%)")
pc2lab      <- paste0("PC2 (",as.character(pc2var),"%)")
plot.table3 <- data.frame(names=colnames(table[,c(17:24,25:ncol(table))]), pca$x, sampleCondition=sampleCondition[c(17:24,25:ncol(table))])
plt3.PCA    <- ggplot(plot.table3, aes(x = PC1, y = PC2, col = (factor(sampleCondition))) ) +
               geom_point(size = 4, alpha=.75 ) + 
               geom_text_repel(aes(label=plot.table3$names), size=3, nudge_x=1, nudge_y=1) +
               xlab(pc1lab) + ylab(pc2lab) + ggtitle(paste("PCA ", comparison,  sep="")) +
               scale_colour_manual(name="Condition", values = c("#1B9E77", "blue")) +
               theme(text = element_text(size=elementTextSize)) 
plt3.hcl    <- ggdendrogram(hclust(dist(t(table[select,c(17:24,25:ncol(table))]))), rotate = FALSE, size = 2) + ggtitle(paste( "Hclust ", comparison,  sep=""))

loadings      <- as.data.frame(pca$rotation)
pca.1         <- loadings[ order(loadings$PC1,decreasing=TRUE), ]
pca.1.25      <- pca.1[c(1:25),]
plt3.pca.1.25 <- ggplot(data=pca.1.25, aes(x=factor(rownames(pca.1.25),levels=unique(rownames(pca.1.25))), y=PC1)) + 
                 geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
pca.2         <- loadings[ order(loadings$PC2,decreasing=TRUE), ]
pca.2.25      <- pca.2[c(1:25),]
plt3.pca.2.25 <- ggplot(data=pca.2.25, aes(x=factor(rownames(pca.2.25),levels=unique(rownames(pca.2.25))), y=PC2)) + 
                 geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1))


#
# Experimental code to test different hclust clustering options
#
# table <- scale(table, center = TRUE, scale = TRUE)
# rv          = rowSums((table[,c(17:24,25:ncol(table))] - rowMeans(table[,c(17:24,25:ncol(table))]))^2)/(dim(table[,c(17:24,25:ncol(table))])[2] - 1)
# select      = order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
# ggdendrogram(hclust(dist(scale(t(table[select,c(17:24,25:ncol(table))])))), rotate = FALSE, size = 2) + ggtitle(paste( "Hclust ", comparison,  sep=""))
# 
# TBL <- table[select,c(17:24,25:ncol(table))]
# 
# sTBL <- t(scale(t(TBL)))
# str(sTBL, max.level = 0, give.attr = FALSE)
# 
# pr.dis <- dist(t(sTBL), method = "euclidean")
# 
# pr.hc.s <- hclust(pr.dis, method = "single")
# pr.hc.c <- hclust(pr.dis, method = "complete")
# pr.hc.a <- hclust(pr.dis, method = "average")
# pr.hc.w <- hclust(pr.dis, method = "ward")
# 
# op <- par(mar = c(0, 4, 4, 2), mfrow = c(2, 2))
# plot(pr.hc.s, labels = pr.hc.s$labels, main = "Single", xlab = "")
# plot(pr.hc.c, labels = pr.hc.c$labels, main = "Complete", xlab = "")
# plot(pr.hc.a, labels = pr.hc.a$labels, main = "Average", xlab = "")
# plot(pr.hc.w, labels = pr.hc.w$labels, main = "Ward", xlab = "")




#
# make a final plot
#

title <- ggdraw() + draw_label(paste(Project, " MeDIP \n Coverage on ", WindowSize, "Kb Genomic Windows, 500 Most Variable", sep=""), fontface='bold')
theme_set(theme_cowplot(font_size=10)) 

#px   <- plot_grid(plt1.PCA, plt1.hcl, plt1.pca.1.25, plt1.pca.2.25, 
#                  plt2.PCA, plt2.hcl, plt2.pca.1.25, plt2.pca.2.25,
#                  plt3.PCA, plt3.hcl, plt3.pca.1.25, plt3.pca.2.25,
#                  labels=c("A.1", "A.2", "A.3", "A.4", "B.1", "B.2", "B.3", "B.4", "C.1", "C.2", "C.3", "C.4"), ncol = 4, nrow = 3)


# Use this version if you only want to see the PCA, without the PC contributors
#px   <- plot_grid(plt1.PCA, plt1.hcl, plt2.PCA, plt2.hcl, plt3.PCA, plt3.hcl, 
#                  labels=c("A", "B", "C", "D", "E", "F"), ncol = 2, nrow = 3)

# Use this version if you only want to see the PCA, without the hclusts
px   <- plot_grid(plt1.PCA, plt2.PCA,  plt3.PCA, labels=c("B", "C", "D" ), ncol = 1, nrow = 3)

px2 <- plot_grid(pltA.PCA , px, labels=c("A", "" ), ncol = 2, nrow = 1, rel_widths = c(1,0.5))

pdf(paste(Project, "_PCA_Replacement_", WindowSize, "KbWindows_", selection, ".pdf", sep=""),width=12.5,height=7.5)
par(bg=NA)
plot_grid(title, px2, ncol=1, rel_heights=c(0.1, 1))
dev.off()
