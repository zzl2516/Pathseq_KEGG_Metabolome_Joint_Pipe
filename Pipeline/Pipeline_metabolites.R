### Pipeline for standard analyses of virome based on pathseq results
### Based on R v4.0.2
setwd("../")
library(vegan)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(ape)
library(picante)
library(dplyr)
library(multcomp)
library(plotrix)
library(RColorBrewer)
library(ggalluvial)
library(pheatmap)
library(VennDiagram)
library(RColorBrewer)
library(gplots)
library(ggrepel)
library(FactoMineR)
library(patchwork)
library(randomForest)
library(forcats)
library(indicspecies)
library(linkET)
library(igraph)
library(ggcor)
library(WGCNA)
library(MASS)
library(pROC)
library(mixOmics)
library(caret)

## Parameters
cbbPalette <- c("#B2182B","#56B4E9","#E69F00","#009E73","#F0E442","#0072B2",
                "#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#99999",
                "#ADD1E5")
group <- read.table(file = "Input/group.txt",header = FALSE,sep = "\t")
colnames(group) <- c("variable","Group")
group$Group <- factor(group$Group)
Group_numb <- length(unique(group[,2]))
Sample_numb <- length(unique(group[,1]))
wid <- ifelse(Group_numb < 4,1,0.8)
species.type <- "Metabolites"
dir.create("Results")
dir.create("Results/Metabolites")

### Input data
metabolites <- read.table("Input/metabolite_abundance.txt",header = FALSE,
                          sep = "\t",quote = "",row.names = 1)

### data processing
colnames(metabolites) <- metabolites[1,]
metabolites <- metabolites[-1,]
metabolites <- metabolites[metabolites$Compound_name != "unknown",]
rownames(metabolites) <- metabolites$Compound_name
metabolites <- metabolites[,-c(1:7)]
for (i in 1:ncol(metabolites)) {
    metabolites[,i] <- as.numeric(metabolites[,i])
}
metabolites[is.na(metabolites)] <- 0

## Beta diveristy
dir.create("Results/Metabolites/Beta_diversity")

### PCoA
dir.create("Results/Metabolites/Beta_diversity/Distance")
dir.create("Results/Metabolites/Beta_diversity/PCoA")
source("Functions/PCoA 3.R")
result <- pcoa.community(t(metabolites),group)
write.table(as.matrix(result[[1]]),
            "Results/Metabolites/Beta_diversity/Distance/distance_bray_curtis.txt",
            sep = "\t")
write.table(result[[2]],"Results/Metabolites/Beta_diversity/PCoA/pcoa.txt",
            sep = "\t",row.names = FALSE)
pdf(file = "Results/Metabolites/Beta_diversity/PCoA/Pcoa_group.pdf",
    width = 7.5,height = 5.4)
result[[3]]
dev.off()
pdf(file = "Results/Metabolites/Beta_diversity/PCoA/Pcoa_ellipse.pdf",
    width = 7.5,height = 5.4)
result[[4]]
dev.off()
pdf(file = "Results/Metabolites/Beta_diversity/PCoA/Pcoa_label.pdf",
    width = 7.5,height = 5.4)
result[[5]]
dev.off()

#### PCA
dir.create("Results/Metabolites/Beta_diversity/PCA")
source("Functions/PCA 3.R")
result <- pca.community(t(metabolites),group)
write.table(result[[1]],"Results/Metabolites/Beta_diversity/PCA/PCA.txt",
            sep = "\t",row.names = FALSE)
pdf(file = "Results/Metabolites/Beta_diversity/PCA/PCA_group.pdf",
    width = 7.5,height = 5.4)
result[[2]]
dev.off()
pdf(file = "Results/Metabolites/Beta_diversity/PCA/PCA_ellipse.pdf",
    width = 7.5,height = 5.4)
result[[3]]
dev.off()
pdf(file = "Results/Metabolites/Beta_diversity/PCA/PCA_label.pdf",
    width = 7.5,height = 5.4)
result[[4]]
dev.off()

#### NMDS
dir.create("Results/Metabolites/Beta_diversity/NMDS")
source("Functions/NMDS 3.R")
result <- nmds.community(t(metabolites),group)
write.table(result[[1]],"Results/Metabolites/Beta_diversity/NMDS/Microbe_stress.txt",
            sep = "\t",row.names = FALSE)
write.table(result[[2]],"Results/Metabolites/Beta_diversity/NMDS/NMDS.txt",
            sep = "\t",row.names = FALSE)
pdf(file = "Results/Metabolites/Beta_diversity/NMDS/NMDS_group.pdf",
    width = 7.5,height = 5.4)
result[[3]]
dev.off()
pdf(file = "Results/Metabolites/Beta_diversity/NMDS/NMDS_ellipse.pdf",
    width = 7.5,height = 5.4)
result[[4]]
dev.off()
pdf(file = "Results/Metabolites/Beta_diversity/NMDS/NMDS_label.pdf",
    width = 7.5,height = 5.4)
result[[5]]
dev.off()

### ANOSIM, Adonis, and MRPP tests
dir.create("Results/Metabolites/Beta_diversity/Beta_tests")
source("Functions/beta.test.R")
result <- beta.test(t(metabolites),group)
sink("Results/Metabolites/Beta_diversity/Beta_tests/Adonis.txt")
result[[1]]
sink("Results/Metabolites/Beta_diversity/Beta_tests/ANOSIM.txt")
result[[2]]
sink("Results/Metabolites/Beta_diversity/Beta_tests/MRPP.txt")
result[[3]]

### Hcluster
dir.create("Results/Metabolites/Beta_diversity/Hcluster")
dist_bray <- vegdist(t(metabolites), method = "bray")
hc <- hclust(dist_bray, method = "average")
hc_tre <- as.phylo(hc)
write.tree(hc_tre, "Results/Metabolites/Beta_diversity/Hcluster/Hcluster.tre")

source("Functions/plot.phylo2.r")
pdf("Results/Metabolites/Beta_diversity/Hcluster/Hcluster.pdf", 
    width = 7, height = 0.22*Sample_numb + 2)
plot.phylo2(hc_tre, type = "phylogram", cex = 1, adj = 1)
title(main = "UPGMA clustering tree", line = 2, cex.main = 1.5, font.main= 2)
dev.off()

### PLS-DA
dir.create("Results/Metabolites/Beta_diversity/PLS-DA")
metabolites.1 <- t(metabolites)
metabolites.1 <- metabolites.1[group$variable,]
plsda.data <- mixOmics::plsda(metabolites.1,group$Group,ncomp = 2)

pdf(file = "Results/Metabolites/Beta_diversity/PLS-DA/PLS-DA.pdf",
    width = 7.5,height = 5.4)
plotIndiv(plsda.data,ind.names = FALSE,ellipse = TRUE,title = "PLS-DA: Metabolites",
          X.label = "Component 1",Y.label =  "Component2",cex = 5,size.xlabel = 20,size.ylabel = 20,
          size.axis = 20,size.title = 24,pch = c(21,22),col.per.group = c("#B2182B","#56B4E9"),)
dev.off()

## Difference analysis
dir.create("Results/Metabolites/Diff_analysis")

## Wilcox sum-rank test
dir.create("Results/Metabolites/Diff_analysis/wilcox")
metabolites$compound <- rownames(metabolites)
metabolites <- metabolites[,c("compound",group$variable)]
source("Functions/wilcox 3.R")
source("Functions/wilcox.abun.R")
set.seed(1111)
result <- wilcox.biomarker(metabolites,group)
write.table(result[[1]],
            "Results/Metabolites/Diff_analysis/wilcox/wilcox_all.txt",
            sep = "\t",row.names = FALSE)
write.table(result[[2]],
            "Results/Metabolites/Diff_analysis/wilcox/wilcox_biomarker.txt",
            sep = "\t",row.names = FALSE)
if (nrow(result[[2]]) > 0) {
    metabolites.wilcox.biomarker <- result[[2]]
    result <- wilcox.abun(metabolites,metabolites.wilcox.biomarker,group)
    pdf("Results/Metabolites/Diff_analysis/wilcox/wilcox_biomarker_abun.pdf",
        width = 14,height = 1.8 + 0.18*nrow(metabolites.wilcox.biomarker))
    print(result[[1]])
    dev.off()
    write.table(result[[2]],"Results/Metabolites/Diff_analysis/wilcox/Group_mean_abun_biomarker.txt",sep = "\t",row.names = FALSE)
    write.table(result[[3]],"Results/Metabolites/Diff_analysis/wilcox/Abun_change_biomarker.txt",sep = "\t",row.names = FALSE)
}

## Random Forest for biomarker identification
dir.create("Results/Metabolites/Diff_analysis/RandomForest")
source("Functions/RandomForest.R")
set.seed(2222)
result <- Random.forest.biomarker(metabolites,group)
pdf("Results/Metabolites/Diff_analysis/RandomForest/classification.pdf",
    height = 3.5,width = 3.5)
result[[1]]
dev.off()
importance.metabolites <- result[[2]]
importance.metabolites$ID <- rownames(importance.metabolites)
importance.metabolites <- importance.metabolites[,c("ID",colnames(importance.metabolites)[1:(ncol(importance.metabolites)-1)])]
write.table(importance.metabolites,"Results/Metabolites/Diff_analysis/RandomForest/importance.txt",
            sep = "\t",row.names = FALSE)
if (nrow(importance.metabolites) > 0) {
    pdf("Results/Metabolites/Diff_analysis/RandomForest/RandomForest_biomarkder.pdf",
        height = nrow(importance.metabolites)*0.3 +1.85,width = 11)
    print(result[[3]])
    dev.off()
    write.table(result[[4]],"Results/Metabolites/Diff_analysis/RandomForest/biomarker_abundance.txt",
                sep = "\t")
}

### LEfSe
dir.create("Results/Metabolites/Diff_analysis/LEfSe")
source("Functions/lefse 2.R")

Lefse.result <- lefse.1(metabolites,group)
write.table(Lefse.result,"Results/Metabolites/Diff_analysis/LEfSe/Lefse_result.txt",
            sep = "\t")
if (nrow(Lefse.result) > 0) {
    source("Functions/lefse_bar.R")
    p <- lefse_bar(Lefse.result)
    pdf("Results/Metabolites/Diff_analysis/LEfSe/Lefse_bar.pdf",
        width = 6,height = 1.15 + 0.3*nrow(Lefse.result))
    print(p)
    dev.off()
}

## Diagnose
dir.create("Results/Metabolites/Diagnose")
source("Functions/diagnose.R")

### Wilcox rank-sum test
data <- read.table("Results/Metabolites/Diff_analysis/Wilcox/wilcox_biomarker.txt",
                   header = TRUE)
if (nrow(data) > 0) {
    dir.create("Results/Metabolites/Diagnose/Wilcox/")
    L1.result <- diagnose(data,metabolites,group)
    write.table(L1.result[[1]],"Results/Metabolites/Diagnose/Wilcox/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Metabolites/Diagnose/Wilcox/index.pdf",
        width = 2.4,height = 3.6)
    print(L1.result[[2]])
    dev.off()
    pdf("Results/Metabolites/Diagnose/Wilcox/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(L1.result[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
         print.thres.best.method = "youden",print.auc.cex = 2.5,print.thres.pch = 19,
         print.thres.pattern.cex = 2.5,auc.polygon = TRUE,auc.polygon.col = "grey90",
         print.auc.x = 0.35,print.auc.y = 0.1,auc.polygon.border = "black",
         max.auc.polygon = TRUE,axes = FALSE,xlab = "",ylab = "",mar = c(6,6,2,2)+.1)
    axis(side = 1,at = c(1.0,0.5,0),labels = c("","",""),lwd = 2,lwd.ticks = 4,line = -1)
    text(x = c(1.0,0.5,0),y = -0.08,labels = c("1.0","0.5","0.0"),cex = 2)
    mtext(side = 1,text = "Speciticity",line = 3,cex = 2.7)
    axis(side = 2,at = c(0,0.2,0.4,0.6,0.8,1.0),labels = c("","","","","",""),lwd = 2,lwd.ticks = 4,line = -2)
    text(x = 1.09,y = c(0,0.2,0.4,0.6,0.8,1.0),labels = c("0.0","0.2","0.4","0.6","0.8","1.0"),cex = 2)
    mtext(side = 2,text = "Sensitivity",line = 2,cex = 2.7)
    segments(x0 = 1.04,y0 = -0.04,x1 = -0.04,y1 = -0.04,lwd = 2)
    segments(x0 = 1.04,y0 = 1.04,x1 = -0.04,y1 = 1.04,lwd = 2)
    segments(x0 = 1.04,y0 = -0.04,x1 = 1.04,y1 = 1.04,lwd = 2)
    segments(x0 = -0.04,y0 = -0.04,x1 = -0.04,y1 = 1.04,lwd = 2)
    dev.off()
}

### RandomForest
data <- read.table("Results/Metabolites/Diff_analysis/RandomForest/importance.txt",
                   header = TRUE)
colnames(data)[1] <- "V1"
if (nrow(data) > 0) {
    dir.create("Results/Metabolites/Diagnose/RandomForest/")
    Phylum <- diagnose(data,metabolites,group)
    write.table(Phylum[[1]],"Results/Metabolites/Diagnose/RandomForest/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Metabolites/Diagnose/RandomForest/index.pdf",
        width = 2.4,height = 3.6)
    print(Phylum[[2]])
    dev.off()
    pdf("Results/Metabolites/Diagnose/RandomForest/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(Phylum[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
         print.thres.best.method = "youden",print.auc.cex = 2.5,print.thres.pch = 19,
         print.thres.pattern.cex = 2.5,auc.polygon = TRUE,auc.polygon.col = "grey90",
         print.auc.x = 0.35,print.auc.y = 0.1,auc.polygon.border = "black",
         max.auc.polygon = TRUE,axes = FALSE,xlab = "",ylab = "",mar = c(6,6,2,2)+.1)
    axis(side = 1,at = c(1.0,0.5,0),labels = c("","",""),lwd = 2,lwd.ticks = 4,line = -1)
    text(x = c(1.0,0.5,0),y = -0.08,labels = c("1.0","0.5","0.0"),cex = 2)
    mtext(side = 1,text = "Speciticity",line = 3,cex = 2.7)
    axis(side = 2,at = c(0,0.2,0.4,0.6,0.8,1.0),labels = c("","","","","",""),lwd = 2,lwd.ticks = 4,line = -2)
    text(x = 1.09,y = c(0,0.2,0.4,0.6,0.8,1.0),labels = c("0.0","0.2","0.4","0.6","0.8","1.0"),cex = 2)
    mtext(side = 2,text = "Sensitivity",line = 2,cex = 2.7)
    segments(x0 = 1.04,y0 = -0.04,x1 = -0.04,y1 = -0.04,lwd = 2)
    segments(x0 = 1.04,y0 = 1.04,x1 = -0.04,y1 = 1.04,lwd = 2)
    segments(x0 = 1.04,y0 = -0.04,x1 = 1.04,y1 = 1.04,lwd = 2)
    segments(x0 = -0.04,y0 = -0.04,x1 = -0.04,y1 = 1.04,lwd = 2)
    dev.off()
}

### LEfse
dir.create("Results/Metabolites/Diagnose/LEfSe")
data <- read.table("Results/Metabolites/Diff_analysis/LEfSe/Lefse_result.txt",
                   header = TRUE,row.names = 1)
data$V1 <- rownames(data)
if (nrow(data) > 0) {
    Species <- diagnose(data,metabolites,group)
    write.table(Species[[1]],"Results/Metabolites/Diagnose/LEfSe/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Metabolites/Diagnose/LEfSe/index.pdf",
        width = 2.4,height = 3.6)
    print(Species[[2]])
    dev.off()
    pdf("Results/Metabolites/Diagnose/LEfSe/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(Species[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
         print.thres.best.method = "youden",print.auc.cex = 2.5,print.thres.pch = 19,
         print.thres.pattern.cex = 2,auc.polygon = TRUE,auc.polygon.col = "grey90",
         print.auc.x = 0.35,print.auc.y = 0.1,auc.polygon.border = "black",
         max.auc.polygon = TRUE,axes = FALSE,xlab = "",ylab = "",mar = c(6,6,2,2)+.1)
    axis(side = 1,at = c(1.0,0.5,0),labels = c("","",""),lwd = 2,lwd.ticks = 4,line = -1)
    text(x = c(1.0,0.5,0),y = -0.08,labels = c("1.0","0.5","0.0"),cex = 2)
    mtext(side = 1,text = "Speciticity",line = 3,cex = 2.7)
    axis(side = 2,at = c(0,0.2,0.4,0.6,0.8,1.0),labels = c("","","","","",""),lwd = 2,lwd.ticks = 4,line = -2)
    text(x = 1.09,y = c(0,0.2,0.4,0.6,0.8,1.0),labels = c("0.0","0.2","0.4","0.6","0.8","1.0"),cex = 2)
    mtext(side = 2,text = "Sensitivity",line = 2,cex = 2.7)
    segments(x0 = 1.04,y0 = -0.04,x1 = -0.04,y1 = -0.04,lwd = 2)
    segments(x0 = 1.04,y0 = 1.04,x1 = -0.04,y1 = 1.04,lwd = 2)
    segments(x0 = 1.04,y0 = -0.04,x1 = 1.04,y1 = 1.04,lwd = 2)
    segments(x0 = -0.04,y0 = -0.04,x1 = -0.04,y1 = 1.04,lwd = 2)
    dev.off()
}

