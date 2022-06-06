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
species.type <- "Gene"
dir.create("Results")
dir.create("Results/Gene")

### Input data
L1 <- read.table(file = "Input/kegg.profile.L1.xls",header = FALSE,row.names = 1,
                 sep = "\t",quote = "")
L2 <- read.table(file = "Input/kegg.profile.L2.xls",header = FALSE,row.names = 1,
                 sep = "\t",quote = "")
L3 <- read.table(file = "Input/kegg.profile.L3.xls",header = FALSE,row.names = 1,
                 sep = "\t",quote = "")
ko <- read.table(file = "Input/kegg.profile.entry.xls",header = FALSE,row.names = 1,
                 sep = "\t",quote = "")

### data processing
colnames(L1) <- L1[1,]
L1 <- L1[-1,]
for (i in 1:ncol(L1)) {
    L1[,i] <- as.numeric(L1[,i])
}
L1 <- L1[,group$variable]
L1$Description <- rownames(L1)
L1[,1:(ncol(L1)-1)] <- t(t(L1[,1:(ncol(L1)-1)])/colSums(L1[,1:(ncol(L1)-1)])*100)
L1 <- L1[,c("Description",colnames(L1)[1:(ncol(L1)-1)])]

colnames(L2) <- L2[1,]
L2 <- L2[-1,]
for (i in 1:ncol(L2)) {
    L2[,i] <- as.numeric(L2[,i])
}
L2 <- L2[,group$variable]
L2$Description <- rownames(L2)
L2[,1:(ncol(L2)-1)] <- t(t(L2[,1:(ncol(L2)-1)])/colSums(L2[,1:(ncol(L2)-1)])*100)
L2 <- L2[,c("Description",colnames(L2)[1:(ncol(L2)-1)])]

colnames(L3) <- L3[1,]
L3 <- L3[-1,]
for (i in 1:ncol(L3)) {
    L3[,i] <- as.numeric(L3[,i])
}
L3 <- L3[,group$variable]
L3$Description <- rownames(L3)
L3[,1:(ncol(L3)-1)] <- t(t(L3[,1:(ncol(L3)-1)])/colSums(L3[,1:(ncol(L3)-1)])*100)
L3 <- L3[,c("Description",colnames(L3)[1:(ncol(L3)-1)])]

colnames(ko) <- ko[1,]
ko <- ko[-1,]
for (i in 1:(ncol(ko)-1)) {
    ko[,i] <- as.numeric(ko[,i])
}
ko <- ko[,c(group$variable,"Description")]
ko[,1:(ncol(ko)-1)] <- t(t(ko[,1:(ncol(ko)-1)])/colSums(ko[,1:(ncol(ko)-1)])*100)
ko <- ko[,c("Description",colnames(ko)[1:(ncol(ko)-1)])]

## Community
dir.create("Results/Gene/Community")

#### relative abundance by sample
dir.create("Results/Gene/Community/summary_r")
write.table(L1,"Results/Gene/Community/summary_r/kegg_L1.xls",
            row.names = FALSE,sep = "\t")
write.table(L2,"Results/Gene/Community/summary_r/kegg_L2.xls",
            row.names = FALSE,sep = "\t")
write.table(L3,"Results/Gene/Community/summary_r/kegg_L3.xls",
            row.names = FALSE,sep = "\t")
write.table(ko,"Results/Gene/Community/summary_r/kegg_ko.xls",
            row.names = FALSE,sep = "\t")

#### relative abundance by group
dir.create("Results/Gene/Community/summary_r_Group")
L1.group <- as.data.frame(t(L1[,-1]))
L1.group$variable <- rownames(L1.group)
L1.group <- merge(L1.group,group)
L1.group <- aggregate(L1.group[,2:(ncol(L1.group)-1)],list(L1.group$Group),mean)
rownames(L1.group) <- L1.group$Group.1
L1.group <- as.data.frame(t(L1.group[,-1]))
L1.group$Description <- rownames(L1.group)
L1.group <- L1.group[,c("Description",levels(group$Group))]
write.table(L1.group,"Results/Gene/Community/summary_r_Group/kegg_L1.xls",
            row.names = FALSE,sep = "\t")

L2.group <- as.data.frame(t(L2[,-1]))
L2.group$variable <- rownames(L2.group)
L2.group <- merge(L2.group,group)
L2.group <- aggregate(L2.group[,2:(ncol(L2.group)-1)],list(L2.group$Group),mean)
rownames(L2.group) <- L2.group$Group.1
L2.group <- as.data.frame(t(L2.group[,-1]))
L2.group$Description <- rownames(L2.group)
L2.group <- L2.group[,c("Description",levels(group$Group))]
write.table(L2.group,"Results/Gene/Community/summary_r_Group/kegg_L2.xls",
            row.names = FALSE,sep = "\t")

L3.group <- as.data.frame(t(L3[,-1]))
L3.group$variable <- rownames(L3.group)
L3.group <- merge(L3.group,group)
L3.group <- aggregate(L3.group[,2:(ncol(L3.group)-1)],list(L3.group$Group),mean)
rownames(L3.group) <- L3.group$Group.1
L3.group <- as.data.frame(t(L3.group[,-1]))
L3.group$Description <- rownames(L3.group)
L3.group <- L3.group[,c("Description",levels(group$Group))]
write.table(L3.group,"Results/Gene/Community/summary_r_Group/kegg_L3.xls",
            row.names = FALSE,sep = "\t")

ko.group <- as.data.frame(t(ko[,-1]))
ko.group$variable <- rownames(ko.group)
ko.group <- merge(ko.group,group)
ko.group <- aggregate(ko.group[,2:(ncol(ko.group)-1)],list(ko.group$Group),mean)
rownames(ko.group) <- ko.group$Group.1
ko.group <- as.data.frame(t(ko.group[,-1]))
ko.group$Description <- rownames(ko.group)
ko.group <- ko.group[,c("Description",levels(group$Group))]
write.table(ko.group,"Results/Gene/Community/summary_r_Group/kegg_ko.xls",
            row.names = FALSE,sep = "\t")

### Community barplot of L1
dir.create("Results/Gene/Community/Barplot")
#### alluvial chart by sample
taxon <- melt(L1)
colnames(taxon) <- c("Taxon","variable","value")
taxon$variable <- factor(taxon$variable,levels = group$variable)

p <- ggplot(data = taxon,aes(x = variable, y = value,
                             alluvium = Taxon, stratum = Taxon)) + 
    geom_alluvium(aes(fill = Taxon),alpha = .5,width = 0.6) + 
    geom_stratum(aes(fill = Taxon),width = 0.6) + 
    ylab(label = "Relative abundance of kegg pathway") + 
    xlab(label = "") + 
    scale_fill_manual(values = cbbPalette,name = "KEGG L1") +
    theme_bw()+ 
    theme(panel.grid=element_blank()) + 
    theme(panel.border = element_blank()) +
    theme(panel.background=element_rect(fill='transparent', color='black'),
          plot.margin = unit(c(3,5,1,1),"mm")) + 
    theme(axis.text.x=element_text(colour="black",size=12,face = "bold",
                                   angle = 45,vjust = 1,hjust = 1)) + 
    theme(axis.text.y=element_text(colour = "black",size = 10)) + 
    theme(axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"))+ 
    theme(axis.title.y = element_text(size = 12,face = "bold",
                                      margin = unit(c(0,1,0,1),"lines"))) + 
    scale_y_continuous(limits = c(0,100.001),expand = c(0,0)) + 
    theme(legend.text = element_text(colour = "black",size = 12)) + 
    theme(legend.title = element_text(size = 14,colour = "black",face = "bold"))

pdf(file = "Results/Gene/Community/Barplot/L1.pdf",
    width = 0.3*Sample_numb + 1.5 + max(str_length(rownames(L1)))*0.1,height = 5)
p
dev.off()

#### alluvial chart by group
taxon <- melt(L1.group)
colnames(taxon) <- c("Taxon","variable","value")

p <- ggplot(data = taxon,aes(x = variable, y = value,
                             alluvium = Taxon, stratum = Taxon)) + 
    geom_alluvium(aes(fill = Taxon),alpha = .5,width = 0.6) + 
    geom_stratum(aes(fill = Taxon),width = 0.6) + 
    ylab(label = "Relative abundance of kegg pathway") + 
    xlab(label = "") + 
    scale_fill_manual(values = cbbPalette,name = "KEGG L1") +
    theme_bw()+ 
    theme(panel.grid=element_blank()) + 
    theme(panel.border = element_blank()) +
    theme(panel.background=element_rect(fill='transparent', color='black'),
          plot.margin = unit(c(3,5,1,1),"mm")) + 
    theme(axis.text.x=element_text(colour="black",size=12,face = "bold",
                                   angle = 45,vjust = 1,hjust = 1)) + 
    theme(axis.text.y=element_text(colour = "black",size = 10)) + 
    theme(axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"))+ 
    theme(axis.title.y = element_text(size = 12,face = "bold",
                                      margin = unit(c(0,1,0,1),"lines"))) + 
    scale_y_continuous(limits = c(0,100.001),expand = c(0,0)) + 
    theme(legend.text = element_text(colour = "black",size = 12)) + 
    theme(legend.title = element_text(size = 14,colour = "black",face = "bold"))

pdf(file = "Results/Gene/Community/Barplot/L1_Group.pdf",
    width = wid*Group_numb + 1.5 + max(str_length(rownames(L1.group)))*0.1,height = 5)
p
dev.off()

source("Functions/top30 2.R")
L2.top <- top30(L2)
L3.top <- top30(L3)

#### bubble chart by sample
dir.create("Results/Gene/Community/Bubble_chart")
source("Functions/Bubblechart.s 2.R")
result <- bubblechart.s(L2,group)
pdf(file = "Results/Gene/Community/Bubble_chart/L2.pdf",
    width = 0.22*Sample_numb + 0.5 + max(str_length(rownames(L2.top)))*0.1,
    height = 12.5)
result
dev.off()

result <- bubblechart.s(L3,group)
pdf(file = "Results/Gene/Community/Bubble_chart/L3.pdf",
    width = 0.22*Sample_numb + 0.5 + max(str_length(rownames(L3.top)))*0.1,
    height = 12.5)
result
dev.off()

result <- bubblechart.s(ko,group)
pdf(file = "Results/Gene/Community/Bubble_chart/ko.pdf",
    width = 0.22*Sample_numb + 0.5 + max(str_length(rownames(ko)))*0.1,
    height = 12.5)
result
dev.off()

#### bubble chart by group
source("Functions/Bubblechart.g 2.R")
result <- bubblechart.g(L2.group,group)
pdf(file = "Results/Gene/Community/Bubble_chart/L2_Group.pdf",
    width = 0.3*Group_numb  + 0.5 + max(str_length(rownames(L2.top)))*0.1,
    height = 12.5)
result
dev.off()

result <- bubblechart.g(L3.group,group)
pdf(file = "Results/Gene/Community/Bubble_chart/L3_Group.pdf",
    width = 0.3*Group_numb  + 0.5 + max(str_length(rownames(L3.top)))*0.1,
    height = 12.5)
result
dev.off()

result <- bubblechart.g(ko.group,group)
pdf(file = "Results/Gene/Community/Bubble_chart/ko_Group.pdf",
    width = 0.3*Group_numb  + 1.8,
    height = 12.5)
result
dev.off()

#### heatmap by sample
dir.create("Results/Gene/Community/Heatmap")
f.abundance <- L2[,-1]
sum <- apply(f.abundance,1,sum) 
f.abundance <- cbind(f.abundance,sum)
f.abundance <- as.data.frame(f.abundance)
f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
f.abundance <- subset(f.abundance, select = -sum)
f.abundance <- f.abundance[rowSums(f.abundance) > 0,]
f.abundance <- f.abundance[1:30,]

p <- pheatmap(f.abundance,color = colorRampPalette(c("navy","white","firebrick3"))(100),
                     fontsize = 12,scale = "row",angle_col = 90)

pdf(file = "Results/Gene/Community/Heatmap/L2.pdf",
    width = 0.22*Sample_numb + 0.5 + max(str_length(rownames(L2.top)))*0.1,
    height = 12.5)
p
dev.off()

f.abundance <- L3[,-1]
sum <- apply(f.abundance,1,sum) 
f.abundance <- cbind(f.abundance,sum)
f.abundance <- as.data.frame(f.abundance)
f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
f.abundance <- subset(f.abundance, select = -sum)
f.abundance <- f.abundance[rowSums(f.abundance) > 0,]
f.abundance <- f.abundance[1:30,]

p <- pheatmap(f.abundance,color = colorRampPalette(c("navy","white","firebrick3"))(100),
              fontsize = 12,scale = "row",angle_col = 90)

pdf(file = "Results/Gene/Community/Heatmap/L3.pdf",
    width = 0.22*Sample_numb + 0.5 + max(str_length(rownames(L3.top)))*0.1,
    height = 12.5)
p
dev.off()

f.abundance <- ko[,-1]
sum <- apply(f.abundance,1,sum) 
f.abundance <- cbind(f.abundance,sum)
f.abundance <- as.data.frame(f.abundance)
f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
f.abundance <- subset(f.abundance, select = -sum)
f.abundance <- f.abundance[rowSums(f.abundance) > 0,]
f.abundance <- f.abundance[1:30,]

p <- pheatmap(f.abundance,color = colorRampPalette(c("navy","white","firebrick3"))(100),
              fontsize = 12,scale = "row",angle_col = 90)

pdf(file = "Results/Gene/Community/Heatmap/ko.pdf",
    width = 0.22*Sample_numb + 0.5 + max(str_length(rownames(ko)))*0.1,
    height = 12.5)
p
dev.off()

#### Heatmap by group
f.abundance <- L2.group[,-1]
sum <- apply(f.abundance,1,sum) 
f.abundance <- cbind(f.abundance,sum)
f.abundance <- as.data.frame(f.abundance)
f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
f.abundance <- subset(f.abundance, select = -sum)
f.abundance <- f.abundance[rowSums(f.abundance) > 0,]
f.abundance <- f.abundance[1:30,]

p <- pheatmap(f.abundance,color = colorRampPalette(c("navy","white","firebrick3"))(100),
              fontsize = 12,scale = "row",angle_col = 90)

pdf(file = "Results/Gene/Community/Heatmap/L2_Group.pdf",
    width = 0.3*Group_numb + 0.5 + max(str_length(rownames(f.abundance)))*0.1,
    height = 12.5)
p
dev.off()

f.abundance <- L3.group[,-1]
sum <- apply(f.abundance,1,sum) 
f.abundance <- cbind(f.abundance,sum)
f.abundance <- as.data.frame(f.abundance)
f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
f.abundance <- subset(f.abundance, select = -sum)
f.abundance <- f.abundance[rowSums(f.abundance) > 0,]
f.abundance <- f.abundance[1:30,]

p <- pheatmap(f.abundance,color = colorRampPalette(c("navy","white","firebrick3"))(100),
              fontsize = 12,scale = "row",angle_col = 90)

pdf(file = "Results/Gene/Community/Heatmap/L3_Group.pdf",
    width = 0.3*Group_numb + 0.5 + max(str_length(rownames(f.abundance)))*0.1,
    height = 12.5)
p
dev.off()

f.abundance <- ko.group[,-1]
sum <- apply(f.abundance,1,sum) 
f.abundance <- cbind(f.abundance,sum)
f.abundance <- as.data.frame(f.abundance)
f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
f.abundance <- subset(f.abundance, select = -sum)
f.abundance <- f.abundance[rowSums(f.abundance) > 0,]
f.abundance <- f.abundance[1:30,]

p <- pheatmap(f.abundance,color = colorRampPalette(c("navy","white","firebrick3"))(100),
              fontsize = 12,scale = "row",angle_col = 90)

pdf(file = "Results/Gene/Community/Heatmap/ko_Group.pdf",
    width = 0.3*Group_numb + 2.2,
    height = 12.5)
p
dev.off()

## Beta diveristy
dir.create("Results/Gene/Beta_diversity")

### PCoA
dir.create("Results/Gene/Beta_diversity/Distance")
dir.create("Results/Gene/Beta_diversity/PCoA")
source("Functions/PCoA 2.R")
result <- pcoa.community(t(ko[,-1]),group)
write.table(as.matrix(result[[1]]),
            "Results/Gene/Beta_diversity/Distance/distance_bray_curtis.txt",
            sep = "\t")
write.table(result[[2]],"Results/Gene/Beta_diversity/PCoA/pcoa.txt",
            sep = "\t",row.names = FALSE)
pdf(file = "Results/Gene/Beta_diversity/PCoA/Pcoa_group.pdf",
    width = 7.5,height = 5.4)
result[[3]]
dev.off()
pdf(file = "Results/Gene/Beta_diversity/PCoA/Pcoa_ellipse.pdf",
    width = 7.5,height = 5.4)
result[[4]]
dev.off()
pdf(file = "Results/Gene/Beta_diversity/PCoA/Pcoa_label.pdf",
    width = 7.5,height = 5.4)
result[[5]]
dev.off()

#### PCA
dir.create("Results/Gene/Beta_diversity/PCA")
source("Functions/PCA 2.R")
result <- pca.community(t(ko[,-1]),group)
write.table(result[[1]],"Results/Gene/Beta_diversity/PCA/PCA.txt",
            sep = "\t",row.names = FALSE)
pdf(file = "Results/Gene/Beta_diversity/PCA/PCA_group.pdf",
    width = 7.5,height = 5.4)
result[[2]]
dev.off()
pdf(file = "Results/Gene/Beta_diversity/PCA/PCA_ellipse.pdf",
    width = 7.5,height = 5.4)
result[[3]]
dev.off()
pdf(file = "Results/Gene/Beta_diversity/PCA/PCA_label.pdf",
    width = 7.5,height = 5.4)
result[[4]]
dev.off()

#### NMDS
dir.create("Results/Gene/Beta_diversity/NMDS")
source("Functions/NMDS 2.R")
result <- nmds.community(t(ko[,-1]),group)
write.table(result[[1]],"Results/Gene/Beta_diversity/NMDS/Microbe_stress.txt",
            sep = "\t",row.names = FALSE)
write.table(result[[2]],"Results/Gene/Beta_diversity/NMDS/NMDS.txt",
            sep = "\t",row.names = FALSE)
pdf(file = "Results/Gene/Beta_diversity/NMDS/NMDS_group.pdf",
    width = 7.5,height = 5.4)
result[[3]]
dev.off()
pdf(file = "Results/Gene/Beta_diversity/NMDS/NMDS_ellipse.pdf",
    width = 7.5,height = 5.4)
result[[4]]
dev.off()
pdf(file = "Results/Gene/Beta_diversity/NMDS/NMDS_label.pdf",
    width = 7.5,height = 5.4)
result[[5]]
dev.off()

### ANOSIM, Adonis, and MRPP tests
dir.create("Results/Gene/Beta_diversity/Beta_tests")
source("Functions/beta.test.R")
result <- beta.test(t(ko[,-1]),group)
sink("Results/Gene/Beta_diversity/Beta_tests/Adonis.txt")
result[[1]]
sink("Results/Gene/Beta_diversity/Beta_tests/ANOSIM.txt")
result[[2]]
sink("Results/Gene/Beta_diversity/Beta_tests/MRPP.txt")
result[[3]]

### Hcluster
dir.create("Results/Gene/Beta_diversity/Hcluster")
dist_bray <- vegdist(t(ko[,-1]), method = "bray")
hc <- hclust(dist_bray, method = "average")
hc_tre <- as.phylo(hc)
write.tree(hc_tre, "Results/Gene/Beta_diversity/Hcluster/Hcluster.tre")

source("Functions/plot.phylo2.r")
pdf("Results/Gene/Beta_diversity/Hcluster/Hcluster.pdf", 
    width = 7, height = 0.22*Sample_numb + 2)
plot.phylo2(hc_tre, type = "phylogram", cex = 1, adj = 1)
title(main = "UPGMA clustering tree", line = 2, cex.main = 1.5, font.main= 2)
dev.off()

## Difference analysis
dir.create("Results/Gene/Diff_analysis")

### Wilcox rank-sum test
dir.create("Results/Gene/Diff_analysis/Wilcox")
source("Functions/wilcox 2.R")
source("Functions/wilcox.abun.R")

#### L1
dir.create("Results/Gene/Diff_analysis/Wilcox/L1")
set.seed(1111)
wilcox.result <- wilcox.biomarker(L1,group)
write.table(wilcox.result[[1]],
            "Results/Gene/Diff_analysis/Wilcox/L1/wilcox_all.txt",
            sep = "\t",row.names = FALSE)
write.table(wilcox.result[[2]],
            "Results/Gene/Diff_analysis/Wilcox/L1/wilcox_biomarker.txt",
            sep = "\t",row.names = FALSE)
if (nrow(wilcox.result[[2]]) > 0) {
    L1.result <- wilcox.abun(L1,wilcox.result[[2]],group)
    pdf("Results/Gene/Diff_analysis/Wilcox/L1/wilcox_biomarker_abun.pdf",
        width = 12,height = 1.8 + 0.18*nrow(wilcox.result[[2]]))
    print(L1.result[[1]])
    dev.off()
    write.table(L1.result[[2]],
                "Results/Gene/Diff_analysis/Wilcox/L1/Group_mean_abun_biomarker.txt",
                sep = "\t",row.names = FALSE)
    write.table(L1.result[[3]],
                "Results/Gene/Diff_analysis/Wilcox/L1/Abun_change_biomarker.txt",
                sep = "\t",row.names = FALSE)
}

#### L2
dir.create("Results/Gene/Diff_analysis/Wilcox/L2")
set.seed(1111)
wilcox.result <- wilcox.biomarker(L2,group)
write.table(wilcox.result[[1]],
            "Results/Gene/Diff_analysis/Wilcox/L2/wilcox_all.txt",
            sep = "\t",row.names = FALSE)
write.table(wilcox.result[[2]],
            "Results/Gene/Diff_analysis/Wilcox/L2/wilcox_biomarker.txt",
            sep = "\t",row.names = FALSE)
if (nrow(wilcox.result[[2]]) > 0) {
    Class <- wilcox.abun(L2,wilcox.result[[2]],group)
    pdf("Results/Gene/Diff_analysis/Wilcox/L2/wilcox_biomarker_abun.pdf",
        width = 12,height = 1.8 + 0.18*nrow(wilcox.result[[2]]))
    print(Class[[1]])
    dev.off()
    write.table(Class[[2]],
                "Results/Gene/Diff_analysis/Wilcox/L2/Group_mean_abun_biomarker.txt",
                sep = "\t",row.names = FALSE)
    write.table(Class[[3]],
                "Results/Gene/Diff_analysis/Wilcox/L2/Abun_change_biomarker.txt",
                sep = "\t",row.names = FALSE)
}

#### L3
dir.create("Results/Gene/Diff_analysis/Wilcox/L3")
set.seed(1111)
wilcox.result <- wilcox.biomarker(L3,group)
write.table(wilcox.result[[1]],
            "Results/Gene/Diff_analysis/Wilcox/L3/wilcox_all.txt",
            sep = "\t",row.names = FALSE)
write.table(wilcox.result[[2]],
            "Results/Gene/Diff_analysis/Wilcox/L3/wilcox_biomarker.txt",
            sep = "\t",row.names = FALSE)
if (nrow(wilcox.result[[2]]) > 0) {
    Order <- wilcox.abun(L3,wilcox.result[[2]],group)
    pdf("Results/Gene/Diff_analysis/Wilcox/L3/wilcox_biomarker_abun.pdf",
        width = 12,height = 1.8 + 0.18*nrow(wilcox.result[[2]]))
    print(Order[[1]])
    dev.off()
    write.table(Order[[2]],
                "Results/Gene/Diff_analysis/Wilcox/L3/Group_mean_abun_biomarker.txt",
                sep = "\t",row.names = FALSE)
    write.table(Order[[3]],
                "Results/Gene/Diff_analysis/Wilcox/L3/Abun_change_biomarker.txt",
                sep = "\t",row.names = FALSE)
}

#### ko
ko1 <- ko
ko1$Description <- rownames(ko)
dir.create("Results/Gene/Diff_analysis/Wilcox/ko")
set.seed(1111)
wilcox.result <- wilcox.biomarker(ko1,group)
write.table(wilcox.result[[1]],
            "Results/Gene/Diff_analysis/Wilcox/ko/wilcox_all.txt",
            sep = "\t",row.names = FALSE)
write.table(wilcox.result[[2]],
            "Results/Gene/Diff_analysis/Wilcox/ko/wilcox_biomarker.txt",
            sep = "\t",row.names = FALSE)
if (nrow(wilcox.result[[2]]) > 0) {
    Order <- wilcox.abun(ko1,wilcox.result[[2]],group)
    pdf("Results/Gene/Diff_analysis/Wilcox/ko/wilcox_biomarker_abun.pdf",
        width = 12,height = 1.8 + 0.18*nrow(wilcox.result[[2]]))
    print(Order[[1]])
    dev.off()
    write.table(Order[[2]],
                "Results/Gene/Diff_analysis/Wilcox/ko/Group_mean_abun_biomarker.txt",
                sep = "\t",row.names = FALSE)
    write.table(Order[[3]],
                "Results/Gene/Diff_analysis/Wilcox/ko/Abun_change_biomarker.txt",
                sep = "\t",row.names = FALSE)
}

### Random Forest for biomarker identification
dir.create("Results/Gene/Diff_analysis/RandomForest")
source("Functions/RandomForest.R")

#### L1
dir.create("Results/gene/Diff_analysis/RandomForest/L1")
set.seed(1111)
RF.result <- Random.forest.biomarker(L1,group)
pdf("Results/Gene/Diff_analysis/RandomForest/L1/classification.pdf",
    height = 3,width = 3)
RF.result[[1]]
dev.off()
importance.otu <- RF.result[[2]]
importance.otu$ID <- rownames(importance.otu)
importance.otu <- importance.otu[,c("ID",colnames(importance.otu)[1:(ncol(importance.otu)-1)])]
write.table(importance.otu,
            "Results/Gene/Diff_analysis/RandomForest/L1/importance.txt",
            sep = "\t",row.names = FALSE)
if (nrow(importance.otu) > 0) {
    pdf("Results/Gene/Diff_analysis/RandomForest/L1/RandomForest_biomarkder.pdf",
        height = nrow(importance.otu)*0.3 + 1.85,width = 11)
    print(RF.result[[3]])
    dev.off()
    write.table(RF.result[[4]],
                "Results/Gene/Diff_analysis/RandomForest/L1/biomarker_abundance.txt",
                sep = "\t") 
}

#### L2
dir.create("Results/Gene/Diff_analysis/RandomForest/L2")
set.seed(1111)
RF.result <- Random.forest.biomarker(L2,group)
pdf("Results/Gene/Diff_analysis/RandomForest/L2/classification.pdf",
    height = 3,width = 3)
RF.result[[1]]
dev.off()
importance.otu <- RF.result[[2]]
importance.otu$ID <- rownames(importance.otu)
importance.otu <- importance.otu[,c("ID",colnames(importance.otu)[1:(ncol(importance.otu)-1)])]
write.table(importance.otu,
            "Results/Gene/Diff_analysis/RandomForest/L2/importance.txt",
            sep = "\t",row.names = FALSE)
if (nrow(importance.otu) > 0) {
    pdf("Results/Gene/Diff_analysis/RandomForest/L2/RandomForest_biomarkder.pdf",
        height = nrow(importance.otu)*0.3 + 1.85,width = 11)
    print(RF.result[[3]])
    dev.off()
    write.table(RF.result[[4]],
                "Results/Gene/Diff_analysis/RandomForest/L2/biomarker_abundance.txt",
                sep = "\t") 
}

#### L3
dir.create("Results/Gene/Diff_analysis/RandomForest/L3")
set.seed(1111)
RF.result <- Random.forest.biomarker(L3,group)
pdf("Results/Gene/Diff_analysis/RandomForest/L3/classification.pdf",
    height = 3,width = 3)
RF.result[[1]]
dev.off()
importance.otu <- RF.result[[2]]
importance.otu$ID <- rownames(importance.otu)
importance.otu <- importance.otu[,c("ID",colnames(importance.otu)[1:(ncol(importance.otu)-1)])]
write.table(importance.otu,
            "Results/Gene/Diff_analysis/RandomForest/L3/importance.txt",
            sep = "\t",row.names = FALSE)
if (nrow(importance.otu) > 0) {
    pdf("Results/Gene/Diff_analysis/RandomForest/L3/RandomForest_biomarkder.pdf",
        height = nrow(importance.otu)*0.3 + 1.85,width = 11)
    print(RF.result[[3]])
    dev.off()
    write.table(RF.result[[4]],
                "Results/Gene/Diff_analysis/RandomForest/L3/biomarker_abundance.txt",
                sep = "\t") 
}

#### ko
dir.create("Results/Gene/Diff_analysis/RandomForest/ko")
set.seed(1111)
RF.result <- Random.forest.biomarker(ko1,group)
pdf("Results/Gene/Diff_analysis/RandomForest/ko/classification.pdf",
    height = 3,width = 3)
RF.result[[1]]
dev.off()
importance.otu <- RF.result[[2]]
importance.otu$ID <- rownames(importance.otu)
importance.otu <- importance.otu[,c("ID",colnames(importance.otu)[1:(ncol(importance.otu)-1)])]
write.table(importance.otu,
            "Results/Gene/Diff_analysis/RandomForest/ko/importance.txt",
            sep = "\t",row.names = FALSE)
if (nrow(importance.otu) > 0) {
    pdf("Results/Gene/Diff_analysis/RandomForest/ko/RandomForest_biomarkder.pdf",
        height = nrow(importance.otu)*0.3 + 1.85,width = 9)
    print(RF.result[[3]])
    dev.off()
    write.table(RF.result[[4]],
                "Results/Gene/Diff_analysis/RandomForest/ko/biomarker_abundance.txt",
                sep = "\t") 
}

### LEfSe
dir.create("Results/Gene/Diff_analysis/LEfSe")
source("Functions/lefse 2.R")

Lefse.result <- lefse.1(ko1,group)
write.table(Lefse.result,"Results/Gene/Diff_analysis/LEfSe/Lefse_result.txt",
            sep = "\t")
if (nrow(Lefse.result) > 0) {
    source("Functions/lefse_bar.R")
    p <- lefse_bar(Lefse.result)
    pdf("Results/Gene/Diff_analysis/LEfSe/Lefse_bar.pdf",
        width = 6,height = 1.15 + 0.3*nrow(Lefse.result))
    print(p)
    dev.off()
}

## Diagnose
dir.create("Results/Gene/Diagnose")

### Wilcox rank-sum test
dir.create("Results/Gene/Diagnose/Wilcox")
source("Functions/diagnose.R")

#### L1
data <- read.table("Results/Gene/Diff_analysis/Wilcox/L1/wilcox_biomarker.txt",
                   header = TRUE,sep = "\t")
if (nrow(data) > 0) {
    dir.create("Results/Gene/Diagnose/Wilcox/L1")
    L1.result <- diagnose(data,L1,group)
    write.table(L1.result[[1]],"Results/Gene/Diagnose/Wilcox/L1/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Gene/Diagnose/Wilcox/L1/index.pdf",
        width = 2.4,height = 3.6)
    print(L1.result[[2]])
    dev.off()
    pdf("Results/Gene/Diagnose/Wilcox/L1/Roc.pdf",width = 10,height = 10)
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

#### L2
data <- read.table("Results/Gene/Diff_analysis/Wilcox/L2/wilcox_biomarker.txt",
                   header = TRUE,sep = "\t")
if (nrow(data) > 0) {
    dir.create("Results/Gene/Diagnose/Wilcox/L2")
    L2.result <- diagnose(data,L2,group)
    write.table(L2.result[[1]],"Results/Gene/Diagnose/Wilcox/L2/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Gene/Diagnose/Wilcox/L2/index.pdf",
        width = 2.4,height = 3.6)
    print(L2.result[[2]])
    dev.off()
    pdf("Results/Gene/Diagnose/Wilcox/L2/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(L2.result[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
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

#### L3
data <- read.table("Results/Gene/Diff_analysis/Wilcox/L3/wilcox_biomarker.txt",
                   header = TRUE,sep = "\t")
if (nrow(data) > 0) {
    dir.create("Results/Gene/Diagnose/Wilcox/L3")
    Order <- diagnose(data,L3,group)
    write.table(Order[[1]],"Results/Gene/Diagnose/Wilcox/L3/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Gene/Diagnose/Wilcox/L3/index.pdf",
        width = 2.4,height = 3.6)
    print(Order[[2]])
    dev.off()
    pdf("Results/Gene/Diagnose/Wilcox/L3/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(Order[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
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

#### ko
data <- read.table("Results/Gene/Diff_analysis/Wilcox/ko/wilcox_biomarker.txt",
                   header = TRUE,sep = "\t")
if (nrow(data) > 0) {
    dir.create("Results/Gene/Diagnose/Wilcox/ko")
    Order <- diagnose(data,ko1,group)
    write.table(Order[[1]],"Results/Gene/Diagnose/Wilcox/ko/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Gene/Diagnose/Wilcox/ko/index.pdf",
        width = 2.4,height = 3.6)
    print(Order[[2]])
    dev.off()
    pdf("Results/Gene/Diagnose/Wilcox/ko/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(Order[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
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
dir.create("Results/Gene/Diagnose/RandomForest")

#### L1
data <- read.table("Results/Gene/Diff_analysis/RandomForest/L1/importance.txt",
                   header = TRUE,sep = "\t")
colnames(data)[1] <- "V1"
if (nrow(data) > 0) {
    dir.create("Results/Gene/Diagnose/RandomForest/L1")
    Phylum <- diagnose(data,L1,group)
    write.table(Phylum[[1]],"Results/Gene/Diagnose/RandomForest/L1/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Gene/Diagnose/RandomForest/L1/index.pdf",
        width = 2.4,height = 3.6)
    print(Phylum[[2]])
    dev.off()
    pdf("Results/Gene/Diagnose/RandomForest/L1/Roc.pdf",width = 10,height = 10)
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

#### L2
data <- read.table("Results/Gene/Diff_analysis/RandomForest/L2/importance.txt",
                   header = TRUE,sep = "\t")
colnames(data)[1] <- "V1"
if (nrow(data) > 0) {
    dir.create("Results/Gene/Diagnose/RandomForest/L2")
    Class <- diagnose(data,L2,group)
    write.table(Class[[1]],"Results/Gene/Diagnose/RandomForest/L2/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Gene/Diagnose/RandomForest/L2/index.pdf",
        width = 2.4,height = 3.6)
    print(Class[[2]])
    dev.off()
    pdf("Results/Gene/Diagnose/RandomForest/L2/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(Class[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
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

#### L3
data <- read.table("Results/Gene/Diff_analysis/RandomForest/L3/importance.txt",
                   header = TRUE,sep = "\t")
colnames(data)[1] <- "V1"
if (nrow(data) > 0) {
    dir.create("Results/Gene/Diagnose/RandomForest/L3")
    Order <- diagnose(data,L3,group)
    write.table(Order[[1]],"Results/Gene/Diagnose/RandomForest/L3/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Gene/Diagnose/RandomForest/L3/index.pdf",
        width = 2.4,height = 3.6)
    print(Order[[2]])
    dev.off()
    pdf("Results/Gene/Diagnose/RandomForest/L3/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(Order[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
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

#### ko
data <- read.table("Results/Gene/Diff_analysis/RandomForest/ko/importance.txt",
                   header = TRUE,sep = "\t")
colnames(data)[1] <- "V1"
if (nrow(data) > 0) {
    dir.create("Results/Gene/Diagnose/RandomForest/ko")
    Order <- diagnose(data,ko1,group)
    write.table(Order[[1]],"Results/Gene/Diagnose/RandomForest/ko/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Gene/Diagnose/RandomForest/ko/index.pdf",
        width = 2.4,height = 3.6)
    print(Order[[2]])
    dev.off()
    pdf("Results/Gene/Diagnose/RandomForest/ko/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(Order[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
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
dir.create("Results/Gene/Diagnose/LEfSe")
data <- read.table("Results/Gene/Diff_analysis/LEfSe/Lefse_result.txt",
                   header = TRUE,sep = "\t",row.names = 1)
data$V1 <- rownames(data)
if (nrow(data) > 0) {
    Species <- diagnose(data,ko1,group)
    write.table(Species[[1]],"Results/Gene/Diagnose/LEfSe/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Gene/Diagnose/LEfSe/index.pdf",
        width = 2.4,height = 3.6)
    print(Species[[2]])
    dev.off()
    pdf("Results/Gene/Diagnose/LEfSe/Roc.pdf",width = 10,height = 10)
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



