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
species.type <- "Virus"
dir.create("Results")
dir.create("Results/Virus")

### Input data
data <- read.table(file = "Input/virus.tpm.xls",header = FALSE,row.names = 1,
                  sep = "\t",quote = "")

### data processing
colnames(data) <- data[1,]
data <- data[-1,]
for (i in 1:(ncol(data)-1)) {
    data[,i] <- as.numeric(data[,i])
}

data$Phylum <- gsub(".*p__","",data$taxonomy)
data$Phylum <- gsub("\\;c__.*","",data$Phylum)
data$Class <- gsub(".*c__","",data$taxonomy)
data$Class <- gsub("\\;o__.*","",data$Class)
data$Order <- gsub(".*o__","",data$taxonomy)
data$Order <- gsub("\\;f__.*","",data$Order)
data$Family <- gsub(".*f__","",data$taxonomy)
data$Family <- gsub("\\;g__.*","",data$Family)
data$Genus <- gsub(".*g__","",data$taxonomy)
data$Genus <- gsub("\\;s__.*","",data$Genus)
data$Species <- gsub(".*s__","",data$taxonomy)
data <- subset(data,select = -taxonomy)
data.tax <- data[,(ncol(data)-5):ncol(data)]
data.tax[data.tax == "norank"] <- "Unclassified"
data.tax[data.tax == "uncultured bacterium"] <- "Unclassified"
data.tax[data.tax == "uncultured"] <- "Unclassified"
for (i in 1:nrow(data.tax)) {
    data.tax[i,grepl("d__Bacteria",data.tax[i,])] <- "Unclassified"
}
abundance <- t(t(data[,1:(ncol(data)-6)])/colSums(data[,1:(ncol(data)-6)])*100)
abundance.tax <- cbind(abundance,data.tax)
data <- cbind(data[,1:(ncol(data)-6)],data.tax)

## Community
dir.create("Results/Virus/Community")

### abundance table
#### read number abundance by sample
dir.create("Results/Virus/Community/tax_summary_a")
source("Functions/tax.summary.R")
result <- tax.summary(data)
write.table(result[[1]],"Results/Virus/Community/tax_summary_a/Phylum.xls",
            row.names = FALSE,sep = "\t")
write.table(result[[2]],"Results/Virus/Community/tax_summary_a/Class.xls",
            row.names = FALSE,sep = "\t")
write.table(result[[3]],"Results/Virus/Community/tax_summary_a/Order.xls",
            row.names = FALSE,sep = "\t")
write.table(result[[4]],"Results/Virus/Community/tax_summary_a/Family.xls",
            row.names = FALSE,sep = "\t")
write.table(result[[5]],"Results/Virus/Community/tax_summary_a/Genus.xls",
            row.names = FALSE,sep = "\t")
## write.table(result[[6]],"Results/Virus/Community/tax_summary_a/Species.xls",
##             row.names = FALSE,sep = "\t")

#### read number abundance by group
dir.create("Results/Virus/Community/tax_summary_a_Group")
data.1 <- as.data.frame(t(data[,1:(ncol(data)-6)]))
aa <- match(rownames(data.1),group$variable)
data.1 <- data.1[aa,]
data.1$Group <- group$Group
data.group <- aggregate(data.1[,1:(ncol(data.1)-1)],
                             list(data.1$Group),mean)
rownames(data.group) <- data.group$Group.1
data.group <- t(data.group[,-1])
data.group <- as.data.frame(cbind(data.group,data.tax))
result <- tax.summary(data.group)
write.table(result[[1]],"Results/Virus/Community/tax_summary_a_Group/Phylum.xls",
            row.names = FALSE,sep = "\t")
write.table(result[[2]],"Results/Virus/Community/tax_summary_a_Group/Class.xls",
            row.names = FALSE,sep = "\t")
write.table(result[[3]],"Results/Virus/Community/tax_summary_a_Group/Order.xls",
            row.names = FALSE,sep = "\t")
write.table(result[[4]],"Results/Virus/Community/tax_summary_a_Group/Family.xls",
            row.names = FALSE,sep = "\t")
write.table(result[[5]],"Results/Virus/Community/tax_summary_a_Group/Genus.xls",
            row.names = FALSE,sep = "\t")
## write.table(result[[6]],"Results/Virus/Community/tax_summary_a_Group/Species.xls",
##             row.names = FALSE,sep = "\t")

#### relative abundance by sample
dir.create("Results/Virus/Community/tax_summary_r")
result <- tax.summary(abundance.tax)
write.table(result[[1]],"Results/Virus/Community/tax_summary_r/Phylum.xls",
            row.names = FALSE,sep = "\t")
write.table(result[[2]],"Results/Virus/Community/tax_summary_r/Class.xls",
            row.names = FALSE,sep = "\t")
write.table(result[[3]],"Results/Virus/Community/tax_summary_r/Order.xls",
            row.names = FALSE,sep = "\t")
write.table(result[[4]],"Results/Virus/Community/tax_summary_r/Family.xls",
            row.names = FALSE,sep = "\t")
write.table(result[[5]],"Results/Virus/Community/tax_summary_r/Genus.xls",
            row.names = FALSE,sep = "\t")
## write.table(result[[6]],"Results/Virus/Community/tax_summary_r/Species.xls",
##             row.names = FALSE,sep = "\t")

#### relative abundance by group
dir.create("Results/Virus/Community/tax_summary_r_Group")
abundance.1 <- as.data.frame(t(abundance))
aa <- match(rownames(abundance.1),group$variable)
abundance.1 <- abundance.1[aa,]
abundance.1$Group <- group$Group
abundance.group <- aggregate(abundance.1[,1:(ncol(abundance.1)-1)],
                             list(abundance.1$Group),mean)
rownames(abundance.group) <- abundance.group$Group.1
abundance.group <- t(abundance.group[,-1])
abundance.group <- as.data.frame(cbind(abundance.group,data.tax))
result <- tax.summary(abundance.group)
write.table(result[[1]],"Results/Virus/Community/tax_summary_r_Group/Phylum.xls",
            row.names = FALSE,sep = "\t")
write.table(result[[2]],"Results/Virus/Community/tax_summary_r_Group/Class.xls",
            row.names = FALSE,sep = "\t")
write.table(result[[3]],"Results/Virus/Community/tax_summary_r_Group/Order.xls",
            row.names = FALSE,sep = "\t")
write.table(result[[4]],"Results/Virus/Community/tax_summary_r_Group/Family.xls",
            row.names = FALSE,sep = "\t")
write.table(result[[5]],"Results/Virus/Community/tax_summary_r_Group/Genus.xls",
            row.names = FALSE,sep = "\t")
## write.table(result[[6]],"Results/Virus/Community/tax_summary_r_Group/Species.xls",
##             row.names = FALSE,sep = "\t")

source("Functions/top10 2.R")
result <- top10(abundance.tax)
aa <- c(max(str_length(colnames(result[[1]]))),max(str_length(colnames(result[[2]]))),
        max(str_length(colnames(result[[3]]))),max(str_length(colnames(result[[4]]))),
        max(str_length(colnames(result[[5]]))),max(str_length(colnames(result[[6]]))))

### Community barplot
#### alluvial chart by sample
dir.create("Results/Virus/Community/Alluvial_chart")
source("Functions/Community.s.R")
result <- community.s(abundance.tax,group)
pdf(file = "Results/Virus/Community/Alluvial_chart/Phylum.pdf",
    width = 1.5 + 0.3*Sample_numb + aa[1]*0.1,height = 5)
result[[1]]
dev.off()
pdf(file = "Results/Virus/Community/Alluvial_chart/Class.pdf",
    width = 1.5 + 0.3*Sample_numb + aa[2]*0.1,height = 5)
result[[2]]
dev.off()
pdf(file = "Results/Virus/Community/Alluvial_chart/Order.pdf",
    width = 1.5 + 0.3*Sample_numb + aa[3]*0.1,height = 5)
result[[3]]
dev.off()
pdf(file = "Results/Virus/Community/Alluvial_chart/Family.pdf",
    width = 1.5 + 0.3*Sample_numb + aa[4]*0.1,height = 5)
result[[4]]
dev.off()
pdf(file = "Results/Virus/Community/Alluvial_chart/Genus.pdf",
    width = 1.5 + 0.3*Sample_numb + aa[5]*0.1,height = 5)
result[[5]]
dev.off()
## pdf(file = "Results/Virus/Community/Alluvial_chart/Species.pdf",
##     width = 1.5 + 0.3*Sample_numb + aa[6]*0.1,height = 5)
## result[[6]]
## dev.off()

#### alluvial chart by group
dir.create("Results/Virus/Community/Alluvial_chart_Group")
source("Functions/Community.g.R")
result <- community.g(abundance.tax,group)
pdf(file = "Results/Virus/Community/Alluvial_chart_Group/Phylum.pdf",
    width = 1.5 + wid*Group_numb + aa[1]*0.1,height = 5)
result[[1]]
dev.off()
pdf(file = "Results/Virus/Community/Alluvial_chart_Group/Class.pdf",
    width = 1.5 + wid*Group_numb + aa[2]*0.1,height = 5)
result[[2]]
dev.off()
pdf(file = "Results/Virus/Community/Alluvial_chart_Group/Order.pdf",
    width = 1.5 + wid*Group_numb + aa[3]*0.1,height = 5)
result[[3]]
dev.off()
pdf(file = "Results/Virus/Community/Alluvial_chart_Group/Family.pdf",
    width = 1.5 + wid*Group_numb + aa[4]*0.1,height = 5)
result[[4]]
dev.off()
pdf(file = "Results/Virus/Community/Alluvial_chart_Group/Genus.pdf",
    width = 1.5 + wid*Group_numb + aa[5]*0.1,height = 5)
result[[5]]
dev.off()
## pdf(file = "Results/Virus/Community/Alluvial_chart_Group/Species.pdf",
##     width = 1.5 + wid*Group_numb + aa[6]*0.1,height = 5)
## result[[6]]
## dev.off()

#### barplot by sample
dir.create("Results/Virus/Community/Barplot")
source("Functions/Barplot.R")
result <- barplot.s(abundance.tax,group)
pdf(file = "Results/Virus/Community/Barplot/Phylum.pdf",
    width = 1.6 + 0.2*Sample_numb + aa[1]*0.1,height = 4.5)
result[[1]]
dev.off()
pdf(file = "Results/Virus/Community/Barplot/Class.pdf",
    width = 1.6 + 0.2*Sample_numb + aa[2]*0.1,height = 4.5)
result[[2]]
dev.off()
pdf(file = "Results/Virus/Community/Barplot/Order.pdf",
    width = 1.6 + 0.2*Sample_numb + aa[3]*0.1,height = 4.5)
result[[3]]
dev.off()
pdf(file = "Results/Virus/Community/Barplot/Family.pdf",
    width = 1.6 + 0.2*Sample_numb + aa[4]*0.1,height = 4.5)
result[[4]]
dev.off()
pdf(file = "Results/Virus/Community/Barplot/Genus.pdf",
    width = 1.6 + 0.2*Sample_numb + aa[5]*0.1,height = 4.5)
result[[5]]
dev.off()
## pdf(file = "Results/Virus/Community/Barplot/Species.pdf",
##     width = 1.6 + 0.2*Sample_numb + aa[6]*0.1,height = 4.5)
## result[[6]]
## dev.off()

#### bubble chart by sample
source("Functions/top30.R")
result <- top30(abundance.tax)
aa <- c(max(str_length(colnames(result[[1]]))),max(str_length(colnames(result[[2]]))),
        max(str_length(colnames(result[[3]]))),max(str_length(colnames(result[[4]]))),
        max(str_length(colnames(result[[5]]))),max(str_length(colnames(result[[6]]))))

dir.create("Results/Virus/Community/Bubble_chart")
source("Functions/Bubblechart.s.R")
result <- bubblechart.s(abundance.tax,group)
pdf(file = "Results/Virus/Community/Bubble_chart/Phylum.pdf",
    width = 0.22*Sample_numb + aa[1]*0.1 + 0.5,
    height = ifelse(length(unique(data.tax$Phylum)) > 30, 12.5, 0.35*length(unique(data.tax$Phylum)) + 2))
result[[1]]
dev.off()
pdf(file = "Results/Virus/Community/Bubble_chart/Class.pdf",
    width = 0.22*Sample_numb + aa[2]*0.1 + 0.5,
    height = ifelse(length(unique(data.tax$Class)) > 30, 12.5, 0.35*length(unique(data.tax$Class)) + 2))
result[[2]]
dev.off()
pdf(file = "Results/Virus/Community/Bubble_chart/Order.pdf",
    width = 0.22*Sample_numb + aa[3]*0.1 + 0.5,
    height = ifelse(length(unique(data.tax$Order)) > 30, 12.5, 0.35*length(unique(data.tax$Order)) + 2))
result[[3]]
dev.off()
pdf(file = "Results/Virus/Community/Bubble_chart/Family.pdf",
    width = 0.22*Sample_numb + aa[4]*0.1 + 0.5,
    height = ifelse(length(unique(data.tax$Family)) > 30, 12.5, 0.35*length(unique(data.tax$Family)) + 2))
result[[4]]
dev.off()
pdf(file = "Results/Virus/Community/Bubble_chart/Genus.pdf",
    width = 0.22*Sample_numb + aa[5]*0.1 + 0.5,
    height = ifelse(length(unique(data.tax$Genus)) > 30, 12.5, 0.35*length(unique(data.tax$Genus)) + 2))
result[[5]]
dev.off()
## pdf(file = "Results/Virus/Community/Bubble_chart/Species.pdf",
##     width = 0.22*Sample_numb + aa[6]*0.1 + 0.5,
##     height = ifelse(length(unique(data.tax$Species)) > 30, 12.5, 0.35*length(unique(data.tax$Species)) + 2))
## result[[6]]
## dev.off()

#### bubble chart by group
dir.create("Results/Virus/Community/Bubble_chart_Group")
source("Functions/Bubblechart.g.R")
result <- bubblechart.g(abundance.tax,group)
pdf(file = "Results/Virus/Community/Bubble_chart_Group/Phylum.pdf",
    width = 0.3*Group_numb + aa[1]*0.1 + 0.5,
    height = ifelse(length(unique(data.tax$Phylum)) > 30, 12.5, 0.35*length(unique(data.tax$Phylum)) + 2))
result[[1]]
dev.off()
pdf(file = "Results/Virus/Community/Bubble_chart_Group/Class.pdf",
    width = 0.3*Group_numb + aa[2]*0.1 + 0.5,
    height = ifelse(length(unique(data.tax$Class)) > 30, 12.5, 0.35*length(unique(data.tax$Class)) + 2))
result[[2]]
dev.off()
pdf(file = "Results/Virus/Community/Bubble_chart_Group/Order.pdf",
    width = 0.3*Group_numb + aa[3]*0.1 + 0.5,
    height = ifelse(length(unique(data.tax$Order)) > 30, 12.5, 0.35*length(unique(data.tax$Order)) + 2))
result[[3]]
dev.off()
pdf(file = "Results/Virus/Community/Bubble_chart_Group/Family.pdf",
    width = 0.3*Group_numb + aa[4]*0.1 + 0.5,
    height = ifelse(length(unique(data.tax$Family)) > 30, 12.5, 0.35*length(unique(data.tax$Family)) + 2))
result[[4]]
dev.off()
pdf(file = "Results/Virus/Community/Bubble_chart_Group/Genus.pdf",
    width = 0.3*Group_numb + aa[5]*0.1 + 0.5,
    height = ifelse(length(unique(data.tax$Genus)) > 30, 12.5, 0.35*length(unique(data.tax$Genus)) + 2))
result[[5]]
dev.off()
## pdf(file = "Results/Virus/Community/Bubble_chart_Group/Species.pdf",
##     width = 0.3*Group_numb + aa[6]*0.1 + 0.5,
##     height = ifelse(length(unique(data.tax$Species)) > 30, 12.5, 0.35*length(unique(data.tax$Species)) + 2))
## result[[6]]
## dev.off()

#### heatmap by sample
dir.create("Results/Virus/Community/Heatmap")
source("Functions/heatmap.s.R")
result <- heatmap.s(abundance.tax,group)
pdf(file = "Results/Virus/Community/Heatmap/Phylum.pdf",
    width = 0.22*Sample_numb + aa[1]*0.1 + 1.5,
    height = ifelse(length(unique(data.tax$Phylum)) > 30, 13, 0.35*length(unique(data.tax$Phylum)) + 3))
result[[1]]
dev.off()
pdf(file = "Results/Virus/Community/Heatmap/Class.pdf",
    width = 0.22*Sample_numb + aa[2]*0.1 + 1.5,
    height = ifelse(length(unique(data.tax$Class)) > 30, 13, 0.35*length(unique(data.tax$Class)) + 3))
result[[2]]
dev.off()
pdf(file = "Results/Virus/Community/Heatmap/Order.pdf",
    width = 0.22*Sample_numb + aa[3]*0.1 + 1.5,
    height = ifelse(length(unique(data.tax$Order)) > 30, 13, 0.35*length(unique(data.tax$Order)) + 3))
result[[3]]
dev.off()
pdf(file = "Results/Virus/Community/Heatmap/Family.pdf",
    width = 0.22*Sample_numb + aa[4]*0.1 + 1.5,
    height = ifelse(length(unique(data.tax$Family)) > 30, 13, 0.35*length(unique(data.tax$Family)) + 3))
result[[4]]
dev.off()
pdf(file = "Results/Virus/Community/Heatmap/Genus.pdf",
    width = 0.22*Sample_numb + aa[5]*0.1 + 1.5,
    height = ifelse(length(unique(data.tax$Genus)) > 30, 13, 0.35*length(unique(data.tax$Genus)) + 3))
result[[5]]
dev.off()
## pdf(file = "Results/Virus/Community/Heatmap/Species.pdf",
##     width = 0.22*Sample_numb + aa[6]*0.1 + 1.5,
##     height = ifelse(length(unique(data.tax$Species)) > 30, 13, 0.35*length(unique(data.tax$Species)) + 3))
## result[[6]]
## dev.off()

#### Heatmap by group
dir.create("Results/Virus/Community/Heatmap_Group")
source("Functions/heatmap.g.R")
result <- heatmap.g(abundance.tax,group)
pdf(file = "Results/Virus/Community/Heatmap_Group/Phylum.pdf",
    width = 0.25*Group_numb + aa[1]*0.1 + 1.5,
    height = ifelse(length(unique(data.tax$Phylum)) > 30, 13, 0.3*length(unique(data.tax$Phylum)) + 0.6))
result[[1]]
dev.off()
pdf(file = "Results/Virus/Community/Heatmap_Group/Class.pdf",
    width = 0.25*Group_numb + aa[2]*0.1 + 1.5,
    height = ifelse(length(unique(data.tax$Class)) > 30, 13, 0.3*length(unique(data.tax$Class)) + 0.6))
result[[2]]
dev.off()
pdf(file = "Results/Virus/Community/Heatmap_Group/Order.pdf",
    width = 0.25*Group_numb + aa[3]*0.1 + 1.5,
    height = ifelse(length(unique(data.tax$Order)) > 30, 13, 0.3*length(unique(data.tax$Order)) + 0.6))
result[[3]]
dev.off()
pdf(file = "Results/Virus/Community/Heatmap_Group/Family.pdf",
    width = 0.25*Group_numb + aa[4]*0.1 + 1.5,
    height = ifelse(length(unique(data.tax$Family)) > 30, 13, 0.3*length(unique(data.tax$Family)) + 0.6))
result[[4]]
dev.off()
pdf(file = "Results/Virus/Community/Heatmap_Group/Genus.pdf",
    width = 0.25*Group_numb + aa[5]*0.1 + 1.5,
    height = ifelse(length(unique(data.tax$Genus)) > 30, 13, 0.3*length(unique(data.tax$Genus)) + 0.6))
result[[5]]
dev.off()
## pdf(file = "Results/Virus/Community/Heatmap_Group/Species.pdf",
##     width = 0.25*Group_numb + aa[6]*0.1 + 1.5,
##     height = ifelse(length(unique(data.tax$Species)) > 30, 13, 0.3*length(unique(data.tax$Species)) + 0.6))
## result[[6]]
## dev.off()

## Alpha diversity
dir.create("Results/Virus/Alpha_diversity")
### Alpha diversity indices calculation
source("Functions/alpha_diversity_indices.R")
alpha.data <- floor(t(data[,1:(ncol(data)-6)]))
alpha <- Alpha_diversity_index(alpha.data)
alpha <- alpha[,-c(3,5)]
colnames(alpha) <- c("Observed_species","Chao1","ACE","Shannon","Simpson","Pielou_J",
                     "Good_coverage")
alpha$ACE[is.na(alpha$ACE)] <- alpha$Chao1[is.na(alpha$ACE)]
alpha1 <- data.frame(ID = rownames(alpha),alpha)
write.table(alpha1,"Results/Virus/Alpha_diversity/alpha_diversity_indices.txt",
            sep = "\t",row.names = FALSE)

### Differences of alpha diversity indices among different groups
dir.create("Results/Virus/Alpha_diversity/alpha_div")
source("Functions/Diff_alpha.R")
result <- diff.alpha(alpha1,group)
pdf(file = "Results/Virus/Alpha_diversity/alpha_div/Observed_species.pdf",
    width = wid*Group_numb,height = 3.6)
result[[1]]
dev.off()
pdf(file = "Results/Virus/Alpha_diversity/alpha_div/Chao1.pdf",
    width = wid*Group_numb,height = 3.6)
result[[2]]
dev.off()
pdf(file = "Results/Virus/Alpha_diversity/alpha_div/ACE.pdf",
    width = wid*Group_numb,height = 3.6)
result[[3]]
dev.off()
pdf(file = "Results/Virus/Alpha_diversity/alpha_div/Shannon.pdf",
    width = wid*Group_numb,height = 3.6)
result[[4]]
dev.off()
pdf(file = "Results/Virus/Alpha_diversity/alpha_div/Simpson.pdf",
    width = wid*Group_numb,height = 3.6)
result[[5]]
dev.off()
pdf(file = "Results/Virus/Alpha_diversity/alpha_div/Pielou_J.pdf",
    width = wid*Group_numb,height = 3.6)
result[[6]]
dev.off()
pdf(file = "Results/Virus/Alpha_diversity/alpha_div/Good_coverage.pdf",
    width = wid*Group_numb,height = 3.6)
result[[7]]
dev.off()

### Venn analysis by sample
dir.create("Results/Virus/Alpha_diversity/Flower")

#### Flower diagram
source("Functions/Flower.s.R")
pdf(file = "Results/Virus/Alpha_diversity/Flower/Flower_diagram.pdf",
    width = 8,height = 8)
flower.s(t(abundance))
dev.off()

#### Shared species list
source("Functions/Shared.species.id.s.R")
core_species_id_s <- shared.species.id.s(t(abundance))
core_species_id_s <- data.tax[core_species_id_s,]
core_species_id_s$Species_ID <- rownames(core_species_id_s)
core_species_id_s <- core_species_id_s[,c("Species_ID",colnames(core_species_id_s)[1:6])]
write.table(core_species_id_s,"Results/Virus/Alpha_diversity/Flower/Shared_species.txt",
            sep = "\t",row.names = FALSE)

#### Total abundance of shared species
source("Functions/Shared.species.abun.s.R")
result <- shared.species.abun.s(t(abundance),core_species_id_s,group)
pdf(file = "Results/Virus/Alpha_diversity/Flower/Shared_sepecies_total_abundance.pdf",
    width = 0.3*Sample_numb,height = 5)
result
dev.off()

#### Total abundance group differences
source("Functions/Shared.species.abun.g.R")
result <- shared.species.abun.g(t(abundance),core_species_id_s,group)
pdf(file = "Results/Virus/Alpha_diversity/Flower/Shared_sepecies_abundance_diff.pdf",
    width = wid*Group_numb,height = 4.5)
result
dev.off()

### Venn analysis by group
dir.create("Results/Virus/Alpha_diversity/Venn")

#### Venn diagram
aa <- levels(group$Group)
species <- list()
for (i in 1:length(aa)) {
    dd <- abundance[,group[group$Group == aa[1],1]]
    dd <- dd[which(rowSums(dd) > 0),]
    species[[aa[i]]] <- rownames(dd)
}

venn.diagram(species,filename = "Results/Virus/Alpha_diversity/Venn/venn.png",
             height = 5400,width = 5400,
             resolution = 600,imagetype = "png",units = "px",
             lwd = 2,lty = 1,fill = cbbPalette[1:length(aa)],cex = 1.5,
             cat.cex = 2,alpha = 0.8,margin = 0.05,fontface = 2,
             cat.fontface = 2, print.mode = c("raw","percent"))

source("Functions//overLapper.new.r")
pdf("Results/Virus/Alpha_diversity/Venn/venn.pdf",width=9,height=9,pointsize=16)
OLlist <- overLapper(setlist=species, sep="", type="vennsets",keepdups=FALSE)
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts,mymain="")
dev.off()

#### Shared species list
inter <- get.venn.partitions(species)
shared.species <- data.tax[unlist(inter$..values..[1]),]
shared.species$Species_ID <- rownames(shared.species)
shared.species <- shared.species[,c("Species_ID",colnames(shared.species)[1:6])]
write.table(shared.species,"Results/Virus/Alpha_diversity/Venn/Shared_species.txt",
            sep = "\t",row.names = FALSE)

#### Total abundance of shared species
result <- shared.species.abun.s(t(abundance),shared.species,group)
pdf(file = "Results/Virus/Alpha_diversity/Venn/Shared_sepecies_total_abundance.pdf",
    width = 0.3*Sample_numb,height = 5)
result
dev.off()

#### Total abundance group differences
result <- shared.species.abun.g(t(abundance),shared.species,group)
pdf(file = "Results/Virus/Alpha_diversity/Venn/Shared_sepecies_abundance_diff.pdf",
    width = wid*Group_numb,height = 4.5)
result
dev.off()

## Beta diveristy
dir.create("Results/Virus/Beta_diversity")

### PCoA
dir.create("Results/Virus/Beta_diversity/Distance")
dir.create("Results/Virus/Beta_diversity/PCoA")
source("Functions/PCoA.R")
result <- pcoa.community(t(abundance),group,species.type)
write.table(as.matrix(result[[1]]),
            "Results/Virus/Beta_diversity/Distance/distance_bray_curtis.txt",
            sep = "\t")
write.table(result[[2]],"Results/Virus/Beta_diversity/PCoA/pcoa.txt",
            sep = "\t",row.names = FALSE)
pdf(file = "Results/Virus/Beta_diversity/PCoA/Pcoa_group.pdf",
    width = 7.5,height = 5.4)
result[[3]]
dev.off()
pdf(file = "Results/Virus/Beta_diversity/PCoA/Pcoa_ellipse.pdf",
    width = 7.5,height = 5.4)
result[[4]]
dev.off()
pdf(file = "Results/Virus/Beta_diversity/PCoA/Pcoa_label.pdf",
    width = 7.5,height = 5.4)
result[[5]]
dev.off()

#### PCA
dir.create("Results/Virus/Beta_diversity/PCA")
source("Functions/PCA.R")
result <- pca.community(t(abundance),group,species.type)
write.table(result[[1]],"Results/Virus/Beta_diversity/PCA/PCA.txt",
            sep = "\t",row.names = FALSE)
pdf(file = "Results/Virus/Beta_diversity/PCA/PCA_group.pdf",
    width = 7.5,height = 5.4)
result[[2]]
dev.off()
pdf(file = "Results/Virus/Beta_diversity/PCA/PCA_ellipse.pdf",
    width = 7.5,height = 5.4)
result[[3]]
dev.off()
pdf(file = "Results/Virus/Beta_diversity/PCA/PCA_label.pdf",
    width = 7.5,height = 5.4)
result[[4]]
dev.off()

#### NMDS
dir.create("Results/Virus/Beta_diversity/NMDS")
source("Functions/NMDS.R")
result <- nmds.community(t(abundance),group,species.type)
write.table(result[[1]],"Results/Virus/Beta_diversity/NMDS/Microbe_stress.txt",
            sep = "\t",row.names = FALSE)
write.table(result[[2]],"Results/Virus/Beta_diversity/NMDS/NMDS.txt",
            sep = "\t",row.names = FALSE)
pdf(file = "Results/Virus/Beta_diversity/NMDS/NMDS_group.pdf",
    width = 7.5,height = 5.4)
result[[3]]
dev.off()
pdf(file = "Results/Virus/Beta_diversity/NMDS/NMDS_ellipse.pdf",
    width = 7.5,height = 5.4)
result[[4]]
dev.off()
pdf(file = "Results/Virus/Beta_diversity/NMDS/NMDS_label.pdf",
    width = 7.5,height = 5.4)
result[[5]]
dev.off()

### ANOSIM, Adonis, and MRPP tests
dir.create("Results/Virus/Beta_diversity/Beta_tests")
source("Functions/beta.test.R")
result <- beta.test(t(abundance),group)
sink("Results/Virus/Beta_diversity/Beta_tests/Adonis.txt")
result[[1]]
sink("Results/Virus/Beta_diversity/Beta_tests/ANOSIM.txt")
result[[2]]
sink("Results/Virus/Beta_diversity/Beta_tests/MRPP.txt")
result[[3]]

### Hcluster
dir.create("Results/Virus/Beta_diversity/Hcluster")
dist_bray <- vegdist(t(abundance), method = "bray")
hc <- hclust(dist_bray, method = "average")
hc_tre <- as.phylo(hc)
write.tree(hc_tre, "Results/Virus/Beta_diversity/Hcluster/Hcluster.tre")

source("Functions/plot.phylo2.r")
pdf("Results/Virus/Beta_diversity/Hcluster/Hcluster.pdf", 
    width = 7, height = 0.22*Sample_numb + 2)
plot.phylo2(hc_tre, type = "phylogram", cex = 1, adj = 1)
title(main = "UPGMA clustering tree", line = 2, cex.main = 1.5, font.main= 2)
dev.off()

### Hcluster barplot
dir.create("Results/Virus/Beta_diversity/Hcluster_bar")
source("Functions/top10 2.R")
result <- top10(abundance.tax)
aa <- c(max(str_length(colnames(result[[1]]))),max(str_length(colnames(result[[2]]))),
        max(str_length(colnames(result[[3]]))),max(str_length(colnames(result[[4]]))),
        max(str_length(colnames(result[[5]]))),max(str_length(colnames(result[[6]]))))

#### Phylum
pdf("Results/Virus/Beta_diversity/Hcluster_bar/Hclusterbar_Phylum.pdf", 
    width = ifelse(aa[1] < 26,14,aa[1]*0.15*3.62), height = 0.2*Sample_numb + 2)
lw <- unlist(strsplit("0.82-1.8-1","-"))
layout(matrix(c(1, 2, 3), 1, 3), widths = lw)

par(mar = c(3.5, 2, 4.5, 0))
plot.phylo2(hc_tre, type = "phylogram", cex = 1, adj = 1)
title(main = "Similarity", line = 2, cex.main = 2, font.main= 4)

par(mar = c(3, 0, 4, 0))
barplot(t(result[[1]]), space =0.3, horiz = T, border = NA, offset = 0.5, 
        col = cbbPalette[1:nrow(t(result[[1]]))], axes = FALSE, axisnames = FALSE)
title(main = "Taxonomic composition", line = 0, cex.main = 2, font.main= 2)

plot.new()
par(mar = c(4, 0, 4, 1))
leng <- rownames(t(result[[1]]))
legend(0.1, 1, legend = leng, fill = cbbPalette[1:nrow(t(result[[1]]))],
       bty = "n", cex = 2)
title(main = "Taxon - Phylum", line = 0, cex.main = 2, font.main= 2)
dev.off()

#### Class
pdf("Results/Virus/Beta_diversity/Hcluster_bar/Hclusterbar_Class.pdf", 
    width = ifelse(aa[2] < 20,14,aa[1]*0.2*3.62), height = 0.2*Sample_numb + 2)
lw <- unlist(strsplit("0.82-1.8-1","-"))
layout(matrix(c(1, 2, 3), 1, 3), widths = lw)

par(mar = c(3.5, 2, 4.5, 0))
plot.phylo2(hc_tre, type = "phylogram", cex = 1, adj = 1)
title(main = "Similarity", line = 2, cex.main = 2, font.main= 4)

par(mar = c(3, 0, 4, 0))
barplot(t(result[[2]]), space =0.3, horiz = T, border = NA, offset = 0.5, 
        col = cbbPalette[1:nrow(t(result[[2]]))], axes = FALSE, axisnames = FALSE)
title(main = "Taxonomic composition", line = 0, cex.main = 2, font.main= 2)

plot.new()
par(mar = c(4, 0, 4, 1))
leng <- rownames(t(result[[2]]))
legend(0.1, 1, legend = leng, fill = cbbPalette[1:nrow(t(result[[2]]))],
       bty = "n", cex = 2)
title(main = "Taxon - Class", line = 0, cex.main = 2, font.main= 2)
dev.off()

#### Order
pdf("Results/Virus/Beta_diversity/Hcluster_bar/Hclusterbar_Order.pdf", 
    width = ifelse(aa[3] < 20,14,aa[1]*0.2*3.62), height = 0.2*Sample_numb + 2)
lw <- unlist(strsplit("0.82-1.8-1","-"))
layout(matrix(c(1, 2, 3), 1, 3), widths = lw)

par(mar = c(3.5, 2, 4.5, 0))
plot.phylo2(hc_tre, type = "phylogram", cex = 1, adj = 1)
title(main = "Similarity", line = 2, cex.main = 2, font.main= 4)

par(mar = c(3, 0, 4, 0))
barplot(t(result[[3]]), space =0.3, horiz = T, border = NA, offset = 0.5, 
        col = cbbPalette[1:nrow(t(result[[3]]))], axes = FALSE, axisnames = FALSE)
title(main = "Taxonomic composition", line = 0, cex.main = 2, font.main= 2)

plot.new()
par(mar = c(4, 0, 4, 1))
leng <- rownames(t(result[[3]]))
legend(0.1, 1, legend = leng, fill = cbbPalette[1:nrow(t(result[[3]]))],
       bty = "n", cex = 2)
title(main = "Taxon - Order", line = 0, cex.main = 2, font.main= 2)
dev.off()

#### Family
pdf("Results/Virus/Beta_diversity/Hcluster_bar/Hclusterbar_Family.pdf", 
    width = ifelse(aa[4] < 20,14,aa[1]*0.2*3.62), height = 0.2*Sample_numb + 2)
lw <- unlist(strsplit("0.82-1.8-1","-"))
layout(matrix(c(1, 2, 3), 1, 3), widths = lw)

par(mar = c(3.5, 2, 4.5, 0))
plot.phylo2(hc_tre, type = "phylogram", cex = 1, adj = 1)
title(main = "Similarity", line = 2, cex.main = 2, font.main= 4)

par(mar = c(3, 0, 4, 0))
barplot(t(result[[4]]), space =0.3, horiz = T, border = NA, offset = 0.5, 
        col = cbbPalette[1:nrow(t(result[[4]]))], axes = FALSE, axisnames = FALSE)
title(main = "Taxonomic composition", line = 0, cex.main = 2, font.main= 2)

plot.new()
par(mar = c(4, 0, 4, 1))
leng <- rownames(t(result[[4]]))
legend(0.1, 1, legend = leng, fill = cbbPalette[1:nrow(t(result[[4]]))],
       bty = "n", cex = 2)
title(main = "Taxon - Family", line = 0, cex.main = 2, font.main= 2)
dev.off()

#### Genus
pdf("Results/Virus/Beta_diversity/Hcluster_bar/Hclusterbar_Genus.pdf", 
    width = ifelse(aa[5] < 20,14,aa[1]*0.2*3.62), height = 0.2*Sample_numb + 2)
lw <- unlist(strsplit("0.82-1.8-1","-"))
layout(matrix(c(1, 2, 3), 1, 3), widths = lw)

par(mar = c(3.5, 2, 4.5, 0))
plot.phylo2(hc_tre, type = "phylogram", cex = 1, adj = 1)
title(main = "Similarity", line = 2, cex.main = 2, font.main= 4)

par(mar = c(3, 0, 4, 0))
barplot(t(result[[5]]), space =0.3, horiz = T, border = NA, offset = 0.5, 
        col = cbbPalette[1:nrow(t(result[[5]]))], axes = FALSE, axisnames = FALSE)
title(main = "Taxonomic composition", line = 0, cex.main = 2, font.main= 2)

plot.new()
par(mar = c(4, 0, 4, 1))
leng <- rownames(t(result[[5]]))
legend(0.1, 1, legend = leng, fill = cbbPalette[1:nrow(t(result[[5]]))],
       bty = "n", cex = 2)
title(main = "Taxon - Genus", line = 0, cex.main = 2, font.main= 2)
dev.off()

#### Species
## pdf("Results/Virus/Beta_diversity/Hcluster_bar/Hclusterbar_Species.pdf", 
##     width = ifelse(aa[6] < 20,14,aa[1]*0.2*3.62), height = 0.2*Sample_numb + 2)
## lw <- unlist(strsplit("0.82-1.8-1","-"))
## layout(matrix(c(1, 2, 3), 1, 3), widths = lw)

## par(mar = c(3.5, 2, 4.5, 0))
## plot.phylo2(hc_tre, type = "phylogram", cex = 1, adj = 1)
## title(main = "Similarity", line = 2, cex.main = 2, font.main= 4)

## par(mar = c(3, 0, 4, 0))
## barplot(t(result[[6]]), space =0.3, horiz = T, border = NA, offset = 0.5, 
##         col = cbbPalette[1:nrow(t(result[[6]]))], axes = FALSE, axisnames = FALSE)
## title(main = "Taxonomic composition", line = 0, cex.main = 2, font.main= 2)

## plot.new()
## par(mar = c(4, 0, 4, 1))
## leng <- rownames(t(result[[6]]))
## legend(0.1, 1, legend = leng, fill = cbbPalette[1:nrow(t(result[[6]]))],
##        bty = "n", cex = 1.5)
## title(main = "Taxon - Species", line = 0, cex.main = 2, font.main= 2)
## dev.off()

## Difference analysis
dir.create("Results/Virus/Diff_analysis")

### Wilcox rank-sum test
dir.create("Results/Virus/Diff_analysis/Wilcox")
source("Functions/wilcox.R")
source("Functions/wilcox.abun.R")
result <- tax.summary(abundance.tax)

#### Phylum
dir.create("Results/Virus/Diff_analysis/Wilcox/Phylum")
set.seed(1111)
wilcox.result <- wilcox.biomarker(result[[1]],group)
write.table(wilcox.result[[1]],
            "Results/Virus/Diff_analysis/Wilcox/Phylum/wilcox_all.txt",
            sep = "\t",row.names = FALSE)
write.table(wilcox.result[[2]],
            "Results/Virus/Diff_analysis/Wilcox/Phylum/wilcox_biomarker.txt",
            sep = "\t",row.names = FALSE)
if (nrow(wilcox.result[[2]]) > 0) {
    Phylum <- wilcox.abun(result[[1]],wilcox.result[[2]],group)
    pdf("Results/Virus/Diff_analysis/Wilcox/Phylum/wilcox_biomarker_abun.pdf",
        width = 10,height = 1.8 + 0.18*nrow(wilcox.result[[2]]))
    print(Phylum[[1]])
    dev.off()
    write.table(Phylum[[2]],
                "Results/Virus/Diff_analysis/Wilcox/Phylum/Group_mean_abun_biomarker.txt",
                sep = "\t",row.names = FALSE)
    write.table(Phylum[[3]],
                "Results/Virus/Diff_analysis/Wilcox/Phylum/Abun_change_biomarker.txt",
                sep = "\t",row.names = FALSE)
}

#### Class
dir.create("Results/Virus/Diff_analysis/Wilcox/Class")
set.seed(1111)
wilcox.result <- wilcox.biomarker(result[[2]],group)
write.table(wilcox.result[[1]],
            "Results/Virus/Diff_analysis/Wilcox/Class/wilcox_all.txt",
            sep = "\t",row.names = FALSE)
write.table(wilcox.result[[2]],
            "Results/Virus/Diff_analysis/Wilcox/Class/wilcox_biomarker.txt",
            sep = "\t",row.names = FALSE)
if (nrow(wilcox.result[[2]]) > 0) {
    Class <- wilcox.abun(result[[2]],wilcox.result[[2]],group)
    pdf("Results/Virus/Diff_analysis/Wilcox/Class/wilcox_biomarker_abun.pdf",
        width = 10,height = 1.8 + 0.18*nrow(wilcox.result[[2]]))
    print(Class[[1]])
    dev.off()
    write.table(Class[[2]],
                "Results/Virus/Diff_analysis/Wilcox/Class/Group_mean_abun_biomarker.txt",
                sep = "\t",row.names = FALSE)
    write.table(Class[[3]],
                "Results/Virus/Diff_analysis/Wilcox/Class/Abun_change_biomarker.txt",
                sep = "\t",row.names = FALSE)
}

#### Order
dir.create("Results/Virus/Diff_analysis/Wilcox/Order")
set.seed(1111)
wilcox.result <- wilcox.biomarker(result[[3]],group)
write.table(wilcox.result[[1]],
            "Results/Virus/Diff_analysis/Wilcox/Order/wilcox_all.txt",
            sep = "\t",row.names = FALSE)
write.table(wilcox.result[[2]],
            "Results/Virus/Diff_analysis/Wilcox/Order/wilcox_biomarker.txt",
            sep = "\t",row.names = FALSE)
if (nrow(wilcox.result[[2]]) > 0) {
    Order <- wilcox.abun(result[[3]],wilcox.result[[2]],group)
    pdf("Results/Virus/Diff_analysis/Wilcox/Order/wilcox_biomarker_abun.pdf",
        width = 10,height = 1.8 + 0.18*nrow(wilcox.result[[2]]))
    print(Order[[1]])
    dev.off()
    write.table(Order[[2]],
                "Results/Virus/Diff_analysis/Wilcox/Order/Group_mean_abun_biomarker.txt",
                sep = "\t",row.names = FALSE)
    write.table(Order[[3]],
                "Results/Virus/Diff_analysis/Wilcox/Order/Abun_change_biomarker.txt",
                sep = "\t",row.names = FALSE)
}

#### Family
dir.create("Results/Virus/Diff_analysis/Wilcox/Family")
set.seed(1111)
wilcox.result <- wilcox.biomarker(result[[4]],group)
write.table(wilcox.result[[1]],
            "Results/Virus/Diff_analysis/Wilcox/Family/wilcox_all.txt",
            sep = "\t",row.names = FALSE)
write.table(wilcox.result[[2]],
            "Results/Virus/Diff_analysis/Wilcox/Family/wilcox_biomarker.txt",
            sep = "\t",row.names = FALSE)
if (nrow(wilcox.result[[2]]) > 0) {
    Family <- wilcox.abun(result[[4]],wilcox.result[[2]],group)
    pdf("Results/Virus/Diff_analysis/Wilcox/Family/wilcox_biomarker_abun.pdf",
        width = 10,height = 1.8 + 0.18*nrow(wilcox.result[[2]]))
    print(Family[[1]])
    dev.off()
    write.table(Family[[2]],
                "Results/Virus/Diff_analysis/Wilcox/Family/Group_mean_abun_biomarker.txt",
                sep = "\t",row.names = FALSE)
    write.table(Family[[3]],
                "Results/Virus/Diff_analysis/Wilcox/Family/Abun_change_biomarker.txt",
                sep = "\t",row.names = FALSE)
}

#### Genus
dir.create("Results/Virus/Diff_analysis/Wilcox/Genus")
set.seed(1111)
wilcox.result <- wilcox.biomarker(result[[5]],group)
write.table(wilcox.result[[1]],
            "Results/Virus/Diff_analysis/Wilcox/Genus/wilcox_all.txt",
            sep = "\t",row.names = FALSE)
write.table(wilcox.result[[2]],
            "Results/Virus/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
            sep = "\t",row.names = FALSE)
if (nrow(wilcox.result[[2]]) > 0) {
    Genus <- wilcox.abun(result[[5]],wilcox.result[[2]],group)
    pdf("Results/Virus/Diff_analysis/Wilcox/Genus/wilcox_biomarker_abun.pdf",
        width = 10,height = 1.8 + 0.18*nrow(wilcox.result[[2]]))
    print(Genus[[1]])
    dev.off()
    write.table(Genus[[2]],
                "Results/Virus/Diff_analysis/Wilcox/Genus/Group_mean_abun_biomarker.txt",
                sep = "\t",row.names = FALSE)
    write.table(Genus[[3]],
                "Results/Virus/Diff_analysis/Wilcox/Genus/Abun_change_biomarker.txt",
                sep = "\t",row.names = FALSE)
}

#### Species
## dir.create("Results/Virus/Diff_analysis/Wilcox/Species")
## set.seed(1111)
## wilcox.result <- wilcox.biomarker(result[[6]],group)
## write.table(wilcox.result[[1]],
##             "Results/Virus/Diff_analysis/Wilcox/Species/wilcox_all.txt",
##             sep = "\t",row.names = FALSE)
## write.table(wilcox.result[[2]],
##             "Results/Virus/Diff_analysis/Wilcox/Species/wilcox_biomarker.txt",
##             sep = "\t",row.names = FALSE)
## if (nrow(wilcox.result[[2]]) > 0) {
##     Species <- wilcox.abun(result[[6]],wilcox.result[[2]],group)
##     pdf("Results/Virus/Diff_analysis/Wilcox/Species/wilcox_biomarker_abun.pdf",
##         width = 10,height = 1.8 + 0.18*nrow(wilcox.result[[2]]))
##     print(Species[[1]])
##     dev.off()
##     write.table(Species[[2]],
##                 "Results/Virus/Diff_analysis/Wilcox/Species/Group_mean_abun_biomarker.txt",
##                 sep = "\t",row.names = FALSE)
##     write.table(Species[[3]],
##                 "Results/Virus/Diff_analysis/Wilcox/Species/Abun_change_biomarker.txt",
##                 sep = "\t",row.names = FALSE)
## }

### Random Forest for biomarker identification
dir.create("Results/Virus/Diff_analysis/RandomForest")
source("Functions/RandomForest.R")

#### Phylum
dir.create("Results/Virus/Diff_analysis/RandomForest/Phylum")
set.seed(1111)
RF.result <- Random.forest.biomarker(result[[1]],group)
pdf("Results/Virus/Diff_analysis/RandomForest/Phylum/classification.pdf",
    height = 3.5,width = 3.5)
RF.result[[1]]
dev.off()
importance.otu <- RF.result[[2]]
importance.otu$ID <- rownames(importance.otu)
importance.otu <- importance.otu[,c("ID",colnames(importance.otu)[1:(ncol(importance.otu)-1)])]
write.table(importance.otu,
            "Results/Virus/Diff_analysis/RandomForest/Phylum/importance.txt",
            sep = "\t",row.names = FALSE)
if (nrow(importance.otu) > 0) {
    pdf("Results/Virus/Diff_analysis/RandomForest/Phylum//RandomForest_biomarkder.pdf",
        height = nrow(importance.otu)*0.2 + 1.85,width = 9)
    print(RF.result[[3]])
    dev.off()
    write.table(RF.result[[4]],
                "Results/Virus/Diff_analysis/RandomForest/Phylum/biomarker_abundance.txt",
                sep = "\t") 
}

#### Class
dir.create("Results/Virus/Diff_analysis/RandomForest/Class")
set.seed(1111)
RF.result <- Random.forest.biomarker(result[[2]],group)
pdf("Results/Virus/Diff_analysis/RandomForest/Class/classification.pdf",
    height = 3.5,width = 3.5)
RF.result[[1]]
dev.off()
importance.otu <- RF.result[[2]]
importance.otu$ID <- rownames(importance.otu)
importance.otu <- importance.otu[,c("ID",colnames(importance.otu)[1:(ncol(importance.otu)-1)])]
write.table(importance.otu,
            "Results/Virus/Diff_analysis/RandomForest/Class/importance.txt",
            sep = "\t",row.names = FALSE)
if (nrow(importance.otu) > 0) {
    pdf("Results/Virus/Diff_analysis/RandomForest/Class//RandomForest_biomarkder.pdf",
        height = nrow(importance.otu)*0.2 + 1.85,width = 9)
    print(RF.result[[3]])
    dev.off()
    write.table(RF.result[[4]],
                "Results/Virus/Diff_analysis/RandomForest/Class/biomarker_abundance.txt",
                sep = "\t") 
}

#### Order
dir.create("Results/Virus/Diff_analysis/RandomForest/Order")
set.seed(1111)
RF.result <- Random.forest.biomarker(result[[3]],group)
pdf("Results/Virus/Diff_analysis/RandomForest/Order/classification.pdf",
    height = 3.5,width = 3.5)
RF.result[[1]]
dev.off()
importance.otu <- RF.result[[2]]
importance.otu$ID <- rownames(importance.otu)
importance.otu <- importance.otu[,c("ID",colnames(importance.otu)[1:(ncol(importance.otu)-1)])]
write.table(importance.otu,
            "Results/Virus/Diff_analysis/RandomForest/Order/importance.txt",
            sep = "\t",row.names = FALSE)
if (nrow(importance.otu) > 0) {
    pdf("Results/Virus/Diff_analysis/RandomForest/Order//RandomForest_biomarkder.pdf",
        height = nrow(importance.otu)*0.2 + 1.85,width = 9)
    print(RF.result[[3]])
    dev.off()
    write.table(RF.result[[4]],
                "Results/Virus/Diff_analysis/RandomForest/Order/biomarker_abundance.txt",
                sep = "\t") 
}

#### Family
dir.create("Results/Virus/Diff_analysis/RandomForest/Family")
set.seed(1111)
RF.result <- Random.forest.biomarker(result[[4]],group)
pdf("Results/Virus/Diff_analysis/RandomForest/Family/classification.pdf",
    height = 3.5,width = 3.5)
RF.result[[1]]
dev.off()
importance.otu <- RF.result[[2]]
importance.otu$ID <- rownames(importance.otu)
importance.otu <- importance.otu[,c("ID",colnames(importance.otu)[1:(ncol(importance.otu)-1)])]
write.table(importance.otu,
            "Results/Virus/Diff_analysis/RandomForest/Family/importance.txt",
            sep = "\t",row.names = FALSE)
if (nrow(importance.otu) > 0) {
    pdf("Results/Virus/Diff_analysis/RandomForest/Family//RandomForest_biomarkder.pdf",
        height = nrow(importance.otu)*0.2 + 1.85,width = 9)
    print(RF.result[[3]])
    dev.off()
    write.table(RF.result[[4]],
                "Results/Virus/Diff_analysis/RandomForest/Family/biomarker_abundance.txt",
                sep = "\t") 
}

#### Genus
dir.create("Results/Virus/Diff_analysis/RandomForest/Genus")
set.seed(1111)
RF.result <- Random.forest.biomarker(result[[5]],group)
pdf("Results/Virus/Diff_analysis/RandomForest/Genus/classification.pdf",
    height = 3.5,width = 3.5)
RF.result[[1]]
dev.off()
importance.otu <- RF.result[[2]]
importance.otu$ID <- rownames(importance.otu)
importance.otu <- importance.otu[,c("ID",colnames(importance.otu)[1:(ncol(importance.otu)-1)])]
write.table(importance.otu,
            "Results/Virus/Diff_analysis/RandomForest/Genus/importance.txt",
            sep = "\t",row.names = FALSE)
if (nrow(importance.otu) > 0) {
    pdf("Results/Virus/Diff_analysis/RandomForest/Genus//RandomForest_biomarkder.pdf",
        height = nrow(importance.otu)*0.2 + 1.85,width = 9)
    print(RF.result[[3]])
    dev.off()
    write.table(RF.result[[4]],
                "Results/Virus/Diff_analysis/RandomForest/Genus/biomarker_abundance.txt",
                sep = "\t") 
}

#### Species
## dir.create("Results/Virus/Diff_analysis/RandomForest/Species")
## set.seed(1111)
## RF.result <- Random.forest.biomarker(result[[6]],group)
## pdf("Results/Virus/Diff_analysis/RandomForest/Species/classification.pdf",
##     height = 3.5,width = 3.5)
## RF.result[[1]]
## dev.off()
## importance.otu <- RF.result[[2]]
## importance.otu$ID <- rownames(importance.otu)
## importance.otu <- importance.otu[,c("ID",colnames(importance.otu)[1:(ncol(importance.otu)-1)])]
## write.table(importance.otu,
##             "Results/Virus/Diff_analysis/RandomForest/Species/importance.txt",
##             sep = "\t",row.names = FALSE)
## if (nrow(importance.otu[importance.otu$MeanDecreaseAccuracy > 1.5,]) > 0) {
##     pdf("Results/Virus/Diff_analysis/RandomForest/Species//RandomForest_biomarkder.pdf",
##         height = nrow(importance.otu[importance.otu$MeanDecreaseAccuracy > 1.5,])*0.2 + 1.85,
##         width = 6 + ncol(result[[1]])*0.1)
##     print(RF.result[[3]])
##     dev.off()
##     write.table(RF.result[[4]],
##                "Results/Virus/Diff_analysis/RandomForest/Species/biomarker_abundance.txt",
##                 sep = "\t") 
## }

### LEfSe
dir.create("Results/Virus/Diff_analysis/LEfSe")
source("Functions/lefse.R")

Lefse.result <- lefse.1(abundance.tax,group)
write.table(Lefse.result,"Results/Virus/Diff_analysis/LEfSe/Lefse_result.txt",
            sep = "\t")
if (nrow(Lefse.result) > 0) {
    source("Functions/lefse_bar.R")
    p <- lefse_bar(Lefse.result)
    pdf("Results/Virus/Diff_analysis/LEfSe/Lefse_bar.pdf",
        width = 6,height = 1.15 + 0.3*nrow(Lefse.result))
    print(p)
    dev.off()
}

## Diagnose
dir.create("Results/Virus/Diagnose")

### Wilcox rank-sum test
dir.create("Results/Virus/Diagnose/Wilcox")
source("Functions/diagnose.R")

#### Phylum
data <- read.table("Results/Virus/Diff_analysis/Wilcox/Phylum/wilcox_biomarker.txt",
                   header = TRUE,sep = "\t")
if (nrow(data) > 0) {
    dir.create("Results/Virus/Diagnose/Wilcox/Phylum")
    Phylum <- diagnose(data,result[[1]],group)
    write.table(Phylum[[1]],"Results/Virus/Diagnose/Wilcox/Phylum/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Virus/Diagnose/Wilcox/Phylum/index.pdf",
        width = 2.4,height = 3.6)
    print(Phylum[[2]])
    dev.off()
    pdf("Results/Virus/Diagnose/Wilcox/Phylum/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(Phylum[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
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

#### Class
data <- read.table("Results/Virus/Diff_analysis/Wilcox/Class/wilcox_biomarker.txt",
                   header = TRUE,sep = "\t")
if (nrow(data) > 0) {
    dir.create("Results/Virus/Diagnose/Wilcox/Class")
    Class <- diagnose(data,result[[2]],group)
    write.table(Class[[1]],"Results/Virus/Diagnose/Wilcox/Class/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Virus/Diagnose/Wilcox/Class/index.pdf",
        width = 2.4,height = 3.6)
    print(Class[[2]])
    dev.off()
    pdf("Results/Virus/Diagnose/Wilcox/Class/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(Class[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
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

#### Order
data <- read.table("Results/Virus/Diff_analysis/Wilcox/Order/wilcox_biomarker.txt",
                   header = TRUE,sep = "\t")
if (nrow(data) > 0) {
    dir.create("Results/Virus/Diagnose/Wilcox/Order")
    Order <- diagnose(data,result[[3]],group)
    write.table(Order[[1]],"Results/Virus/Diagnose/Wilcox/Order/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Virus/Diagnose/Wilcox/Order/index.pdf",
        width = 2.4,height = 3.6)
    print(Order[[2]])
    dev.off()
    pdf("Results/Virus/Diagnose/Wilcox/Order/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(Order[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
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

#### Family
data <- read.table("Results/Virus/Diff_analysis/Wilcox/Family/wilcox_biomarker.txt",
                   header = TRUE,sep = "\t")
if (nrow(data) > 0) {
    dir.create("Results/Virus/Diagnose/Wilcox/Family")
    Family <- diagnose(data,result[[4]],group)
    write.table(Family[[1]],"Results/Virus/Diagnose/Wilcox/Family/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Virus/Diagnose/Wilcox/Family/index.pdf",
        width = 2.4,height = 3.6)
    print(Family[[2]])
    dev.off()
    pdf("Results/Virus/Diagnose/Wilcox/Family/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(Family[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
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

#### Genus
data <- read.table("Results/Virus/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                   header = TRUE,sep = "\t")
if (nrow(data) > 0) {
    dir.create("Results/Virus/Diagnose/Wilcox/Genus")
    Genus <- diagnose(data,result[[5]],group)
    write.table(Genus[[1]],"Results/Virus/Diagnose/Wilcox/Genus/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Virus/Diagnose/Wilcox/Genus/index.pdf",
        width = 2.4,height = 3.6)
    print(Genus[[2]])
    dev.off()
    pdf("Results/Virus/Diagnose/Wilcox/Genus/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(Genus[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
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

#### Species
## data <- read.table("Results/Virus/Diff_analysis/Wilcox/Species/wilcox_biomarker.txt",
##                    header = TRUE,sep = "\t")
## if (nrow(data) > 0) {
##     dir.create("Results/Virus/Diagnose/Wilcox/Species")
##     Species <- diagnose(data,result[[6]],group)
##     write.table(Species[[1]],"Results/Virus/Diagnose/Wilcox/Species/index.txt",
##                 sep = "\t",row.names = FALSE)
##     pdf("Results/Virus/Diagnose/Wilcox/Species/index.pdf",
##         width = 2.4,height = 3.6)
##     print(Species[[2]])
##     dev.off()
##     pdf("Results/Virus/Diagnose/Wilcox/Species/Roc.pdf",width = 10,height = 10)
##     par(xpd = TRUE)
##     plot(Species[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
##          print.thres.best.method = "youden",print.auc.cex = 2.5,print.thres.pch = 19,
##          print.thres.pattern.cex = 2.5,auc.polygon = TRUE,auc.polygon.col = "grey90",
##          print.auc.x = 0.35,print.auc.y = 0.1,auc.polygon.border = "black",
##          max.auc.polygon = TRUE,axes = FALSE,xlab = "",ylab = "",mar = c(6,6,2,2)+.1)
##     axis(side = 1,at = c(1.0,0.5,0),labels = c("","",""),lwd = 2,lwd.ticks = 4,line = -1)
##     text(x = c(1.0,0.5,0),y = -0.08,labels = c("1.0","0.5","0.0"),cex = 2)
##     mtext(side = 1,text = "Speciticity",line = 3,cex = 2.7)
##     axis(side = 2,at = c(0,0.2,0.4,0.6,0.8,1.0),labels = c("","","","","",""),lwd = 2,lwd.ticks = 4,line = -2)
##     text(x = 1.09,y = c(0,0.2,0.4,0.6,0.8,1.0),labels = c("0.0","0.2","0.4","0.6","0.8","1.0"),cex = 2)
##     mtext(side = 2,text = "Sensitivity",line = 2,cex = 2.7)
##     segments(x0 = 1.04,y0 = -0.04,x1 = -0.04,y1 = -0.04,lwd = 2)
##     segments(x0 = 1.04,y0 = 1.04,x1 = -0.04,y1 = 1.04,lwd = 2)
##     segments(x0 = 1.04,y0 = -0.04,x1 = 1.04,y1 = 1.04,lwd = 2)
##     segments(x0 = -0.04,y0 = -0.04,x1 = -0.04,y1 = 1.04,lwd = 2)
##     dev.off()
## }

### RandomForest
dir.create("Results/Virus/Diagnose/RandomForest")

#### Phylum
data <- read.table("Results/Virus/Diff_analysis/RandomForest/Phylum/importance.txt",
                   header = TRUE,sep = "\t")
colnames(data)[1] <- "V1"
if (nrow(data) > 0) {
    dir.create("Results/Virus/Diagnose/RandomForest/Phylum")
    Phylum <- diagnose(data,result[[1]],group)
    write.table(Phylum[[1]],"Results/Virus/Diagnose/RandomForest/Phylum/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Virus/Diagnose/RandomForest/Phylum/index.pdf",
        width = 2.4,height = 3.6)
    print(Phylum[[2]])
    dev.off()
    pdf("Results/Virus/Diagnose/RandomForest/Phylum/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(Phylum[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
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

#### Class
data <- read.table("Results/Virus/Diff_analysis/RandomForest/Class/importance.txt",
                   header = TRUE,sep = "\t")
colnames(data)[1] <- "V1"
if (nrow(data) > 0) {
    dir.create("Results/Virus/Diagnose/RandomForest/Class")
    Class <- diagnose(data,result[[2]],group)
    write.table(Class[[1]],"Results/Virus/Diagnose/RandomForest/Class/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Virus/Diagnose/RandomForest/Class/index.pdf",
        width = 2.4,height = 3.6)
    print(Class[[2]])
    dev.off()
    pdf("Results/Virus/Diagnose/RandomForest/Class/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(Class[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
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

#### Order
data <- read.table("Results/Virus/Diff_analysis/RandomForest/Order/importance.txt",
                   header = TRUE,sep = "\t")
colnames(data)[1] <- "V1"
if (nrow(data) > 0) {
    dir.create("Results/Virus/Diagnose/RandomForest/Order")
    Order <- diagnose(data,result[[3]],group)
    write.table(Order[[1]],"Results/Virus/Diagnose/RandomForest/Order/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Virus/Diagnose/RandomForest/Order/index.pdf",
        width = 2.4,height = 3.6)
    print(Order[[2]])
    dev.off()
    pdf("Results/Virus/Diagnose/RandomForest/Order/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(Order[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
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

#### Family
data <- read.table("Results/Virus/Diff_analysis/RandomForest/Family/importance.txt",
                   header = TRUE,sep = "\t")
colnames(data)[1] <- "V1"
if (nrow(data) > 0) {
    dir.create("Results/Virus/Diagnose/RandomForest/Family")
    Family <- diagnose(data,result[[4]],group)
    write.table(Family[[1]],"Results/Virus/Diagnose/RandomForest/Family/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Virus/Diagnose/RandomForest/Family/index.pdf",
        width = 2.4,height = 3.6)
    print(Family[[2]])
    dev.off()
    pdf("Results/Virus/Diagnose/RandomForest/Family/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(Family[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
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

#### Genus
data <- read.table("Results/Virus/Diff_analysis/RandomForest/Genus/importance.txt",
                   header = TRUE,sep = "\t")
colnames(data)[1] <- "V1"
if (nrow(data) > 0) {
    dir.create("Results/Virus/Diagnose/RandomForest/Genus")
    Genus <- diagnose(data,result[[5]],group)
    write.table(Genus[[1]],"Results/Virus/Diagnose/RandomForest/Genus/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Virus/Diagnose/RandomForest/Genus/index.pdf",
        width = 2.4,height = 3.6)
    print(Genus[[2]])
    dev.off()
    pdf("Results/Virus/Diagnose/RandomForest/Genus/Roc.pdf",width = 10,height = 10)
    par(xpd = TRUE)
    plot(Genus[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
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

#### Species
## data <- read.table("Results/Virus/Diff_analysis/RandomForest/Species/importance.txt",
##                    header = TRUE,sep = "\t")
## colnames(data)[1] <- "V1"
## if (nrow(data) > 0) {
##     dir.create("Results/Virus/Diagnose/RandomForest/Species")
##     Species <- diagnose(data,result[[6]],group)
##     write.table(Species[[1]],"Results/Virus/Diagnose/RandomForest/Species/index.txt",
##                 sep = "\t",row.names = FALSE)
##     pdf("Results/Virus/Diagnose/RandomForest/Species/index.pdf",
##         width = 2.4,height = 3.6)
##     print(Species[[2]])
##     dev.off()
##     pdf("Results/Virus/Diagnose/RandomForest/Species/Roc.pdf",width = 10,height = 10)
##     par(xpd = TRUE)
##     plot(Species[[3]],print.auc=T,print.thres=T,lwd = 5,identity = FALSE,print.thres.cex = 4,
##          print.thres.best.method = "youden",print.auc.cex = 2.5,print.thres.pch = 19,
##          print.thres.pattern.cex = 2.5,auc.polygon = TRUE,auc.polygon.col = "grey90",
##          print.auc.x = 0.35,print.auc.y = 0.1,auc.polygon.border = "black",
##          max.auc.polygon = TRUE,axes = FALSE,xlab = "",ylab = "",mar = c(6,6,2,2)+.1)
##     axis(side = 1,at = c(1.0,0.5,0),labels = c("","",""),lwd = 2,lwd.ticks = 4,line = -1)
##     text(x = c(1.0,0.5,0),y = -0.08,labels = c("1.0","0.5","0.0"),cex = 2)
##     mtext(side = 1,text = "Speciticity",line = 3,cex = 2.7)
##     axis(side = 2,at = c(0,0.2,0.4,0.6,0.8,1.0),labels = c("","","","","",""),lwd = 2,lwd.ticks = 4,line = -2)
##     text(x = 1.09,y = c(0,0.2,0.4,0.6,0.8,1.0),labels = c("0.0","0.2","0.4","0.6","0.8","1.0"),cex = 2)
##     mtext(side = 2,text = "Sensitivity",line = 2,cex = 2.7)
##     segments(x0 = 1.04,y0 = -0.04,x1 = -0.04,y1 = -0.04,lwd = 2)
##     segments(x0 = 1.04,y0 = 1.04,x1 = -0.04,y1 = 1.04,lwd = 2)
##     segments(x0 = 1.04,y0 = -0.04,x1 = 1.04,y1 = 1.04,lwd = 2)
##     segments(x0 = -0.04,y0 = -0.04,x1 = -0.04,y1 = 1.04,lwd = 2)
##     dev.off()
## }

### LEfse
dir.create("Results/Virus/Diagnose/LEfSe")
data <- read.table("Results/Virus/Diff_analysis/LEfSe/Lefse_result.txt",
                   header = TRUE,sep = "\t",row.names = 1)
data$V1 <- rownames(data)
data <- data[grep("g__",data$V1),]
data$V1 <- gsub("g__","",data$V1)
if (nrow(data) > 0) {
    Species <- diagnose(data,result[[5]],group)
    write.table(Species[[1]],"Results/Virus/Diagnose/LEfSe/index.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Virus/Diagnose/LEfSe/index.pdf",
        width = 2.4,height = 3.6)
    print(Species[[2]])
    dev.off()
    pdf("Results/Virus/Diagnose/LEfSe/Roc.pdf",width = 10,height = 10)
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






