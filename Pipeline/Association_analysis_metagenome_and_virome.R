### Pipeline for standard analyses of virome based on pathseq results
### Based on R v4.0.2
setwd("../")
library(vegan)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(ape)
library(dplyr)
library(RColorBrewer)
library(patchwork)
library(linkET)
library(ggraph)
library(tidygraph)
library(pheatmap)
library(psych)
library(VennDiagram)

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
dir.create("Results/")
dir.create("Results/Association_analysis_virome")

### Input bacteria
bacteria <- read.table(file = "Input/bacteria.tpm.xls",header = FALSE,row.names = 1,
                   sep = "\t",quote = "")

### bacteria processing
colnames(bacteria) <- bacteria[1,]
bacteria <- bacteria[-1,]
for (i in 1:(ncol(bacteria)-1)) {
    bacteria[,i] <- as.numeric(bacteria[,i])
}

bacteria$Phylum <- gsub(".*p__","",bacteria$taxonomy)
bacteria$Phylum <- gsub("\\;c__.*","",bacteria$Phylum)
bacteria$Class <- gsub(".*c__","",bacteria$taxonomy)
bacteria$Class <- gsub("\\;o__.*","",bacteria$Class)
bacteria$Order <- gsub(".*o__","",bacteria$taxonomy)
bacteria$Order <- gsub("\\;f__.*","",bacteria$Order)
bacteria$Family <- gsub(".*f__","",bacteria$taxonomy)
bacteria$Family <- gsub("\\;g__.*","",bacteria$Family)
bacteria$Genus <- gsub(".*g__","",bacteria$taxonomy)
bacteria$Genus <- gsub("\\;s__.*","",bacteria$Genus)
bacteria$Species <- gsub(".*s__","",bacteria$taxonomy)
bacteria <- subset(bacteria,select = -taxonomy)
bacteria.tax <- bacteria[,(ncol(bacteria)-5):ncol(bacteria)]
bacteria.tax[bacteria.tax == "norank"] <- "Unclassified"
bacteria.tax[bacteria.tax == "uncultured bacterium"] <- "Unclassified"
bacteria.tax[bacteria.tax == "uncultured"] <- "Unclassified"
for (i in 1:nrow(bacteria.tax)) {
    bacteria.tax[i,grepl("d__Bacteria",bacteria.tax[i,])] <- "Unclassified"
}
abundance.bacteria <- t(t(bacteria[,1:(ncol(bacteria)-6)])/colSums(bacteria[,1:(ncol(bacteria)-6)])*100)
tax.bacteria <- cbind(abundance.bacteria,bacteria.tax)

### Input virus
virus <- read.table(file = "Input/virus.tpm.xls",header = FALSE,row.names = 1,
                   sep = "\t",quote = "")

### virus processing
colnames(virus) <- virus[1,]
virus <- virus[-1,]
for (i in 1:(ncol(virus)-1)) {
    virus[,i] <- as.numeric(virus[,i])
}

virus$Phylum <- gsub(".*p__","",virus$taxonomy)
virus$Phylum <- gsub("\\;c__.*","",virus$Phylum)
virus$Class <- gsub(".*c__","",virus$taxonomy)
virus$Class <- gsub("\\;o__.*","",virus$Class)
virus$Order <- gsub(".*o__","",virus$taxonomy)
virus$Order <- gsub("\\;f__.*","",virus$Order)
virus$Family <- gsub(".*f__","",virus$taxonomy)
virus$Family <- gsub("\\;g__.*","",virus$Family)
virus$Genus <- gsub(".*g__","",virus$taxonomy)
virus$Genus <- gsub("\\;s__.*","",virus$Genus)
virus$Species <- gsub(".*s__","",virus$taxonomy)
virus <- subset(virus,select = -taxonomy)
virus.tax <- virus[,(ncol(virus)-5):ncol(virus)]
virus.tax[virus.tax == "norank"] <- "Unclassified"
virus.tax[virus.tax == "uncultured bacterium"] <- "Unclassified"
virus.tax[virus.tax == "uncultured"] <- "Unclassified"
for (i in 1:nrow(virus.tax)) {
    virus.tax[i,grepl("d__Bacteria",virus.tax[i,])] <- "Unclassified"
}
abundance.virus <- t(t(virus[,1:(ncol(virus)-6)])/colSums(virus[,1:(ncol(virus)-6)])*100)
aa <- match(colnames(abundance.virus),colnames(abundance.bacteria))
abundance.virus <- abundance.virus[,aa]
tax.virus <- cbind(abundance.virus,virus.tax)

### Input kegg
L1 <- read.table(file = "Input/kegg.profile.L1.xls",header = FALSE,row.names = 1,
                 sep = "\t",quote = "")
L2 <- read.table(file = "Input/kegg.profile.L2.xls",header = FALSE,row.names = 1,
                 sep = "\t",quote = "")
L3 <- read.table(file = "Input/kegg.profile.L3.xls",header = FALSE,row.names = 1,
                 sep = "\t",quote = "")
ko <- read.table(file = "Input/kegg.profile.entry.xls",header = FALSE,row.names = 1,
                 sep = "\t",quote = "")

### kegg processing
colnames(L1) <- L1[1,]
L1 <- L1[-1,]
for (i in 1:ncol(L1)) {
    L1[,i] <- as.numeric(L1[,i])
}
aa <- match(colnames(L1),colnames(abundance.bacteria))
L1 <- L1[,aa]
L1$Description <- rownames(L1)
L1[,1:(ncol(L1)-1)] <- t(t(L1[,1:(ncol(L1)-1)])/colSums(L1[,1:(ncol(L1)-1)])*100)
L1 <- L1[,c("Description",colnames(L1)[1:(ncol(L1)-1)])]

colnames(L2) <- L2[1,]
L2 <- L2[-1,]
for (i in 1:ncol(L2)) {
    L2[,i] <- as.numeric(L2[,i])
}
aa <- match(colnames(L2),colnames(abundance.bacteria))
L2 <- L2[,aa]
L2$Description <- rownames(L2)
L2[,1:(ncol(L2)-1)] <- t(t(L2[,1:(ncol(L2)-1)])/colSums(L2[,1:(ncol(L2)-1)])*100)
L2 <- L2[,c("Description",colnames(L2)[1:(ncol(L2)-1)])]

colnames(L3) <- L3[1,]
L3 <- L3[-1,]
for (i in 1:ncol(L3)) {
    L3[,i] <- as.numeric(L3[,i])
}
aa <- match(colnames(L3),colnames(abundance.bacteria))
L3 <- L3[,aa]
L3$Description <- rownames(L3)
L3[,1:(ncol(L3)-1)] <- t(t(L3[,1:(ncol(L3)-1)])/colSums(L3[,1:(ncol(L3)-1)])*100)
L3 <- L3[,c("Description",colnames(L3)[1:(ncol(L3)-1)])]

colnames(ko) <- ko[1,]
ko <- ko[-1,]
ko <- ko[,-ncol(ko)]
for (i in 1:ncol(ko)) {
    ko[,i] <- as.numeric(ko[,i])
}
aa <- match(colnames(ko),colnames(abundance.bacteria))
ko <- ko[,aa]
ko$Description <- rownames(ko)
ko[,1:(ncol(ko)-1)] <- t(t(ko[,1:(ncol(ko)-1)])/colSums(ko[,1:(ncol(ko)-1)])*100)
ko <- ko[,c("Description",colnames(ko)[1:(ncol(ko)-1)])]

### Mantel test and Prco analysis
dir.create("Results/Association_analysis_virome//Mantel")
source("Functions/Mantel.proc 2.R")
d <- mantel.proc.t(abundance.bacteria,abundance.virus,ko[,-1])
write.table(d,file = "Results/Association_analysis_virome//Mantel/Mantel_Proc.txt",sep = "\t",row.names = FALSE)
result <- mantel.proc(abundance.bacteria,abundance.virus,ko[,-1])
pdf(file = "Results/Association_analysis_virome//Mantel/Bac_Vir.pdf",width = 4,height = 3)
result[[1]]
dev.off()
pdf(file = "Results/Association_analysis_virome//Mantel/Bac_Fun.pdf",width = 4,height = 3)
result[[2]]
dev.off()
pdf(file = "Results/Association_analysis_virome//Mantel/Vir_Fun.pdf",width = 4,height = 3)
result[[3]]
dev.off()
pdf(file = "Results/Association_analysis_virome//Mantel/All.pdf",width = 12,height = 3)
result[[4]]
dev.off()

## Network
dir.create("Results/Association_analysis_virome/Network")
source("Functions/tax.summary.R")
bacteria<- tax.summary(tax.bacteria)
virus <- tax.summary(tax.virus)

### Phylum-L2
dir.create("Results/Association_analysis_virome/Network/Phylum-L2")
network_data <- cbind(t(bacteria[[1]][,2:ncol(bacteria[[1]])]),
                      t(virus[[1]][,2:ncol(virus[[1]])]),
                      t(L2[,2:ncol(L2)]))
colnames(network_data) <- c(bacteria[[1]]$Phylum,virus[[1]]$Phylum,L2$Description)
network_data <- network_data[,colMeans(network_data) >= 0.1]
f <- function(x) sum(x==0)
d <- apply(network_data,2,f)
network_data1 <- network_data[,d < nrow(network_data)*0.4]

nodes <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    as_tibble(what = "vertices")
nodes$Type <- ifelse(nodes$name %in% bacteria[[1]]$Phylum,"Bacteria",
                     ifelse(nodes$name %in% virus[[1]]$Phylum,"Virus","Function"))
network_data1 <- as.data.frame(network_data1)
Total.abun <- data.frame(name = colnames(network_data1),
                         Abun = colMeans(network_data1))
nodes1 <- merge(nodes,Total.abun)
rownames(nodes1) <- nodes1$name
nodes1 <- nodes1[nodes$name,]

net <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    activate("nodes") %>%
    mutate(group = nodes1$Type,
           Abundance = nodes1$Abun)

p <- ggraph(net, "circular") +
    geom_edge_fan(color = "grey70", show.legend = FALSE, alpha = 0.2) +
    geom_node_point(aes(colour = group, size = Abundance, shape = group),
                    show.legend = TRUE) +
    geom_node_text(aes(label = name),size = 1.5) +
    scale_edge_color_gradient2() +
    scale_size_continuous(name = "Abundance") +
    coord_fixed() +
    ggplot2::theme(panel.background = ggplot2::element_blank())

pdf("Results/Association_analysis_virome/Network/Phylum-L2/network_abun.pdf",
    width = 5,height = 5)
p
dev.off()

p2 <- ggraph(net, "circular") +
    geom_edge_fan(color = "grey70", show.legend = FALSE, alpha = 0.2) +
    geom_node_point(aes(colour = group, size = Abundance, shape = group),
                    show.legend = TRUE) +
    scale_edge_color_gradient2() +
    scale_size_continuous(name = "Abundance") +
    coord_fixed() +
    ggplot2::theme(panel.background = ggplot2::element_blank())

pdf("Results/Association_analysis_virome/Network/Phylum-L2/network_abun_nolabel.pdf",
    width = 5,height = 5)
p2
dev.off()

net2 <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    activate("nodes") %>%
    mutate(group = nodes1$Type,
           Degree = centrality_degree())

p <- ggraph(net2, "circular") +
    geom_edge_fan(color = "grey70", show.legend = FALSE, alpha = 0.2) +
    geom_node_point(aes(colour = group, size = Degree, shape = group),
                    show.legend = TRUE) +
    geom_node_text(aes(label = name),size = 1.5) +
    scale_edge_color_gradient2() +
    scale_size_continuous(name = "Degree") +
    coord_fixed() +
    ggplot2::theme(panel.background = ggplot2::element_blank())

pdf("Results/Association_analysis_virome/Network/Phylum-L2/network_degree.pdf",
    width = 5,height = 5)
p
dev.off()

p2 <- ggraph(net2, "circular") +
    geom_edge_fan(color = "grey70", show.legend = FALSE, alpha = 0.2) +
    geom_node_point(aes(colour = group, size = Degree, shape = group),
                    show.legend = TRUE) +
    scale_edge_color_gradient2() +
    scale_size_continuous(name = "Degree") +
    coord_fixed() +
    ggplot2::theme(panel.background = ggplot2::element_blank())

pdf("Results/Association_analysis_virome/Network/Phylum-L2/network_degree_nolabel.pdf",
    width = 5,height = 5)
p2
dev.off()

res <- linkET::fast_correlate2(network_data1,method = "spearman",adjust = TRUE,
                             adjust_method = "BH")

res$r[lower.tri(res$r)] = 0
res$p[lower.tri(res$p)] = 0
data.r <- melt(res$r)
data.p <- melt(res$p)
result <- cbind(data.r,data.p)
result <- result[,c(1,2,3,6)]
colnames(result) <- c("Var1","Var2","R","P")

result.1 <- result[result$R > 0.8,]
result.2 <- result[result$R < -0.8,]
result <- rbind(result.1,result.2)
result <- result[result$P < 0.05,]
result<- result[result$R < 1,]
result <- result[result$Var1 != result$Var2,]
result$R[result$R > 0.8] = 1
result$R[result$R < -0.8] = -1
colnames(result) <- c("Source","Target","R","P")
write.table(result,"Results/Association_analysis_virome/Network/Phylum-L2/edge.csv",
            sep=",",quote=F,row.names = FALSE)

### Genus-Pathway
dir.create("Results/Association_analysis_virome/Network/Genus-Pathway")
network_data <- cbind(t(bacteria[[5]][,2:ncol(bacteria[[5]])]),
                      t(virus[[5]][,2:ncol(virus[[5]])]),
                      t(L3[,2:ncol(L3)]))
colnames(network_data) <- c(bacteria[[5]]$Genus,virus[[5]]$Genus,L3$Description)
network_data <- network_data[,colMeans(network_data) >= 0.1]
f <- function(x) sum(x==0)
d <- apply(network_data,2,f)
network_data1 <- network_data[,d < nrow(network_data)*0.4]
network_data1 <- network_data[,unique(colnames(network_data))]

nodes <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    as_tibble(what = "vertices")
nodes$Type <- ifelse(nodes$name %in% bacteria[[5]]$Genus,"Bacteria",
                     ifelse(nodes$name %in% virus[[5]]$Genus,"Virus","Function"))
network_data1 <- as.data.frame(network_data1)
Total.abun <- data.frame(name = colnames(network_data1),
                         Abun = colMeans(network_data1))
nodes1 <- merge(nodes,Total.abun)
rownames(nodes1) <- nodes1$name
nodes1 <- nodes1[nodes$name,]

net <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    activate("nodes") %>%
    mutate(group = nodes1$Type,
           Abundance = nodes1$Abun)

p <- ggraph(net, "circular") +
    geom_edge_fan(color = "grey70", show.legend = FALSE, alpha = 0.2) +
    geom_node_point(aes(colour = group, size = Abundance, shape = group),
                    show.legend = TRUE) +
    geom_node_text(aes(label = name),size = 1.5) +
    scale_edge_color_gradient2() +
    scale_size_continuous(name = "Abundance") +
    coord_fixed() +
    ggplot2::theme(panel.background = ggplot2::element_blank())

pdf("Results/Association_analysis_virome/Network/Genus-Pathway/network_abun.pdf",
    width = 5,height = 5)
p
dev.off()

p2 <- ggraph(net, "circular") +
    geom_edge_fan(color = "grey70", show.legend = FALSE, alpha = 0.2) +
    geom_node_point(aes(colour = group, size = Abundance, shape = group),
                    show.legend = TRUE) +
    scale_edge_color_gradient2() +
    scale_size_continuous(name = "Abundance") +
    coord_fixed() +
    ggplot2::theme(panel.background = ggplot2::element_blank())

pdf("Results/Association_analysis_virome/Network/Genus-Pathway/network_abun_nolabel.pdf",
    width = 5,height = 5)
p2
dev.off()

net2 <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    activate("nodes") %>%
    mutate(group = nodes1$Type,
           Degree = centrality_degree())

p <- ggraph(net2, "circular") +
    geom_edge_fan(color = "grey70", show.legend = FALSE, alpha = 0.2) +
    geom_node_point(aes(colour = group, size = Degree, shape = group),
                    show.legend = TRUE) +
    geom_node_text(aes(label = name),size = 1.5) +
    scale_edge_color_gradient2() +
    scale_size_continuous(name = "Degree") +
    coord_fixed() +
    ggplot2::theme(panel.background = ggplot2::element_blank())

pdf("Results/Association_analysis_virome/Network/Genus-Pathway/network_degree.pdf",
    width = 5,height = 5)
p
dev.off()

p2 <- ggraph(net2, "circular") +
    geom_edge_fan(color = "grey70", show.legend = FALSE, alpha = 0.2) +
    geom_node_point(aes(colour = group, size = Degree, shape = group),
                    show.legend = TRUE) +
    scale_edge_color_gradient2() +
    scale_size_continuous(name = "Degree") +
    coord_fixed() +
    ggplot2::theme(panel.background = ggplot2::element_blank())

pdf("Results/Association_analysis_virome/Network/Genus-Pathway/network_degree_nolabel.pdf",
    width = 5,height = 5)
p2
dev.off()

res <- linkET::fast_correlate2(network_data1,method = "spearman",adjust = TRUE)

res$r[lower.tri(res$r)] = 0
res$p[lower.tri(res$p)] = 0
data.r <- melt(res$r)
data.p <- melt(res$p)
result <- cbind(data.r,data.p)
result <- result[,c(1,2,3,6)]
colnames(result) <- c("Var1","Var2","R","P")

result.1 <- result[result$R > 0.8,]
result.2 <- result[result$R < -0.8,]
result <- rbind(result.1,result.2)
result <- result[result$P < 0.05,]
result<- result[result$R < 1,]
result <- result[result$Var1 != result$Var2,]
result$R[result$R > 0.8] = 1
result$R[result$R < -0.8] = -1
colnames(result) <- c("Source","Target","R","P")
write.table(result,"Results/Association_analysis_virome/Network/Genus-Pathway/edge.csv",
            sep=",",quote=F,row.names = FALSE)

### Genus-ko
dir.create("Results/Association_analysis_virome/Network/Genus-ko")
network_data <- cbind(t(bacteria[[5]][,2:ncol(bacteria[[5]])]),
                      t(virus[[5]][,2:ncol(virus[[5]])]),
                      t(ko[,2:ncol(ko)]))
colnames(network_data) <- c(bacteria[[5]]$Genus,virus[[5]]$Genus,rownames(ko))
network_data <- network_data[,colMeans(network_data) >= 0.1]
f <- function(x) sum(x==0)
d <- apply(network_data,2,f)
network_data1 <- network_data[,d < nrow(network_data)*0.4]
network_data1 <- network_data[,unique(colnames(network_data))]

nodes <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    as_tibble(what = "vertices")
nodes$Type <- ifelse(nodes$name %in% bacteria[[5]]$Genus,"Bacteria",
                     ifelse(nodes$name %in% virus[[5]]$Genus,"Virus","Function"))
network_data1 <- as.data.frame(network_data1)
Total.abun <- data.frame(name = colnames(network_data1),
                         Abun = colMeans(network_data1))
nodes1 <- merge(nodes,Total.abun)
rownames(nodes1) <- nodes1$name
nodes1 <- nodes1[nodes$name,]

net <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    activate("nodes") %>%
    mutate(group = nodes1$Type,
           Abundance = nodes1$Abun)

p <- ggraph(net, "circular") +
    geom_edge_fan(color = "grey70", show.legend = FALSE, alpha = 0.2) +
    geom_node_point(aes(colour = group, size = Abundance, shape = group),
                    show.legend = TRUE) +
    geom_node_text(aes(label = name),size = 1.5) +
    scale_edge_color_gradient2() +
    scale_size_continuous(name = "Abundance") +
    coord_fixed() +
    ggplot2::theme(panel.background = ggplot2::element_blank())

pdf("Results/Association_analysis_virome/Network/Genus-ko/network_abun.pdf",
    width = 5,height = 5)
p
dev.off()

p2 <- ggraph(net, "circular") +
    geom_edge_fan(color = "grey70", show.legend = FALSE, alpha = 0.2) +
    geom_node_point(aes(colour = group, size = Abundance, shape = group),
                    show.legend = TRUE) +
    scale_edge_color_gradient2() +
    scale_size_continuous(name = "Abundance") +
    coord_fixed() +
    ggplot2::theme(panel.background = ggplot2::element_blank())

pdf("Results/Association_analysis_virome/Network/Genus-ko/network_abun_nolabel.pdf",
    width = 5,height = 5)
p2
dev.off()

net2 <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    activate("nodes") %>%
    mutate(group = nodes1$Type,
           Degree = centrality_degree())

p <- ggraph(net2, "circular") +
    geom_edge_fan(color = "grey70", show.legend = FALSE, alpha = 0.2) +
    geom_node_point(aes(colour = group, size = Degree, shape = group),
                    show.legend = TRUE) +
    geom_node_text(aes(label = name),size = 1.5) +
    scale_edge_color_gradient2() +
    scale_size_continuous(name = "Degree") +
    coord_fixed() +
    ggplot2::theme(panel.background = ggplot2::element_blank())

pdf("Results/Association_analysis_virome/Network/Genus-ko/network_degree.pdf",
    width = 5,height = 5)
p
dev.off()

p2 <- ggraph(net2, "circular") +
    geom_edge_fan(color = "grey70", show.legend = FALSE, alpha = 0.2) +
    geom_node_point(aes(colour = group, size = Degree, shape = group),
                    show.legend = TRUE) +
    scale_edge_color_gradient2() +
    scale_size_continuous(name = "Degree") +
    coord_fixed() +
    ggplot2::theme(panel.background = ggplot2::element_blank())

pdf("Results/Association_analysis_virome/Network/Genus-ko/network_degree_nolabel.pdf",
    width = 5,height = 5)
p2
dev.off()

res <- linkET::fast_correlate2(network_data1,method = "spearman",adjust = TRUE)

res$r[lower.tri(res$r)] = 0
res$p[lower.tri(res$p)] = 0
data.r <- melt(res$r)
data.p <- melt(res$p)
result <- cbind(data.r,data.p)
result <- result[,c(1,2,3,6)]
colnames(result) <- c("Var1","Var2","R","P")

result.1 <- result[result$R > 0.8,]
result.2 <- result[result$R < -0.8,]
result <- rbind(result.1,result.2)
result <- result[result$P < 0.05,]
result<- result[result$R < 1,]
result <- result[result$Var1 != result$Var2,]
result$R[result$R > 0.8] = 1
result$R[result$R < -0.8] = -1
colnames(result) <- c("Source","Target","R","P")
write.table(result,"Results/Association_analysis_virome/Network/Genus-ko/edge.csv",
            sep=",",quote=F,row.names = FALSE)

## Heatmap
dir.create("Results/Association_analysis_virome/Heatmap")
source("Functions/Correlation.heatmap.R")

### Bacteria-Virus
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Virus")
#### Wilcox
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox")

#### Phylum
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Phylum/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
marker2 <- read.table("Results/Virus/Diff_analysis/Wilcox/Phylum/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[1]],virus[[1]])
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Class/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
marker2 <- read.table("Results/Virus/Diff_analysis/Wilcox/Class/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[2]],virus[[2]])
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Order/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
marker2 <- read.table("Results/Virus/Diff_analysis/Wilcox/Order/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[3]],virus[[3]])
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Family/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
marker2 <- read.table("Results/Virus/Diff_analysis/Wilcox/Family/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[4]],virus[[4]])
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
marker2 <- read.table("Results/Virus/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[5]],virus[[5]])
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Wilcox/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### RandomForest
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest")

#### Phylum
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Phylum/importance.txt",
                      header = TRUE,sep = "\t")
marker2 <- read.table("Results/Virus/Diff_analysis/RandomForest/Phylum/importance.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[1]],virus[[1]])
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Class/importance.txt",
                      header = TRUE,sep = "\t")
marker2 <- read.table("Results/Virus/Diff_analysis/RandomForest/Class/importance.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[2]],virus[[2]])
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Order/importance.txt",
                      header = TRUE,sep = "\t")
marker2 <- read.table("Results/Virus/Diff_analysis/RandomForest/Order/importance.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[3]],virus[[3]])
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Family/importance.txt",
                      header = TRUE,sep = "\t")
marker2 <- read.table("Results/Virus/Diff_analysis/RandomForest/Family/importance.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[4]],virus[[4]])
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Genus/importance.txt",
                      header = TRUE,sep = "\t")
marker2 <- read.table("Results/Virus/Diff_analysis/RandomForest/Genus/importance.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[5]],virus[[5]])
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/RandomForest/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Lefse
marker1 <- read.table("Results/Bacteria/Diff_analysis/LEfSe/Lefse_result.txt",
                      header = TRUE,sep = "\t")
marker2 <- read.table("Results/Virus/Diff_analysis/LEfSe/Lefse_result.txt",
                      header = TRUE,sep = "\t")
marker1$V1 <- rownames(marker1)
marker1 <- marker1[grep("g__",marker1$V1),]
marker1$V1 <- gsub("g__","",marker1$V1)
marker1 <- marker1[,c("V1",colnames(marker1)[1:(ncol(marker1)-1)])]
marker2$V1 <- rownames(marker2)
marker2 <- marker2[grep("g__",marker2$V1),]
marker2$V1 <- gsub("g__","",marker2$V1)
marker2 <- marker2[,c("V1",colnames(marker2)[1:(ncol(marker2)-1)])]

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Lefse")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[5]],virus[[5]])
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Lefse/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Lefse/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Lefse/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Share
source("Functions/overLapper.new.r")
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Shared")
marker1.w <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                        header = TRUE,sep = "\t")
marker1.r <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Genus/importance.txt",
                        header = TRUE,sep = "\t")

species.v <- list(wilcox = marker1.w$V1,RandomForest = marker1.r$ID,
                  Lefse = marker1$V1)

pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Shared/venn_bacteria.pdf",
    width=9,height=9,pointsize=16)
OLlist <- overLapper(setlist=species.v, sep="", type="vennsets",keepdups=FALSE)
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts,mymain="")
dev.off()

inter.v <- get.venn.partitions(species.v)
shared.species.v <- marker1[marker1$V1 %in% unlist(inter.v$..values..[1]),]

marker2.w <- read.table("Results/Virus/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                        header = TRUE,sep = "\t")
marker2.r <- read.table("Results/Virus/Diff_analysis/RandomForest/Genus/importance.txt",
                        header = TRUE,sep = "\t")

species.f <- list(wilcox = marker2.w$V1,RandomForest = marker2.r$ID,
                  Lefse = marker2$V1)

pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Shared/venn_virus.pdf",
    width=9,height=9,pointsize=16)
OLlist <- overLapper(setlist=species.f, sep="", type="vennsets",keepdups=FALSE)
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts,mymain="")
dev.off()

inter.f <- get.venn.partitions(species.f)
shared.species.f <- marker2[unlist(inter.f$..values..[1]),]

if (nrow(shared.species.v) > 0 & nrow(shared.species.f) > 0) {
    result <- correlation.heatmap(shared.species.v,shared.species.f,
                                  virus[[5]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Shared/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Shared/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Virus/Shared/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

### Bacteria-Function
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function")

#### L1
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1")

#### Wilcox
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox")
marker2 <- read.table("Results/Gene/Diff_analysis/Wilcox/L1/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
#### Phylum
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Phylum/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[1]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Class/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[2]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Order/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[3]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Family/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[4]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[5]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/Wilcox/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### RandomForest
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest")
marker2 <- read.table("Results/Gene/Diff_analysis/RandomForest/L1/importance.txt",
                      header = TRUE,sep = "\t")

#### Phylum
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Phylum/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[1]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Class/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[2]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Order/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[3]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Family/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[4]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Genus/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[5]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L1/RandomForest/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### L2
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2")

#### Wilcox
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox")
marker2 <- read.table("Results/Gene/Diff_analysis/Wilcox/L2/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
#### Phylum
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Phylum/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[1]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Class/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[2]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Order/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[3]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Family/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[4]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[5]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/Wilcox/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### RandomForest
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest")
marker2 <- read.table("Results/Gene/Diff_analysis/RandomForest/L2/importance.txt",
                      header = TRUE,sep = "\t")

#### Phylum
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Phylum/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[1]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Class/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[2]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Order/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[3]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Family/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[4]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Genus/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[5]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L2/RandomForest/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### L3
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3")

#### Wilcox
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox")
marker2 <- read.table("Results/Gene/Diff_analysis/Wilcox/L3/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
#### Phylum
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Phylum/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[1]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Class/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[2]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Order/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[3]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Family/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[4]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[5]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/Wilcox/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### RandomForest
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest")
marker2 <- read.table("Results/Gene/Diff_analysis/RandomForest/L3/importance.txt",
                      header = TRUE,sep = "\t")

#### Phylum
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Phylum/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[1]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Class/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[2]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Order/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[3]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Family/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[4]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Genus/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[5]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/L3/RandomForest/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### ko
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko")

#### Wilcox
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox")
marker2 <- read.table("Results/Gene/Diff_analysis/Wilcox/ko/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
#### Phylum
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Phylum/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[1]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Class/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[2]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Order/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[3]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Family/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[4]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[5]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/Wilcox/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### RandomForest
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest")
marker2 <- read.table("Results/Gene/Diff_analysis/RandomForest/ko/importance.txt",
                      header = TRUE,sep = "\t")

#### Phylum
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Phylum/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[1]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Class/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[2]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Order/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[3]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Family/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[4]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Genus/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[5]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/ko/RandomForest/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Lefse
marker1 <- read.table("Results/Bacteria/Diff_analysis/LEfSe/Lefse_result.txt",
                      header = TRUE,sep = "\t")
marker2 <- read.table("Results/Gene/Diff_analysis/LEfSe/Lefse_result.txt",
                      header = TRUE,sep = "\t")
marker1$V1 <- rownames(marker1)
marker1 <- marker1[grep("g__",marker1$V1),]
marker1$V1 <- gsub("g__","",marker1$V1)
marker1 <- marker1[,c("V1",colnames(marker1)[1:(ncol(marker1)-1)])]
marker2$V1 <- rownames(marker2)
marker2 <- marker2[,c("V1",colnames(marker2)[1:(ncol(marker2)-1)])]

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/Lefse")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[5]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/Lefse/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/Lefse/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/Lefse/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Share
source("Functions/overLapper.new.r")
dir.create("Results/Association_analysis_virome/Heatmap/Bacteria_Function/Shared")
marker1.w <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                        header = TRUE,sep = "\t")
marker1.r <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Genus/importance.txt",
                        header = TRUE,sep = "\t")

species.v <- list(wilcox = marker1.w$V1,RandomForest = marker1.r$ID,
                  Lefse = marker1$V1)

pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/Shared/venn_bacteria.pdf",
    width=9,height=9,pointsize=16)
OLlist <- overLapper(setlist=species.v, sep="", type="vennsets",keepdups=FALSE)
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts,mymain="")
dev.off()

inter.v <- get.venn.partitions(species.v)
shared.species.v <- marker1[marker1$V1 %in% unlist(inter.v$..values..[1]),]

marker2.w <- read.table("Results/Gene/Diff_analysis/Wilcox/ko/wilcox_biomarker.txt",
                        header = TRUE,sep = "\t")
marker2.r <- read.table("Results/Gene/Diff_analysis/RandomForest/ko/importance.txt",
                        header = TRUE,sep = "\t")

species.f <- list(wilcox = marker2.w$V1,RandomForest = marker2.r$ID,
                  Lefse = marker2$V1)

pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/Shared/venn_function.pdf",
    width=9,height=9,pointsize=16)
OLlist <- overLapper(setlist=species.f, sep="", type="vennsets",keepdups=FALSE)
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts,mymain="")
dev.off()

inter.f <- get.venn.partitions(species.f)
shared.species.f <- marker2[unlist(inter.f$..values..[1]),]

if (nrow(shared.species.v) > 0 & nrow(shared.species.f) > 0) {
    result <- correlation.heatmap(shared.species.v,shared.species.f,
                                  bacteria[[5]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/Shared/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Bacteria_Function/Shared/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Bacteria_Function/Shared/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

### Virus-Function
dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function")

#### L1
dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L1")

#### Wilcox
dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox")
marker2 <- read.table("Results/Gene/Diff_analysis/Wilcox/L1/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
#### Phylum
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Phylum/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[1]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Class/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[2]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Order/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[3]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Family/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[4]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[5]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/Wilcox/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### RandomForest
dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest")
marker2 <- read.table("Results/Gene/Diff_analysis/RandomForest/L1/importance.txt",
                      header = TRUE,sep = "\t")

#### Phylum
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Phylum/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[1]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Class/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[2]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Order/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[3]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Family/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[4]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Genus/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[5]],L1)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L1/RandomForest/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### L2
dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L2")

#### Wilcox
dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox")
marker2 <- read.table("Results/Gene/Diff_analysis/Wilcox/L2/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
#### Phylum
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Phylum/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[1]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Class/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[2]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Order/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[3]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Family/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[4]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[5]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/Wilcox/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### RandomForest
dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest")
marker2 <- read.table("Results/Gene/Diff_analysis/RandomForest/L2/importance.txt",
                      header = TRUE,sep = "\t")

#### Phylum
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Phylum/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[1]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Class/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[2]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Order/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[3]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Family/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[4]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Genus/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[5]],L2)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L2/RandomForest/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### L3
dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L3")

#### Wilcox
dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox")
marker2 <- read.table("Results/Gene/Diff_analysis/Wilcox/L3/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
#### Phylum
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Phylum/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[1]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Class/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[2]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Order/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[3]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Family/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[4]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[5]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/Wilcox/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### RandomForest
dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest")
marker2 <- read.table("Results/Gene/Diff_analysis/RandomForest/L3/importance.txt",
                      header = TRUE,sep = "\t")

#### Phylum
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Phylum/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[1]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Class/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[2]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Order/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[3]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Family/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[4]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Genus/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[5]],L3)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/L3/RandomForest/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### ko
dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/ko")

#### Wilcox
dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox")
marker2 <- read.table("Results/Gene/Diff_analysis/Wilcox/ko/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
#### Phylum
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Phylum/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[1]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Class/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[2]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Order/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[3]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Family/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[4]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[5]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/Wilcox/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### RandomForest
dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest")
marker2 <- read.table("Results/Gene/Diff_analysis/RandomForest/ko/importance.txt",
                      header = TRUE,sep = "\t")

#### Phylum
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Phylum/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[1]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Class/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[2]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Order/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[3]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Family/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[4]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Genus/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[5]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/ko/RandomForest/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Lefse
marker1 <- read.table("Results/Virus/Diff_analysis/LEfSe/Lefse_result.txt",
                      header = TRUE,sep = "\t")
marker2 <- read.table("Results/Gene/Diff_analysis/LEfSe/Lefse_result.txt",
                      header = TRUE,sep = "\t")
marker1$V1 <- rownames(marker1)
marker1 <- marker1[grep("g__",marker1$V1),]
marker1$V1 <- gsub("g__","",marker1$V1)
marker1 <- marker1[,c("V1",colnames(marker1)[1:(ncol(marker1)-1)])]
marker2$V1 <- rownames(marker2)
marker2 <- marker2[,c("V1",colnames(marker2)[1:(ncol(marker2)-1)])]

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/Lefse")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[5]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/Lefse/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/Lefse/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/Lefse/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Share
source("Functions/overLapper.new.r")
dir.create("Results/Association_analysis_virome/Heatmap/Virus_Function/Shared")
marker1.w <- read.table("Results/Virus/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
marker1.r <- read.table("Results/Virus/Diff_analysis/RandomForest/Genus/importance.txt",
                      header = TRUE,sep = "\t")

species.v <- list(wilcox = marker1.w$V1,RandomForest = marker1.r$ID,
                Lefse = marker1$V1)

pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/Shared/venn_virus.pdf",
    width=9,height=9,pointsize=16)
OLlist <- overLapper(setlist=species.v, sep="", type="vennsets",keepdups=FALSE)
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts,mymain="")
dev.off()

inter.v <- get.venn.partitions(species.v)
shared.species.v <- marker1[marker1$V1 %in% unlist(inter.v$..values..[1]),]

marker2.w <- read.table("Results/Gene/Diff_analysis/Wilcox/ko/wilcox_biomarker.txt",
                        header = TRUE,sep = "\t")
marker2.r <- read.table("Results/Gene/Diff_analysis/RandomForest/ko/importance.txt",
                        header = TRUE,sep = "\t")

species.f <- list(wilcox = marker2.w$V1,RandomForest = marker2.r$ID,
                  Lefse = marker2$V1)

pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/Shared/venn_function.pdf",
    width=9,height=9,pointsize=16)
OLlist <- overLapper(setlist=species.f, sep="", type="vennsets",keepdups=FALSE)
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts,mymain="")
dev.off()

inter.f <- get.venn.partitions(species.f)
shared.species.f <- marker2[unlist(inter.f$..values..[1]),]

if (nrow(shared.species.v) > 0 & nrow(shared.species.f) > 0) {
    result <- correlation.heatmap(shared.species.v,shared.species.f,
                                  virus[[5]],ko)
    write.table(result[[1]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/Shared/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_virome/Heatmap/Virus_Function/Shared/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_virome/Heatmap/Virus_Function/Shared/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}











