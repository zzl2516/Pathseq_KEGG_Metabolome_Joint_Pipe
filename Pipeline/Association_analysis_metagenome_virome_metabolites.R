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
dir.create("Results/Association_analysis_metabolites")

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

### Input metabolites
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
aa <- match(colnames(metabolites),colnames(abundance.bacteria))
metabolites <- metabolites[,aa]
metabolites <- metabolites[rowSums(metabolites) > 0,]
metabolites <- scale(metabolites)
metabolites <- as.data.frame(metabolites)
metabolites$compound <- rownames(metabolites)
metabolites <- metabolites[,c("compound",colnames(metabolites)[1:(ncol(metabolites)-1)])]

### Mantel test and Prco analysis
dir.create("Results/Association_analysis_metabolites/Mantel")
source("Functions/Mantel.proc 3.R")
d <- mantel.proc.t(abundance.bacteria,abundance.virus,ko[,-1],metabolites[,-1])
write.table(d,file = "Results/Association_analysis_metabolites/Mantel/Mantel_Proc.txt",sep = "\t",row.names = FALSE)
result <- mantel.proc(abundance.bacteria,abundance.virus,ko[,-1],metabolites[,-1])
pdf(file = "Results/Association_analysis_metabolites/Mantel/Bac_Metabolites.pdf",width = 4,height = 3)
result[[1]]
dev.off()
pdf(file = "Results/Association_analysis_metabolites/Mantel/Vir_Metabolites.pdf",width = 4,height = 3)
result[[2]]
dev.off()
pdf(file = "Results/Association_analysis_metabolites/Mantel/Fun_Metabolites.pdf",width = 4,height = 3)
result[[3]]
dev.off()
pdf(file = "Results/Association_analysis_metabolites/Mantel/All.pdf",width = 12,height = 3)
result[[4]]
dev.off()

### Network
dir.create("Results/Association_analysis_metabolites/Network")
source("Functions/tax.summary.R")
bacteria<- tax.summary(tax.bacteria)
virus <- tax.summary(tax.virus)

### Phylum-L2-DAMs
dir.create("Results/Association_analysis_metabolites/Network/Phylum-L2-DAMs")

#### Wilcox
dir.create("Results/Association_analysis_metabolites/Network/Phylum-L2-DAMs/Wilcox")
marker <- read.table("Results/Metabolites/Diff_analysis/wilcox/wilcox_biomarker.txt",
                     header = TRUE)
metabolites1 <- metabolites[metabolites$compound %in% marker$V1,]
network_data <- cbind(t(bacteria[[1]][,2:ncol(bacteria[[1]])]),
                      t(virus[[1]][,2:ncol(virus[[1]])]),
                      t(L2[,2:ncol(L2)]))
colnames(network_data) <- c(bacteria[[1]]$Phylum,virus[[1]]$Phylum,L2$Description)
network_data <- network_data[,colMeans(network_data) >= 0.1]
f <- function(x) sum(x==0)
d <- apply(network_data,2,f)
network_data1 <- network_data[,d < nrow(network_data)*0.4]
network_data1 <- cbind(network_data1,t(metabolites1[,2:ncol(metabolites1)]))

nodes <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    as_tibble(what = "vertices")
nodes$Type <- ifelse(nodes$name %in% bacteria[[1]]$Phylum,"Bacteria",
                     ifelse(nodes$name %in% virus[[1]]$Phylum,"Virus",
                            ifelse(nodes$name %in% L2$Description,"Function",
                                   "Metabolites")))
network_data1 <- as.data.frame(network_data1)
Total.abun <- data.frame(name = colnames(network_data1),
                         Abun = colMeans(network_data1))
nodes1 <- merge(nodes,Total.abun)
rownames(nodes1) <- nodes1$name
nodes1 <- nodes1[nodes$name,]

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

pdf("Results/Association_analysis_metabolites/Network/Phylum-L2-DAMs/Wilcox/network_degree.pdf",
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

pdf("Results/Association_analysis_metabolites/Network/Phylum-L2-DAMs/Wilcox/network_degree_nolabel.pdf",
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
write.table(result,"Results/Association_analysis_metabolites/Network/Phylum-L2-DAMs/Wilcox/edge.csv",
            sep=",",quote=F,row.names = FALSE)

#### RandomForest
dir.create("Results/Association_analysis_metabolites/Network/Phylum-L2-DAMs/RandomForest")
marker <- read.table("Results/Metabolites/Diff_analysis/RandomForest/importance.txt",
                     header = TRUE)
metabolites1 <- metabolites[metabolites$compound %in% marker$ID,]
network_data <- cbind(t(bacteria[[1]][,2:ncol(bacteria[[1]])]),
                      t(virus[[1]][,2:ncol(virus[[1]])]),
                      t(L2[,2:ncol(L2)]))
colnames(network_data) <- c(bacteria[[1]]$Phylum,virus[[1]]$Phylum,L2$Description)
network_data <- network_data[,colMeans(network_data) >= 0.1]
f <- function(x) sum(x==0)
d <- apply(network_data,2,f)
network_data1 <- network_data[,d < nrow(network_data)*0.4]
network_data1 <- cbind(network_data1,t(metabolites1[,2:ncol(metabolites1)]))

nodes <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    as_tibble(what = "vertices")
nodes$Type <- ifelse(nodes$name %in% bacteria[[1]]$Phylum,"Bacteria",
                     ifelse(nodes$name %in% virus[[1]]$Phylum,"Virus",
                            ifelse(nodes$name %in% L2$Description,"Function",
                                   "Metabolites")))
network_data1 <- as.data.frame(network_data1)
Total.abun <- data.frame(name = colnames(network_data1),
                         Abun = colMeans(network_data1))
nodes1 <- merge(nodes,Total.abun)
rownames(nodes1) <- nodes1$name
nodes1 <- nodes1[nodes$name,]

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

pdf("Results/Association_analysis_metabolites/Network/Phylum-L2-DAMs/RandomForest/network_degree.pdf",
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

pdf("Results/Association_analysis_metabolites/Network/Phylum-L2-DAMs/RandomForest/network_degree_nolabel.pdf",
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
write.table(result,"Results/Association_analysis_metabolites/Network/Phylum-L2-DAMs/RandomForest/edge.csv",
            sep=",",quote=F,row.names = FALSE)

#### Lefse
dir.create("Results/Association_analysis_metabolites/Network/Phylum-L2-DAMs/Lefse")
marker <- read.table("Results/Metabolites/Diff_analysis/LEfSe/Lefse_result.txt",
                     header = TRUE)
metabolites1 <- metabolites[metabolites$compound %in% rownames(marker),]
network_data <- cbind(t(bacteria[[1]][,2:ncol(bacteria[[1]])]),
                      t(virus[[1]][,2:ncol(virus[[1]])]),
                      t(L2[,2:ncol(L2)]))
colnames(network_data) <- c(bacteria[[1]]$Phylum,virus[[1]]$Phylum,L2$Description)
network_data <- network_data[,colMeans(network_data) >= 0.1]
f <- function(x) sum(x==0)
d <- apply(network_data,2,f)
network_data1 <- network_data[,d < nrow(network_data)*0.4]
network_data1 <- cbind(network_data1,t(metabolites1[,2:ncol(metabolites1)]))

nodes <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    as_tibble(what = "vertices")
nodes$Type <- ifelse(nodes$name %in% bacteria[[1]]$Phylum,"Bacteria",
                     ifelse(nodes$name %in% virus[[1]]$Phylum,"Virus",
                            ifelse(nodes$name %in% L2$Description,"Function",
                                   "Metabolites")))
network_data1 <- as.data.frame(network_data1)
Total.abun <- data.frame(name = colnames(network_data1),
                         Abun = colMeans(network_data1))
nodes1 <- merge(nodes,Total.abun)
rownames(nodes1) <- nodes1$name
nodes1 <- nodes1[nodes$name,]

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

pdf("Results/Association_analysis_metabolites/Network/Phylum-L2-DAMs/Lefse/network_degree.pdf",
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

pdf("Results/Association_analysis_metabolites/Network/Phylum-L2-DAMs/Lefse/network_degree_nolabel.pdf",
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
write.table(result,"Results/Association_analysis_metabolites/Network/Phylum-L2-DAMs/Lefse/edge.csv",
            sep=",",quote=F,row.names = FALSE)

### Genus-Pathway-DAMs
dir.create("Results/Association_analysis_metabolites/Network/Genus-Pathway-DAMs")

#### Wilcox
dir.create("Results/Association_analysis_metabolites/Network/Genus-Pathway-DAMs/Wilcox")
marker <- read.table("Results/Metabolites/Diff_analysis/wilcox/wilcox_biomarker.txt",
                     header = TRUE)
metabolites1 <- metabolites[metabolites$compound %in% marker$V1,]
network_data <- cbind(t(bacteria[[5]][,2:ncol(bacteria[[5]])]),
                      t(virus[[5]][,2:ncol(virus[[5]])]),
                      t(L3[,2:ncol(L3)]))
colnames(network_data) <- c(bacteria[[5]]$Genus,virus[[5]]$Genus,L3$Description)
network_data <- network_data[,colMeans(network_data) >= 0.1]
f <- function(x) sum(x==0)
d <- apply(network_data,2,f)
network_data1 <- network_data[,d < nrow(network_data)*0.4]
network_data1 <- cbind(network_data1,t(metabolites1[,2:ncol(metabolites1)]))
network_data1 <- network_data1[,unique(colnames(network_data1))]

nodes <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    as_tibble(what = "vertices")
nodes$Type <- ifelse(nodes$name %in% bacteria[[5]]$Genus,"Bacteria",
                     ifelse(nodes$name %in% virus[[5]]$Genus,"Virus",
                            ifelse(nodes$name %in% L3$Description,"Function",
                                   "Metabolites")))
network_data1 <- as.data.frame(network_data1)
Total.abun <- data.frame(name = colnames(network_data1),
                         Abun = colMeans(network_data1))
nodes1 <- merge(nodes,Total.abun)
rownames(nodes1) <- nodes1$name
nodes1 <- nodes1[nodes$name,]

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

pdf("Results/Association_analysis_metabolites/Network/Genus-Pathway-DAMs/Wilcox/network_degree.pdf",
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

pdf("Results/Association_analysis_metabolites/Network/Genus-Pathway-DAMs/Wilcox/network_degree_nolabel.pdf",
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
write.table(result,"Results/Association_analysis_metabolites/Network/Genus-Pathway-DAMs/Wilcox/edge.csv",
            sep=",",quote=F,row.names = FALSE)

#### RandomForest
dir.create("Results/Association_analysis_metabolites/Network/Genus-Pathway-DAMs/RandomForest")
marker <- read.table("Results/Metabolites/Diff_analysis/RandomForest/importance.txt",
                     header = TRUE)
metabolites1 <- metabolites[metabolites$compound %in% marker$ID,]
network_data <- cbind(t(bacteria[[5]][,2:ncol(bacteria[[5]])]),
                      t(virus[[5]][,2:ncol(virus[[5]])]),
                      t(L3[,2:ncol(L3)]))
colnames(network_data) <- c(bacteria[[5]]$Genus,virus[[5]]$Genus,L3$Description)
network_data <- network_data[,colMeans(network_data) >= 0.1]
f <- function(x) sum(x==0)
d <- apply(network_data,2,f)
network_data1 <- network_data[,d < nrow(network_data)*0.4]
network_data1 <- cbind(network_data1,t(metabolites1[,2:ncol(metabolites1)]))
network_data1 <- network_data1[,unique(colnames(network_data1))]

nodes <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    as_tibble(what = "vertices")
nodes$Type <- ifelse(nodes$name %in% bacteria[[5]]$Genus,"Bacteria",
                     ifelse(nodes$name %in% virus[[5]]$Genus,"Virus",
                            ifelse(nodes$name %in% L3$Description,"Function",
                                   "Metabolites")))
network_data1 <- as.data.frame(network_data1)
Total.abun <- data.frame(name = colnames(network_data1),
                         Abun = colMeans(network_data1))
nodes1 <- merge(nodes,Total.abun)
rownames(nodes1) <- nodes1$name
nodes1 <- nodes1[nodes$name,]

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

pdf("Results/Association_analysis_metabolites/Network/Genus-Pathway-DAMs/RandomForest/network_degree.pdf",
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

pdf("Results/Association_analysis_metabolites/Network/Genus-Pathway-DAMs/RandomForest/network_degree_nolabel.pdf",
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
write.table(result,"Results/Association_analysis_metabolites/Network/Genus-Pathway-DAMs/RandomForest/edge.csv",
            sep=",",quote=F,row.names = FALSE)

#### Lefse
dir.create("Results/Association_analysis_metabolites/Network/Genus-Pathway-DAMs/Lefse")
marker <- read.table("Results/Metabolites/Diff_analysis/LEfSe/Lefse_result.txt",
                     header = TRUE)
metabolites1 <- metabolites[metabolites$compound %in% rownames(marker),]
network_data <- cbind(t(bacteria[[5]][,2:ncol(bacteria[[5]])]),
                      t(virus[[5]][,2:ncol(virus[[5]])]),
                      t(L3[,2:ncol(L3)]))
colnames(network_data) <- c(bacteria[[5]]$Genus,virus[[5]]$Genus,L3$Description)
network_data <- network_data[,colMeans(network_data) >= 0.1]
f <- function(x) sum(x==0)
d <- apply(network_data,2,f)
network_data1 <- network_data[,d < nrow(network_data)*0.4]
network_data1 <- cbind(network_data1,t(metabolites1[,2:ncol(metabolites1)]))
network_data1 <- network_data1[,unique(colnames(network_data1))]

nodes <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    as_tibble(what = "vertices")
nodes$Type <- ifelse(nodes$name %in% bacteria[[5]]$Genus,"Bacteria",
                     ifelse(nodes$name %in% virus[[5]]$Genus,"Virus",
                            ifelse(nodes$name %in% L3$Description,"Function",
                                   "Metabolites")))
network_data1 <- as.data.frame(network_data1)
Total.abun <- data.frame(name = colnames(network_data1),
                         Abun = colMeans(network_data1))
nodes1 <- merge(nodes,Total.abun)
rownames(nodes1) <- nodes1$name
nodes1 <- nodes1[nodes$name,]

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

pdf("Results/Association_analysis_metabolites/Network/Genus-Pathway-DAMs/Lefse/network_degree.pdf",
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

pdf("Results/Association_analysis_metabolites/Network/Genus-Pathway-DAMs/Lefse/network_degree_nolabel.pdf",
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
write.table(result,"Results/Association_analysis_metabolites/Network/Genus-Pathway-DAMs/Lefse/edge.csv",
            sep=",",quote=F,row.names = FALSE)

### Genus-ko-DAMs
dir.create("Results/Association_analysis_metabolites/Network/Genus-ko-DAMs")

#### Wilcox
dir.create("Results/Association_analysis_metabolites/Network/Genus-ko-DAMs/Wilcox")
marker <- read.table("Results/Metabolites/Diff_analysis/wilcox/wilcox_biomarker.txt",
                     header = TRUE)
metabolites1 <- metabolites[metabolites$compound %in% marker$V1,]
network_data <- cbind(t(bacteria[[5]][,2:ncol(bacteria[[5]])]),
                      t(virus[[5]][,2:ncol(virus[[5]])]),
                      t(ko[,2:ncol(ko)]))
colnames(network_data) <- c(bacteria[[5]]$Genus,virus[[5]]$Genus,ko$Description)
network_data <- network_data[,colMeans(network_data) >= 0.1]
f <- function(x) sum(x==0)
d <- apply(network_data,2,f)
network_data1 <- network_data[,d < nrow(network_data)*0.4]
network_data1 <- cbind(network_data1,t(metabolites1[,2:ncol(metabolites1)]))
network_data1 <- network_data1[,unique(colnames(network_data1))]

nodes <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    as_tibble(what = "vertices")
nodes$Type <- ifelse(nodes$name %in% bacteria[[5]]$Genus,"Bacteria",
                     ifelse(nodes$name %in% virus[[5]]$Genus,"Virus",
                            ifelse(nodes$name %in% ko$Description,"Function",
                                   "Metabolites")))
network_data1 <- as.data.frame(network_data1)
Total.abun <- data.frame(name = colnames(network_data1),
                         Abun = colMeans(network_data1))
nodes1 <- merge(nodes,Total.abun)
rownames(nodes1) <- nodes1$name
nodes1 <- nodes1[nodes$name,]

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

pdf("Results/Association_analysis_metabolites/Network/Genus-ko-DAMs/Wilcox/network_degree.pdf",
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

pdf("Results/Association_analysis_metabolites/Network/Genus-ko-DAMs/Wilcox/network_degree_nolabel.pdf",
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
write.table(result,"Results/Association_analysis_metabolites/Network/Genus-ko-DAMs/Wilcox/edge.csv",
            sep=",",quote=F,row.names = FALSE)

#### RandomForest
dir.create("Results/Association_analysis_metabolites/Network/Genus-ko-DAMs/RandomForest")
marker <- read.table("Results/Metabolites/Diff_analysis/RandomForest/importance.txt",
                     header = TRUE)
metabolites1 <- metabolites[metabolites$compound %in% marker$ID,]
network_data <- cbind(t(bacteria[[5]][,2:ncol(bacteria[[5]])]),
                      t(virus[[5]][,2:ncol(virus[[5]])]),
                      t(ko[,2:ncol(ko)]))
colnames(network_data) <- c(bacteria[[5]]$Genus,virus[[5]]$Genus,ko$Description)
network_data <- network_data[,colMeans(network_data) >= 0.1]
f <- function(x) sum(x==0)
d <- apply(network_data,2,f)
network_data1 <- network_data[,d < nrow(network_data)*0.4]
network_data1 <- cbind(network_data1,t(metabolites1[,2:ncol(metabolites1)]))
network_data1 <- network_data1[,unique(colnames(network_data1))]

nodes <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    as_tibble(what = "vertices")
nodes$Type <- ifelse(nodes$name %in% bacteria[[5]]$Genus,"Bacteria",
                     ifelse(nodes$name %in% virus[[5]]$Genus,"Virus",
                            ifelse(nodes$name %in% ko$Description,"Function",
                                   "Metabolites")))
network_data1 <- as.data.frame(network_data1)
Total.abun <- data.frame(name = colnames(network_data1),
                         Abun = colMeans(network_data1))
nodes1 <- merge(nodes,Total.abun)
rownames(nodes1) <- nodes1$name
nodes1 <- nodes1[nodes$name,]

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

pdf("Results/Association_analysis_metabolites/Network/Genus-ko-DAMs/RandomForest/network_degree.pdf",
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

pdf("Results/Association_analysis_metabolites/Network/Genus-ko-DAMs/RandomForest/network_degree_nolabel.pdf",
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
write.table(result,"Results/Association_analysis_metabolites/Network/Genus-ko-DAMs/RandomForest/edge.csv",
            sep=",",quote=F,row.names = FALSE)

#### Lefse
dir.create("Results/Association_analysis_metabolites/Network/Genus-ko-DAMs/Lefse")
marker <- read.table("Results/Metabolites/Diff_analysis/LEfSe/Lefse_result.txt",
                     header = TRUE)
metabolites1 <- metabolites[metabolites$compound %in% rownames(marker),]
network_data <- cbind(t(bacteria[[5]][,2:ncol(bacteria[[5]])]),
                      t(virus[[5]][,2:ncol(virus[[5]])]),
                      t(ko[,2:ncol(ko)]))
colnames(network_data) <- c(bacteria[[5]]$Genus,virus[[5]]$Genus,ko$Description)
network_data <- network_data[,colMeans(network_data) >= 0.1]
f <- function(x) sum(x==0)
d <- apply(network_data,2,f)
network_data1 <- network_data[,d < nrow(network_data)*0.4]
network_data1 <- cbind(network_data1,t(metabolites1[,2:ncol(metabolites1)]))
network_data1 <- network_data1[,unique(colnames(network_data1))]

nodes <- linkET::correlate(network_data1) %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    as_tibble(what = "vertices")
nodes$Type <- ifelse(nodes$name %in% bacteria[[5]]$Genus,"Bacteria",
                     ifelse(nodes$name %in% virus[[5]]$Genus,"Virus",
                            ifelse(nodes$name %in% ko$Description,"Function",
                                   "Metabolites")))
network_data1 <- as.data.frame(network_data1)
Total.abun <- data.frame(name = colnames(network_data1),
                         Abun = colMeans(network_data1))
nodes1 <- merge(nodes,Total.abun)
rownames(nodes1) <- nodes1$name
nodes1 <- nodes1[nodes$name,]

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

pdf("Results/Association_analysis_metabolites/Network/Genus-ko-DAMs/Lefse/network_degree.pdf",
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

pdf("Results/Association_analysis_metabolites/Network/Genus-ko-DAMs/Lefse/network_degree_nolabel.pdf",
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
write.table(result,"Results/Association_analysis_metabolites/Network/Genus-ko-DAMs/Lefse/edge.csv",
            sep=",",quote=F,row.names = FALSE)

## Heatmap
dir.create("Results/Association_analysis_metabolites/Heatmap")
source("Functions/Correlation.heatmap.R")

### Bacteria-DAM
dir.create("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM")

#### Wilcox
dir.create("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox")
marker2 <- read.table("Results/Metabolites//Diff_analysis/Wilcox/wilcox_biomarker.txt",
                      header = TRUE)
#### Phylum
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Phylum/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[1]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Class/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[2]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Order/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[3]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Family/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[4]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[5]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Wilcox/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### RandomForest
dir.create("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest")
marker2 <- read.table("Results/Metabolites/Diff_analysis/RandomForest/importance.txt",
                      header = TRUE)

#### Phylum
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Phylum/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[1]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Class/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[2]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Order/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[3]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Family/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[4]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Genus/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[5]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/RandomForest/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Lefse
marker1 <- read.table("Results/Bacteria/Diff_analysis/LEfSe/Lefse_result.txt",
                      header = TRUE,sep = "\t")
marker2 <- read.table("Results/Metabolites/Diff_analysis/LEfSe/Lefse_result.txt",
                      header = TRUE)
marker1$V1 <- rownames(marker1)
marker1 <- marker1[grep("g__",marker1$V1),]
marker1$V1 <- gsub("g__","",marker1$V1)
marker1 <- marker1[,c("V1",colnames(marker1)[1:(ncol(marker1)-1)])]
marker2$V1 <- rownames(marker2)
marker2 <- marker2[,c("V1",colnames(marker2)[1:(ncol(marker2)-1)])]

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Lefse")
    result <- correlation.heatmap(marker1,marker2,
                                  bacteria[[5]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Lefse/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Lefse/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Lefse/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Share
source("Functions/overLapper.new.r")
dir.create("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Shared")
marker1.w <- read.table("Results/Bacteria/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                        header = TRUE,sep = "\t")
marker1.r <- read.table("Results/Bacteria/Diff_analysis/RandomForest/Genus/importance.txt",
                        header = TRUE,sep = "\t")

species.v <- list(wilcox = marker1.w$V1,RandomForest = marker1.r$ID,
                  Lefse = marker1$V1)

pdf("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Shared/venn_bacteria.pdf",
    width=9,height=9,pointsize=16)
OLlist <- overLapper(setlist=species.v, sep="", type="vennsets",keepdups=FALSE)
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts,mymain="")
dev.off()

inter.v <- get.venn.partitions(species.v)
shared.species.v <- marker1[marker1$V1 %in% unlist(inter.v$..values..[1]),]

marker2.w <- read.table("Results/Metabolites/Diff_analysis/Wilcox/wilcox_biomarker.txt",
                        header = TRUE)
marker2.r <- read.table("Results/Metabolites/Diff_analysis/RandomForest/importance.txt",
                        header = TRUE)

species.f <- list(wilcox = marker2.w$V1,RandomForest = marker2.r$ID,
                  Lefse = marker2$V1)

pdf("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Shared/venn_DAM.pdf",
    width=9,height=9,pointsize=16)
OLlist <- overLapper(setlist=species.f, sep="", type="vennsets",keepdups=FALSE)
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts,mymain="")
dev.off()

inter.f <- get.venn.partitions(species.f)
shared.species.f <- marker2[unlist(inter.f$..values..[1]),]

if (nrow(shared.species.v) > 0 & nrow(shared.species.f) > 0) {
    result <- correlation.heatmap(shared.species.v,shared.species.f,
                                  bacteria[[5]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Shared/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Shared/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Bacteria_DAM/Shared/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

### Virus-DAM
dir.create("Results/Association_analysis_metabolites/Heatmap/Virus_DAM")

#### Wilcox
dir.create("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox")
marker2 <- read.table("Results/Metabolites//Diff_analysis/Wilcox/wilcox_biomarker.txt",
                      header = TRUE)
#### Phylum
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Phylum/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[1]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Class/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[2]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Order/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[3]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Family/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[4]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Virus/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[5]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Wilcox/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### RandomForest
dir.create("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest")
marker2 <- read.table("Results/Metabolites/Diff_analysis/RandomForest/importance.txt",
                      header = TRUE)

#### Phylum
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Phylum/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Phylum")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[1]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Phylum/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Phylum/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Phylum/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Class
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Class/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Class")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[2]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Class/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Class/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Class/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Order
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Order/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Order")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[3]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Order/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Order/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Order/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Family
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Family/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Family")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[4]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Family/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Family/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Family/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Genus
marker1 <- read.table("Results/Virus/Diff_analysis/RandomForest/Genus/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Genus")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[5]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Genus/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Genus/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/RandomForest/Genus/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Lefse
marker1 <- read.table("Results/Virus/Diff_analysis/LEfSe/Lefse_result.txt",
                      header = TRUE,sep = "\t")
marker2 <- read.table("Results/Metabolites/Diff_analysis/LEfSe/Lefse_result.txt",
                      header = TRUE)
marker1$V1 <- rownames(marker1)
marker1 <- marker1[grep("g__",marker1$V1),]
marker1$V1 <- gsub("g__","",marker1$V1)
marker1 <- marker1[,c("V1",colnames(marker1)[1:(ncol(marker1)-1)])]
marker2$V1 <- rownames(marker2)
marker2 <- marker2[,c("V1",colnames(marker2)[1:(ncol(marker2)-1)])]

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Lefse")
    result <- correlation.heatmap(marker1,marker2,
                                  virus[[5]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Lefse/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Lefse/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Lefse/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Share
source("Functions/overLapper.new.r")
dir.create("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Shared")
marker1.w <- read.table("Results/Virus/Diff_analysis/Wilcox/Genus/wilcox_biomarker.txt",
                        header = TRUE,sep = "\t")
marker1.r <- read.table("Results/Virus/Diff_analysis/RandomForest/Genus/importance.txt",
                        header = TRUE,sep = "\t")

species.v <- list(wilcox = marker1.w$V1,RandomForest = marker1.r$ID,
                  Lefse = marker1$V1)

pdf("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Shared/venn_virus.pdf",
    width=9,height=9,pointsize=16)
OLlist <- overLapper(setlist=species.v, sep="", type="vennsets",keepdups=FALSE)
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts,mymain="")
dev.off()

inter.v <- get.venn.partitions(species.v)
shared.species.v <- marker1[marker1$V1 %in% unlist(inter.v$..values..[1]),]

marker2.w <- read.table("Results/Metabolites/Diff_analysis/Wilcox/wilcox_biomarker.txt",
                        header = TRUE)
marker2.r <- read.table("Results/Metabolites/Diff_analysis/RandomForest/importance.txt",
                        header = TRUE)

species.f <- list(wilcox = marker2.w$V1,RandomForest = marker2.r$ID,
                  Lefse = marker2$V1)

pdf("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Shared/venn_DAM.pdf",
    width=9,height=9,pointsize=16)
OLlist <- overLapper(setlist=species.f, sep="", type="vennsets",keepdups=FALSE)
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts,mymain="")
dev.off()

inter.f <- get.venn.partitions(species.f)
shared.species.f <- marker2[unlist(inter.f$..values..[1]),]

if (nrow(shared.species.v) > 0 & nrow(shared.species.f) > 0) {
    result <- correlation.heatmap(shared.species.v,shared.species.f,
                                  virus[[5]],metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Shared/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Shared/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Virus_DAM/Shared/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

### Function-DAM
dir.create("Results/Association_analysis_metabolites/Heatmap/Function_DAM")

#### Wilcox
dir.create("Results/Association_analysis_metabolites/Heatmap/Function_DAM/Wilcox")
marker2 <- read.table("Results/Metabolites//Diff_analysis/Wilcox/wilcox_biomarker.txt",
                      header = TRUE)
#### L1
marker1 <- read.table("Results/Gene/Diff_analysis/Wilcox/L1/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Function_DAM/Wilcox/L1")
    result <- correlation.heatmap(marker1,marker2,
                                  L1,metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/Wilcox/L1/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/Wilcox/L1/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Function_DAM/Wilcox/L1/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### L2
marker1 <- read.table("Results/Gene/Diff_analysis/Wilcox/L2/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Function_DAM/Wilcox/L2")
    result <- correlation.heatmap(marker1,marker2,
                                  L2,metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/Wilcox/L2/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/Wilcox/L2/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Function_DAM/Wilcox/L2/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### L3
marker1 <- read.table("Results/Gene/Diff_analysis/Wilcox/L3/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Function_DAM/Wilcox/L3")
    result <- correlation.heatmap(marker1,marker2,
                                  L3,metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/Wilcox/L3/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/Wilcox/L3/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Function_DAM/Wilcox/L3/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### ko
marker1 <- read.table("Results/Gene/Diff_analysis/Wilcox/ko/wilcox_biomarker.txt",
                      header = TRUE,sep = "\t")
if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Function_DAM/Wilcox/ko")
    result <- correlation.heatmap(marker1,marker2,
                                  ko,metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/Wilcox/ko/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/Wilcox/ko/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Function_DAM/Wilcox/ko/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}


#### RandomForest
dir.create("Results/Association_analysis_metabolites/Heatmap/Function_DAM/RandomForest")
marker2 <- read.table("Results/Metabolites/Diff_analysis/RandomForest/importance.txt",
                      header = TRUE,sep = "\t")

#### L1
marker1 <- read.table("Results/Gene/Diff_analysis/RandomForest/L1/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Function_DAM/RandomForest/L1")
    result <- correlation.heatmap(marker1,marker2,
                                  L1,metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/RandomForest/L1/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/RandomForest/L1/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Function_DAM/RandomForest/L1/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### L2
marker1 <- read.table("Results/Gene/Diff_analysis/RandomForest/L2/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Function_DAM/RandomForest/L2")
    result <- correlation.heatmap(marker1,marker2,
                                  L2,metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/RandomForest/L2/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/RandomForest/L2/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Function_DAM/RandomForest/L2/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### L3
marker1 <- read.table("Results/Gene/Diff_analysis/RandomForest/L3/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Function_DAM/RandomForest/L3")
    result <- correlation.heatmap(marker1,marker2,
                                  L3,metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/RandomForest/L3/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/RandomForest/L3/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Function_DAM/RandomForest/L3/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### ko
marker1 <- read.table("Results/Gene/Diff_analysis/RandomForest/ko/importance.txt",
                      header = TRUE,sep = "\t")

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Function_DAM/RandomForest/ko")
    result <- correlation.heatmap(marker1,marker2,
                                  ko,metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/RandomForest/ko/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/RandomForest/ko/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Function_DAM/RandomForest/ko/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Lefse
marker1 <- read.table("Results/Gene/Diff_analysis/LEfSe/Lefse_result.txt",
                      header = TRUE,sep = "\t")
marker2 <- read.table("Results/Metabolites/Diff_analysis/LEfSe/Lefse_result.txt",
                      header = TRUE,sep = "\t")
marker1$V1 <- rownames(marker1)
marker1 <- marker1[,c("V1",colnames(marker1)[1:(ncol(marker1)-1)])]
marker2$V1 <- rownames(marker2)
marker2 <- marker2[,c("V1",colnames(marker2)[1:(ncol(marker2)-1)])]

if (nrow(marker1) > 0 & nrow(marker2) > 0) {
    dir.create("Results/Association_analysis_metabolites/Heatmap/Function_DAM/Lefse")
    result <- correlation.heatmap(marker1,marker2,
                                  ko,metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/Lefse/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/Lefse/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Function_DAM/Lefse/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}

#### Share
source("Functions/overLapper.new.r")
dir.create("Results/Association_analysis_metabolites/Heatmap/Function_DAM/Shared")
marker1.w <- read.table("Results/Gene/Diff_analysis/Wilcox/ko/wilcox_biomarker.txt",
                        header = TRUE,sep = "\t")
marker1.r <- read.table("Results/Gene/Diff_analysis/RandomForest/ko/importance.txt",
                        header = TRUE,sep = "\t")

species.v <- list(wilcox = marker1.w$V1,RandomForest = marker1.r$ID,
                  Lefse = marker1$V1)

pdf("Results/Association_analysis_metabolites/Heatmap/Function_DAM/Shared/venn_function.pdf",
    width=9,height=9,pointsize=16)
OLlist <- overLapper(setlist=species.v, sep="", type="vennsets",keepdups=FALSE)
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts,mymain="")
dev.off()

inter.v <- get.venn.partitions(species.v)
shared.species.v <- marker1[marker1$V1 %in% unlist(inter.v$..values..[1]),]

marker2.w <- read.table("Results/Metabolites/Diff_analysis/Wilcox/wilcox_biomarker.txt",
                        header = TRUE)
marker2.r <- read.table("Results/Metabolites/Diff_analysis/RandomForest/importance.txt",
                        header = TRUE,sep = "\t")

species.f <- list(wilcox = marker2.w$V1,RandomForest = marker2.r$ID,
                  Lefse = marker2$V1)

pdf("Results/Association_analysis_metabolites/Heatmap/Function_DAM/Shared/venn_DAM.pdf",
    width=9,height=9,pointsize=16)
OLlist <- overLapper(setlist=species.f, sep="", type="vennsets",keepdups=FALSE)
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts,mymain="")
dev.off()

inter.f <- get.venn.partitions(species.f)
shared.species.f <- marker2[unlist(inter.f$..values..[1]),]

if (nrow(shared.species.v) > 0 & nrow(shared.species.f) > 0) {
    result <- correlation.heatmap(shared.species.v,shared.species.f,
                                  ko,metabolites)
    write.table(result[[1]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/Shared/correlation.xls",
                sep = "\t",quote = F,col.names = NA)
    write.table(result[[2]],
                "Results/Association_analysis_metabolites/Heatmap/Function_DAM/Shared/pvalue.xls",
                sep = "\t",quote = F,col.names = NA)
    pdf("Results/Association_analysis_metabolites/Heatmap/Function_DAM/Shared/heatmap.pdf",
        height = 3 + 0.2*nrow(marker1),width = 5 + 0.2*nrow(marker2))
    print(result[[3]])
    dev.off()
}





