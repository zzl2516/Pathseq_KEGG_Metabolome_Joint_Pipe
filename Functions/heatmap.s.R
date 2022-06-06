heatmap.s <- function(abundance.tax,group){
    abundance.1 <- aggregate(abundance.tax[,1:(ncol(abundance.tax)-6)],
                             list(abundance.tax$Phylum),sum)
    colnames(abundance.1)[1] <- "Phylum"
    abundance.2 <- as.data.frame(t(abundance.1))
    colnames(abundance.2) <- abundance.2[1,]
    abundance.2 <- abundance.2[-1,]
    abundance.2 <- as.matrix(abundance.2)
    f.abundance <- matrix(as.numeric(abundance.2),nrow = nrow(abundance.2))
    rownames(f.abundance) <- rownames(abundance.2)
    colnames(f.abundance) <- colnames(abundance.2)
    if ("Unclassified" %in% colnames(f.abundance)) {
        f.abundance <- subset(f.abundance,select = -Unclassified)
    }
    f.abundance <- t(f.abundance)
    sum <- apply(f.abundance,1,sum) 
    f.abundance <- cbind(f.abundance,sum)
    f.abundance <- as.data.frame(f.abundance)
    f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
    f.abundance <- subset(f.abundance, select = -sum)
    f.abundance <- f.abundance[rowSums(f.abundance) > 0,]
    if (nrow(f.abundance) > 30) {
        f.abundance <- f.abundance[1:30,]
    }
    
    p.Phylum <- pheatmap(f.abundance,color = colorRampPalette(c("navy","white","firebrick3"))(100),
                         fontsize = 12,scale = "row",angle_col = 90)
    
    abundance.1 <- aggregate(abundance.tax[,1:(ncol(abundance.tax)-6)],
                             list(abundance.tax$Class),sum)
    colnames(abundance.1)[1] <- "Class"
    abundance.2 <- as.data.frame(t(abundance.1))
    colnames(abundance.2) <- abundance.2[1,]
    abundance.2 <- abundance.2[-1,]
    abundance.2 <- as.matrix(abundance.2)
    f.abundance <- matrix(as.numeric(abundance.2),nrow = nrow(abundance.2))
    rownames(f.abundance) <- rownames(abundance.2)
    colnames(f.abundance) <- colnames(abundance.2)
    if ("Unclassified" %in% colnames(f.abundance)) {
        f.abundance <- subset(f.abundance,select = -Unclassified)
    }
    f.abundance <- t(f.abundance)
    sum <- apply(f.abundance,1,sum) 
    f.abundance <- cbind(f.abundance,sum)
    f.abundance <- as.data.frame(f.abundance)
    f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
    f.abundance <- subset(f.abundance, select = -sum)
    f.abundance <- f.abundance[rowSums(f.abundance) > 0,]
    if (nrow(f.abundance) > 30) {
        f.abundance <- f.abundance[1:30,]
    }
    
    p.Class <- pheatmap(f.abundance,color = colorRampPalette(c("navy","white","firebrick3"))(100),
                         fontsize = 12,scale = "row",angle_col = 90)
    
    abundance.1 <- aggregate(abundance.tax[,1:(ncol(abundance.tax)-6)],
                             list(abundance.tax$Order),sum)
    colnames(abundance.1)[1] <- "Order"
    abundance.2 <- as.data.frame(t(abundance.1))
    colnames(abundance.2) <- abundance.2[1,]
    abundance.2 <- abundance.2[-1,]
    abundance.2 <- as.matrix(abundance.2)
    f.abundance <- matrix(as.numeric(abundance.2),nrow = nrow(abundance.2))
    rownames(f.abundance) <- rownames(abundance.2)
    colnames(f.abundance) <- colnames(abundance.2)
    if ("Unclassified" %in% colnames(f.abundance)) {
        f.abundance <- subset(f.abundance,select = -Unclassified)
    }
    f.abundance <- t(f.abundance)
    sum <- apply(f.abundance,1,sum) 
    f.abundance <- cbind(f.abundance,sum)
    f.abundance <- as.data.frame(f.abundance)
    f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
    f.abundance <- subset(f.abundance, select = -sum)
    f.abundance <- f.abundance[rowSums(f.abundance) > 0,]
    if (nrow(f.abundance) > 30) {
        f.abundance <- f.abundance[1:30,]
    }
    
    p.Order <- pheatmap(f.abundance,color = colorRampPalette(c("navy","white","firebrick3"))(100),
                         fontsize = 12,scale = "row",angle_col = 90)
    
    abundance.1 <- aggregate(abundance.tax[,1:(ncol(abundance.tax)-6)],
                             list(abundance.tax$Family),sum)
    colnames(abundance.1)[1] <- "Family"
    abundance.2 <- as.data.frame(t(abundance.1))
    colnames(abundance.2) <- abundance.2[1,]
    abundance.2 <- abundance.2[-1,]
    abundance.2 <- as.matrix(abundance.2)
    f.abundance <- matrix(as.numeric(abundance.2),nrow = nrow(abundance.2))
    rownames(f.abundance) <- rownames(abundance.2)
    colnames(f.abundance) <- colnames(abundance.2)
    if ("Unclassified" %in% colnames(f.abundance)) {
        f.abundance <- subset(f.abundance,select = -Unclassified)
    }
    f.abundance <- t(f.abundance)
    sum <- apply(f.abundance,1,sum) 
    f.abundance <- cbind(f.abundance,sum)
    f.abundance <- as.data.frame(f.abundance)
    f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
    f.abundance <- subset(f.abundance, select = -sum)
    f.abundance <- f.abundance[rowSums(f.abundance) > 0,]
    if (nrow(f.abundance) > 30) {
        f.abundance <- f.abundance[1:30,]
    }
    
    p.Family <- pheatmap(f.abundance,color = colorRampPalette(c("navy","white","firebrick3"))(100),
                         fontsize = 12,scale = "row",angle_col = 90)
    
    abundance.1 <- aggregate(abundance.tax[,1:(ncol(abundance.tax)-6)],
                             list(abundance.tax$Genus),sum)
    colnames(abundance.1)[1] <- "Genus"
    abundance.2 <- as.data.frame(t(abundance.1))
    colnames(abundance.2) <- abundance.2[1,]
    abundance.2 <- abundance.2[-1,]
    abundance.2 <- as.matrix(abundance.2)
    f.abundance <- matrix(as.numeric(abundance.2),nrow = nrow(abundance.2))
    rownames(f.abundance) <- rownames(abundance.2)
    colnames(f.abundance) <- colnames(abundance.2)
    if ("Unclassified" %in% colnames(f.abundance)) {
        f.abundance <- subset(f.abundance,select = -Unclassified)
    }
    f.abundance <- t(f.abundance)
    sum <- apply(f.abundance,1,sum) 
    f.abundance <- cbind(f.abundance,sum)
    f.abundance <- as.data.frame(f.abundance)
    f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
    f.abundance <- subset(f.abundance, select = -sum)
    f.abundance <- f.abundance[rowSums(f.abundance) > 0,]
    if (nrow(f.abundance) > 30) {
        f.abundance <- f.abundance[1:30,]
    }
    
    p.Genus <- pheatmap(f.abundance,color = colorRampPalette(c("navy","white","firebrick3"))(100),
                         fontsize = 12,scale = "row",angle_col = 90)
    
    abundance.1 <- aggregate(abundance.tax[,1:(ncol(abundance.tax)-6)],
                             list(abundance.tax$Species),sum)
    colnames(abundance.1)[1] <- "Species"
    abundance.2 <- as.data.frame(t(abundance.1))
    colnames(abundance.2) <- abundance.2[1,]
    abundance.2 <- abundance.2[-1,]
    abundance.2 <- as.matrix(abundance.2)
    f.abundance <- matrix(as.numeric(abundance.2),nrow = nrow(abundance.2))
    rownames(f.abundance) <- rownames(abundance.2)
    colnames(f.abundance) <- colnames(abundance.2)
    if ("Unclassified" %in% colnames(f.abundance)) {
        f.abundance <- subset(f.abundance,select = -Unclassified)
    }
    f.abundance <- t(f.abundance)
    sum <- apply(f.abundance,1,sum) 
    f.abundance <- cbind(f.abundance,sum)
    f.abundance <- as.data.frame(f.abundance)
    f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
    f.abundance <- subset(f.abundance, select = -sum)
    f.abundance <- f.abundance[rowSums(f.abundance) > 0,]
    if (nrow(f.abundance) > 30) {
        f.abundance <- f.abundance[1:30,]
    }
    
    p.Species <- pheatmap(f.abundance,color = colorRampPalette(c("navy","white","firebrick3"))(100),
                         fontsize = 12,scale = "row",angle_col = 90)
    
    result <- list(p.Phylum,p.Class,p.Order,p.Family,p.Genus,p.Species)
    return(result)
}