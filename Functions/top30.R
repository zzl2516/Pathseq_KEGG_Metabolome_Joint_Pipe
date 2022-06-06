top30 <- function(abundance){
    abundance.1 <- aggregate(abundance[,1:(ncol(abundance)-6)],
                             list(abundance$Phylum),sum)
    colnames(abundance.1)[1] <- "Phylum"
    abundance.2 <- as.data.frame(t(abundance.1))
    colnames(abundance.2) <- abundance.2[1,]
    abundance.2 <- abundance.2[-1,]
    abundance.2 <- as.matrix(abundance.2)
    f.abundance <- matrix(as.numeric(abundance.2),nrow = nrow(abundance.2))
    rownames(f.abundance) <- rownames(abundance.2)
    colnames(f.abundance) <- colnames(abundance.2)
    f.abundance <- t(f.abundance)
    sum <- apply(f.abundance,1,sum) 
    f.abundance <- cbind(f.abundance,sum)
    f.abundance <- as.data.frame(f.abundance)
    f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
    f.abundance <- subset(f.abundance, select = -sum)
    f.abundance.1 <- f.abundance
    sum1 <- apply(f.abundance.1,2,sum)
    if (nrow(f.abundance.1) > 30) {
        f.abundance.1 <- f.abundance.1[1:30,]
        f.abundance.1 <- t(f.abundance.1)
        sum2 <- apply(f.abundance.1,1,sum) 
        Others <- sum1-sum2
        f.abundance.1 <- as.data.frame(cbind(f.abundance.1,Others))
    }else{
        f.abundance.1 <- as.data.frame(t(f.abundance.1))
    }
    Phylum <- f.abundance.1
    
    abundance.1 <- aggregate(abundance[,1:(ncol(abundance)-6)],
                             list(abundance$Class),sum)
    colnames(abundance.1)[1] <- "Class"
    abundance.2 <- as.data.frame(t(abundance.1))
    colnames(abundance.2) <- abundance.2[1,]
    abundance.2 <- abundance.2[-1,]
    abundance.2 <- as.matrix(abundance.2)
    f.abundance <- matrix(as.numeric(abundance.2),nrow = nrow(abundance.2))
    rownames(f.abundance) <- rownames(abundance.2)
    colnames(f.abundance) <- colnames(abundance.2)
    f.abundance <- t(f.abundance)
    sum <- apply(f.abundance,1,sum) 
    f.abundance <- cbind(f.abundance,sum)
    f.abundance <- as.data.frame(f.abundance)
    f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
    f.abundance <- subset(f.abundance, select = -sum)
    f.abundance.1 <- f.abundance
    sum1 <- apply(f.abundance.1,2,sum)
    if (nrow(f.abundance.1) > 30) {
        f.abundance.1 <- f.abundance.1[1:30,]
        f.abundance.1 <- t(f.abundance.1)
        sum2 <- apply(f.abundance.1,1,sum) 
        Others <- sum1-sum2
        f.abundance.1 <- as.data.frame(cbind(f.abundance.1,Others))
    }else{
        f.abundance.1 <- as.data.frame(t(f.abundance.1))
    }
    Class <- f.abundance.1
    
    abundance.1 <- aggregate(abundance[,1:(ncol(abundance)-6)],
                             list(abundance$Order),sum)
    colnames(abundance.1)[1] <- "Order"
    abundance.2 <- as.data.frame(t(abundance.1))
    colnames(abundance.2) <- abundance.2[1,]
    abundance.2 <- abundance.2[-1,]
    abundance.2 <- as.matrix(abundance.2)
    f.abundance <- matrix(as.numeric(abundance.2),nrow = nrow(abundance.2))
    rownames(f.abundance) <- rownames(abundance.2)
    colnames(f.abundance) <- colnames(abundance.2)
    f.abundance <- t(f.abundance)
    sum <- apply(f.abundance,1,sum) 
    f.abundance <- cbind(f.abundance,sum)
    f.abundance <- as.data.frame(f.abundance)
    f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
    f.abundance <- subset(f.abundance, select = -sum)
    f.abundance.1 <- f.abundance
    sum1 <- apply(f.abundance.1,2,sum)
    if (nrow(f.abundance.1) > 30) {
        f.abundance.1 <- f.abundance.1[1:30,]
        f.abundance.1 <- t(f.abundance.1)
        sum2 <- apply(f.abundance.1,1,sum) 
        Others <- sum1-sum2
        f.abundance.1 <- as.data.frame(cbind(f.abundance.1,Others))
    }else{
        f.abundance.1 <- as.data.frame(t(f.abundance.1))
    }
    Order <- f.abundance.1
    
    abundance.1 <- aggregate(abundance[,1:(ncol(abundance)-6)],
                             list(abundance$Family),sum)
    colnames(abundance.1)[1] <- "Family"
    abundance.2 <- as.data.frame(t(abundance.1))
    colnames(abundance.2) <- abundance.2[1,]
    abundance.2 <- abundance.2[-1,]
    abundance.2 <- as.matrix(abundance.2)
    f.abundance <- matrix(as.numeric(abundance.2),nrow = nrow(abundance.2))
    rownames(f.abundance) <- rownames(abundance.2)
    colnames(f.abundance) <- colnames(abundance.2)
    f.abundance <- t(f.abundance)
    sum <- apply(f.abundance,1,sum) 
    f.abundance <- cbind(f.abundance,sum)
    f.abundance <- as.data.frame(f.abundance)
    f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
    f.abundance <- subset(f.abundance, select = -sum)
    f.abundance.1 <- f.abundance
    sum1 <- apply(f.abundance.1,2,sum)
    if (nrow(f.abundance.1) > 30) {
        f.abundance.1 <- f.abundance.1[1:30,]
        f.abundance.1 <- t(f.abundance.1)
        sum2 <- apply(f.abundance.1,1,sum) 
        Others <- sum1-sum2
        f.abundance.1 <- as.data.frame(cbind(f.abundance.1,Others))
    }else{
        f.abundance.1 <- as.data.frame(t(f.abundance.1))
    }
    Family <- f.abundance.1
    
    abundance.1 <- aggregate(abundance[,1:(ncol(abundance)-6)],
                             list(abundance$Genus),sum)
    colnames(abundance.1)[1] <- "Genus"
    abundance.2 <- as.data.frame(t(abundance.1))
    colnames(abundance.2) <- abundance.2[1,]
    abundance.2 <- abundance.2[-1,]
    abundance.2 <- as.matrix(abundance.2)
    f.abundance <- matrix(as.numeric(abundance.2),nrow = nrow(abundance.2))
    rownames(f.abundance) <- rownames(abundance.2)
    colnames(f.abundance) <- colnames(abundance.2)
    f.abundance <- t(f.abundance)
    sum <- apply(f.abundance,1,sum) 
    f.abundance <- cbind(f.abundance,sum)
    f.abundance <- as.data.frame(f.abundance)
    f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
    f.abundance <- subset(f.abundance, select = -sum)
    f.abundance.1 <- f.abundance
    sum1 <- apply(f.abundance.1,2,sum)
    if (nrow(f.abundance.1) > 30) {
        f.abundance.1 <- f.abundance.1[1:30,]
        f.abundance.1 <- t(f.abundance.1)
        sum2 <- apply(f.abundance.1,1,sum) 
        Others <- sum1-sum2
        f.abundance.1 <- as.data.frame(cbind(f.abundance.1,Others))
    }else{
        f.abundance.1 <- as.data.frame(t(f.abundance.1))
    }
    Genus <- f.abundance.1
    
    abundance.1 <- aggregate(abundance[,1:(ncol(abundance)-6)],
                             list(abundance$Species),sum)
    colnames(abundance.1)[1] <- "Species"
    abundance.2 <- as.data.frame(t(abundance.1))
    colnames(abundance.2) <- abundance.2[1,]
    abundance.2 <- abundance.2[-1,]
    abundance.2 <- as.matrix(abundance.2)
    f.abundance <- matrix(as.numeric(abundance.2),nrow = nrow(abundance.2))
    rownames(f.abundance) <- rownames(abundance.2)
    colnames(f.abundance) <- colnames(abundance.2)
    f.abundance <- t(f.abundance)
    sum <- apply(f.abundance,1,sum) 
    f.abundance <- cbind(f.abundance,sum)
    f.abundance <- as.data.frame(f.abundance)
    f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
    f.abundance <- subset(f.abundance, select = -sum)
    f.abundance.1 <- f.abundance
    sum1 <- apply(f.abundance.1,2,sum)
    if (nrow(f.abundance.1) > 30) {
        f.abundance.1 <- f.abundance.1[1:30,]
        f.abundance.1 <- t(f.abundance.1)
        sum2 <- apply(f.abundance.1,1,sum) 
        Others <- sum1-sum2
        f.abundance.1 <- as.data.frame(cbind(f.abundance.1,Others))
    }else{
        f.abundance.1 <- as.data.frame(t(f.abundance.1))
    }
    Species <- f.abundance.1
    
    result <- list(Phylum,Class,Order,Family,Genus,Species)
    return(result)
    
}