top30 <- function(abundance){
    f.abundance <- abundance[,-1]
    sum <- apply(f.abundance,1,sum) 
    f.abundance <- cbind(f.abundance,sum)
    f.abundance <- as.data.frame(f.abundance)
    f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
    f.abundance <- subset(f.abundance, select = -sum)
    f.abundance <- f.abundance[rowSums(f.abundance) > 0,]
    if (nrow(f.abundance) > 30) {
        f.abundance <- f.abundance[1:30,]
    }
    
    return(f.abundance)
    
}