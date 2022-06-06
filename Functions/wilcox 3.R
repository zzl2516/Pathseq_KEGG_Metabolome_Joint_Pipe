wilcox.biomarker <- function(otu,group){
    f <- function(x)sum(x==0)
    a <- apply(otu[,2:ncol(otu)],1,f)
    a <- as.data.frame(a)
    otu1 <- cbind(otu[,2:ncol(otu)],a)
    rownames(otu1) <- otu[,1]
    otu2 <- dplyr::filter(otu1,a < (ncol(otu1)-1)*0.4)
    otu2 <- otu2[,1:(ncol(otu2)-1)]
    
    n <- c(1:nrow(otu2))
    p.value <- rep(NA,nrow(otu2))
    for(i in n){
        p.value[i] <- wilcox.test(as.numeric(otu2[i,group[group$Group == unique(group$Group)[1],"variable"]]),
                                  as.numeric(otu2[i,group[group$Group == unique(group$Group)[2],"variable"]]))$p.value
    }
    wilcox <- as.data.frame(cbind(rownames(otu2),p.value))
    wilcox.biomarker <- wilcox[wilcox$p.value < 0.05,]
    wilcox.biomarker <- wilcox.biomarker[order(wilcox.biomarker$p.value),]
    result <- list(wilcox,wilcox.biomarker)
    return(result)
}