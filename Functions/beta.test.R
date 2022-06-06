beta.test <- function(ARG_sub,group){
    colnames(group) <- c("sample","Group")
    ARG_beta <- ARG_sub[group[,1],]
    data <- vegdist(ARG_beta)
    
    pcoa<- pcoa(data, correction = "none", rn = NULL)
    PC1 = pcoa$vectors[,1]
    PC2 = pcoa$vectors[,2]
    plotdata <- data.frame(rownames(pcoa$vectors),PC1,PC2)
    colnames(plotdata) <-c("sample","PC1","PC2")
    plotdata <- merge(plotdata,group)
    
    ARG.adonis <- adonis(data~Group,data = plotdata,distance = "bray")
    ARG.anosim <- with(plotdata,anosim(data,Group))
    ARG.mrpp <- with(plotdata,mrpp(data,Group))
    result <- list(ARG.adonis,ARG.anosim,ARG.mrpp)
    return(result)
}