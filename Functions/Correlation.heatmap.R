correlation.heatmap <- function(otu.wilcox.biomarker,metabolites.wilcox.biomarker,
                                otu,metabolites){
    colnames(otu.wilcox.biomarker)[1] <- "V1"
    colnames(metabolites.wilcox.biomarker)[1] <- "V1"
    rownames(otu) <- otu[,1]
    otu <- otu[,-1]
    rownames(metabolites) <- metabolites[,1]
    metabolites <- metabolites[,-1]
    otu1 <- otu[rownames(otu) %in% otu.wilcox.biomarker$V1,]
    metabolites1 <- metabolites[rownames(metabolites) %in% metabolites.wilcox.biomarker$V1,]
    aa <- match(colnames(otu1),colnames(metabolites1))
    metabolites1 <- metabolites1[,aa]
    res=ggcor::fast_correlate(t(metabolites1),t(otu1),method = "pearson",
                    p.adjust = TRUE,p.adjust.method = "holm")
    if ((nrow(otu1) + nrow(metabolites1)) > 2) {
        p1 <- pheatmap(res$r,
                       display_numbers = matrix(ifelse(res$p <= 0.01, "++", 
                                                       ifelse(res$p <= 0.05 ,"+"," ")), 
                                                nrow(res$p)),
                       fontsize_number=8, number_color = "white", fontsize=8,
                       cluster_rows = ifelse(nrow(metabolites1) > 1,TRUE,FALSE),
                       cluster_cols = ifelse(nrow(otu1) > 1,TRUE,FALSE))
    }
    if ((nrow(otu1) + nrow(metabolites1)) > 2) {
        result <- list(res$r,res$p,p1)
    }else{
        result <- list(res$r,res$p,NA)
    }
    return(result)
}


