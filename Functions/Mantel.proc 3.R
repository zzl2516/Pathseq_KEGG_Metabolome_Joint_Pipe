mantel.proc <- function(ARG_sub,MRG_sub,MGE_sub,BRG_sub){
    ARG_beta <- t(ARG_sub)
    ARG_dis <- vegdist(ARG_beta)
    MRG_beta <- t(MRG_sub)
    aa <- match(rownames(ARG_beta),rownames(MRG_beta))
    MRG_beta <- MRG_beta[aa,]
    MRG_dis <- vegdist(MRG_beta)
    BRG_beta <- t(BRG_sub)
    aa <- match(rownames(ARG_beta),rownames(BRG_beta))
    BRG_beta <- BRG_beta[aa,]
    BRG_dis <- vegdist(BRG_beta,method = "euclidean")
    MGE_beta <- t(MGE_sub)
    aa <- match(rownames(ARG_beta),rownames(MGE_beta))
    MGE_beta <- MGE_beta[aa,]
    MGE_dis <- vegdist(MGE_beta)
    
    Mantel_ARG_BRG <- mantel(ARG_dis,BRG_dis)
    Mantel_MRG_BRG <- mantel(MRG_dis,BRG_dis)
    Mantel_MGE_BRG <- mantel(MGE_dis,BRG_dis)
    
    ARG_mds <- monoMDS(ARG_dis)
    MRG_mds <- monoMDS(MRG_dis)
    BRG_mds <- monoMDS(BRG_dis)
    MGE_mds <- monoMDS(MGE_dis)
    
    pro_ARG_BRG <- procrustes(ARG_mds,BRG_mds)
    pro_MRG_BRG <- procrustes(MRG_mds,BRG_mds)
    pro_MGE_BRG <- procrustes(MGE_mds,BRG_mds)
    
    protest_ARG_BRG <- protest(ARG_mds,BRG_mds)
    protest_MRG_BRG <- protest(MRG_mds,BRG_mds)
    protest_MGE_BRG <- protest(MGE_mds,BRG_mds)
    
    Y <- cbind(data.frame(pro_ARG_BRG$Yrot), data.frame(pro_ARG_BRG$X))
    X <- data.frame(pro_ARG_BRG$rotation)
    Y$ID <- rownames(Y)
    ARG_BRG <- ggplot(Y) +
        geom_segment(aes(x = X1, y = X2, xend = (X1 + MDS1)/2, yend = (X2 + MDS2)/2), 
                     arrow = arrow(length = unit(0, 'cm')),
                     color = "#B2182B", size = 1) +
        geom_segment(aes(x = (X1 + MDS1)/2, y = (X2 + MDS2)/2, xend = MDS1, yend = MDS2), 
                     arrow = arrow(length = unit(0, 'cm')),
                     color = "#56B4E9", size = 1) +
        geom_point(aes(X1, X2), fill = "#B2182B", size = 3, shape = 21) +
        geom_point(aes(MDS1, MDS2), fill = "#56B4E9", size = 3, shape = 21) +
        theme(panel.grid = element_blank(), 
              panel.background = element_rect(color = 'black', fill = 'transparent'),
              legend.key = element_rect(fill = 'transparent'),
              axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_text(colour='black', size=10),
              axis.title.y=element_text(colour='black', size=10),
              axis.text=element_text(colour='black',size=8)) +
        labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
        labs(title="Correlation between bacteria and metabolome") + 
        geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
        geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
        geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
        geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
        theme(plot.title = element_text(size=12,colour = "black",hjust = 0.5,face = "bold"))
    
    Y <- cbind(data.frame(pro_MRG_BRG$Yrot), data.frame(pro_MRG_BRG$X))
    X <- data.frame(pro_MRG_BRG$rotation)
    Y$ID <- rownames(Y)
    MRG_BRG <- ggplot(Y) +
        geom_segment(aes(x = X1, y = X2, xend = (X1 + MDS1)/2, yend = (X2 + MDS2)/2), 
                     arrow = arrow(length = unit(0, 'cm')),
                     color = "#B2182B", size = 1) +
        geom_segment(aes(x = (X1 + MDS1)/2, y = (X2 + MDS2)/2, xend = MDS1, yend = MDS2), 
                     arrow = arrow(length = unit(0, 'cm')),
                     color = "#56B4E9", size = 1) +
        geom_point(aes(X1, X2), fill = "#B2182B", size = 3, shape = 21) +
        geom_point(aes(MDS1, MDS2), fill = "#56B4E9", size = 3, shape = 21) +
        theme(panel.grid = element_blank(), 
              panel.background = element_rect(color = 'black', fill = 'transparent'),
              legend.key = element_rect(fill = 'transparent'),
              axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_text(colour='black', size=10),
              axis.title.y=element_text(colour='black', size=10),
              axis.text=element_text(colour='black',size=8)) +
        labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
        labs(title="Correlation between viruses and metabolome") + 
        geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
        geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
        geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
        geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
        theme(plot.title = element_text(size=12,colour = "black",hjust = 0.5,face = "bold"))
    
    Y <- cbind(data.frame(pro_MGE_BRG$Yrot), data.frame(pro_MGE_BRG$X))
    X <- data.frame(pro_MGE_BRG$rotation)
    Y$ID <- rownames(Y)
    MGE_BRG <- ggplot(Y) +
        geom_segment(aes(x = X1, y = X2, xend = (X1 + MDS1)/2, yend = (X2 + MDS2)/2), 
                     arrow = arrow(length = unit(0, 'cm')),
                     color = "#B2182B", size = 1) +
        geom_segment(aes(x = (X1 + MDS1)/2, y = (X2 + MDS2)/2, xend = MDS1, yend = MDS2), 
                     arrow = arrow(length = unit(0, 'cm')),
                     color = "#56B4E9", size = 1) +
        geom_point(aes(X1, X2), fill = "#B2182B", size = 3, shape = 21) +
        geom_point(aes(MDS1, MDS2), fill = "#56B4E9", size = 3, shape = 21) +
        theme(panel.grid = element_blank(), 
              panel.background = element_rect(color = 'black', fill = 'transparent'),
              legend.key = element_rect(fill = 'transparent'),
              axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_text(colour='black', size=10),
              axis.title.y=element_text(colour='black', size=10),
              axis.text=element_text(colour='black',size=8)) +
        labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
        labs(title="Correlation between functions and metabolome") + 
        geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
        geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
        geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
        geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
        theme(plot.title = element_text(size=12,colour = "black",hjust = 0.5,face = "bold"))
    
    
    
    Proc <-  ARG_BRG + MRG_BRG + MGE_BRG +
        plot_layout(ncol = 3)
    result <- list(ARG_BRG,MRG_BRG,MGE_BRG,Proc)
    return(result)
}


mantel.proc.t <- function(ARG_sub,MRG_sub,MGE_sub,BRG_sub){
    ARG_beta <- t(ARG_sub)
    ARG_dis <- vegdist(ARG_beta)
    MRG_beta <- t(MRG_sub)
    aa <- match(rownames(ARG_beta),rownames(MRG_beta))
    MRG_beta <- MRG_beta[aa,]
    MRG_dis <- vegdist(MRG_beta)
    BRG_beta <- t(BRG_sub)
    aa <- match(rownames(ARG_beta),rownames(BRG_beta))
    BRG_beta <- BRG_beta[aa,]
    BRG_dis <- vegdist(BRG_beta,method = "euclidean")
    MGE_beta <- t(MGE_sub)
    aa <- match(rownames(ARG_beta),rownames(MGE_beta))
    MGE_beta <- MGE_beta[aa,]
    MGE_dis <- vegdist(MGE_beta)
    
    Mantel_ARG_BRG <- mantel(ARG_dis,BRG_dis)
    Mantel_MRG_BRG <- mantel(MRG_dis,BRG_dis)
    Mantel_MGE_BRG <- mantel(MGE_dis,BRG_dis)
    
    ARG_mds <- monoMDS(ARG_dis)
    MRG_mds <- monoMDS(MRG_dis)
    BRG_mds <- monoMDS(BRG_dis)
    MGE_mds <- monoMDS(MGE_dis)
    
    pro_ARG_BRG <- procrustes(ARG_mds,BRG_mds)
    pro_MRG_BRG <- procrustes(MRG_mds,BRG_mds)
    pro_MGE_BRG <- procrustes(MGE_mds,BRG_mds)
    
    protest_ARG_BRG <- protest(ARG_mds,BRG_mds)
    protest_MRG_BRG <- protest(MRG_mds,BRG_mds)
    protest_MGE_BRG <- protest(MGE_mds,BRG_mds)
    
    Mantel_r = c(Mantel_ARG_BRG$statistic,Mantel_MRG_BRG$statistic,Mantel_MGE_BRG$statistic)
    Mantel_p = c(Mantel_ARG_BRG$signif,Mantel_MRG_BRG$signif,Mantel_MGE_BRG$signif)
    Proc_r = c(protest_ARG_BRG$ss,protest_MRG_BRG$ss,protest_MGE_BRG$ss)
    Proc_p = c(protest_ARG_BRG$signif,protest_MRG_BRG$signif,protest_MGE_BRG$signif)
    Mantel_Proc <- data.frame(Comparison = c("Bacteria vs Metabolites",
                                             "Virus vs Metabolites","Function vs Metabolites"),
                              Mantel_r,Mantel_p,Proc_r,Proc_p)
    return(Mantel_Proc)
}