mantel.proc <- function(ARG_sub,MRG_sub,MGE_sub){
    ARG_beta <- t(ARG_sub)
    ARG_dis <- vegdist(ARG_beta)
    MRG_beta <- t(MRG_sub)
    aa <- match(rownames(ARG_beta),rownames(MRG_beta))
    MRG_beta <- MRG_beta[aa,]
    MRG_dis <- vegdist(MRG_beta)
    MGE_beta <- t(MGE_sub)
    aa <- match(rownames(ARG_beta),rownames(MGE_beta))
    MGE_beta <- MGE_beta[aa,]
    MGE_dis <- vegdist(MGE_beta)
    
    Mantel_ARG_MRG <- mantel(ARG_dis,MRG_dis)
    Mantel_ARG_MGE <- mantel(ARG_dis,MGE_dis)
    Mantel_MRG_MGE <- mantel(MRG_dis,MGE_dis)
    
    ARG_mds <- monoMDS(ARG_dis)
    MRG_mds <- monoMDS(MRG_dis)
    MGE_mds <- monoMDS(MGE_dis)
    
    pro_ARG_MRG <- procrustes(ARG_mds,MRG_mds)
    pro_ARG_MGE <- procrustes(ARG_mds,MGE_mds)
    pro_MRG_MGE <- procrustes(MRG_mds,MGE_mds)
    
    protest_ARG_MRG <- protest(ARG_mds,MRG_mds)
    protest_ARG_MGE <- protest(ARG_mds,MGE_mds)
    protest_MRG_MGE <- protest(MRG_mds,MGE_mds)
    
    Y <- cbind(data.frame(pro_ARG_MRG$Yrot), data.frame(pro_ARG_MRG$X))
    X <- data.frame(pro_ARG_MRG$rotation)
    Y$ID <- rownames(Y)
    ARG_MRG <- ggplot(Y) +
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
        labs(title="Correlation between bacteria and viruses") + 
        geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
        geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
        geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
        geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
        theme(plot.title = element_text(size=12,colour = "black",hjust = 0.5,face = "bold"))
    
    Y <- cbind(data.frame(pro_ARG_MGE$Yrot), data.frame(pro_ARG_MGE$X))
    X <- data.frame(pro_ARG_MGE$rotation)
    Y$ID <- rownames(Y)
    ARG_MGE <- ggplot(Y) +
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
        labs(title="Correlation between bacteria and functions") + 
        geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
        geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
        geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
        geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
        theme(plot.title = element_text(size=12,colour = "black",hjust = 0.5,face = "bold"))
    
    Y <- cbind(data.frame(pro_MRG_MGE$Yrot), data.frame(pro_MRG_MGE$X))
    X <- data.frame(pro_MRG_MGE$rotation)
    Y$ID <- rownames(Y)
    MRG_MGE <- ggplot(Y) +
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
        labs(title="Correlation between viruses and functions") + 
        geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
        geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
        geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
        geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
        theme(plot.title = element_text(size=12,colour = "black",hjust = 0.5,face = "bold"))
    
    Proc <-  ARG_MRG + ARG_MGE + MRG_MGE + 
        plot_layout(ncol = 3)
    result <- list(ARG_MRG,ARG_MGE,MRG_MGE,Proc)
    return(result)
}


mantel.proc.t <- function(ARG_sub,MRG_sub,MGE_sub){
    ARG_beta <- t(ARG_sub)
    ARG_dis <- vegdist(ARG_beta)
    MRG_beta <- t(MRG_sub)
    aa <- match(rownames(ARG_beta),rownames(MRG_beta))
    MRG_beta <- MRG_beta[aa,]
    MRG_dis <- vegdist(MRG_beta)
    MGE_beta <- t(MGE_sub)
    aa <- match(rownames(ARG_beta),rownames(MGE_beta))
    MGE_beta <- MGE_beta[aa,]
    MGE_dis <- vegdist(MGE_beta)
    
    Mantel_ARG_MRG <- mantel(ARG_dis,MRG_dis)
    Mantel_ARG_MGE <- mantel(ARG_dis,MGE_dis)
    Mantel_MRG_MGE <- mantel(MRG_dis,MGE_dis)
    
    ARG_mds <- monoMDS(ARG_dis)
    MRG_mds <- monoMDS(MRG_dis)
    MGE_mds <- monoMDS(MGE_dis)
    
    pro_ARG_MRG <- procrustes(ARG_mds,MRG_mds)
    pro_ARG_MGE <- procrustes(ARG_mds,MGE_mds)
    pro_MRG_MGE <- procrustes(MRG_mds,MGE_mds)
    
    protest_ARG_MRG <- protest(ARG_mds,MRG_mds)
    protest_ARG_MGE <- protest(ARG_mds,MGE_mds)
    protest_MRG_MGE <- protest(MRG_mds,MGE_mds)
    
    Mantel_r = c(Mantel_ARG_MRG$statistic,Mantel_ARG_MGE$statistic,
                 Mantel_MRG_MGE$statistic)
    Mantel_p = c(Mantel_ARG_MRG$signif,Mantel_ARG_MGE$signif,
                 Mantel_MRG_MGE$signif)
    Proc_r = c(protest_ARG_MRG$ss,protest_ARG_MGE$ss,protest_MRG_MGE$ss)
    Proc_p = c(protest_ARG_MRG$signif,protest_ARG_MGE$signif,
               protest_MRG_MGE$signif)
    Mantel_Proc <- data.frame(Comparison = c("Bacteria vs Virus","Bacteria vs Function",
                                             "Virus vs Function"),
                              Mantel_r,Mantel_p,Proc_r,Proc_p)
    return(Mantel_Proc)
}