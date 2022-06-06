bubblechart.g <- function(abundance,group){
    f.abundance <- abundance[,-1]
    sum <- apply(f.abundance,1,sum) 
    f.abundance <- cbind(f.abundance,sum)
    f.abundance <- as.data.frame(f.abundance)
    f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
    f.abundance <- subset(f.abundance, select = -sum)
    f.abundance.1 <- f.abundance
    if (nrow(f.abundance.1) > 30) {
        f.abundance.1 <- f.abundance.1[1:30,]
    }
    f.abundance.1$Taxon <- rownames(f.abundance.1)
    taxon <- melt(f.abundance.1)
    colnames(taxon) <- c("Taxon","variable","value")
    
    p.Phylum <- ggplot(data = taxon,aes(x = variable, y = Taxon)) + 
        geom_point(aes(color = variable,size = value)) + 
        ylab(label = "Relative abundance of dominant KEGG terms") + 
        xlab(label = "") + 
        scale_color_manual(values = cbbPalette) +
        theme_bw()+ 
        theme(panel.grid=element_blank()) + 
        theme(panel.border = element_blank()) +
        theme(panel.background=element_rect(fill='transparent', color='black'),
              plot.margin = unit(c(3,5,1,1),"mm")) + 
        theme(axis.text.x=element_text(colour="black",size=10,face = "bold",
                                       angle = 45,hjust = 1,vjust = 1)) + 
        theme(axis.text.y=element_text(colour = "black",size = 10)) + 
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"))+ 
        theme(axis.title.y = element_text(size = 18,face = "bold",
                                          margin = unit(c(0,1,0,1),"lines"))) + 
        theme(legend.position = "none")
    
    return(p.Phylum)
}