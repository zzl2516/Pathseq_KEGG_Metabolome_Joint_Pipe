bubblechart.s <- function(abundance.tax,group){
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
    f.abundance.1 <- f.abundance
    if (nrow(f.abundance.1) > 30) {
        f.abundance.1 <- f.abundance.1[1:30,]
    }
    f.abundance.1$Taxon <- rownames(f.abundance.1)
    taxon <- melt(f.abundance.1)
    colnames(taxon) <- c("Taxon","variable","value")
    taxon$variable <- factor(taxon$variable,levels = group$variable)
    taxon <- merge(taxon,group)
    
    p.Phylum <- ggplot(data = taxon,aes(x = variable, y = Taxon)) + 
        geom_point(aes(color = Group,size = value)) + 
        ylab(label = "Relative abundance of dominant phyla") + 
        xlab(label = "") + 
        scale_color_manual(values = cbbPalette) +
        facet_grid(.~Group,scales = "free",space = "free") + 
        theme_bw()+ 
        theme(panel.grid=element_blank()) + 
        theme(panel.border = element_blank()) +
        theme(panel.background=element_rect(fill='transparent', color='black'),
              plot.margin = unit(c(3,5,1,1),"mm")) + 
        theme(axis.text.x=element_text(colour="black",size=10,face = "bold",
                                       angle = 90,vjust = 0.5,hjust = 1)) + 
        theme(axis.text.y=element_text(colour = "black",size = 10)) + 
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"))+ 
        theme(axis.title.y = element_text(size = 18,face = "bold",
                                          margin = unit(c(0,1,0,1),"lines"))) + 
        theme(legend.position = "none",
              strip.text = element_text(size = 14,colour = "black",face = "bold"))
    
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
    f.abundance.1 <- f.abundance
    if (nrow(f.abundance.1) > 30) {
        f.abundance.1 <- f.abundance.1[1:30,]
    }
    f.abundance.1$Taxon <- rownames(f.abundance.1)
    taxon <- melt(f.abundance.1)
    colnames(taxon) <- c("Taxon","variable","value")
    taxon$variable <- factor(taxon$variable,levels = group$variable)
    taxon <- merge(taxon,group)
    
    p.Class <- ggplot(data = taxon,aes(x = variable, y = Taxon)) + 
        geom_point(aes(color = Group,size = value)) + 
        ylab(label = "Relative abundance of dominant classes") + 
        xlab(label = "") + 
        scale_color_manual(values = cbbPalette) +
        facet_grid(.~Group,scales = "free",space = "free") + 
        theme_bw()+ 
        theme(panel.grid=element_blank()) + 
        theme(panel.border = element_blank()) +
        theme(panel.background=element_rect(fill='transparent', color='black'),
              plot.margin = unit(c(3,5,1,1),"mm")) + 
        theme(axis.text.x=element_text(colour="black",size=10,face = "bold",
                                       angle = 90,vjust = 0.5,hjust = 1)) + 
        theme(axis.text.y=element_text(colour = "black",size = 10)) + 
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"))+ 
        theme(axis.title.y = element_text(size = 18,face = "bold",
                                          margin = unit(c(0,1,0,1),"lines"))) + 
        theme(legend.position = "none",
              strip.text = element_text(size = 14,colour = "black",face = "bold"))
    
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
    f.abundance.1 <- f.abundance
    if (nrow(f.abundance.1) > 30) {
        f.abundance.1 <- f.abundance.1[1:30,]
    }
    f.abundance.1$Taxon <- rownames(f.abundance.1)
    taxon <- melt(f.abundance.1)
    colnames(taxon) <- c("Taxon","variable","value")
    taxon$variable <- factor(taxon$variable,levels = group$variable)
    taxon <- merge(taxon,group)
    
    p.Order <- ggplot(data = taxon,aes(x = variable, y = Taxon)) + 
        geom_point(aes(color = Group,size = value)) + 
        ylab(label = "Relative abundance of dominant orders") + 
        xlab(label = "") + 
        scale_color_manual(values = cbbPalette) +
        facet_grid(.~Group,scales = "free",space = "free") + 
        theme_bw()+ 
        theme(panel.grid=element_blank()) + 
        theme(panel.border = element_blank()) +
        theme(panel.background=element_rect(fill='transparent', color='black'),
              plot.margin = unit(c(3,5,1,1),"mm")) + 
        theme(axis.text.x=element_text(colour="black",size=10,face = "bold",
                                       angle = 90,vjust = 0.5,hjust = 1)) + 
        theme(axis.text.y=element_text(colour = "black",size = 10)) + 
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"))+ 
        theme(axis.title.y = element_text(size = 18,face = "bold",
                                          margin = unit(c(0,1,0,1),"lines"))) + 
        theme(legend.position = "none",
              strip.text = element_text(size = 14,colour = "black",face = "bold"))
    
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
    f.abundance.1 <- f.abundance
    if (nrow(f.abundance.1) > 30) {
        f.abundance.1 <- f.abundance.1[1:30,]
    }
    f.abundance.1$Taxon <- rownames(f.abundance.1)
    taxon <- melt(f.abundance.1)
    colnames(taxon) <- c("Taxon","variable","value")
    taxon$variable <- factor(taxon$variable,levels = group$variable)
    taxon <- merge(taxon,group)
    
    p.Family <- ggplot(data = taxon,aes(x = variable, y = Taxon)) + 
        geom_point(aes(color = Group,size = value)) + 
        ylab(label = "Relative abundance of dominant families") + 
        xlab(label = "") + 
        scale_color_manual(values = cbbPalette) +
        facet_grid(.~Group,scales = "free",space = "free") + 
        theme_bw()+ 
        theme(panel.grid=element_blank()) + 
        theme(panel.border = element_blank()) +
        theme(panel.background=element_rect(fill='transparent', color='black'),
              plot.margin = unit(c(3,5,1,1),"mm")) + 
        theme(axis.text.x=element_text(colour="black",size=10,face = "bold",
                                       angle = 90,vjust = 0.5,hjust = 1)) + 
        theme(axis.text.y=element_text(colour = "black",size = 10)) + 
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"))+ 
        theme(axis.title.y = element_text(size = 18,face = "bold",
                                          margin = unit(c(0,1,0,1),"lines"))) + 
        theme(legend.position = "none",
              strip.text = element_text(size = 14,colour = "black",face = "bold"))
    
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
    f.abundance.1 <- f.abundance
    if (nrow(f.abundance.1) > 30) {
        f.abundance.1 <- f.abundance.1[1:30,]
    }
    f.abundance.1$Taxon <- rownames(f.abundance.1)
    taxon <- melt(f.abundance.1)
    colnames(taxon) <- c("Taxon","variable","value")
    taxon$variable <- factor(taxon$variable,levels = group$variable)
    taxon <- merge(taxon,group)
    
    p.Genus <- ggplot(data = taxon,aes(x = variable, y = Taxon)) + 
        geom_point(aes(color = Group,size = value)) + 
        ylab(label = "Relative abundance of dominant genera") + 
        xlab(label = "") + 
        scale_color_manual(values = cbbPalette) +
        facet_grid(.~Group,scales = "free",space = "free") + 
        theme_bw()+ 
        theme(panel.grid=element_blank()) + 
        theme(panel.border = element_blank()) +
        theme(panel.background=element_rect(fill='transparent', color='black'),
              plot.margin = unit(c(3,5,1,1),"mm")) + 
        theme(axis.text.x=element_text(colour="black",size=10,face = "bold",
                                       angle = 90,vjust = 0.5,hjust = 1)) + 
        theme(axis.text.y=element_text(colour = "black",size = 10)) + 
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"))+ 
        theme(axis.title.y = element_text(size = 18,face = "bold",
                                          margin = unit(c(0,1,0,1),"lines"))) + 
        theme(legend.position = "none",
              strip.text = element_text(size = 14,colour = "black",face = "bold"))
    
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
    f.abundance.1 <- f.abundance
    if (nrow(f.abundance.1) > 30) {
        f.abundance.1 <- f.abundance.1[1:30,]
    }
    f.abundance.1$Taxon <- rownames(f.abundance.1)
    taxon <- melt(f.abundance.1)
    colnames(taxon) <- c("Taxon","variable","value")
    taxon$variable <- factor(taxon$variable,levels = group$variable)
    taxon <- merge(taxon,group)
    
    p.Species <- ggplot(data = taxon,aes(x = variable, y = Taxon)) + 
        geom_point(aes(color = Group,size = value)) + 
        ylab(label = "Relative abundance of dominant species") + 
        xlab(label = "") + 
        scale_color_manual(values = cbbPalette) +
        facet_grid(.~Group,scales = "free",space = "free") + 
        theme_bw()+ 
        theme(panel.grid=element_blank()) + 
        theme(panel.border = element_blank()) +
        theme(panel.background=element_rect(fill='transparent', color='black'),
              plot.margin = unit(c(3,5,1,1),"mm")) + 
        theme(axis.text.x=element_text(colour="black",size=10,face = "bold",
                                       angle = 90,vjust = 0.5,hjust = 1)) + 
        theme(axis.text.y=element_text(colour = "black",size = 10)) + 
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"))+ 
        theme(axis.title.y = element_text(size = 18,face = "bold",
                                          margin = unit(c(0,1,0,1),"lines"))) + 
        theme(legend.position = "none",
              strip.text = element_text(size = 14,colour = "black",face = "bold"))

    result <- list(p.Phylum,p.Class,p.Order,p.Family,p.Genus,p.Species)
    return(result)
}