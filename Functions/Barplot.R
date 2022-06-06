barplot.s <- function(abundance,group){
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
    if (nrow(f.abundance.1) > 10) {
        f.abundance.1 <- f.abundance.1[1:10,]
        f.abundance.1 <- t(f.abundance.1)
        sum2 <- apply(f.abundance.1,1,sum) 
        Others <- 100.00001-sum2
        f.abundance.1 <- as.data.frame(cbind(f.abundance.1,Others))
    }
    f.abundance.1$Taxon <- rownames(f.abundance.1)
    taxon <- melt(f.abundance.1)
    if (nrow(f.abundance) > 10) {
        colnames(taxon) <- c("variable","Taxon","value")
        }
    taxon <- merge(taxon,group)
    aa <- f.abundance[,order(f.abundance[1,],decreasing = T)]
    taxon$variable <- factor(taxon$variable,levels = colnames(aa))
    
    p.Phylum <- ggplot(data = taxon,aes(x = variable, y = value)) + 
        geom_bar(aes(fill = Taxon),stat = "identity",
                 position = position_stack(reverse = TRUE)) + 
        ylab(label = "Relative abundance of dominant phyla") + 
        xlab(label = "") + 
        scale_fill_manual(values = cbbPalette,name = "Phylum") +
        facet_grid(.~Group,scales = "free",space = "free") + 
        theme_bw()+ 
        theme(panel.grid=element_blank()) + 
        theme(panel.border = element_blank()) +
        theme(panel.background=element_rect(fill='transparent', color='black'),
              plot.margin = unit(c(3,5,1,1),"mm")) + 
        theme(axis.text.x=element_blank()) + 
        theme(axis.text.y=element_text(colour = "black",size = 10)) + 
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"),
              axis.ticks.x = element_blank())+ 
        theme(axis.title.y = element_text(size = 12,face = "bold",
                                          margin = unit(c(0,1,0,1),"lines"))) + 
        scale_y_continuous(limits = c(0,100.001),expand = c(0,0)) + 
        theme(legend.text = element_text(colour = "black",size = 12)) + 
        theme(legend.title = element_text(size = 14,colour = "black",face = "bold"),
              strip.text = element_text(size = 14,colour = "black",face = "bold"))
    
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
    if (nrow(f.abundance.1) > 10) {
        f.abundance.1 <- f.abundance.1[1:10,]
        f.abundance.1 <- t(f.abundance.1)
        sum2 <- apply(f.abundance.1,1,sum) 
        Others <- 100.00001-sum2
        f.abundance.1 <- as.data.frame(cbind(f.abundance.1,Others))
    }
    f.abundance.1$Taxon <- rownames(f.abundance.1)
    taxon <- melt(f.abundance.1)
    if (nrow(f.abundance) > 10) {
        colnames(taxon) <- c("variable","Taxon","value")
        }
    taxon <- merge(taxon,group)
    aa <- f.abundance[,order(f.abundance[1,],decreasing = T)]
    taxon$variable <- factor(taxon$variable,levels = colnames(aa))
    
    p.Class <- ggplot(data = taxon,aes(x = variable, y = value)) + 
        geom_bar(aes(fill = Taxon),stat = "identity",
                 position = position_stack(reverse = TRUE)) + 
        ylab(label = "Relative abundance of dominant classes") + 
        xlab(label = "") + 
        scale_fill_manual(values = cbbPalette,name = "Class") +
        facet_grid(.~Group,scales = "free",space = "free") + 
        theme_bw()+ 
        theme(panel.grid=element_blank()) + 
        theme(panel.border = element_blank()) +
        theme(panel.background=element_rect(fill='transparent', color='black'),
              plot.margin = unit(c(3,5,1,1),"mm")) + 
        theme(axis.text.x=element_blank()) + 
        theme(axis.text.y=element_text(colour = "black",size = 10)) + 
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"),
              axis.ticks.x = element_blank())+ 
        theme(axis.title.y = element_text(size = 12,face = "bold",
                                          margin = unit(c(0,1,0,1),"lines"))) + 
        scale_y_continuous(limits = c(0,100.001),expand = c(0,0)) + 
        theme(legend.text = element_text(colour = "black",size = 12)) + 
        theme(legend.title = element_text(size = 14,colour = "black",face = "bold"),
              strip.text = element_text(size = 14,colour = "black",face = "bold"))
    
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
    if (nrow(f.abundance.1) > 10) {
        f.abundance.1 <- f.abundance.1[1:10,]
        f.abundance.1 <- t(f.abundance.1)
        sum2 <- apply(f.abundance.1,1,sum) 
        Others <- 100.00001-sum2
        f.abundance.1 <- as.data.frame(cbind(f.abundance.1,Others))
    }
    f.abundance.1$Taxon <- rownames(f.abundance.1)
    taxon <- melt(f.abundance.1)
    if (nrow(f.abundance) > 10) {
        colnames(taxon) <- c("variable","Taxon","value")
        }
    taxon <- merge(taxon,group)
    aa <- f.abundance[,order(f.abundance[1,],decreasing = T)]
    taxon$variable <- factor(taxon$variable,levels = colnames(aa))
    
    p.Order <- ggplot(data = taxon,aes(x = variable, y = value)) + 
        geom_bar(aes(fill = Taxon),stat = "identity",
                 position = position_stack(reverse = TRUE)) + 
        ylab(label = "Relative abundance of dominant orders") + 
        xlab(label = "") + 
        scale_fill_manual(values = cbbPalette,name = "Order") +
        facet_grid(.~Group,scales = "free",space = "free") + 
        theme_bw()+ 
        theme(panel.grid=element_blank()) + 
        theme(panel.border = element_blank()) +
        theme(panel.background=element_rect(fill='transparent', color='black'),
              plot.margin = unit(c(3,5,1,1),"mm")) + 
        theme(axis.text.x=element_blank()) + 
        theme(axis.text.y=element_text(colour = "black",size = 10)) + 
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"),
              axis.ticks.x = element_blank())+ 
        theme(axis.title.y = element_text(size = 12,face = "bold",
                                          margin = unit(c(0,1,0,1),"lines"))) + 
        scale_y_continuous(limits = c(0,100.001),expand = c(0,0)) + 
        theme(legend.text = element_text(colour = "black",size = 12)) + 
        theme(legend.title = element_text(size = 14,colour = "black",face = "bold"),
              strip.text = element_text(size = 14,colour = "black",face = "bold"))
    
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
    if (nrow(f.abundance.1) > 10) {
        f.abundance.1 <- f.abundance.1[1:10,]
        f.abundance.1 <- t(f.abundance.1)
        sum2 <- apply(f.abundance.1,1,sum) 
        Others <- 100.00001-sum2
        f.abundance.1 <- as.data.frame(cbind(f.abundance.1,Others))
    }
    f.abundance.1$Taxon <- rownames(f.abundance.1)
    taxon <- melt(f.abundance.1)
    if (nrow(f.abundance) > 10) {
        colnames(taxon) <- c("variable","Taxon","value")
        }
    taxon <- merge(taxon,group)
    aa <- f.abundance[,order(f.abundance[1,],decreasing = T)]
    taxon$variable <- factor(taxon$variable,levels = colnames(aa))
    
    p.Family <- ggplot(data = taxon,aes(x = variable, y = value)) + 
        geom_bar(aes(fill = Taxon),stat = "identity",
                 position = position_stack(reverse = TRUE)) + 
        ylab(label = "Relative abundance of dominant families") + 
        xlab(label = "") + 
        scale_fill_manual(values = cbbPalette,name = "Family") +
        facet_grid(.~Group,scales = "free",space = "free") + 
        theme_bw()+ 
        theme(panel.grid=element_blank()) + 
        theme(panel.border = element_blank()) +
        theme(panel.background=element_rect(fill='transparent', color='black'),
              plot.margin = unit(c(3,5,1,1),"mm")) + 
        theme(axis.text.x=element_blank()) + 
        theme(axis.text.y=element_text(colour = "black",size = 10)) + 
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"),
              axis.ticks.x = element_blank())+ 
        theme(axis.title.y = element_text(size = 12,face = "bold",
                                          margin = unit(c(0,1,0,1),"lines"))) + 
        scale_y_continuous(limits = c(0,100.001),expand = c(0,0)) + 
        theme(legend.text = element_text(colour = "black",size = 12)) + 
        theme(legend.title = element_text(size = 14,colour = "black",face = "bold"),
              strip.text = element_text(size = 14,colour = "black",face = "bold"))
    
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
    if (nrow(f.abundance.1) > 10) {
        f.abundance.1 <- f.abundance.1[1:10,]
        f.abundance.1 <- t(f.abundance.1)
        sum2 <- apply(f.abundance.1,1,sum) 
        Others <- 100.00001-sum2
        f.abundance.1 <- as.data.frame(cbind(f.abundance.1,Others))
    }
    f.abundance.1$Taxon <- rownames(f.abundance.1)
    taxon <- melt(f.abundance.1)
    if (nrow(f.abundance) > 10) {
        colnames(taxon) <- c("variable","Taxon","value")
        }
    taxon <- merge(taxon,group)
    aa <- f.abundance[,order(f.abundance[1,],decreasing = T)]
    taxon$variable <- factor(taxon$variable,levels = colnames(aa))
    
    p.Genus <- ggplot(data = taxon,aes(x = variable, y = value)) + 
        geom_bar(aes(fill = Taxon),stat = "identity",
                 position = position_stack(reverse = TRUE)) + 
        ylab(label = "Relative abundance of dominant genera") + 
        xlab(label = "") + 
        scale_fill_manual(values = cbbPalette,name = "Genus") +
        facet_grid(.~Group,scales = "free",space = "free") + 
        theme_bw()+ 
        theme(panel.grid=element_blank()) + 
        theme(panel.border = element_blank()) +
        theme(panel.background=element_rect(fill='transparent', color='black'),
              plot.margin = unit(c(3,5,1,1),"mm")) + 
        theme(axis.text.x=element_blank()) + 
        theme(axis.text.y=element_text(colour = "black",size = 10)) + 
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"),
              axis.ticks.x = element_blank())+ 
        theme(axis.title.y = element_text(size = 12,face = "bold",
                                          margin = unit(c(0,1,0,1),"lines"))) + 
        scale_y_continuous(limits = c(0,100.001),expand = c(0,0)) + 
        theme(legend.text = element_text(colour = "black",size = 12)) + 
        theme(legend.title = element_text(size = 14,colour = "black",face = "bold"),
              strip.text = element_text(size = 14,colour = "black",face = "bold"))
    
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
    if (nrow(f.abundance.1) > 10) {
        f.abundance.1 <- f.abundance.1[1:10,]
        f.abundance.1 <- t(f.abundance.1)
        sum2 <- apply(f.abundance.1,1,sum) 
        Others <- 100.00001-sum2
        f.abundance.1 <- as.data.frame(cbind(f.abundance.1,Others))
    }
    f.abundance.1$Taxon <- rownames(f.abundance.1)
    taxon <- melt(f.abundance.1)
    if (nrow(f.abundance) > 10) {
        colnames(taxon) <- c("variable","Taxon","value")
        }
    taxon <- merge(taxon,group)
    aa <- f.abundance[,order(f.abundance[1,],decreasing = T)]
    taxon$variable <- factor(taxon$variable,levels = colnames(aa))
    
    p.Species <- ggplot(data = taxon,aes(x = variable, y = value)) + 
        geom_bar(aes(fill = Taxon),stat = "identity",
                 position = position_stack(reverse = TRUE)) + 
        ylab(label = "Relative abundance of dominant species") + 
        xlab(label = "") + 
        scale_fill_manual(values = cbbPalette,name = "Species") +
        facet_grid(.~Group,scales = "free",space = "free") + 
        theme_bw()+ 
        theme(panel.grid=element_blank()) + 
        theme(panel.border = element_blank()) +
        theme(panel.background=element_rect(fill='transparent', color='black'),
              plot.margin = unit(c(3,5,1,1),"mm")) + 
        theme(axis.text.x=element_blank()) + 
        theme(axis.text.y=element_text(colour = "black",size = 10)) + 
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"),
              axis.ticks.x = element_blank())+ 
        theme(axis.title.y = element_text(size = 12,face = "bold",
                                          margin = unit(c(0,1,0,1),"lines"))) + 
        scale_y_continuous(limits = c(0,100.001),expand = c(0,0)) + 
        theme(legend.text = element_text(colour = "black",size = 12)) + 
        theme(legend.title = element_text(size = 14,colour = "black",face = "bold"),
              strip.text = element_text(size = 14,colour = "black",face = "bold"))

    result <- list(p.Phylum,p.Class,p.Order,p.Family,p.Genus,p.Species)
    return(result)
}