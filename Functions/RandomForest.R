Random.forest.biomarker <- function(otu,group){
    f <- function(x)sum(x==0)
    a <- apply(otu[,2:ncol(otu)],1,f)
    a <- as.data.frame(a)
    otu1 <- cbind(otu[,2:ncol(otu)],a)
    rownames(otu1) <- otu[,1]
    otu2 <- filter(otu1,a < (ncol(otu1)-1)*0.4)
    otu2 <- otu2[,1:(ncol(otu2)-1)]
    otu2 <- otu2[rowMeans(otu2) >= 0.1,]
    if (nrow(otu2) == 1) {
        otu2 <- otu1[,1:(ncol(otu1)-1)]
    }
    data1 <- t(otu2)
    aa <- match(rownames(data1),group$variable)
    group <- group[aa,]
    data1 <- data.frame(data1,Group = group$Group)
    colnames(data1) <- c(rownames(otu2),"Group")
    data1$Group <- factor(data1$Group)
    
    phy.rf <- randomForest(data1[,1:(ncol(data1)-1)],data1$Group,ntree = 1000,
                           importance = TRUE,proximity = TRUE)
    
    p1 <- pheatmap(phy.rf$confusion[,1:length(unique(data1$Group))],cluster_rows = FALSE,
             cluster_cols = FALSE,legend = FALSE,display_numbers = TRUE,
             number_format = "%i",fontsize = 10,number_color = "black",
             main = paste("RF classifier (",(1-round(phy.rf$err.rate[1000,1],4))*100,"% correct)",sep = ""),
             fontsize_number = 26)
    
    imp <- as.data.frame(round(importance(phy.rf), 4))
    imp <- imp[order(imp$MeanDecreaseAccuracy,decreasing = TRUE),]
    imp1 <- imp[1:phy.rf$mtry,]
    if (nrow(imp1) > 0) {
        imp1$ID <- rownames(imp1)
        
        p2<- imp1 %>%
            mutate(ID = fct_reorder(ID,MeanDecreaseAccuracy)) %>%
            ggplot(aes(ID,MeanDecreaseAccuracy)) +
            geom_errorbar(aes(ymin = MeanDecreaseAccuracy - MeanDecreaseGini,
                              ymax = MeanDecreaseAccuracy + MeanDecreaseGini),width = 0.3) + 
            geom_vline(xintercept = 1:nrow(imp),linetype = 3) +
            geom_point(shape = 21,fill = "grey80",size = 5) +
            coord_flip() +
            theme_bw()+ theme(panel.grid=element_blank()) + 
            theme(panel.border = element_blank()) +
            theme(panel.background=element_rect(fill='transparent', color='black'),
                  plot.margin = unit(c(3,5,1,1),"mm")) +
            scale_y_continuous(limits = c(ifelse(min(imp1$MeanDecreaseAccuracy - imp1$MeanDecreaseGini) > 0,
                                                 min((imp1$MeanDecreaseAccuracy - imp1$MeanDecreaseGini)*0.9),
                                                 min((imp1$MeanDecreaseAccuracy - imp1$MeanDecreaseGini)*1.1)),
                                          max((imp1$MeanDecreaseAccuracy + imp1$MeanDecreaseGini)*1.1)),
                               expand = c(0,0)) +
            theme(axis.text.y = element_text(size = 10,face = "bold.italic",
                                             colour = "black"),
                  axis.text.x = element_text(size = 12,face = "bold",colour = "black"),
                  axis.title.x = element_text(size = 16,face = "bold",
                                              colour = "black")) +
            xlab("") +ylab("Importance")
        
        abundance <- otu2[rownames(otu2) %in% rownames(imp1),]
        abundance$ID <- rownames(abundance)
        abundance1 <- melt(abundance)
        colnames(group)[1] <- "variable"
        abundance1 <- merge(abundance1,group)
        abundance1$ID <- factor(abundance1$ID,levels = rev(imp1$ID))
        p3 <- ggplot(abundance1) +
            geom_tile(aes(variable,ID,fill = value)) +
            theme_classic() +
            theme(axis.ticks = element_blank(),axis.line = element_blank()) +
            xlab("") + ylab("") +
            facet_grid(.~Group,scales = "free",space = "free") + 
            scale_fill_gradient2('Relative abundance (%)',low = 'black',high = 'red') +
            theme(panel.background = element_rect(color = "black"),
                  legend.position = "bottom",
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  strip.text.x = element_text(colour = "black",face = "bold",size = 12))
        p3 <- p2 + p3 + plot_layout(widths = c(1,2))
    }
    
    if (nrow(imp1) > 0) {
        result <- list(p1,imp1,p3,abundance)
    }else{
        result <- list(p1,imp1)
    }
    
    return(result)
}