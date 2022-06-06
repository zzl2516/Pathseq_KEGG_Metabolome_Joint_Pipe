nmds.community <- function(ARG_sub,group){
    colnames(group) <- c("sample","Group")
    pcoa <- metaMDS(ARG_sub)
    stress <- pcoa$stress
    nmds_scores <- scores(pcoa, choices = c(1,2))
    nmds_scores <- as.data.frame(nmds_scores)
    nmds_scores$sample <- rownames(nmds_scores)
    nmds_scores <- nmds_scores[,c("sample","NMDS1","NMDS2")]
    colnames(nmds_scores) <-c("sample","PC1","PC2")
    plotdata <- merge(nmds_scores,group)
    
    ARG_pcoa1<-ggplot(plotdata, aes(PC1, PC2)) +
        geom_point(aes(fill=Group),size=4.5,color = "black",shape = 21,alpha = 0.8)+
        geom_vline(aes(xintercept = 0),linetype="dotted")+
        geom_hline(aes(yintercept = 0),linetype="dotted")+
        scale_fill_manual(values=cbbPalette)+
        labs(title="NMDS - Metabolites") +
        xlab("NMDS1") + 
        ylab("NMDS2")+
        theme(text=element_text(size=18))+
        theme(panel.background = element_rect(fill='white', colour='black'),
              panel.grid=element_blank(), 
              axis.title = element_text(color='black',size=18),
              axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_text(colour='black', size=18),
              axis.title.y=element_text(colour='black', size=18),
              axis.text=element_text(colour='black',size=16),
              legend.title=element_text(size = 14,face = "bold"),
              legend.text=element_text(size=12),
              legend.key=element_blank(),legend.position = "right",
              legend.background = element_rect(colour = "black"))+
        theme(plot.title = element_text(size=20,colour = "black",hjust = 0.5,face = "bold"))
    
    ARG_pcoa2<-ggplot(plotdata, aes(PC1, PC2)) +
        geom_point(aes(fill=Group),size=4.5,color = "black",shape = 21,alpha = 0.8)+
        stat_ellipse(aes(fill = Group),geom = "polygon",level = 0.95,alpha = 0.3)+
        geom_vline(aes(xintercept = 0),linetype="dotted")+
        geom_hline(aes(yintercept = 0),linetype="dotted")+
        scale_fill_manual(values=cbbPalette)+
        labs(title="NMDS - Metabolites") +
        xlab("NMDS1") + 
        ylab("NMDS2")+
        theme(text=element_text(size=18))+
        theme(panel.background = element_rect(fill='white', colour='black'),
              panel.grid=element_blank(), 
              axis.title = element_text(color='black',size=18),
              axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_text(colour='black', size=18),
              axis.title.y=element_text(colour='black', size=18),
              axis.text=element_text(colour='black',size=16),
              legend.title=element_text(size = 14,face = "bold"),
              legend.text=element_text(size=12),
              legend.key=element_blank(),legend.position = "right",
              legend.background = element_rect(colour = "black"))+
        theme(plot.title = element_text(size=20,colour = "black",hjust = 0.5,face = "bold"))
    
    ARG_pcoa3<-ggplot(plotdata, aes(PC1, PC2)) +
        geom_point(aes(fill=Group),size=4.5,color = "black",shape = 21,alpha = 0.8)+
        geom_label_repel(aes(PC1,PC2,label = sample),fill = "white",color = "black",
                         box.padding = unit(0.3,"lines"),segment.colour = "grey50",
                         label.padding = unit(0.15,"lines"),size = 2) +
        geom_vline(aes(xintercept = 0),linetype="dotted")+
        geom_hline(aes(yintercept = 0),linetype="dotted")+
        scale_fill_manual(values=cbbPalette)+
        labs(title="NMDS - Metabolites") + 
        xlab("NMDS1") + 
        ylab("NMDS2")+
        theme(text=element_text(size=18))+
        theme(panel.background = element_rect(fill='white', colour='black'),
              panel.grid=element_blank(), 
              axis.title = element_text(color='black',size=18),
              axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_text(colour='black', size=18),
              axis.title.y=element_text(colour='black', size=18),
              axis.text=element_text(colour='black',size=16),
              legend.title=element_text(size = 14,face = "bold"),
              legend.text=element_text(size=12),
              legend.key=element_blank(),legend.position = "right",
              legend.background = element_rect(colour = "black"))+
        theme(plot.title = element_text(size=20,colour = "black",hjust = 0.5,face = "bold"))
    
    colnames(plotdata) <- c("sample","NMDS1","NMDS2","Group")
    result <- list(stress,plotdata,ARG_pcoa1,ARG_pcoa2,ARG_pcoa3)
    return(result)
}