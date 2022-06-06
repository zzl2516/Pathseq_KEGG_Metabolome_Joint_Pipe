pcoa.community <- function(ARG_sub,group,species.type){
    colnames(group) <- c("sample","Group")
    ARG_beta <- ARG_sub[group[,1],]
    ARG_dis <- vegdist(ARG_beta)
    
    pcoa<- pcoa(ARG_dis, correction = "none", rn = NULL)
    PC1 = pcoa$vectors[,1]
    PC2 = pcoa$vectors[,2]
    plotdata <- data.frame(rownames(pcoa$vectors),PC1,PC2)
    colnames(plotdata) <-c("sample","PC1","PC2")
    plotdata <- merge(plotdata,group)
    pc1 <-floor(pcoa$values$Relative_eig[1]*100)
    pc2 <-floor(pcoa$values$Relative_eig[2]*100)
    
    ARG_pcoa1<-ggplot(plotdata, aes(PC1, PC2)) +
        geom_point(aes(fill=Group),size=4.5,color = "black",shape = 21,alpha = 0.8)+
        geom_vline(aes(xintercept = 0),linetype="dotted")+
        geom_hline(aes(yintercept = 0),linetype="dotted")+
        scale_fill_manual(values=cbbPalette)+
        labs(title=paste("PCoA - ",species.type," communities",sep = )) + 
        xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
        ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
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
        labs(title=paste("PCoA - ",species.type," communities",sep = )) + 
        xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
        ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
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
                         label.padding = unit(0.2,"lines"),size = 2) +
        geom_vline(aes(xintercept = 0),linetype="dotted")+
        geom_hline(aes(yintercept = 0),linetype="dotted")+
        scale_fill_manual(values=cbbPalette)+
        labs(title=paste("PCoA - ",species.type," communities",sep = )) + 
        xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
        ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
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
    
    
    result <- list(ARG_dis,plotdata,ARG_pcoa1,ARG_pcoa2,ARG_pcoa3)
    return(result)
}