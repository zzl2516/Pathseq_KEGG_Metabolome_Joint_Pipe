shared.species.abun.s <- function(ARG_sub,core_ARG_id,group){
    ARG_sub <- as.data.frame(t(ARG_sub))
    ARG_sub <- t(t(ARG_sub)/colSums(ARG_sub)*100)
    ARG_sub <- as.data.frame(ARG_sub)
    ARG_sub$Species_ID <- rownames(ARG_sub)
    ARG_sub <- ARG_sub[,c("Species_ID",colnames(ARG_sub)[1:(ncol(ARG_sub)-1)])]
    ARG_sub1 <- melt(ARG_sub)
    colnames(group) <- c("variable","Group")
    colnames(ARG_sub1)[3] <- "Abun"
    core_ARG_abun <- ARG_sub1 %>%
        filter(Species_ID %in% core_ARG_id$Species_ID)
    core_ARG_ta <- core_ARG_abun %>%
        group_by(variable) %>%
        summarise(Abun = sum(Abun))
    core_ARG_ta <- merge(core_ARG_ta,group)

    core_ARG_ta1 <- ggplot(core_ARG_ta,aes(variable,Abun,fill = Group)) + 
        geom_bar(stat = "identity",position = "stack",width = 0.6) +
        labs(y = "Relative abundance of shared species",x = "") +
        scale_fill_manual(values = cbbPalette) +
        theme_bw()+
        theme(panel.grid=element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_blank(),
              axis.title.y = element_text(face = "bold",color = "black",size = 12),
              axis.text.y=element_text(colour='black',size=10),
              axis.text.x=element_text(colour = "black",size = 10,
                                       angle = 45,hjust = 1,vjust = 1),
              legend.position = "top",
              legend.title = element_text(face = "bold",color = "black",size = 12),
              legend.text = element_text(face = "bold",color = "black",size = 12))
    return(core_ARG_ta1)
}