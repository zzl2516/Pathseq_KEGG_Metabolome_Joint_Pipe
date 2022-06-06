lefse_bar = function(taxtree){
    taxtree$ID = row.names(taxtree)
    taxtree <- taxtree %>%
        arrange(class,LDAscore)
    taxtree$ID = factor(taxtree$ID,levels=taxtree$ID)
    taxtree$class = factor(taxtree$class,levels = unique(taxtree$class))
    taxtree$LDAscore <- ifelse(taxtree$class == levels(taxtree$class)[2],
                               taxtree$LDAscore,-taxtree$LDAscore)
    
    pbar <- ggplot(taxtree,aes(y =ID, x = LDAscore,fill = class)) + 
        geom_bar(stat = "identity") +
        xlab(label = "LDAscore (Log10)") + 
        ylab(label = "") + 
        scale_fill_manual(values = cbbPalette,name = "Group") +
        theme_bw()+ 
        theme(panel.grid=element_blank()) + 
        theme(panel.border = element_blank()) +
        theme(panel.background=element_rect(fill='transparent', color='black')) + 
        theme(axis.text.x=element_text(colour = "black",size = 10)) + 
        theme(axis.text.y=element_text(colour = "black",size = 10)) + 
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"))+ 
        theme(axis.title.x = element_text(size = 12,face = "bold")) + 
        theme(legend.text = element_text(colour = "black",size = 12)) + 
        theme(legend.title = element_text(size = 14,colour = "black",face = "bold"),
              legend.position = "top")
        
    return(pbar)
}