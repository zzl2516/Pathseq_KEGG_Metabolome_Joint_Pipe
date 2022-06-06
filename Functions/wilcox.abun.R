wilcox.abun <- function(otu,otu.wilcox.biomarker,group){
    data1 <- t(otu[,2:ncol(otu)])
    rownames(data1) <- colnames(otu)[2:ncol(otu)]
    aa <- match(rownames(data1),group$variable)
    group <- group[aa,]
    data1 <- data.frame(data1,Group = group$Group)
    colnames(data1) <- c(otu[,1],"Group")
    data1$Group <- factor(data1$Group)
    abun.bar <- data1[,c(otu.wilcox.biomarker$V1,"Group")] %>% 
        gather(variable,value,-Group) %>% 
        group_by(variable,Group) %>% 
        summarise(Mean = mean(value))
    
    diff1 <- data1[,c(otu.wilcox.biomarker$V1,"Group")] %>% 
        select_if(is.numeric) %>%
        map_df(~ broom::tidy(t.test(. ~ Group,data = data1)), .id = 'var')

    diff.mean <- diff1[,c("var","estimate","conf.low","conf.high")]
    aa <- match(diff.mean$var,otu.wilcox.biomarker$V1)
    otu.wilcox.biomarker <- otu.wilcox.biomarker[aa,]
    diff.mean <- data.frame(diff.mean,p.value = otu.wilcox.biomarker$p.value)
    diff.mean$Group <- c(ifelse(diff.mean$estimate >0,levels(data1$Group)[1],
                                levels(data1$Group)[2]))
    diff.mean <- diff.mean[order(diff.mean$estimate,decreasing = TRUE),]
    
    abun.bar$variable <- factor(abun.bar$variable,levels = rev(diff.mean$var))
    p1 <- ggplot(abun.bar,aes(variable,Mean,fill = Group)) +
        scale_x_discrete(limits = levels(diff.mean$var)) +
        coord_flip() +
        xlab("") +
        ylab("Mean proportion") +
        theme(panel.background = element_rect(fill = 'transparent'),
              panel.grid = element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"),
              axis.title.x=element_text(colour='black', size=12,face = "bold"),
              axis.text=element_text(colour='black',size=10,face = "bold"),
              legend.title=element_blank(),
              legend.text=element_text(size=12,face = "bold",colour = "black"),
              legend.position = "top",
              legend.direction = "horizontal",
              legend.key.width = unit(0.8,"cm"),
              legend.key.height = unit(0.5,"cm"))
    
    
    for (i in 1:(nrow(diff.mean) - 1)) 
        p1 <- p1 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, 
                            fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
    
    p1 <- p1 + 
        geom_bar(stat = "identity",position = "dodge",width = 0.7,colour = "black") +
        scale_fill_manual(values=cbbPalette)
    
    
    diff.mean$var <- factor(diff.mean$var,levels = levels(abun.bar$variable))
    diff.mean$p.value <- as.character(diff.mean$p.value)
    p2 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
        theme(panel.background = element_rect(fill = 'transparent'),
              panel.grid = element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"),
              axis.title.x=element_text(colour='black', size=12,face = "bold"),
              axis.text=element_text(colour='black',size=10,face = "bold"),
              axis.text.y = element_blank(),
              legend.position = "none",
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              plot.title = element_text(size = 15,face = "bold",colour = "black",hjust = 0.5)) +
        scale_x_discrete(limits = levels(diff.mean$var)) +
        coord_flip() +
        xlab("") +
        ylab("Difference in mean proportions") +
        labs(title="95% confidence intervals") 
    
    for (i in 1:(nrow(diff.mean) - 1)) 
        p2 <- p2 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, 
                            fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
    
    p2 <- p2 +
        geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                      position = position_dodge(0.8), width = 0.5, size = 0.5) +
        geom_point(shape = 21,size = 3) +
        scale_fill_manual(values = if(length(unique(diff.mean$Group)) == 1 &
                                      unique(diff.mean$Group) == levels(group$Group)[2]){
            cbbPalette[2]
        }else{
            cbbPalette
        }) +
        geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black')
    
    
    p3 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
        geom_text(aes(y = 0,x = var),label = substr(diff.mean$p.value,1,6),
                  hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
        geom_text(aes(x = nrow(diff.mean)/2 +0.5,y = 0.85),label = "P-value",
                  srt = 90,fontface = "bold",size = 3) +
        coord_flip() +
        ylim(c(0,1)) +
        theme(panel.background = element_blank(),
              panel.grid = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank())
    
    ## 图像拼接
    p <- p1 + p2 + p3 + plot_layout(widths = c(4,6,2))
    result <- list(p,abun.bar,diff.mean)
    return(result)
}