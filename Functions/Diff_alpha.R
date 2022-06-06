diff.alpha <- function(alpha1,group){
    colnames(group) <- c("ID","Group")
    alpha <- merge(alpha1,group)
    alpha.test <- alpha[,c(2,9)]
    colnames(alpha.test) <- c("Num","Group")
    alpha.test$Group <- factor(alpha.test$Group)
    x <- c("a","b")
    y <- c("a","a")
    if (length(levels(alpha.test$Group)) == 2) {
        fit1 <- t.test(Num~Group,data = alpha.test)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = if(fit1$p.value < 0.05){
                               x
                           }else{
                               y
                           },
                           value.y = dd$Max*1.03)
    }else{
        fit1 <- aov(Num~Group,data = alpha.test)
        tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
        res1 <- cld(tuk1,alpah=0.05)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = res1$mcletters$Letters,
                           value.y = dd$Max*1.03)
    }
    
    alpha.test1 <- ggplot(alpha.test,aes(Group,Num,color = Group)) + 
        geom_boxplot(width = 0.6,outlier.color = "transparent") +
        geom_jitter(width = 0.3,size = 1.5,alpha = 0.5) +
        geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
                  size = 5,color = "black",fontface = "bold") +
        labs(y = "Observed species",x = "") +
        scale_color_manual(values = cbbPalette) +
        theme_bw()+
        theme(panel.grid=element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_blank(),
              axis.title.y = element_text(face = "bold",color = "black",size = 16),
              axis.text.y=element_text(colour='black',size=10),
              axis.text.x=element_text(colour = "black",size = 10,face = "bold",
                                       angle = 45,hjust = 1,vjust = 1),
              legend.position = "none")
    alpha.test <- alpha[,c(3,9)]
    colnames(alpha.test) <- c("Num","Group")
    alpha.test$Group <- factor(alpha.test$Group)
    x <- c("a","b")
    y <- c("a","a")
    if (length(levels(alpha.test$Group)) == 2) {
        fit1 <- t.test(Num~Group,data = alpha.test)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = if(fit1$p.value < 0.05){
                               x
                           }else{
                               y
                           },
                           value.y = dd$Max*1.03)
    }else{
        fit1 <- aov(Num~Group,data = alpha.test)
        tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
        res1 <- cld(tuk1,alpah=0.05)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = res1$mcletters$Letters,
                           value.y = dd$Max*1.03)
    }
    
    alpha.test2 <- ggplot(alpha.test,aes(Group,Num,color = Group)) + 
        geom_boxplot(width = 0.6,outlier.color = "transparent") +
        geom_jitter(width = 0.3,size = 1.5,alpha = 0.5) +
        geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
                  size = 5,color = "black",fontface = "bold") +
        labs(y = "Chao1 index",x = "") +
        scale_color_manual(values = cbbPalette) +
        theme_bw()+
        theme(panel.grid=element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_blank(),
              axis.title.y = element_text(face = "bold",color = "black",size = 16),
              axis.text.y=element_text(colour='black',size=10),
              axis.text.x=element_text(colour = "black",size = 10,face = "bold",
                                       angle = 45,hjust = 1,vjust = 1),
              legend.position = "none")
    alpha.test <- alpha[,c(4,9)]
    colnames(alpha.test) <- c("Num","Group")
    alpha.test$Group <- factor(alpha.test$Group)
    x <- c("a","b")
    y <- c("a","a")
    if (length(levels(alpha.test$Group)) == 2) {
        fit1 <- t.test(Num~Group,data = alpha.test)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = if(fit1$p.value < 0.05){
                               x
                           }else{
                               y
                           },
                           value.y = dd$Max*1.03)
    }else{
        fit1 <- aov(Num~Group,data = alpha.test)
        tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
        res1 <- cld(tuk1,alpah=0.05)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = res1$mcletters$Letters,
                           value.y = dd$Max*1.03)
    }
    
    alpha.test3 <- ggplot(alpha.test,aes(Group,Num,color = Group)) + 
        geom_boxplot(width = 0.6,outlier.color = "transparent") +
        geom_jitter(width = 0.3,size = 1.5,alpha = 0.5) +
        geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
                  size = 5,color = "black",fontface = "bold") +
        labs(y = "ACE index",x = "") +
        scale_color_manual(values = cbbPalette) +
        theme_bw()+
        theme(panel.grid=element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_blank(),
              axis.title.y = element_text(face = "bold",color = "black",size = 16),
              axis.text.y=element_text(colour='black',size=10),
              axis.text.x=element_text(colour = "black",size = 10,face = "bold",
                                       angle = 45,hjust = 1,vjust = 1),
              legend.position = "none")
    alpha.test <- alpha[,c(5,9)]
    colnames(alpha.test) <- c("Num","Group")
    alpha.test$Group <- factor(alpha.test$Group)
    x <- c("a","b")
    y <- c("a","a")
    if (length(levels(alpha.test$Group)) == 2) {
        fit1 <- t.test(Num~Group,data = alpha.test)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = if(fit1$p.value < 0.05){
                               x
                           }else{
                               y
                           },
                           value.y = dd$Max*1.005)
    }else{
        fit1 <- aov(Num~Group,data = alpha.test)
        tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
        res1 <- cld(tuk1,alpah=0.05)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = res1$mcletters$Letters,
                           value.y = dd$Max*1.005)
    }
    
    alpha.test4 <- ggplot(alpha.test,aes(Group,Num,color = Group)) + 
        geom_boxplot(width = 0.6,outlier.color = "transparent") +
        geom_jitter(width = 0.3,size = 1.5,alpha = 0.5) +
        geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
                  size = 5,color = "black",fontface = "bold") +
        labs(y = "Shannon index",x = "") +
        scale_color_manual(values = cbbPalette) +
        theme_bw()+
        theme(panel.grid=element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_blank(),
              axis.title.y = element_text(face = "bold",color = "black",size = 16),
              axis.text.y=element_text(colour='black',size=10),
              axis.text.x=element_text(colour = "black",size = 10,face = "bold",
                                       angle = 45,hjust = 1,vjust = 1),
              legend.position = "none")
    alpha.test <- alpha[,c(6,9)]
    colnames(alpha.test) <- c("Num","Group")
    alpha.test$Group <- factor(alpha.test$Group)
    x <- c("a","b")
    y <- c("a","a")
    if (length(levels(alpha.test$Group)) == 2) {
        fit1 <- t.test(Num~Group,data = alpha.test)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = if(fit1$p.value < 0.05){
                               x
                           }else{
                               y
                           },
                           value.y = dd$Max*1.001)
    }else{
        fit1 <- aov(Num~Group,data = alpha.test)
        tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
        res1 <- cld(tuk1,alpah=0.05)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = res1$mcletters$Letters,
                           value.y = dd$Max*1.001)
    }
    
    alpha.test5 <- ggplot(alpha.test,aes(Group,Num,color = Group)) + 
        geom_boxplot(width = 0.6,outlier.color = "transparent") +
        geom_jitter(width = 0.3,size = 1.5,alpha = 0.5) +
        geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
                  size = 5,color = "black",fontface = "bold") +
        labs(y = "Simpson index",x = "") +
        scale_color_manual(values = cbbPalette) +
        theme_bw()+
        theme(panel.grid=element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_blank(),
              axis.title.y = element_text(face = "bold",color = "black",size = 16),
              axis.text.y=element_text(colour='black',size=10),
              axis.text.x=element_text(colour = "black",size = 10,face = "bold",
                                       angle = 45,hjust = 1,vjust = 1),
              legend.position = "none")
    alpha.test <- alpha[,c(7,9)]
    colnames(alpha.test) <- c("Num","Group")
    alpha.test$Group <- factor(alpha.test$Group)
    x <- c("a","b")
    y <- c("a","a")
    if (length(levels(alpha.test$Group)) == 2) {
        fit1 <- t.test(Num~Group,data = alpha.test)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = if(fit1$p.value < 0.05){
                               x
                           }else{
                               y
                           },
                           value.y = dd$Max*1.005)
    }else{
        fit1 <- aov(Num~Group,data = alpha.test)
        tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
        res1 <- cld(tuk1,alpah=0.05)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = res1$mcletters$Letters,
                           value.y = dd$Max*1.005)
    }
    
    alpha.test6 <- ggplot(alpha.test,aes(Group,Num,color = Group)) + 
        geom_boxplot(width = 0.6,outlier.color = "transparent") +
        geom_jitter(width = 0.3,size = 1.5,alpha = 0.5) +
        geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
                  size = 5,color = "black",fontface = "bold") +
        labs(y = "Pielou_J index",x = "") +
        scale_color_manual(values = cbbPalette) +
        theme_bw()+
        theme(panel.grid=element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_blank(),
              axis.title.y = element_text(face = "bold",color = "black",size = 16),
              axis.text.y=element_text(colour='black',size=10),
              axis.text.x=element_text(colour = "black",size = 10,face = "bold",
                                       angle = 45,hjust = 1,vjust = 1),
              legend.position = "none")
    alpha.test <- alpha[,c(8,9)]
    colnames(alpha.test) <- c("Num","Group")
    alpha.test$Group <- factor(alpha.test$Group)
    x <- c("a","b")
    y <- c("a","a")
    if (length(levels(alpha.test$Group)) == 2) {
        fit1 <- t.test(Num~Group,data = alpha.test)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = if(fit1$p.value < 0.05){
                               x
                           }else{
                               y
                           },
                           value.y = dd$Max*1.001)
    }else{
        fit1 <- aov(Num~Group,data = alpha.test)
        tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
        res1 <- cld(tuk1,alpah=0.05)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = res1$mcletters$Letters,
                           value.y = dd$Max*1.001)
    }
    
    alpha.test7 <- ggplot(alpha.test,aes(Group,Num,color = Group)) + 
        geom_boxplot(width = 0.6,outlier.color = "transparent") +
        geom_jitter(width = 0.3,size = 1.5,alpha = 0.5) +
        geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
                  size = 5,color = "black",fontface = "bold") +
        labs(y = "Good_coverage",x = "") +
        scale_color_manual(values = cbbPalette) +
        theme_bw()+
        theme(panel.grid=element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_blank(),
              axis.title.y = element_text(face = "bold",color = "black",size = 16),
              axis.text.y=element_text(colour='black',size=10),
              axis.text.x=element_text(colour = "black",size = 10,face = "bold",
                                       angle = 45,hjust = 1,vjust = 1),
              legend.position = "none")
    result <- list(alpha.test1,alpha.test2,alpha.test3,alpha.test4,alpha.test5,
                   alpha.test6,alpha.test7)
    return(result)
}