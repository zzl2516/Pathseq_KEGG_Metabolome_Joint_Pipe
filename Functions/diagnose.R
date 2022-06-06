diagnose <- function(data,abundance.biomarker,group){
    rownames(abundance.biomarker) <- abundance.biomarker[,1]
    abundance.biomarker <- abundance.biomarker[,-1]
    abundance.biomarker <- abundance.biomarker[rownames(abundance.biomarker) %in% data$V1,]
    abundance.1 <- as.data.frame(t(abundance.biomarker))
    abundance.1$variable <- rownames(abundance.1)
    abundance.1 <- merge(abundance.1,group)
    abundance.1 <- melt(abundance.1)
    abundance.1 <- abundance.1[,-1] %>%
        group_by(Group,variable) %>%
        summarise(value = mean(value))
    biomarker <- data.frame(V1 = rownames(abundance.biomarker),
                            V2 = rep(NA,nrow(abundance.biomarker)))
    for (i in rownames(abundance.biomarker)) {
        a <- abundance.1[abundance.1$variable == i,]
        if (a[a$Group == levels(group$Group)[1],"value"] > 
            a[a$Group == levels(group$Group)[2],"value"]) {
            biomarker[biomarker$V1 == i,"V2"] <- levels(group$Group)[1]
        }else{
            biomarker[biomarker$V1 == i,"V2"] <- levels(group$Group)[2]
        }
    }
    colnames(biomarker) <- c("biomarker","Group")
    abundance.biomarker$biomarker <- rownames(abundance.biomarker)
    abundance.biomarker <- merge(abundance.biomarker,biomarker)
    abundance.2 <- melt(abundance.biomarker)
    abundance.2 <- abundance.2 %>%
        group_by(biomarker,Group,variable) %>%
        summarise(value = sum(value))
    index <- abundance.2 %>%
        group_by(Group,variable) %>%
        summarise(value = sum(value))
    index <- spread(index,variable,value)
    index <- t(index)
    if (ncol(index) == 1) {
        a <- c(ifelse(index[1,1] == unique(group$Group)[1],unique(group$Group)[2],
                      unique(group$Group)[1]),rep(0,(nrow(index)-1)))
        index.1 <- cbind(index,a)
    }else{
        index.1 <- index
    }
    colnames(index.1) <- index.1[1,]
    index.1 <- as.data.frame(index.1[-1,])
    index.1 <- index.1[,unique(group$Group)]
    for (i in 1:2) {
        index.1[,i] <- as.numeric(index.1[,i])
    }
    index.2 <- data.frame(variable = rownames(index.1),
                          index = index.1[,1] - index.1[,2])
    index.3 <- merge(index.2,group)
    
    alpha.test <- index.3[,2:3]
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
    
    p1 <- ggplot(alpha.test,aes(Group,Num,color = Group)) + 
        geom_boxplot(width = 0.6,outlier.color = "transparent") +
        geom_jitter(width = 0.3,size = 1.5,alpha = 0.5) +
        geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
                  size = 5,color = "black",fontface = "bold") +
        labs(y = "Disease index",x = "") +
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
    
    roc1 <- roc(index.3$Group,index.3$index)
    auc1 <- auc(roc1)
    
    result <- list(index.2,p1,roc1)
    return(result)
}