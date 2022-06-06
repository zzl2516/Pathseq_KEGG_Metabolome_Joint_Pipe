lefse.1 <- function(lefse,group,p.lvl = 0.05,lda.lvl = 2,seed = 11,adjust.p = F){
    lefse$Phylum <- paste("p__",lefse$Phylum,sep = "")
    lefse$Class <- paste("c__",lefse$Class,sep = "")
    lefse$Order <- paste("o__",lefse$Order,sep = "")
    lefse$Family <- paste("f__",lefse$Family,sep = "")
    lefse$Genus <- paste("g__",lefse$Genus,sep = "")
    phylum <- aggregate(lefse[,1:(ncol(lefse)-6)],list(lefse$Phylum),sum)
    class <- aggregate(lefse[,1:(ncol(lefse)-6)],list(lefse$Class),sum)
    order <- aggregate(lefse[,1:(ncol(lefse)-6)],list(lefse$Order),sum)
    family <- aggregate(lefse[,1:(ncol(lefse)-6)],list(lefse$Family),sum)
    genus <- aggregate(lefse[,1:(ncol(lefse)-6)],list(lefse$Genus),sum)
    lefse <- rbind(phylum,class,order,family,genus)
    colnames(lefse)[1] <- "variable"
    rownames(lefse) <- lefse[,1]
    lefse <- lefse[,-1]
    rownames(lefse) <- gsub("-","_",rownames(lefse))
    lefse <- as.data.frame(t(lefse))
    lefse$variable <- rownames(lefse)
    lefse <- merge(lefse,group)
    lefse <- lefse[,c("Group",colnames(lefse)[1:(ncol(lefse)-1)])]
    lefse <- as.data.frame(t(lefse))
    lefse.test <- lefse[3:nrow(lefse),]
    for (i in 1:ncol(lefse.test)) {
        lefse.test[,i] <- as.numeric(lefse.test[,i])
    }
    otu <- lefse.test
    f <- function(x)sum(x==0)
    a <- apply(otu,1,f)
    otu1 <- cbind(otu,a)
    otu2 <- filter(otu1,a < ncol(otu1)*0.4)
    lefse.test <- otu2[,1:(ncol(otu2)-1)]
    lefse.test <- lefse.test[rowMeans(lefse.test) > 0.1,]
    colnames(lefse.test) <- lefse[2,]
    rawpvalues <- apply(lefse.test,1, 
                        function(x) kruskal.test(x,as.vector(t(lefse[1,])))$p.value)
    ord.inx <- order(rawpvalues)
    rawpvalues <- rawpvalues[ord.inx]
    clapvalues <- p.adjust(rawpvalues, method ="fdr")
    lefse.test <- lefse.test[ord.inx,]
    lefse.test <- data.frame(t(lefse.test),Group = as.vector(t(lefse[1,])))
    ldares <- lda(Group~ .,data = lefse.test)
    class_no <- length(unique(group$Group))
    ldamean <- as.data.frame(t(ldares$means))
    ldamean$max <- apply(ldamean[,1:class_no],1,max)
    ldamean$min <- apply(ldamean[,1:class_no],1,min)
    ldamean$LDAscore <- as.vector(abs(ldamean$max-ldamean$min)/2)
    dd <- as.data.frame(ldamean$LDAscore)
    dd1 <- preProcess(dd,method = c("range"),rangeBounds = c(1,1000000))
    dd2 <- predict(dd1,dd)
    ldamean$LDAscore <- signif(log10(dd2[,1]),digits = 3)
    colnames(ldamean)[5] <- "LDAscore"
    ldamean$Pvalues <- signif(rawpvalues,digits=5)
    ldamean$FDR <- signif(clapvalues,digits=5)
    
    a = rep("A",length(ldamean$max))
    for (i in 1:length(ldamean$max)) {
        name =colnames(ldamean[,1:class_no])
        a[i] = name[ldamean[,1:class_no][i,] %in% ldamean$max[i]]
    }
    ldamean$class = a
    rownames(ldamean) <- gsub("\\.", ' ', rownames(ldamean))
    rownames(ldamean) <- gsub("sp", 'sp.', rownames(ldamean))
    rownames(ldamean) <- gsub("`", '', rownames(ldamean))
    ldamean$Pvalues <- signif(rawpvalues[match(row.names(ldamean),names(rawpvalues))],digits=5)
    ldamean$FDR <- signif(clapvalues,digits=5)
    resTable <- ldamean
    resTable <- resTable[complete.cases(resTable$Pvalues),]
    
    if (adjust.p) {
        de.Num <- sum(resTable$FDR <= p.lvl & resTable$LDAscore>=lda.lvl)
    } else {
        de.Num <- sum(resTable$Pvalues <= p.lvl & resTable$LDAscore >= lda.lvl)
    }
    
    if(de.Num == 0){
        current.msg <<- "No significant features were identified with given criteria.";
    }else{
        current.msg <<- paste("A total of", de.Num, "significant features with given criteria.")
    }
    print(current.msg)
    # sort by p value
    ord.inx <- order(resTable$Pvalues, resTable$LDAscore)
    resTable <- resTable[ord.inx, ,drop=FALSE]
    resTable <- resTable[,c(ncol(resTable),1:(ncol(resTable)-1))]
    resTable <- resTable[,c(ncol(resTable),1:(ncol(resTable)-1))]
    
    # resTable %>% tail()
    if (adjust.p) {
        taxtree = resTable[resTable$FDR <= p.lvl & resTable$LDAscore >= lda.lvl,]
    } else {
        taxtree = resTable[resTable$Pvalues <= p.lvl & resTable$LDAscore >= lda.lvl,]
    }
    
    #-提取所需要的颜色
    if (de.Num > 0) {
        colour = c('darkgreen','red',"blue","#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
        selececol = colour[1:length(levels(as.factor(taxtree$class)))]
        names(selececol) = levels(as.factor(taxtree$class))
        A = rep("a",length(row.names(taxtree)))
        
        for (i in 1:length(row.names(taxtree))) {
            A[i] = selececol [taxtree$class[i]]
        }
        
        taxtree$color = A
    }
    return(taxtree)
}