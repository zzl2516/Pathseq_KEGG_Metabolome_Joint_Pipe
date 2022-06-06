tax.summary <- function(data){
    Phylum <- aggregate(data[,1:(ncol(data)-6)],
                             list(data$Phylum),sum)
    colnames(Phylum)[1] <- "Phylum"
    Class <- aggregate(data[,1:(ncol(data)-6)],
                        list(data$Class),sum)
    colnames(Class)[1] <- "Class"
    Order <- aggregate(data[,1:(ncol(data)-6)],
                        list(data$Order),sum)
    colnames(Order)[1] <- "Order"
    Family <- aggregate(data[,1:(ncol(data)-6)],
                        list(data$Family),sum)
    colnames(Family)[1] <- "Family"
    Genus <- aggregate(data[,1:(ncol(data)-6)],
                        list(data$Genus),sum)
    colnames(Genus)[1] <- "Genus"
    Species <- aggregate(data[,1:(ncol(data)-6)],
                        list(data$Species),sum)
    colnames(Species)[1] <- "Species"
    result <- list(Phylum,Class,Order,Family,Genus,Species)
    return(result)
}