#!/usr/bin/env Rscript
args = commandArgs(TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<=0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==2) {
  k <- 3 # default cluster value
} else if (length(args)==3) {
  k <- as.numeric(commandArgs(TRUE)[3])
}

setwd(args[2]) #The file working directory location
library("Ckmeans.1d.dp")
resTable <- read.table( args[1], sep="\t", header=TRUE)
name <- as.vector(resTable$globalScore) # Currently 13 is the score column -- globalScore is the header name
name
result <- Ckmeans.1d.dp(name, k)
kmeanCluster <- as.numeric(result$cluster)
tmpVal <- read.table( args[1], sep="\t", header=TRUE)
tmpVal <- cbind(tmpVal, kmeanCluster)

#Create a new filename by concatinating clustered.csv in old file args[0]
newFilename <- paste(args[1], "clustered.csv", sep = "_")
write.table(tmpVal, file = newFilename, append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

plot(result, main="Optimal univariate k-median given cluster value")
plot(name, col = result$cluster, pch = result$cluster, cex = 1,  main="Optimal univariate k-means clustering given cluster value", sub=paste("Number of clusters given:", k))
abline(h = result$centers, col = 1:2, lty="dashed", lwd=1, )
legend("bottomleft", paste("Cluster", 1:k), col=1:k, pch=1:k,
cex=1, bty="n")


####------------------------------------------------------------------
# Plot globalScore versus GC plot

library("ggplot2")
library("ggrepel")
Alien <- read.table( newFilename, sep="\t", header=TRUE)
Alien <- Alien[,c("spsClass","len","alienCnt","gc","globalScore")]

#Remove all rows with any NAs to make this simple
#Alien <- data.frame(na.omit(Alien))

ggplot(Alien, aes(gc, globalScore)) + geom_point(aes(colour = spsClass), alpha = 0.4, size=2)
#+ geom_text_repel(aes(label=Len), size = 3)
#theme(legend.position = "bottom")

dev.off()
