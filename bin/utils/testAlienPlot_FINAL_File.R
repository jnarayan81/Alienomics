#!/usr/bin/env Rscript

# 20-06-2016
# works with R 2.15.2 and ggplot 0.9.3.1
# Check ggplot2 help forums or contact Jitendra Narayan jnarayan81@gmail.com if something doesn't run
# because of updated programs/packages

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)<= 2) {
  # default output file
  args[3] = "finalPlot.pdf"
}

setwd(args[2]) #The file working directory location
pdf (file=paste (Sys.time(), ".pdf", sep=""))
plot (rnorm (100))
dev.off()

Alien <- read.table( args[1], sep="\t", header=TRUE)
Alien <- Alien[,c("spsClass","len","alienCnt","gc","globalScore")]

#Remove all rows with any NAs to make this simple
Alien <- data.frame(na.omit(Alien))

#png(file=args[2],width=950,height=600)

library("ggplot2")
library("ggrepel")
#ggplot(Alien, aes(GC, Score)) + geom_point(aes(label=Len, shape = Name, colour = Len), alpha = 0.4)
ggplot(Alien, aes(gc, globalScore)) + geom_point(aes(colour = spsClass), alpha = 0.4, size=2) #+ geom_text_repel(aes(label=Len), size = 3)

dev.off()


