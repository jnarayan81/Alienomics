####------------------------------------------------------------------
#Rscript plotGCvsCpv.R final_results.out

# Plot globalScore versus GC plot
#!/usr/bin/env Rscript
args = commandArgs(TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<=0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 


library("ggplot2")
library("ggrepel")
Alien <- read.table( args[1], sep="\t", header=TRUE)
Alien <- Alien[,c("spsClass","len","alienCnt","gc","globalScore")]

#Remove all rows with any NAs to make this simple
#Alien <- data.frame(na.omit(Alien))

ggplot(Alien, aes(gc, globalScore)) + geom_point(aes(colour = spsClass, size = len), alpha = 0.4)
#+ geom_text_repel(aes(label=Len), size = 3)
#theme(legend.position = "bottom")

dev.off()
