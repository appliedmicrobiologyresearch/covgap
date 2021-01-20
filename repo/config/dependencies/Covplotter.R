#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Sushi")
library("Sushi")
#setwd("~/Desktop/bedgrapph/")
args <- commandArgs(trailingOnly = T)
file <- args[1]
out <- args[2]
name <- args[3]
init <-read.csv(file, sep = "\t", header = F, row.names = NULL, stringsAsFactors = F)
plottable <- init#[,c(1,2,3,5)]
colnames(plottable) <- c("chrom","start","end","value")
head(plottable)
chrom = "MN908947.3"
chromstart = min(plottable$start)
chromend = max(plottable$end)
#new_plot <- paste(name,".pdf", sep = "")
pdf(out, width = 20, height = 10)
range <- c(min(plottable$value),max(plottable$value))
plotBedgraph(plottable,chrom, chromstart,chromend,
               transparency=.50, colorbycol= SushiColors(5))
labelgenome(chrom,chromstart,chromend,n=20,scale="Kb")
mtext("Read Depth",side=12,line=1.75,cex=1,font=2)
mtext(name, outer=TRUE,  cex=3, line=-3)
axis(side=2,las=2,tcl=.2)
dev.off()
