library(stringi)
library(seqinr)
args <- commandArgs(trailingOnly=TRUE)
file <- args[1]
output <- args[2] 
chr <- read.fasta(file, as.string=T)
chromosome <- chr[[1]][1]
# find all position in one chromosome
n.pos <- stri_locate_all_fixed(chromosome,"n")[[1]][,1]
# group consecutive Ns in elements of a list
n.pos.list <- split(n.pos, cumsum(c(1, diff(n.pos) != 1)))
# extract for each group of Ns the length (thus the length of the gap)
n.length <-sapply(n.pos.list,length)
# construct results
n.pos.res <- cbind(do.call(rbind,lapply(n.pos.list,function(x){return(c(min(x),max(x)))})),length=n.length)
# order by length
n.pos.res <- n.pos.res[order(n.pos.res[,2],decreasing = F),]
# name column
colnames(n.pos.res)<-c("start","end","length")

write.table(n.pos.res,output, quote=F, sep = "\t", col.names=T, row.names=F, append=TRUE)
