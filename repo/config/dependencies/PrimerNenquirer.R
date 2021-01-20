library(stringi)
library(seqinr)
library(stringr)
bed <- read.csv("~/Desktop/V3.primers.bed", stringsAsFactors = F, sep = "\t", row.names = NULL, header=F)
head(bed)
args <- commandArgs(trailingOnly = T)
samplename <- args[1]
dir <- args[2]
consfile <- paste(dir,"/", samplename,".consensus.fasta", sep="\t")
colnames(bed) <- c("chromosome","start","end","ID","version","strand")
consensus <- read.fasta(consfile, as.string=F)#read.fasta("~/Desktop/42180271.consensus.fasta", as.string=F)
head(consensus)
seq <- consensus$MN908947.3
width <- c()
sequence_in_consensus <- c()
Ns_count <- c()
Ns_percentage <- c()
nrbed <- dim(bed)[1]
for (i in 1:nrbed){
list1 <- seq[bed$start[i]:bed$end[i]]
seq1 <- toString(list1)
seq2 <- gsub(", ","",seq1)
sequence_in_consensus <- rbind(sequence_in_consensus,seq2)
}
head(bed)

bed$sequence_in_consensus <- sequence_in_consensus
bed$Ns_count <- str_count(bed$sequence_in_consensus, "n")
bed$width <- bed$end-bed$start+1
bed$Npercent <- bed$Ns_count/bed$width*100
bed$Sample<- samplename
head(bed)
bedN <- bed[which(bed$Npercent>0),]
bedname <- paste(samplename,"_primer_enquiry.tab", sep ="")
bedname1 <- paste(samplename,"_primer_enquiryN.tab", sep ="")
write.table(bed,bedname, quote = F, col.names = T, row.names = F, sep = "\t")
write.table(bedN,bedname1, quote = F, col.names = T, row.names = F, sep = "\t")




