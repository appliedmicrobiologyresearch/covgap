library(rtracklayer)
library(GenomicFeatures) 
library(GenomicRanges )
library(Gviz)
library(RMariaDB)
library(reshape2)

#inputs and outputs
lista <- commandArgs(trailingOnly =  TRUE)
samplename <- as.character(lista[1])
dir <- lista[2] #directory in which to find inputs, and to ut outputs
tag <- lista[3] #coverage tag, in the pipe is HighCov or LowCov

N_file <- paste(dir,"/",samplename,".NPosition.tab", sep = "") ##XXXX.NPosition.tab
annvcf <- paste(dir,"/", samplename,".annotated.variants.vcf", sep = "") ##XXX.annotated variants.vcf
sndannvcf <- paste(dir,"/",samplename,".annotated.minorityvariants.vcf", sep = "") ##XXX.annotated.minority variants.vcf
techannot <- paste(dir,"/",samplename,".minority_alleles_report.tech", sep = "") ##XXX.minority variants report.vcf --> includes data like vcf annotation and perc of support, stuff that couldn't be added in the vcf

#first let's call in the bed-like nfile --> output from Ner.R
Ns <- read.csv(N_file, sep = "\t", header = T, stringsAsFactors = F, skipNul=TRUE)
n <- data.frame(chromosome="NC_045512.2", start=Ns$start, end=Ns$end, width=Ns$length, strand="+", feature= "N", 
                gene="_", exon="_", transcript="_", symbol="_")

#now lwt's call the variants, we will have to structure them before making a track out of them/draft them in the report
preopen1 <- read.csv(annvcf, sep = "\t", header=F,row.names=NULL)
skipcoord1 <- which(grepl("#CHROM",preopen1$V1)==TRUE)-1
vcfR_test <- read.csv(annvcf, sep = "\t", skip= skipcoord1, header=T,row.names=NULL, skipNul=TRUE)
add_info <- colsplit(vcfR_test$INFO, "\\|", c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", 
                                              "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", 
                                              "Rank", "HGVS.c", "HGVS.p", "cDNA.pos_cDNA.length", "CDS.pos_CDS.length",
                                              "AA.pos","AA.length", "Distance", "Notes"))
preopen2 <- read.csv(sndannvcf, sep = "\t", header=F,row.names=NULL)
skipcoord2 <- which(grepl("#CHROM",preopen2$V1)==TRUE)-1
vcfR2_test <- read.csv(sndannvcf, sep = "\t", skip=skipcoord2, header=T,row.names=NULL, skipNul=TRUE)
add_info2 <- colsplit(vcfR2_test$INFO, "\\|", c("Allele", "Annotation", "Annotation_Impact", "Gene_Name",
                                              "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType",
                                              "Rank", "HGVS.c", "HGVS.p", "cDNA.pos_cDNA.length", "CDS.pos_CDS.length",
                                              "AA.pos","AA.length", "Distance", "Notes"))
techinfo <- read.csv(techannot, sep = "\t", header=T, row.names=NULL, skipNul=TRUE)
techadds <- techinfo[,-c(1,2,3)]
#Read in all the ranges -binding pockets of farmaceuticals- useful both for the plot and for the report
#Bcells epitopes
epi <- read.csv("/scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/ranges_epitopes.txt", sep = "\t", header = T, row.names = NULL,stringsAsFactors = T )
wid <- epi$end-epi$start
#in good shape
epitope <- data.frame(chromosome="NC_045512.2", start=epi$start, end=epi$end, width=wid, strand="+", feature= "B_cell_epitope",
                      gene="_", exon="_", transcript="_", symbol="_", stringsAsFactors = F)

#Lopinavir
lopi <- read.csv("/scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/Mproranges.txt", sep = "\t", header = T, row.names = NULL,stringsAsFactors = T )
lwid <- lopi$end-lopi$start
#in good shape
lopinavir <- data.frame(chromosome="NC_045512.2", start=lopi$start, end=lopi$end, width=lwid, strand="+", feature= "_",
                      gene="_", exon="_", transcript=lopi$residues, symbol="_", stringsAsFactors = F)
#Chloroquin
chlor <- read.csv("/scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/ranges_chloroquin.txt", sep = "\t", header = T, row.names = NULL,stringsAsFactors = T )
cwid <- chlor$end-chlor$start
#in good shape
chloroquin <- data.frame(chromosome="NC_045512.2", start=chlor$start, end=chlor$end, width=cwid, strand="+", feature= "Chloroquin-ACE2_binding_site",
                      gene="_", exon="_", transcript="_", symbol="_", stringsAsFactors = F)

#Remdesivir
remde <- read.csv("/scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/ranges_remdesivir1.txt", sep = "\t", header = T, row.names = NULL,stringsAsFactors = T )
rwid <- remde$end-remde$start
#in good shape
remdesivir <- data.frame(chromosome="NC_045512.2", start=remde$start, end=remde$end, width=rwid, strand="+", feature= "Remdesivir binding site",
                        gene="_", exon="_", transcript=remde$residues, symbol="_", stringsAsFactors = F)


#####DRAFTING THE REPORTS#####

#dominant alternative alles
final_df <- vcfR_test[,c(1,2,3,4,5,6)]
final_df$GENE <- add_info$Gene_Name
final_df$ANNOTATION <- add_info$Annotation
final_df$AA_INVOLVED <- add_info$HGVS.p
final_df$SAMPLEID <- samplename
final_df$PCR <- "No"
final_df$PCR[which(final_df$POS>=22364 & final_df$POS<=22475)] <- "Yes"
final_df$PCR[which(final_df$POS>=28119 & final_df$POS<=28215)] <- "Yes"
#more outputs for the pharmaceuticals

final_df$Mpro_pockets <- "No"
final_df$Mpro_pockets[which(final_df$POS>=lopi$start[1] & final_df$POS<=lopi$end[1])] <- "Yes_Pocket_Lopinavir1"
final_df$Mpro_pockets[which(final_df$POS>=lopi$start[2] & final_df$POS<=lopi$end[2])] <- "Yes_Pocket_Lopinavir1_Thr24"
final_df$Mpro_pockets[which(final_df$POS>=lopi$start[3] & final_df$POS<=lopi$end[3])] <- "Yes_Pocket_Lopinavir1_Thr26"
final_df$Mpro_pockets[which(final_df$POS>=lopi$start[4] & final_df$POS<=lopi$end[4])] <- "Yes_Pocket_Lopinavir2"
final_df$Mpro_pockets[which(final_df$POS>=lopi$start[5] & final_df$POS<=lopi$end[5])] <- "Yes_Pocket_Lopinavir2_Asn119"
final_df$Mpro_pockets[which(final_df$POS>=lopi$start[6] & final_df$POS<=lopi$end[6])] <- "Yes_PocketGC373"
final_df$Mpro_pockets[which(final_df$POS>=lopi$start[7] & final_df$POS<=lopi$end[7])] <- "Yes_PocketGC373_Cys145"
final_df$Mpro_pockets[which(final_df$POS>=lopi$start[8] & final_df$POS<=lopi$end[8])] <- "Yes_PocketGC373_Gly143"
final_df$Mpro_pockets[which(final_df$POS>=lopi$start[9] & final_df$POS<=lopi$end[9])] <- "Yes_PocketGC373_Ser144"
final_df$Mpro_pockets[which(final_df$POS>=lopi$start[10] & final_df$POS<=lopi$end[10])] <- "Yes_PocketGC373_His163"
final_df$Mpro_pockets[which(final_df$POS>=lopi$start[11] & final_df$POS<=lopi$end[11])] <- "Yes_PocketGC373_Glu166"
final_df$Mpro_pockets[which(final_df$POS>=lopi$start[12] & final_df$POS<=lopi$end[12])] <- "Yes_Virulence_Ala285"

final_df$Mpro <- "No"
final_df$Mpro[which(final_df$POS>=lopi$start[13] & final_df$POS<=lopi$end[13])] <- "Yes"

final_df$Chloroquine_binding_affecting <- "No"
final_df$Chloroquine_binding_affecting[which(final_df$POS>=chlor$start[1] & final_df$POS<=chlor$end[1])] <- "Yes"

final_df$Remdesivir_binding_affecting <- "No"
final_df$Remdesivir_binding_affecting[which(final_df$POS>=remde$start[1] & final_df$POS<=remde$end[1])] <- "Yes_pocket_F480"
final_df$Remdesivir_binding_affecting[which(final_df$POS>=remde$start[2] & final_df$POS<=remde$end[2])] <- "Yes_pocket_codon_F480"
final_df$Remdesivir_binding_affecting[which(final_df$POS>=remde$start[3] & final_df$POS<=remde$end[3])] <- "Yes_pocket_V557"
final_df$Remdesivir_binding_affecting[which(final_df$POS>=remde$start[4] & final_df$POS<=remde$end[4])] <- "Yes_pocket_codon_V557"
final_df$FiveFU_binding_affecting <- "No"
final_df$FiveFU_binding_affecting[which(final_df$POS>=remde$start[5] & final_df$POS<=remde$end[5])] <- "Yes_pocket_M615"
final_df$FiveFU_binding_affecting[which(final_df$POS>=remde$start[6] & final_df$POS<=remde$end[6])] <- "Yes_pocket_codon_M615"
final_df$RdRp_NSP12 <- "No"
final_df$RdRp_NSP12[which(final_df$POS>=remde$start[7] & final_df$POS<=remde$end[7])] <- "Yes"
#now we output it, as a improved version of primer_flagged version
outdf <- paste(dir,"/", samplename,".variant.regions.annotated.primary.alleles.report.",tag,".tab", sep ="")
write.table(final_df,outdf,quote=F, row.names=F, col.names=T, sep = "\t")
if(dim(vcfR2_test)[1]!=0){
#secondary alternative alleles
final_df2 <- vcfR2_test[,c(1,2,3,4,5,6)]
final_df2$GENE <- add_info2$Gene_Name
final_df2$ANNOTATION <- add_info2$Annotation
final_df2$AA_INVOLVED <- add_info2$HGVS.p
final_df2$SAMPLEID <- samplename
final_df2 <- cbind(final_df2,techadds)
final_df2$Label <- "None"
final_df2$AlleleSupport_AF <- as.numeric(as.character(final_df2$AlleleSupport_AF))
final_df2$Label[which(final_df2$AlleleSupport_AF < 0.70)] <- "secondary_allele"
final_df2$Label[which(final_df2$AlleleSupport_AF > 0.70 & final_df2$AlleleSupport_AF < 0.99)] <- "alternative_allele" 
final_df2$PCR <- "No"
final_df2$PCR[which(final_df2$POS>=22364 & final_df2$POS<=22475)] <- "Yes"
final_df2$PCR[which(final_df2$POS>=28119 & final_df2$POS<=28215)] <- "Yes"

final_df2$PCR[which(final_df2$POS>=22364 & final_df2$POS<=22475)] <- "Yes"
final_df2$PCR[which(final_df2$POS>=28119 & final_df2$POS<=28215)] <- "Yes"

#now re outputting the same things for the secondary alleles if any
final_df2$Mpro_pockets <- "No"
final_df2$Mpro_pockets[which(final_df2$POS>=lopi$start[1] & final_df2$POS<=lopi$end[1])] <- "Yes_Pocket_Lopinavir1"
final_df2$Mpro_pockets[which(final_df2$POS>=lopi$start[2] & final_df2$POS<=lopi$end[2])] <- "Yes_Pocket_Lopinavir1_Thr24"
final_df2$Mpro_pockets[which(final_df2$POS>=lopi$start[3] & final_df2$POS<=lopi$end[3])] <- "Yes_Pocket_Lopinavir1_Thr26"
final_df2$Mpro_pockets[which(final_df2$POS>=lopi$start[4] & final_df2$POS<=lopi$end[4])] <- "Yes_Pocket_Lopinavir2"
final_df2$Mpro_pockets[which(final_df2$POS>=lopi$start[5] & final_df2$POS<=lopi$end[5])] <- "Yes_Pocket_Lopinavir2_Asn119"
final_df2$Mpro_pockets[which(final_df2$POS>=lopi$start[6] & final_df2$POS<=lopi$end[6])] <- "Yes_PocketGC373"
final_df2$Mpro_pockets[which(final_df2$POS>=lopi$start[7] & final_df2$POS<=lopi$end[7])] <- "Yes_PocketGC373_Cys145"
final_df2$Mpro_pockets[which(final_df2$POS>=lopi$start[8] & final_df2$POS<=lopi$end[8])] <- "Yes_PocketGC373_Gly143"
final_df2$Mpro_pockets[which(final_df2$POS>=lopi$start[9] & final_df2$POS<=lopi$end[9])] <- "Yes_PocketGC373_Ser144"
final_df2$Mpro_pockets[which(final_df2$POS>=lopi$start[10] & final_df2$POS<=lopi$end[10])] <- "Yes_PocketGC373_His163"
final_df2$Mpro_pockets[which(final_df2$POS>=lopi$start[11] & final_df2$POS<=lopi$end[11])] <- "Yes_PocketGC373_Glu166"
final_df2$Mpro_pockets[which(final_df2$POS>=lopi$start[12] & final_df2$POS<=lopi$end[12])] <- "Yes_Virulence_Ala285"

final_df2$Mpro <- "No"
final_df2$Mpro[which(final_df2$POS>=lopi$start[13] & final_df2$POS<=lopi$end[13])] <- "Yes"

final_df2$Chloroquine_binding_affecting <- "No"
final_df2$Chloroquine_binding_affecting[which(final_df2$POS>=chlor$start[1] & final_df2$POS<=chlor$end[1])] <- "Yes"

final_df2$Remdesivir_binding_affecting <- "No"
final_df2$Remdesivir_binding_affecting[which(final_df2$POS>=remde$start[1] & final_df2$POS<=remde$end[1])] <- "Yes_pocket_F480"
final_df2$Remdesivir_binding_affecting[which(final_df2$POS>=remde$start[2] & final_df2$POS<=remde$end[2])] <- "Yes_pocket_codon_F480"
final_df2$Remdesivir_binding_affecting[which(final_df2$POS>=remde$start[3] & final_df2$POS<=remde$end[3])] <- "Yes_pocket_V557"
final_df2$Remdesivir_binding_affecting[which(final_df2$POS>=remde$start[4] & final_df2$POS<=remde$end[4])] <- "Yes_pocket_codon_V557"
final_df2$FiveFU_binding_affecting <- "No"
final_df2$FiveFU_binding_affecting[which(final_df2$POS>=remde$start[5] & final_df2$POS<=remde$end[5])] <- "Yes_pocket_M615"
final_df2$FiveFU_binding_affecting[which(final_df2$POS>=remde$start[6] & final_df2$POS<=remde$end[6])] <- "Yes_pocket_codon_M615"
final_df2$RdRp_NSP12 <- "No"
final_df2$RdRp_NSP12[which(final_df2$POS>=remde$start[7] & final_df2$POS<=remde$end[7])] <- "Yes"
#now we output it, as a improved version of primer_flagged version
outdf2 <- paste(dir,"/", samplename,".variant.regions.annotated.minority.alleles.report.",tag,".tab", sep ="")
write.table(final_df2,outdf2,quote=F, row.names=F, col.names=T, sep = "\t")
}else{
print("No secondary alleles found, moving forward")
}
############DRAFTING THE PLOT#############
#-we need different tracks that will be later on stacked together, unfortunately, has to be done one by one
#NOTE: only the dominat alleles are ending up in the plot

#1:SNP(variant) formatting
snp <- data.frame(chromosome="NC_045512.2", start=final_df$POS, end=final_df$POS, width=1, strand="+", feature= final_df$PCR, 
                  gene=final_df$GENE, exon="_", transcript=final_df$ANNOTATION, symbol="_", stringsAsFactors = F)

#renaming some snps for better comprehension
snp$transcript <- gsub("upstream_gene_variant", "non coding variants", snp$transcript)
snp$transcript <- gsub("synonymous_variant", "variants Not causing amino acid change", snp$transcript)
snp$transcript <- gsub("missense_variant", "variants Causing amino acid change", snp$transcript)

#2:Reference track
po1<-import.gff("/scicore/home/egliadr/GROUP/projects/COVID-19/Illumina/Softwares/dependencies/GCA_009858895.3_ASM985889v3_genomic1.gff")
options(ucscChromosomeNames=FALSE)
gen<- (genome(po1))
chr <- as.character(unique(seqnames(po1)))
gtrack <- GeneRegionTrack(po1, genome = "wuhCor1", name= "Encoded Genes (NCBI-Reference)", transcriptAnnotation="gene", 
                          collapseTranscripts="longest", cex.group=1.2 , fontsize=16, background.panel="transparent", 
                          fill="deepskyblue3", shape="arrow") 
#3: we take all data.frames and have them in GenomeRanges stackable tracks
#3.1= the genome scale size
ntrack <- GenomeAxisTrack(background.title="white", fontsize=13, col="darkgrey",name = "Genome size", 
			showTitle=T, ticksAt=c(0,5000,10000,15000,20000,25000,30000))
#3.2= the variants
Snptrack <- GeneRegionTrack(snp, genome = "wuhCor1", just.group="below", col.line="transparent",
                          chromosome = chr, name = "Variants found (SNPs)", shape = "box", col="transparent",
                           fill = "black",background.panel="peachpuff", fontsize=16, cex.group=1.1,
				transcriptAnnotation="transcript", feature=as.vector(snp$feature), No="black", Yes="yellow")
#3.3= The Ns along the genomes, here we use the output from Ner.R
Ntrack <- GeneRegionTrack(n, genome = "wuhCor1", chromosome = chr, col.line="transparent", 
			name = "Regions not sequenced", background.panel="ghostwhite", cex=1.4,fill="black",
                          shape = "box",fontsize=13)
#3.5 The Bcell epitope range
epitrack<- GeneRegionTrack(epitope, genome = "wuhCor1", chromosome = chr, name = "Antigen binding sites - B Cell", 
				col.line="transparent",background.panel="transparent", cex=1.4,fill="black",
                           shape = "box",fontsize=13)
#3.6: Lopinavir binding sites
lopitrack<- GeneRegionTrack(lopinavir, genome = "wuhCor1", chromosome = chr, name = "Lopinavir binding sites",  
				col.line="transparent",background.panel="transparent", cex=1.4,fill="black",
                           shape = "box",fontsize=16,background.panel="ghostwhite")
#3.7: Chloroquine (ACE2) binding sites
chlorotrack<- GeneRegionTrack(chloroquin, genome = "wuhCor1", chromosome = chr, name = "Chloroquine -ACE2- binding site", 
				col.line="transparent",background.panel="transparent", cex=1.4,fill="black",
                           shape = "box",fontsize=13)
#3.8: Remdesivir binding sites
remdetrack<- GeneRegionTrack(remdesivir, genome = "wuhCor1", chromosome = chr, name = "Remdesivir binding sites",
                                col.line="transparent",background.panel="transparent", cex=1.4,fill="black",
                           shape = "box",fontsize=16,background.panel="ghostwhite")

#3.9: the downsampled coverage plot, to give a flavour of the coverage along the genome
bedfile <- paste(dir,"/",samplename,".rd.1000uptrim.bedgraph", sep = "")
bgFile <- read.csv(bedfile, sep = "\t", stringsAsFactors = F, header=F)
colnames(bgFile) <- c("chromosome", "start", "end", "coverage")
Covtrack <- DataTrack(range = bgFile, genome = "wuhCor1", type = "histogram", fontsize=16, 
                     chromosome = chr, name = "Coverage", window=-1, col.histogram="cornflowerblue", windowSize = 50)

#Plots, one with coverage and one without.
pdf_out <- paste(dir,"/", samplename,".complete.trackplot.",tag,".pdf", sep = "") # output
graphname <- paste("Sample number: ",samplename, sep = "")
pdf(pdf_out, height = 12, width=11)
plotTracks(list(Covtrack,Snptrack,Ntrack,epitrack,lopitrack,chlorotrack,remdetrack,gtrack,ntrack), background.title = "darkgrey", col="grey", sizes=c(3,1.5,1,1,1.5,1,1.5,2,1), panel.only = F, main=graphname)
dev.off()

pdf_out1 <- paste(dir,"/", samplename,".nocoverage.trackplot.",tag,".pdf", sep = "") # output
graphname1 <- paste("Sample number: ",samplename, sep = "")
pdf(pdf_out1, height = 10, width=11)
plotTracks(list(Snptrack,Ntrack,epitrack,lopitrack,chlorotrack,remdetrack,gtrack,ntrack), background.title = "darkgrey", col="grey", sizes=c(1.5,1,1,1.5,1,1.5,2,1), panel.only = F, main=graphname1)
dev.off()
