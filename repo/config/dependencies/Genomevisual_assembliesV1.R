library(reshape2)

lista <- commandArgs(trailingOnly =  TRUE)
samplename <- as.character(lista[1])
indir <- lista[2] #directory in which to find inputs, and to ut outputs
outdir <- lista[3]
annvcf <- paste(indir,"/", samplename,"_anno.vcf", sep = "") ##XXX.annotated variants.vcf

#now lwt's call the variants, we will have to structure them before making a track out of them/draft them in the report
preopen1 <- read.csv(annvcf, sep = "\t", header=F,row.names=NULL)
skipcoord1 <- which(grepl("#CHROM",preopen1$V1)==TRUE)-1
vcfR_test <- read.csv(annvcf, sep = "\t", skip= skipcoord1, header=T,row.names=NULL)
add_info <- colsplit(vcfR_test$INFO, "\\|", c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", 
                                              "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", 
                                              "Rank", "HGVS.c", "HGVS.p", "cDNA.pos_cDNA.length", "CDS.pos_CDS.length",
                                              "AA.pos","AA.length", "Distance", "Notes"))
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
outdf <- paste(outdir,"/", samplename,".variant.regions.annotated.primary.alleles.report.tab", sep ="")
write.table(final_df,outdf,quote=F, row.names=F, col.names=T, sep = "\t")

