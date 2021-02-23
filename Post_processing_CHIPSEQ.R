

library(IRanges)
library(GenomicRanges)
require(ChIPpeakAnno)
library(stringi)

source("/home/nostromo/data/pepe/ZCL/library_functions_ZCL_chip.R")
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesVik.R")
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesEli.R")
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesAlba.R")


MACS<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/MACS2_EGR1_peaks.narrowPeak")
BAYES<-read.table("/home/nostromo/data/pepe/peak_calling_EGR1/bayes_output1.txt")
SICER<-read.table("/home/nostromo/data/pepe/peak_calling_EGR1/EGR1_chip_T-W200-G200-FDR0.01-island.bed")
#ZCL<-read.table("/home/nostromo/data/pepe/peak_calling_EGR1/ZCL_peaks.final.bed")
ZCL<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/ZCL_peaks_g9_a1.2_k0_modified.bed")
JAMM<-read.table("/home/nostromo/data/pepe/peak_calling_EGR1/selected_JAMM_EGR1.narrowPeak")
RANGER<-read.table("/home/nostromo/data/pepe/peak_calling_EGR1/RANGER_EGR1.bed_region.bed")

MACS_bed<-MACS[,c("V1","V2","V3")]
SICER_bed<-SICER[,c("V1","V2","V3")]
JAMM_bed<-JAMM[,c("V1","V2","V3")]
BAYES_bed<-BAYES[,c("space","start","end")]
ZCL_bed<-ZCL[,c("V1","V2","V3")]
RANGER_bed<-RANGER[,c("V1","V2","V3")]

MACS_bed$diff<-(MACS_bed$V3 - MACS_bed$V2)
SICER_bed$diff<-(SICER_bed$V3 - SICER_bed$V2)
JAMM_bed$diff<-(JAMM_bed$V3 - JAMM_bed$V2)
BAYES_bed$diff<-(BAYES_bed$end - BAYES_bed$start)
ZCL_bed$diff<-(ZCL_bed$V3 - ZCL_bed$V2)
RANGER_bed$diff<-(RANGER_bed$V3 - RANGER_bed$V2)

MACS_bed_lessthan500<-MACS_bed[(MACS_bed$diff < 500),]
MACS_bed_lessthan500$diff<-NULL

SICER_bed_lessthan500<-SICER_bed[(SICER_bed$diff < 500),]
SICER_bed_lessthan500$diff<-NULL

JAMM_bed_lessthan500<-JAMM_bed[(JAMM_bed$diff < 500),]
JAMM_bed_lessthan500$diff<-NULL

BAYES_bed_lessthan500<-BAYES_bed[(BAYES_bed$diff < 500),]
BAYES_bed_lessthan500$diff<-NULL

RANGER_bed_lessthan500<-RANGER_bed[(RANGER_bed$diff < 500),]
RANGER_bed_lessthan500$diff<-NULL

ZCL_bed_lessthan500<-ZCL_bed[(ZCL_bed$diff < 500),]
ZCL_bed_lessthan500$diff<-NULL

genecodev25_all <- read.table(file ="/home/nostromo/data/pepe/Gencodev25SinExones.gtf", skip = 5, header = FALSE, sep = "\t")

genecodev25 <- genecodev25_all[genecodev25_all$V3 == "gene", ]
genecodev25_tmp <- apply(as.data.frame(genecodev25[,9]), 1, parseENCODE)
genecodev25_tmp <- t(genecodev25_tmp)

colnames(genecodev25_tmp) <- c("gene_id", "gene_type", "gene_status", "gene_name", "level", "havana_gene")
genecodev25_Annot <- data.frame(genecodev25[,1:8], genecodev25_tmp[, c(1, 2,3, 4, 5,6)])

genecodev25_Annot_Gene <- genecodev25_Annot

G25AnnotChIP <- unique(genecodev25_Annot[genecodev25_Annot[, "V3"] == "gene", c("V1", "V4", "V5", "V7", "gene_id")])
colnames(G25AnnotChIP) <- c("chr", "start", "end", "strand", "gene_id")
G25AnnotChIP_RL <- RangedData(IRanges(start = G25AnnotChIP$start, end = G25AnnotChIP$end, names = paste(G25AnnotChIP$gene_id)), space = paste(G25AnnotChIP$chr), strand = paste(G25AnnotChIP$strand))

MACS_bed<-MACS_bed_lessthan500
SICER_bed<-SICER_bed_lessthan500
JAMM_bed<-JAMM_bed_lessthan500
BAYES_bed<-BAYES_bed_lessthan500
ZCL_bed<-ZCL_bed_lessthan500
RANGER_bed<-RANGER_bed_lessthan500

fileName_chip<-basename("/home/nostromo/data/pepe/peak_calling_EGR1/MACS2_EGR1_peaks.narrowPeak")
bed<-MACS_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)
bed_data_annotated_unique<-unique(stri_split_fixed(paste(bed_data_annotated$gene_name), ".", simplify = TRUE)[,1])
out_bed<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
list_ID_unique<-unique(stri_split_fixed(paste(out_bed$feature), ".", simplify = TRUE)[,1])
list_genomes_uniques<-unique(stri_split_fixed(paste(out_bed$gene_name), ".", simplify = TRUE)[,1])

write.table(bed_data_annotated, file =paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_ALL_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(out_bed, file =paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_proteincoding_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_genomes_uniques, file=paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_geneslist_pcoding_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(bed_data_annotated_unique, file=paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_geneslist_ALL_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique, file=paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSMBL_ID_pcoding_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")

fileName_chip<-basename("/home/nostromo/data/pepe/peak_calling_EGR1/EGR1_chip_T-W200-G200-FDR0.01-island.bed")
bed<-SICER_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)
bed_data_annotated_unique<-unique(stri_split_fixed(paste(bed_data_annotated$gene_name), ".", simplify = TRUE)[,1])
out_bed<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
list_ID_unique<-unique(stri_split_fixed(paste(out_bed$feature), ".", simplify = TRUE)[,1])
list_genomes_uniques<-unique(stri_split_fixed(paste(out_bed$gene_name), ".", simplify = TRUE)[,1])

write.table(bed_data_annotated, file =paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_ALL_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(out_bed, file =paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_proteincoding_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_genomes_uniques, file=paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_geneslist_pcoding_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(bed_data_annotated_unique, file=paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_geneslist_ALL_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique, file=paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSMBL_ID_pcoding_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")

fileName_chip<-basename("/home/nostromo/data/pepe/peak_calling_EGR1/selected_JAMM_EGR1.narrowPeak")
bed<-JAMM_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)
bed_data_annotated_unique<-unique(stri_split_fixed(paste(bed_data_annotated$gene_name), ".", simplify = TRUE)[,1])
out_bed<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
list_ID_unique<-unique(stri_split_fixed(paste(out_bed$feature), ".", simplify = TRUE)[,1])
list_genomes_uniques<-unique(stri_split_fixed(paste(out_bed$gene_name), ".", simplify = TRUE)[,1])

write.table(bed_data_annotated, file =paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_ALL_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(out_bed, file =paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_proteincoding_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_genomes_uniques, file=paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_geneslist_pcoding_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(bed_data_annotated_unique, file=paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_geneslist_ALL_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique, file=paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSMBL_ID_pcoding_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")

fileName_chip<-basename("/home/nostromo/data/pepe/peak_calling_EGR1/bayes_output1.txt")
bed<-BAYES_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)
bed_data_annotated_unique<-unique(stri_split_fixed(paste(bed_data_annotated$gene_name), ".", simplify = TRUE)[,1])
out_bed<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
list_ID_unique<-unique(stri_split_fixed(paste(out_bed$feature), ".", simplify = TRUE)[,1])
list_genomes_uniques<-unique(stri_split_fixed(paste(out_bed$gene_name), ".", simplify = TRUE)[,1])

write.table(bed_data_annotated, file =paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_ALL_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(out_bed, file =paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_proteincoding_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_genomes_uniques, file=paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_geneslist_pcoding_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(bed_data_annotated_unique, file=paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_geneslist_ALL_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique, file=paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSMBL_ID_pcoding_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")

fileName_chip<-basename("/home/nostromo/data/pepe/peak_calling_EGR1/RANGER_EGR1.bed_region.bed")
bed<-RANGER_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)
bed_data_annotated_unique<-unique(stri_split_fixed(paste(bed_data_annotated$gene_name), ".", simplify = TRUE)[,1])
out_bed<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
list_ID_unique<-unique(stri_split_fixed(paste(out_bed$feature), ".", simplify = TRUE)[,1])
list_genomes_uniques<-unique(stri_split_fixed(paste(out_bed$gene_name), ".", simplify = TRUE)[,1])

write.table(bed_data_annotated, file =paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_ALL_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(out_bed, file =paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_proteincoding_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_genomes_uniques, file=paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_geneslist_pcoding_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(bed_data_annotated_unique, file=paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_geneslist_ALL_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique, file=paste0("/home/nostromo/data/pepe/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSMBL_ID_pcoding_less500.txt"), quote = FALSE, row.names = FALSE, sep="\t")

fileName_chip<-basename("/home/nostromo/data/pepe/peak_calling_EGR1/ZCL_peaks.final.bed")
dirname_chip<-dirname("/home/nostromo/data/pepe/01_CHIPSEQ/probando_ZCL/CHIP/ZCL_peaks.final.bed")

ZCL_bed<-ZCL[,c("V1","V2","V3")]
bed<-ZCL_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)
bed_data_annotated_unique<-unique(stri_split_fixed(paste(bed_data_annotated$gene_name), ".", simplify = TRUE)[,1])
out_bed<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
list_ID_unique<-unique(stri_split_fixed(paste(out_bed$feature), ".", simplify = TRUE)[,1])
list_genomes_uniques<-unique(stri_split_fixed(paste(out_bed$gene_name), ".", simplify = TRUE)[,1])
list_ID_unique_ALL<-unique(stri_split_fixed(paste(bed_data_annotated$feature), ".", simplify = TRUE)[,1])

write.table(bed_data_annotated, file =paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_ALL.txt"), quote = FALSE, col.names=FALSE, row.names = FALSE, sep="\t")
write.table(out_bed, file =paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_proteincoding.txt"), quote = FALSE, col.names=FALSE, row.names = FALSE, sep="\t")
write.table(list_genomes_uniques, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_geneslist_pcoding.txt"), quote = FALSE,col.names=FALSE, row.names = FALSE, sep="\t")
write.table(bed_data_annotated_unique, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_geneslist_ALL.txt"), quote = FALSE, col.names=FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSMBL_ID_pcoding.txt"), quote = FALSE, col,names=FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique_ALL, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSMBL_ID_ALL.txt"), quote = FALSE, col.names=FALSE, row.names = FALSE, sep="\t")
