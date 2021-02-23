library(IRanges)
library(GenomicRanges)
require(ChIPpeakAnno)
library(stringi)

source("/home/nostromo/data/pepe/ZCL/library_functions_ZCL_chip.R")
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesVik.R")
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesEli.R")
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesAlba.R")


MACS<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/MACS2_EGR1_peaks.narrowPeak")
BAYES<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/bayes_output1.txt")
SICER<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/EGR1_chip_T-W200-G200-FDR0.01-island.bed")
#ZCL<-read.table("/home/nostromo/data/pepe/peak_calling_EGR1/ZCL_peaks.final.bed")
ZCL<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/ZCL_peaks_g8_a1.2_k0_modified.bed")
JAMM<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/selected_JAMM_EGR1.narrowPeak")
RANGER<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/RANGER_EGR1.bed_region.bed")

MACS_bed<-MACS[,c("V1","V2","V3")]
SICER_bed<-SICER[,c("V1","V2","V3")]
JAMM_bed<-JAMM[,c("V1","V2","V3")]
BAYES_bed<-BAYES[,c("space","start","end")]
ZCL_bed<-ZCL[,c("V1","V2","V3")]
RANGER_bed<-RANGER[,c("V1","V2","V3")]

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

####### MACS


fileName_chip<-basename("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/MACS2_EGR1_peaks.narrowPeak")
bed<-MACS_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)
bed_data_annotated_less_PC<-bed_data_annotated[(bed_data_annotated$gene_type!="protein_coding"),]
bed_protein_coding<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
list_ID_unique_less_PC<-unique(stri_split_fixed(paste(bed_data_annotated_less_PC$feature), ".", simplify = TRUE)[,1])
list_ID_unique_protein_coding<-unique(stri_split_fixed(paste(bed_protein_coding$feature), ".", simplify = TRUE)[,1])

MACS_FOR_BED_LESS_PC<-bed_data_annotated_less_PC[,c("seqnames","start","end")]
MACS_FOR_BED_PROTCODING<-bed_protein_coding[,c("seqnames","start","end")]

write.table(MACS_FOR_BED_LESS_PC,file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/BED_TO_FASTAS/", fileName_chip,"_LESS_PC.bed"),quote = FALSE, row.names = FALSE, sep="\t")
write.table(MACS_FOR_BED_PROTCODING,file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/BED_TO_FASTAS/", fileName_chip,"_PROTCODING.bed"),quote = FALSE, row.names = FALSE, sep="\t")

write.table(bed_data_annotated, file =paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_ALL.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(bed_data_annotated_less_PC, file =paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_less_PC.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(bed_protein_coding, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_pcoding.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique_less_PC, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSEMBLID_less_PC.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique_protein_coding, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSMBL_ID_pcoding.txt"), quote = FALSE, row.names = FALSE, sep="\t")


#BAYES

fileName_chip<-basename("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/bayes_output1.txt")
bed<-BAYES_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)
bed_data_annotated_less_PC<-bed_data_annotated[(bed_data_annotated$gene_type!="protein_coding"),]
bed_protein_coding<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
list_ID_unique_less_PC<-unique(stri_split_fixed(paste(bed_data_annotated_less_PC$feature), ".", simplify = TRUE)[,1])
list_ID_unique_protein_coding<-unique(stri_split_fixed(paste(bed_protein_coding$feature), ".", simplify = TRUE)[,1])

BAYES_FOR_BED_LESS_PC<-bed_data_annotated_less_PC[,c("seqnames","start","end")]
BAYES_FOR_BED_PROTCODING<-bed_protein_coding[,c("seqnames","start","end")]

write.table(BAYES_FOR_BED_LESS_PC,file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/BED_TO_FASTAS/", fileName_chip,"_LESS_PC.bed"),quote = FALSE, row.names = FALSE, sep="\t")
write.table(BAYES_FOR_BED_PROTCODING,file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/BED_TO_FASTAS/", fileName_chip,"_PROTCODING.bed"),quote = FALSE, row.names = FALSE, sep="\t")

write.table(bed_data_annotated, file =paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_ALL.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(bed_data_annotated_less_PC, file =paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_less_PC.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(bed_protein_coding, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_pcoding.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique_less_PC, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSEMBLID_less_PC.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique_protein_coding, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSMBL_ID_pcoding.txt"), quote = FALSE, row.names = FALSE, sep="\t")


####ZCL

fileName_chip<-basename("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/ZCL_peaks_g8_a1.2_k0_modified.bed")

bed<-ZCL_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)
bed_data_annotated_less_PC<-bed_data_annotated[(bed_data_annotated$gene_type!="protein_coding"),]
bed_protein_coding<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
list_ID_unique_less_PC<-unique(stri_split_fixed(paste(bed_data_annotated_less_PC$feature), ".", simplify = TRUE)[,1])
list_ID_unique_protein_coding<-unique(stri_split_fixed(paste(bed_protein_coding$feature), ".", simplify = TRUE)[,1])

ZCL_FOR_BED_LESS_PC<-bed_data_annotated_less_PC[,c("seqnames","start","end")]
ZCL_FOR_BED_PROTCODING<-bed_protein_coding[,c("seqnames","start","end")]

write.table(ZCL_FOR_BED_LESS_PC,file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/BED_TO_FASTAS/", fileName_chip,"_LESS_PC.bed"),quote = FALSE, row.names = FALSE, sep="\t")
write.table(ZCL_FOR_BED_PROTCODING,file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/BED_TO_FASTAS/", fileName_chip,"_PROTCODING.bed"),quote = FALSE, row.names = FALSE, sep="\t")

write.table(bed_data_annotated, file =paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_ALL.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(bed_data_annotated_less_PC, file =paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_less_PC.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(bed_protein_coding, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_pcoding.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique_less_PC, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSEMBLID_less_PC.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique_protein_coding, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSMBL_ID_pcoding.txt"), quote = FALSE, row.names = FALSE, sep="\t")

####SICER

fileName_chip<-basename("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/EGR1_chip_T-W200-G200-FDR0.01-island.bed")

bed<-SICER_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)
bed_data_annotated_less_PC<-bed_data_annotated[(bed_data_annotated$gene_type!="protein_coding"),]
bed_protein_coding<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
list_ID_unique_less_PC<-unique(stri_split_fixed(paste(bed_data_annotated_less_PC$feature), ".", simplify = TRUE)[,1])
list_ID_unique_protein_coding<-unique(stri_split_fixed(paste(bed_protein_coding$feature), ".", simplify = TRUE)[,1])

SICER_FOR_BED_LESS_PC<-bed_data_annotated_less_PC[,c("seqnames","start","end")]
SICER_FOR_BED_PROTCODING<-bed_protein_coding[,c("seqnames","start","end")]

write.table(SICER_FOR_BED_LESS_PC,file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/BED_TO_FASTAS/", fileName_chip,"_LESS_PC.bed"),quote = FALSE, row.names = FALSE, sep="\t")
write.table(SICER_FOR_BED_PROTCODING,file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/BED_TO_FASTAS/", fileName_chip,"_PROTCODING.bed"),quote = FALSE, row.names = FALSE, sep="\t")

write.table(bed_data_annotated, file =paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_ALL.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(bed_data_annotated_less_PC, file =paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_less_PC.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(bed_protein_coding, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_pcoding.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique_less_PC, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSEMBLID_less_PC.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique_protein_coding, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSMBL_ID_pcoding.txt"), quote = FALSE, row.names = FALSE, sep="\t")


###JAMM

fileName_chip<-basename("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/selected_JAMM_EGR1.narrowPeak")

bed<-JAMM_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)
bed_data_annotated_less_PC<-bed_data_annotated[(bed_data_annotated$gene_type!="protein_coding"),]
bed_protein_coding<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
list_ID_unique_less_PC<-unique(stri_split_fixed(paste(bed_data_annotated_less_PC$feature), ".", simplify = TRUE)[,1])
list_ID_unique_protein_coding<-unique(stri_split_fixed(paste(bed_protein_coding$feature), ".", simplify = TRUE)[,1])

JAMM_FOR_BED_LESS_PC<-bed_data_annotated_less_PC[,c("seqnames","start","end")]
JAMM_FOR_BED_PROTCODING<-bed_protein_coding[,c("seqnames","start","end")]

write.table(JAMM_FOR_BED_LESS_PC,file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/BED_TO_FASTAS/", fileName_chip,"_LESS_PC.bed"),quote = FALSE, row.names = FALSE, sep="\t")
write.table(JAMM_FOR_BED_PROTCODING,file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/BED_TO_FASTAS/", fileName_chip,"_PROTCODING.bed"),quote = FALSE, row.names = FALSE, sep="\t")

write.table(bed_data_annotated, file =paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_ALL.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(bed_data_annotated_less_PC, file =paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_less_PC.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(bed_protein_coding, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_pcoding.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique_less_PC, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSEMBLID_less_PC.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique_protein_coding, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSMBL_ID_pcoding.txt"), quote = FALSE, row.names = FALSE, sep="\t")

###RANGER

fileName_chip<-basename("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/RANGER_EGR1.bed_region.bed")

bed<-RANGER_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)
bed_data_annotated_less_PC<-bed_data_annotated[(bed_data_annotated$gene_type!="protein_coding"),]
bed_protein_coding<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
list_ID_unique_less_PC<-unique(stri_split_fixed(paste(bed_data_annotated_less_PC$feature), ".", simplify = TRUE)[,1])
list_ID_unique_protein_coding<-unique(stri_split_fixed(paste(bed_protein_coding$feature), ".", simplify = TRUE)[,1])

RANGER_FOR_BED_LESS_PC<-bed_data_annotated_less_PC[,c("seqnames","start","end")]
RANGER_FOR_BED_PROTCODING<-bed_protein_coding[,c("seqnames","start","end")]
write.table(bed_data_annotated, file =paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_ALL.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(bed_data_annotated_less_PC, file =paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_less_PC.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(bed_protein_coding, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_annot_pcoding.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique_less_PC, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSEMBLID_less_PC.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(list_ID_unique_protein_coding, file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/",fileName_chip,"_unique_ENSMBL_ID_pcoding.txt"), quote = FALSE, row.names = FALSE, sep="\t")
write.table(RANGER_FOR_BED_LESS_PC,file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/BED_TO_FASTAS/", fileName_chip,"_LESS_PC.bed"),quote = FALSE, row.names = FALSE, sep="\t")
write.table(RANGER_FOR_BED_PROTCODING,file=paste0("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/BED_TO_FASTAS/", fileName_chip,"_PROTCODING.bed"),quote = FALSE, row.names = FALSE, sep="\t")
