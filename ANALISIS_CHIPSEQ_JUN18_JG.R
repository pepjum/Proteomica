source("/home/nostromo/data/pepe/ZCL/library_functions_ZCL_chip.R")
source("/home/nostromo/data/pepe/ZCL/ZCL_CHIP.R")

#ZCL_list<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/probando_ZCL/CHIP/OUTPUT_g5_a1.2/ZCL_peaks.final.bed_unique_EMSEMBL_ID_ALL.txt")
ZCL_list<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/ZCL_peaks_g8_a1.2_k0_modified.bed_annot_ALL.txt")
bayes_list<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/anotaciones_EGR1/bayes_output1.txt_annot.txt", header=T)
MSIG_DB<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/mac_resultados/dianas_GSEA_ensemble")
ENSEMBL<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/ALL_PEAKS/Gencodev25SinExones.gtf_ENSEMBL_ID_ALL.txt")

library(dplyr)

###LISTA COMPLETA TODO G5 A 1.2
#msigdb<-269
#ensembl<-57992
TP<-semi_join(ZCL_list,MSIG_DB)   #185
FP<-nrow(ZCL_list)-nrow(MSIG_DB)   #30909
TN<-57992-nrow(ZCL_list)
FN<-269-nrow(TP) #84


write.table(TP,file="/home/nostromo/data/pepe/01_CHIPSEQ/probando_ZCL/CHIP/ZCL_list__g8_a1.2_k0_NEW_ingenuity.txt", quote = FALSE, col.names=FALSE, row.names = FALSE, sep="\t")
###LISTA COMPLETA <500 tamaño G5 a 1.2

ZCL_list_less500<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/probando_ZCL/CHIP/OUTPUT_g5_a1.2/ZCL_peaks.final.bed_unique_EMSEMBL_ID_ALL_less500.txt")
TP<-semi_join(ZCL_list_less500,MSIG_DB)
FP<-nrow(ZCL_list_less500)-nrow(MSIG_DB)
FN<-269-nrow(TP)
TN<-57992-nrow(ZCL_list_less500)

#PCODING LISTA TODO g5 a 1.2

ENSEMBL_pcoding<-19932
ENSEMBL_PCODING<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/ALL_PEAKS/Gencodev25SinExones.gtf_genecodev25_ENSMBL_ID_filtered_pcoding.txt")
ZCL_list_pcoding_All<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/probando_ZCL/CHIP/OUTPUT_g5_a1.2/ZCL_peaks.final.bed_unique_EMSEMBL_ID_pcoding.txt")
TP<-semi_join(ZCL_list_pcoding_All,MSIG_DB)
TN<-19932-nrow(ZCL_list_pcoding_All)
FN<-269-nrow(TP)
FP<-nrow(ZCL_list_pcoding_All)-nrow(MSIG_DB)

#PCODING LESS 500 g5 a 1.2
ZCL_list_pcoding_less500<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/probando_ZCL/CHIP/OUTPUT_g5_a1.2/ZCL_peaks.final.bed_unique_EMSEMBL_ID_pcoding_less500.txt")
#ZCL_list_pcoding_less500<-1820
TP<-semi_join(ZCL_list_pcoding_less500,MSIG_DB)
TN<-19932-nrow(ZCL_list_pcoding_less500)
FP<-nrow(ZCL_list_pcoding_less500)-nrow(MSIG_DB)

#####LISTA COMPLETA G6 a 1.2
#nrow(ZCL_list)<-24842
ZCL_list<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/probando_ZCL/CHIP/OUTPUT_g6_a1.2/ZCL_peaks.final.bed_unique_ENSEMBL_ID_ALL.txt")
TP<-semi_join(ZCL_list,MSIG_DB)   #185
TN<-57992-nrow(ZCL_list)-FN

###LISTA COMPLETA pcoding

ZCL_list_pcoding_All<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/probando_ZCL/CHIP/OUTPUT_g6_a1.2/ZCL_peaks.final.bed_unique_ENSEMBL_ID_pcoding_ALL.txt")
#nrow(ZCL_list_pcoding_All)<-10393
TP<-semi_join(ZCL_list_pcoding_All,MSIG_DB)   #185
TN<-19932-nrow(ZCL_list_pcoding_All)
FP<-nrow(ZCL_list_pcoding_All)-nrow(MSIG_DB)

#####ZCL PRIMERA

ZCL<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/ALL_PEAKS/ZCL_peaks.final.bed_annot_ALL.txt")
bed_data_annotated<-ZCL
bed_data_annotated_unique<-unique(stri_split_fixed(paste(bed_data_annotated$gene_name), ".", simplify = TRUE)[,1])
out_bed<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
list_ID_unique<-unique(stri_split_fixed(paste(out_bed$feature), ".", simplify = TRUE)[,1])
list_genomes_uniques<-unique(stri_split_fixed(paste(out_bed$gene_name), ".", simplify = TRUE)[,1])
list_ID_unique_ALL<-unique(stri_split_fixed(paste(bed_data_annotated$feature), ".", simplify = TRUE)[,1])

ZCL_list<-list_ID_unique_ALL
ZCL_list<-as.data.frame(ZCL_list)
names(ZCL_list)<-"V1"
TP<-semi_join(ZCL_list,MSIG_DB)   #185
FP<-nrow(ZCL_list)-nrow(MSIG_DB)
TN<-57992-nrow(ZCL_list)-FN

ZCL_list<-out_bed
list_ID_unique<-unique(stri_split_fixed(paste(out_bed$feature), ".", simplify = TRUE)[,1])
ZCL_list<-list_ID_unique
ZCL_list<-as.data.frame(ZCL_list)
names(ZCL_list)<-"V1"
TP<-semi_join(ZCL_list,MSIG_DB)   #185
FP<-nrow(ZCL_list)-nrow(MSIG_DB)
TN<-19932-nrow(ZCL_list)-FN

RANGER<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/ALL_PEAKS/RANGER_EGR1.bed_region.bed_annot_ALL.txt",header=T)
bed_data_annotated<-RANGER
bed_data_annotated_unique<-unique(stri_split_fixed(paste(bed_data_annotated$gene_name), ".", simplify = TRUE)[,1])
out_bed<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
list_ID_unique<-unique(stri_split_fixed(paste(out_bed$feature), ".", simplify = TRUE)[,1])
list_genomes_uniques<-unique(stri_split_fixed(paste(out_bed$gene_name), ".", simplify = TRUE)[,1])
list_ID_unique_ALL<-unique(stri_split_fixed(paste(bed_data_annotated$feature), ".", simplify = TRUE)[,1])

RANGER_list<-as.data.frame(list_ID_unique)
names(RANGER_list)<-"V1"
TP<-semi_join(RANGER_list,MSIG_DB)
RANGER_ALL<-as.data.frame(list_ID_unique_ALL)
names(RANGER_ALL)<-"V1"
TP<-semi_join(RANGER_ALL,MSIG_DB)

write.table(TP,file="/home/nostromo/data/pepe/01_CHIPSEQ/probando_ZCL/CHIP/OUTPUT_g8_a1.2_k0.1/RANGER_list_ingenuity.txt", quote = FALSE, col.names=FALSE, row.names = FALSE, sep="\t")



JAMM<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/ALL_PEAKS/selected_JAMM_EGR1.narrowPeak_annot_ALL.txt",header=T)
bed_data_annotated<-JAMM
bed_data_annotated_unique<-unique(stri_split_fixed(paste(bed_data_annotated$gene_name), ".", simplify = TRUE)[,1])
out_bed<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
list_ID_unique<-unique(stri_split_fixed(paste(out_bed$feature), ".", simplify = TRUE)[,1])
list_genomes_uniques<-unique(stri_split_fixed(paste(out_bed$gene_name), ".", simplify = TRUE)[,1])
list_ID_unique_ALL<-unique(stri_split_fixed(paste(bed_data_annotated$feature), ".", simplify = TRUE)[,1])

JAMM_list<-as.data.frame(list_ID_unique)
names(JAMM_list)<-"V1"
TP<-semi_join(JAMM_list,MSIG_DB)
JAMM_ALL<-as.data.frame(list_ID_unique_ALL)
names(JAMM_ALL)<-"V1"
TP<-semi_join(JAMM_ALL,MSIG_DB)

write.table(TP,file="/home/nostromo/data/pepe/01_CHIPSEQ/probando_ZCL/CHIP/OUTPUT_g8_a1.2_k0.1/JAMM_list_ingenuity.txt", quote = FALSE, col.names=FALSE, row.names = FALSE, sep="\t")


MACS<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/ALL_PEAKS/MACS2_EGR1_peaks.narrowPeak_annot_ALL.txt",header=T)
bed_data_annotated<-MACS
bed_data_annotated_unique<-unique(stri_split_fixed(paste(bed_data_annotated$gene_name), ".", simplify = TRUE)[,1])
out_bed<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
list_ID_unique<-unique(stri_split_fixed(paste(out_bed$feature), ".", simplify = TRUE)[,1])
list_genomes_uniques<-unique(stri_split_fixed(paste(out_bed$gene_name), ".", simplify = TRUE)[,1])
list_ID_unique_ALL<-unique(stri_split_fixed(paste(bed_data_annotated$feature), ".", simplify = TRUE)[,1])

MACS_list<-as.data.frame(list_ID_unique)
names(MACS_list)<-"V1"
TP<-semi_join(MACS_list,MSIG_DB)
MACS_ALL<-as.data.frame(list_ID_unique_ALL)
names(MACS_ALL)<-"V1"
TP<-semi_join(MACS_ALL,MSIG_DB)
write.table(TP,file="/home/nostromo/data/pepe/01_CHIPSEQ/probando_ZCL/CHIP/OUTPUT_g8_a1.2_k0.1/MACS_list_ingenuity.txt", quote = FALSE, col.names=FALSE, row.names = FALSE, sep="\t")




SICER<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/ALL_PEAKS/EGR1_chip_T-W200-G200-FDR0.01-island.bed_annot_ALL.txt",header=T)
bed_data_annotated<-SICER
bed_data_annotated_unique<-unique(stri_split_fixed(paste(bed_data_annotated$gene_name), ".", simplify = TRUE)[,1])
out_bed<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
list_ID_unique<-unique(stri_split_fixed(paste(out_bed$feature), ".", simplify = TRUE)[,1])
list_genomes_uniques<-unique(stri_split_fixed(paste(out_bed$gene_name), ".", simplify = TRUE)[,1])
list_ID_unique_ALL<-unique(stri_split_fixed(paste(bed_data_annotated$feature), ".", simplify = TRUE)[,1])

SICER_list<-as.data.frame(list_ID_unique)
names(SICER_list)<-"V1"
TP<-semi_join(SICER_list,MSIG_DB)
SICER_ALL<-as.data.frame(list_ID_unique_ALL)
names(SICER_ALL)<-"V1"
TP<-semi_join(SICER_ALL,MSIG_DB)

write.table(TP,file="/home/nostromo/data/pepe/01_CHIPSEQ/probando_ZCL/CHIP/OUTPUT_g8_a1.2_k0.1/SICER_list_ingenuity.txt", quote = FALSE, col.names=FALSE, row.names = FALSE, sep="\t")

library(stringi)
BAYES<-read.table("/home/nostromo/data/pepe/01_CHIPSEQ/peak_calling_EGR1/anotaciones/ALL_PEAKS/bayes_output1.txt_annot_ALL.txt",header=T)
bed_data_annotated<-BAYES
bed_data_annotated_unique<-unique(stri_split_fixed(paste(bed_data_annotated$gene_name), ".", simplify = TRUE)[,1])
out_bed<-bed_data_annotated[(bed_data_annotated$gene_type=="protein_coding"),]
#out_bed<-bed_data_annotated
list_ID_unique<-unique(stri_split_fixed(paste(out_bed$feature), ".", simplify = TRUE)[,1])
list_genomes_uniques<-unique(stri_split_fixed(paste(out_bed$gene_name), ".", simplify = TRUE)[,1])
list_ID_unique_ALL<-unique(stri_split_fixed(paste(bed_data_annotated$feature), ".", simplify = TRUE)[,1])

BAYES_list<-as.data.frame(list_ID_unique)
names(BAYES_list)<-"V1"
TP<-semi_join(BAYES_list,MSIG_DB)
BAYES_ALL<-as.data.frame(list_ID_unique_ALL)
names(BAYES_ALL)<-"V1"
TP<-semi_join(BAYES_ALL,MSIG_DB)

write.table(TP,file="/home/nostromo/data/pepe/01_CHIPSEQ/probando_ZCL/CHIP/OUTPUT_g8_a1.2_k0.1/BAYES_list_ingenuity.txt", quote = FALSE, col.names=FALSE, row.names = FALSE, sep="\t")
