library(ape)
library(IRanges)
library(GenomicRanges)
require(ChIPpeakAnno)
library(stringi)

source("/home/margaret/data/pepe/scripts/ZCL/library_functions_ZCL_chip.R")
source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesVikv2.R")
source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesEli.R")
#source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesAlba.R")
#source("")

genecodev25_all <- read.table(file ="/home/margaret/data/pepe/Gencodev25SinExones.gtf", skip = 5, header = FALSE, sep = "\t")

genecodev25 <- genecodev25_all[genecodev25_all$V3 == "gene", ]
genecodev25_tmp <- apply(as.data.frame(genecodev25[,9]), 1, parseENCODE)
genecodev25_tmp <- t(genecodev25_tmp)

colnames(genecodev25_tmp) <- c("gene_id", "gene_type", "gene_status", "gene_name", "level", "havana_gene")
genecodev25_Annot <- data.frame(genecodev25[,1:8], genecodev25_tmp[, c(1, 2,3, 4, 5,6)])

genecodev25_Annot_Gene <- genecodev25_Annot

G25AnnotChIP <- unique(genecodev25_Annot[genecodev25_Annot[, "V3"] == "gene", c("V1", "V4", "V5", "V7", "gene_id")])
colnames(G25AnnotChIP) <- c("chr", "start", "end", "strand", "gene_id")
G25AnnotChIP_RL <- RangedData(IRanges(start = G25AnnotChIP$start, end = G25AnnotChIP$end, names = paste(G25AnnotChIP$gene_id)), space = paste(G25AnnotChIP$chr), strand = paste(G25AnnotChIP$strand))

#####MACS

file_MACS<-read.gff("/home/margaret/data/pepe/01_FASTAS_A_MEME_COMPLETAS/MACS_EGR1/fimo_out_1/fimo.gff", na.strings = c(".", "?"), GFF3 = TRUE)
MACS_bed<-file_MACS[,c("seqid","start","end")]
bed<-MACS_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated_MACS<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)

MACS_df<-data.frame("gene_type"=paste(bed_data_annotated_MACS$gene_type),"METHOD"="MACS")
MACS_df_end<-as.data.frame(table(MACS_df$gene_type))

####SICER

file_SICER<-read.gff("/home/margaret/data/pepe/01_FASTAS_A_MEME_COMPLETAS/SICER_MEME/fimo_out_3/fimo.gff", na.strings = c(".", "?"), GFF3 = TRUE)
SICER_bed<-file_SICER[,c("seqid","start","end")]
bed<-SICER_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated_SICER<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)

SICER_df<-data.frame("gene_type"=paste(bed_data_annotated_SICER$gene_type),"METHOD"="SICER")
SICER_df_end<-as.data.frame(table(SICER_df$gene_type))

####RANGER

file_RANGER<-read.gff("/home/margaret/data/pepe/01_FASTAS_A_MEME_COMPLETAS/RANGER_MEME/RANGER/fimo_out_1/fimo.gff", na.strings = c(".", "?"), GFF3 = TRUE)
RANGER_bed<-file_RANGER[,c("seqid","start","end")]
bed<-RANGER_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated_RANGER<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)

RANGER_df<-data.frame("gene_type"=paste(bed_data_annotated_RANGER$gene_type),"METHOD"="RANGER")
RANGER_df_end<-as.data.frame(table(RANGER_df$gene_type))

####ZCL

file_ZCL<-read.gff("/home/margaret/data/pepe/01_FASTAS_A_MEME_COMPLETAS/ZCL_NEW_MEME/fimo_out_2/fimo.gff", na.strings = c(".", "?"), GFF3 = TRUE)
ZCL_bed<-file_ZCL[,c("seqid","start","end")]
ZCL_bed$seqid<-sapply(strsplit(paste(ZCL_bed$seqid),"_"),"[",2)
bed<-ZCL_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated_ZCL<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)

ZCL_df<-data.frame("gene_type"=paste(bed_data_annotated_ZCL$gene_type),"METHOD"="ZCL")
ZCL_df_end<-as.data.frame(table(ZCL_df$gene_type))
####BAYES

file_BAYES<-read.gff("/home/margaret/data/pepe/01_FASTAS_A_MEME_COMPLETAS/BAYES_MEME/fimo_out_1/fimo.gff", na.strings = c(".", "?"), GFF3 = TRUE)
BAYES_bed<-file_BAYES[,c("seqid","start","end")]
bed<-BAYES_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated_BAYES<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)

BAYES_df<-data.frame("gene_type"=paste(bed_data_annotated_BAYES$gene_type),"METHOD"="BAYES")
BAYES_df_end<-as.data.frame(table(BAYES_df$gene_type))

####JAMM

file_JAMM<-read.gff("/home/margaret/data/pepe/01_FASTAS_A_MEME_COMPLETAS/JAMM_MEME/fimo_out_1/fimo.gff", na.strings = c(".", "?"), GFF3 = TRUE)
JAMM_bed<-file_JAMM[,c("seqid","start","end")]
bed<-JAMM_bed
names(bed)<-c("V1","V2","V3")                                       ###change this
bed_data_annotated_JAMM<-annotChIP(bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)

JAMM_df<-data.frame("gene_type"=paste(bed_data_annotated_JAMM$gene_type),"METHOD"="JAMM")
JAMM_df_end<-as.data.frame(table(JAMM_df$gene_type))

#plotter<-rbind(ZCL_df,MACS_df,SICER_df,RANGER_df,JAMM_df,BAYES_df)

library(ggplot2)

pdf("ZCL_piechart.pdf",width = 22, height=15)
ggplot(ZCL_df_end, aes(x="", y=Freq, fill=Var1))+ geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_color_brewer(type="div", palette="Spectral") + theme_bw() + theme(legend.justification=c(0.5,0), legend.position=c(1,0))
dev.off()

pdf("MACS_piechart.pdf",width = 22, height=15)
ggplot(MACS_df_end, aes(x="", y=Freq, fill=Var1))+ geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_color_brewer(type="div", palette="Spectral") + theme_bw() + theme(legend.justification=c(0.5,0), legend.position=c(1,0))
dev.off()

pdf("SICER_piechart.pdf",width = 22, height=15)
ggplot(SICER_df_end, aes(x="", y=Freq, fill=Var1))+ geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_color_brewer(type="div", palette="Spectral") + theme_bw() + theme(legend.justification=c(0.5,0), legend.position=c(1,0))
dev.off()

pdf("RANGER_piechart.pdf",width = 22, height=15)
ggplot(RANGER_df_end, aes(x="", y=Freq, fill=Var1))+ geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_color_brewer(type="div", palette="Spectral") + theme_bw() + theme(legend.justification=c(0.5,0), legend.position=c(1,0))
dev.off()

pdf("BAYES_piechart.pdf",width = 22, height=15)
ggplot(BAYES_df_end, aes(x="", y=Freq, fill=Var1))+ geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_color_brewer(type="div", palette="Spectral") + theme_bw() + theme(legend.justification=c(0.5,0), legend.position=c(1,0))
dev.off()

pdf("JAMM_piechart.pdf",width = 22, height=15)
ggplot(JAMM_df_end, aes(x="", y=Freq, fill=Var1))+ geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_color_brewer(type="div", palette="Spectral") + theme_bw() + theme(legend.justification=c(0.5,0), legend.position=c(1,0))
dev.off()
