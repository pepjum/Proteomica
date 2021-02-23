source("/home/nostromo/data/pepe/ZCL/library_functions_ZCL_chip.R") # cambia por analogos en el cluster
source("/home/nostromo/data/pepe/ZCL/ZCL_CHIP.R")

library(dplyr)
library(IRanges)
library(ChIPpeakAnno)

##### READING FILES FROM ELI, MACARENA AND ME

exp1<-read.table("/mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/06_SEÑALES_CHIPSEQ_ZCL_h3k4me3_ENE18/regions_H3K4me3_PEPE.bed")
exp2<-read.table("/mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/06_SEÑALES_CHIPSEQ_ZCL_h3k4me3_ENE18/regions_eli_todo.bed")
exp3<-read.table("/mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/06_SEÑALES_CHIPSEQ_ZCL_h3k4me3_ENE18/CHIP_h3k4me3_T_regions_macarena.bed")

#merging in one dataframe

groundtruth_raw <- rbind(exp1,exp2, exp3)

# order BY COLUMN V2, increasing

groundtruth_sorted<- groundtruth_raw[order(groundtruth_raw$V2),]

#merging peaks less than 500 Da

groundtruth_sorted_clean<-groundtruth_sorted[c("V2","V3")]
exp1_clean<-exp1[c("V2","V3")]
exp2_clean<-exp2[c("V2","V3")]
exp3_clean<-exp3[c("V2","V3")]

unique<-unique(groundtruth_sorted_clean)

bining<-function(dataframe,Dalton){
    dataframe<-as.data.frame(dataframe)
    #dataframe<-dataframe[!(dataframe$M=="NA"),]
	new = 0
	old = nrow(dataframe)

	dataframe <- setNames(dataframe,c("M","Z"))
	while( old != new){
        print(paste0(new, "  Vs  ", old))

		tmp = data.frame('M' = NULL, 'Z' = NULL)
		SKIP_LINE = F
		iteration = nrow(dataframe)- 1
		for (i in 1:(iteration)){
			if(SKIP_LINE){
				SKIP_LINE = F
				next
			}
			dataframe<-dataframe[order(as.numeric(paste(dataframe$M))),]
			if(as.numeric(paste(dataframe[i+1, 'M'])) - as.numeric(paste(dataframe[i, 'M']))  <= as.numeric(Dalton)){
				#print('shrinking')
				maxM = max(as.numeric(paste(dataframe[i+1, 'M'])), as.numeric(paste(dataframe[i, 'M'] )))
				maxZ = max(as.numeric(paste(dataframe[i+1, 'Z'])), as.numeric(paste(dataframe[i, 'Z'] )))
				tmp = rbind(tmp, data.frame('M' = as.numeric(maxM), 'Z' = as.numeric(maxZ)))
				SKIP_LINE = T
			}else{
				tmp = rbind(tmp, dataframe[i,])
				SKIP_LINE = F
			}
		}
		tmp <- rbind(tmp, dataframe[i+1,])
		new <- nrow(tmp)
		old <- nrow(dataframe)
		dataframe <- tmp
	}
	if(as.numeric(paste(dataframe[nrow(dataframe), 'M'])) - as.numeric(paste(dataframe[nrow(dataframe)-1, 'M'])) <= as.numeric(Dalton)){
		maxM = max(as.numeric(paste(dataframe[nrow(dataframe), 'M'])), as.numeric(paste(dataframe[nrow(dataframe)-1 , 'M'] )))
		maxZ = max(as.numeric(paste(dataframe[nrow(dataframe), 'Z'])), as.numeric(paste(dataframe[nrow(dataframe)-1, 'Z'] )))
		dataframe2 <- dataframe[c(seq(1,nrow(dataframe)-2)), ]
		dataframe2 <- rbind(dataframe2, data.frame('M' = as.numeric(maxM), 'Z' = as.numeric(maxZ)))
	}else{
        dataframe2 <- dataframe

    }
    dataframe2<-unique(dataframe2)
	return(dataframe2)

}

groundtruth_clustered<-bining(unique,500)
exp1_clustered<-bining(exp1_clean,500)
exp2_clustered<-bining(exp2_clean,500)
exp3_clustered<-bining(exp3_clean,500)

write.table(groundtruth_clustered,"/mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/06_SEÑALES_CHIPSEQ_ZCL_h3k4me3_ENE18/GROUND_TRUTH.bed",sep="\t",row.names=FALSE, col.names=FALSE)
#1418
write.table(exp1_clustered,"/mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/06_SEÑALES_CHIPSEQ_ZCL_h3k4me3_ENE18/exp1_clustered.bed",sep="\t",row.names=FALSE, col.names=FALSE)
#1293
write.table(exp2_clustered,"/mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/06_SEÑALES_CHIPSEQ_ZCL_h3k4me3_ENE18/exp2_clustered.bed",sep="\t",row.names=FALSE, col.names=FALSE)
#853
write.table(exp3_clustered,"/mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/06_SEÑALES_CHIPSEQ_ZCL_h3k4me3_ENE18/exp3_clustered.bed",sep="\t",row.names=FALSE, col.names=FALSE)
#927
###### LOAD MACS RESULTS TO COMPARE WITH exp1

MACS<-read.table("/mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/10_CHIPSEQ_ANOTACIONES_METODOS_2018/MACS/CHIP_h3k4me3_peaks.narrowPeak", header=T)

#only chr16
MACS_chr16<-MACS[which(MACS$V1=="chr16"),]
#1011
####ANOTAR REGIONES PARA PODER CREAR UN OBJETO GRANGES

# parseENCODE <- function(x) {
#
# 	tmp <- unlist(strsplit(unlist(strsplit(x, ";")), " "))
# 	tmp <- tmp[tmp != ""]
# 	tmp <- tmp[seq(2, 12, 2)]
# 	return(tmp)
# }

# annotChIP<- function(file_bed, file_annot, file_annot_gene) {
#
# 	beddata_RL <- BED2RangedData(file_bed, header=FALSE)
#
# 	beddata_Annot <- annotatePeakInBatch(beddata_RL, AnnotationData =file_annot, output = "both", multiple = F, maxgap = 0)
# 	beddata_Annot.df <- as.data.frame(beddata_Annot)
# 	#write.table(beddata_Annot.df, file = paste(nameOut, "_Annot.txt", sep = ""), quote = FALSE, row.names = FALSE, sep="\t")
# 	#beddata_Annot_G15.df <- merge(beddata_Annot.df, file_annot_gene, by.x = 7, by.y = 9)
# 	#write.table(beddata_Annot_G15.df, file = paste(nameOut, "_Annot_G15.txt", sep = ""), quote = FALSE, row.names = FALSE, sep="\t")
#
# 	#beddata.Filter10 <- beddata_Annot_G15.df[((beddata_Annot_G15.df$fromOverlappingOrNearest == "NearestLocation") & (beddata_Annot_G15.df$insideFeature == "upstream") & (abs(beddata_Annot_G15.df$distancetoFeature)<10000)) | ((beddata_Annot_G15.df$fromOverlappingOrNearest == "NearestLocation") & (beddata_Annot_G15.df$insideFeature == "inside") & (abs(beddata_Annot_G15.df$distancetoFeature)<10000)) | ((beddata_Annot_G15.df$fromOverlappingOrNearest == "NearestLocation") & (beddata_Annot_G15.df$insideFeature == "overlapStart")) | ((beddata_Annot_G15.df$fromOverlappingOrNearest == "NearestLocation") & (beddata_Annot_G15.df$insideFeature == "includeFeature")),]
#
# 	#write.table(beddata.Filter10, file = paste(nameOut, "_Filter10.txt", sep = ""), quote = FALSE, row.names = FALSE, sep="\t")
# 	return(beddata_Annot.df)
# }

# cat("annotation of peaks...",fileName_chip,"\n")
# genecodev25_all <- read.table(file ="/mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/Gencodev25SinExones.gtf" , skip = 5, header = FALSE, sep = "\t")
# genecodev25 <- genecodev25_all[genecodev25_all$V3 == "gene", ]
# genecodev25_tmp <- apply(as.data.frame(genecodev25[,9]), 1, parseENCODE)
# genecodev25_tmp <- t(genecodev25_tmp)
#
# colnames(genecodev25_tmp) <- c("gene_id", "gene_type", "gene_status", "gene_name", "level", "havana_gene")
# genecodev25_Annot <- data.frame(genecodev25[,1:8], genecodev25_tmp[, c(1, 3, 4, 5)])
#
# genecodev25_Annot_Gene <- genecodev25_Annot
#
# G25AnnotChIP <- unique(genecodev25_Annot[genecodev25_Annot[, "V3"] == "gene", c("V1", "V4", "V5", "V7", "gene_id")])
# colnames(G25AnnotChIP) <- c("chr", "start", "end", "strand", "gene_id")
# G25AnnotChIP_RL <- RangedData(IRanges(start = G25AnnotChIP$start, end = G25AnnotChIP$end, names = paste(G25AnnotChIP$gene_id)), space = paste(G25AnnotChIP$chr), strand = paste(G25AnnotChIP$strand))
#
# chr16_names_exp1<-as.data.frame(rep("chr16",nrow(exp1_clustered)))

BED_MACS_prepared<-MACS_chr16[,c(1:3)]
BED_MACS_prepared$name<-paste0("MACS_peak_",1:nrow(BED_MACS_prepared))
chr16_names_exp1<-as.data.frame(rep("chr16",nrow(exp1_clustered)))
exp1_clustered<-cbind(chr16_names_exp1,exp1_clustered)
exp1_clustered$name<-paste0("Exp1_peak_",1:nrow(exp1_clustered))

exp1_clustered<-setNames(exp1_clustered,c("chr","start","end","names"))
BED_MACS_prepared<-setNames(BED_MACS_prepared,c("chr","start","end","names"))

MACS_ranges<-GRanges(BED_MACS_prepared)
exp1_ranges<-GRanges(exp1_clustered)

# n_peaks_overlapping_connectedPeaks_merge_EXP1_MACS<-findOverlapsOfPeaks(exp1_ranges,MACS_ranges, connectedPeaks="merge")[1]$venn_cnt
##860!!

# peaks_overlapped_EXP1_MACS<-unname(countOverlaps(exp1_ranges,MACS_ranges))
# peaks_overlapped_EXP1_MACS<-replace(peaks_overlapped_EXP1_MACS, peaks_overlapped_EXP1_MACS>1,1)
#1157
peaks_overlapped_EXP1_MACS<-annotatePeakInBatch(MACS_ranges, exp1_ranges, output = "overlapping")
peaks_from_EXP1_matched_MACS<-peaks_overlapped_EXP1_MACS$peaks
location_EXP1_matched_MACS<-peaks_overlapped_EXP1_MACS$fromOverlappingOrNearest
MACS_df<-data.frame("peaks_ref"=peaks_from_EXP1_matched_MACS,"location"=location)
MACS_df_res<-MACS_df[which(MACS_df$location=="Overlapping"),]
#nrow(MACS_df_res)==1223
###LOAD BAYESPEAK

BAYES<-read.table("/mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/10_CHIPSEQ_ANOTACIONES_METODOS_2018/BAYESPEAK/output1.txt")
BAYES_chr16<-BAYES[which(BAYES$space=="chr16"),]
#1539
BAYES_chr16<-BAYES_chr16[,c(1:3)]
names(BAYES_chr16)<-c("chr","start","end")
BAYES_chr16$name<-paste0("BAYES_peak_",1:nrow(BAYES_chr16))
BAYES_ranges<-GRanges(BAYES_chr16)

#n_peaks_overlapping_connectedPeaks_merge_EXP1_BAYES<-findOverlapsOfPeaks(exp1_ranges,BAYES_ranges, connectedPeaks="merge")[1]$venn_cnt
#854 !!!!

#peaks_overlapped_EXP1_BAYES<-unname(countOverlaps(exp1_ranges,BAYES_ranges))
#peaks_overlapped_EXP1_BAYES<-replace(peaks_overlapped_EXP1_BAYES, peaks_overlapped_EXP1_BAYES>1,1)
#1078

peaks_overlapped_EXP1_BAYES<-annotatePeakInBatch(BAYES_ranges, exp1_ranges, output = "overlapping")
peaks_from_EXP1_matched_BAYES<-peaks_overlapped_EXP1_BAYES$peak
location_EXP1_matched_BAYES<-peaks_overlapped_EXP1_BAYES$fromOverlappingOrNearest
BAYES_df<-data.frame("peaks_ref"=peaks_from_EXP1_matched_BAYES,"location"=location_EXP1_matched_BAYES)
BAYES_df_res<-BAYES_df[which(BAYES_df$location=="Overlapping"),]
#nrow(BAYES_df_res)==1200

JAMM<-read.table("/mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/10_CHIPSEQ_ANOTACIONES_METODOS_2018/JAMM/peaks/selected_JAMM.narrowPeak")
JAMM_chr16<-JAMM[which(JAMM$V1=="chr16"),]

BED_JAMM_prepared<-JAMM_chr16[,c(1:3)]

BED_JAMM_prepared$name<-paste0("JAMM_peak_",1:nrow(JAMM_chr16))
names(BED_JAMM_prepared)<-c("chr","start","end","name")
JAMM_ranges<-GRanges(BED_JAMM_prepared)

peaks_overlapped_EXP1_JAMM<-annotatePeakInBatch(JAMM_ranges, exp1_ranges, output = "overlapping")
peaks_from_EXP1_matched_JAMM<-peaks_overlapped_EXP1_JAMM$peak
location_EXP1_matched_JAMM<-peaks_overlapped_EXP1_JAMM$fromOverlappingOrNearest
JAMM_df<-data.frame("peaks_ref"=peaks_from_EXP1_matched_JAMM,"location"=location_EXP1_matched_JAMM)
JAMM_df_res<-JAMM_df[which(JAMM_df$location=="Overlapping"),]
#nrow(JAMM_df_res)==991

SICER<-read.table("/mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/10_CHIPSEQ_ANOTACIONES_METODOS_2018/SICER/CHIP_h3k4me3-W200-G200-FDR0.01-island.bed")
SICER_chr16<-SICER[which(SICER$V1=="chr16"),]
#1067
SICER_chr16$V4<-NULL
SICER_chr16$name<-paste0("SICER_peak_",1:nrow(SICER_chr16))
names(SICER_chr16)<-c("chr","start","end","name")
SICER_ranges<-GRanges(SICER_chr16)

peaks_overlapped_EXP1_SICER<-annotatePeakInBatch(SICER_ranges, exp1_ranges, output = "overlapping")
peaks_from_EXP1_matched_SICER<-peaks_overlapped_EXP1_SICER$peak
location_EXP1_matched_SICER<-peaks_overlapped_EXP1_SICER$fromOverlappingOrNearest
SICER_df<-data.frame("peaks_ref"=peaks_from_EXP1_matched_SICER,"location"=location_EXP1_matched_SICER)
SICER_df_res<-SICER_df[which(SICER_df$location=="Overlapping"),]
#nrow(SICER_df_res)==1262
