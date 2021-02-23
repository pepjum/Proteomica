#!/usr/bin/Rscript

args=(commandArgs(TRUE))
projectDir <- args[1]		# /home/nostromo/data/03_Analysis/agarin/28_PFortes_Shotgun_lncRNA_Feb18/
fastafileroot <- args[2]	# /home/nostromo/data/03_Analysis/agarin/28_PFortes_Shotgun_lncRNA_Feb18/nextProtDB20170801
nextprotfolder <- args[3]	# /home/nostromo/data/03_Analysis/agarin/23_neXtprot_20170801_Nov17/

library(Biostrings)
library(doBy)
library(reshape2)
library(ggplot2)
library(stringr)
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesShotgun.R")
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesVikv2.R")

load(paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_percolator.rda"))
load(paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_percolator.rda"))
load(paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_percolator.rda"))
load(paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_percolator.rda"))

results_comet_NX_annot_disc_missing_peptides_uniqueness <- read.csv2(paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_peptides_uniqueness_percolator.csv"), header = TRUE, sep = ",")
if(dim(results_comet_NX_annot_disc_missing_peptides_uniqueness)[1] != 0) {
	results_comet_NX_annot_disc_missing_uniqueness <- results_comet_NX_annot_disc_missing[paste(results_comet_NX_annot_disc_missing$PeptideSeq) %in% paste(results_comet_NX_annot_disc_missing_peptides_uniqueness$peptide), ]
	save(results_comet_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_uniqueness_percolator.rda"))
	#tmp <- unique(results_comet_NX_annot_disc_missing_uniqueness[, c("sample", "NextprotID", "PeptideSeq")])
	tmp <- unique(results_comet_NX_annot_disc_missing_uniqueness[, c("sample", "ProteinID", "PeptideSeq")])

	#tmp2 <- summaryBy(PeptideSeq~NextprotID + sample, data=tmp, FUN=length, keep.names=TRUE)
	tmp2 <- summaryBy(PeptideSeq~ProteinID + sample, data=tmp, FUN=length, keep.names=TRUE)

	colnames(tmp2)[3] <- "NofDiscPeptides"
	results_comet_NX_annot_disc_missing_uniqueness <- merge(results_comet_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("NextprotID", "sample"), by.y = c("NextprotID", "sample"), all.x = TRUE)
#	results_comet_NX_annot_disc_missing_uniqueness <- merge(results_comet_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("ProteinID", "sample"), by.y = c("ProteinID", "sample"), all.x = TRUE)
	results_comet_NX_annot_disc_missing_uniqueness[is.na(results_comet_NX_annot_disc_missing_uniqueness$NofDiscPeptides), "NofDiscPeptides"] <- 0
	results_comet_NX_annot_disc_missing_uniqueness_2UniqPep <- results_comet_NX_annot_disc_missing_uniqueness[results_comet_NX_annot_disc_missing_uniqueness$NofDiscPeptides > 1,]
	save(results_comet_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_uniqueness_2UniqPep_percolator.rda"))
	results_comet_NX_annot_disc_missing_uniqueness_f <- unique(results_comet_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "NextprotID", "search_engine", "Chr", "sample", "NofDiscPeptides")])
#	results_comet_NX_annot_disc_missing_uniqueness_f <- unique(results_comet_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "ProteinID", "search_engine", "Chr", "sample", "NofDiscPeptides")])

	results_comet_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + NextprotID + search_engine + Chr + NofDiscPeptides, data = results_comet_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")
#	results_comet_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + ProteinID + search_engine + Chr + NofDiscPeptides, data = results_comet_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")

} else {
	cat("Comet has not identified any missing proteins!:")
	results_comet_NX_annot_disc_missing_uniqueness <- results_comet_NX_annot_disc_missing[0,]
	save(results_comet_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_uniqueness_percolator.rda"))
	results_comet_NX_annot_disc_missing_uniqueness_2UniqPep <- results_comet_NX_annot_disc_missing_uniqueness
	save(results_comet_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_uniqueness_2UniqPep_percolator.rda"))
	results_comet_NX_annot_disc_missing_uniqueness_2UniqPep_f <- results_comet_NX_annot_disc_missing_uniqueness_2UniqPep
	results_comet_NX_annot_disc_missing_uniqueness_2UniqPep_ff <- results_comet_NX_annot_disc_missing_uniqueness_2UniqPep
}

results_mascot_NX_annot_disc_missing_peptides_uniqueness <- read.csv2(paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_peptides_uniqueness_percolator.csv"), header = TRUE, sep = ",")
if(dim(results_mascot_NX_annot_disc_missing_peptides_uniqueness)[1] != 0) {
	results_mascot_NX_annot_disc_missing_uniqueness <- results_mascot_NX_annot_disc_missing[paste(results_mascot_NX_annot_disc_missing$PeptideSeq) %in% paste(results_mascot_NX_annot_disc_missing_peptides_uniqueness$peptide), ]
	save(results_mascot_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_uniqueness_percolator.rda"))
	tmp <- unique(results_mascot_NX_annot_disc_missing_uniqueness[, c("sample", "NextprotID", "PeptideSeq")])
	#tmp <- unique(results_mascot_NX_annot_disc_missing_uniqueness[, c("sample", "ProteinID", "PeptideSeq")])

	tmp2 <- summaryBy(PeptideSeq~NextprotID + sample, data=tmp, FUN=length, keep.names=TRUE)
	#tmp2 <- summaryBy(PeptideSeq~ProteinID + sample, data=tmp, FUN=length, keep.names=TRUE)

	colnames(tmp2)[3] <- "NofDiscPeptides"
	results_mascot_NX_annot_disc_missing_uniqueness <- merge(results_mascot_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("NextprotID", "sample"), by.y = c("NextprotID", "sample"), all.x = TRUE)
#	results_mascot_NX_annot_disc_missing_uniqueness <- merge(results_mascot_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("ProteinID", "sample"), by.y = c("ProteinID", "sample"), all.x = TRUE)

	results_mascot_NX_annot_disc_missing_uniqueness[is.na(results_mascot_NX_annot_disc_missing_uniqueness$NofDiscPeptides), "NofDiscPeptides"] <- 0
	results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep <- results_mascot_NX_annot_disc_missing_uniqueness[results_mascot_NX_annot_disc_missing_uniqueness$NofDiscPeptides > 1,]
	save(results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep_percolator.rda"))
#	results_mascot_NX_annot_disc_missing_uniqueness_f <- unique(results_mascot_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "ProteinID", "search_engine", "Chr", "sample", "NofDiscPeptides")])
	results_mascot_NX_annot_disc_missing_uniqueness_f <- unique(results_mascot_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "NextprotID", "search_engine", "Chr", "sample", "NofDiscPeptides")])

#	results_mascot_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + ProteinID + search_engine + Chr + NofDiscPeptides, data = results_mascot_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")
	results_mascot_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + NextprotID + search_engine + Chr + NofDiscPeptides, data = results_mascot_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")

} else  {
	cat("Mascot has not identified any missing proteins!:")
	results_mascot_NX_annot_disc_missing_uniqueness <- results_mascot_NX_annot_disc_missing[0,]
	save(results_mascot_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_uniqueness_percolator.rda"))
	results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep <- results_mascot_NX_annot_disc_missing_uniqueness
	save(results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep_percolator.rda"))
	results_mascot_NX_annot_disc_missing_uniqueness_f <- results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep
	results_mascot_NX_annot_disc_missing_uniqueness_ff <- results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep
}

results_omssa_NX_annot_disc_missing_peptides_uniqueness <- read.csv2(paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_peptides_uniqueness_percolator.csv"), header = TRUE, sep = ",")
if(dim(results_omssa_NX_annot_disc_missing_peptides_uniqueness)[1] != 0) {
	results_omssa_NX_annot_disc_missing_uniqueness <- results_omssa_NX_annot_disc_missing[paste(results_omssa_NX_annot_disc_missing$PeptideSeq) %in% paste(results_omssa_NX_annot_disc_missing_peptides_uniqueness$peptide), ]
	save(results_omssa_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_uniqueness_percolator.rda"))
	tmp <- unique(results_omssa_NX_annot_disc_missing_uniqueness[, c("sample", "NextprotID", "PeptideSeq")])
	#tmp <- unique(results_omssa_NX_annot_disc_missing_uniqueness[, c("sample", "ProteinID", "PeptideSeq")])

	tmp2 <- summaryBy(PeptideSeq~NextprotID + sample, data=tmp, FUN=length, keep.names=TRUE)
	#tmp2 <- summaryBy(PeptideSeq~ProteinID + sample, data=tmp, FUN=length, keep.names=TRUE)

	colnames(tmp2)[3] <- "NofDiscPeptides"
	results_omssa_NX_annot_disc_missing_uniqueness <- merge(results_omssa_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("NextprotID", "sample"), by.y = c("NextprotID", "sample"), all.x = TRUE)
	#results_omssa_NX_annot_disc_missing_uniqueness <- merge(results_omssa_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("ProteinID", "sample"), by.y = c("ProteinID", "sample"), all.x = TRUE)

	results_omssa_NX_annot_disc_missing_uniqueness[is.na(results_omssa_NX_annot_disc_missing_uniqueness$NofDiscPeptides), "NofDiscPeptides"] <- 0
	results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep <- results_omssa_NX_annot_disc_missing_uniqueness[results_omssa_NX_annot_disc_missing_uniqueness$NofDiscPeptides > 1,]
	save(results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep_percolator.rda"))
	results_omssa_NX_annot_disc_missing_uniqueness_f <- unique(results_omssa_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "NextprotID", "search_engine", "Chr", "sample", "NofDiscPeptides")])
#	results_omssa_NX_annot_disc_missing_uniqueness_f <- unique(results_omssa_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "ProteinID", "search_engine", "Chr", "sample", "NofDiscPeptides")])

	results_omssa_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + NextprotID + search_engine + Chr + NofDiscPeptides, data = results_omssa_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")
#	results_omssa_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + ProteinID + search_engine + Chr + NofDiscPeptides, data = results_omssa_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")

} else  {
	cat("Omssa has not identified any missing proteins!:")
	results_omssa_NX_annot_disc_missing_uniqueness <- results_omssa_NX_annot_disc_missing[0,]
	save(results_omssa_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_uniqueness_percolator.rda"))
	results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep <- results_omssa_NX_annot_disc_missing_uniqueness
	save(results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep_percolator.rda"))
	results_omssa_NX_annot_disc_missing_uniqueness_f <- results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep
	results_omssa_NX_annot_disc_missing_uniqueness_ff <- results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep
}

results_tandem_NX_annot_disc_missing_peptides_uniqueness <- read.csv2(paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_peptides_uniqueness_percolator.csv"), header = TRUE, sep = ",")
if(dim(results_tandem_NX_annot_disc_missing_peptides_uniqueness)[1] != 0) {
	results_tandem_NX_annot_disc_missing_uniqueness <- results_tandem_NX_annot_disc_missing[paste(results_tandem_NX_annot_disc_missing$PeptideSeq) %in% paste(results_tandem_NX_annot_disc_missing_peptides_uniqueness$peptide), ]
	results_tandem_NX_annot_disc_missing_uniqueness$sample<-"PXD001381"
	save(results_tandem_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_uniqueness_percolator.rda"))
	tmp <- unique(results_tandem_NX_annot_disc_missing_uniqueness[, c("sample", "ProteinID", "PeptideSeq")])
	tmp2 <- summaryBy(PeptideSeq~ProteinID + sample, data=tmp, FUN=length, keep.names=TRUE)
#	tmp <- unique(results_omssa_NX_annot_disc_missing_uniqueness[, c("sample", "ProteinID", "PeptideSeq")])
#	tmp2 <- summaryBy(PeptideSeq~ProteinID + sample, data=tmp, FUN=length, keep.names=TRUE)

	colnames(tmp2)[3] <- "NofDiscPeptides"
	results_tandem_NX_annot_disc_missing_uniqueness <- merge(results_tandem_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("ProteinID", "sample"), by.y = c("ProteinID", "sample"), all.x = TRUE)
	#results_tandem_NX_annot_disc_missing_uniqueness <- merge(results_tandem_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("ProteinID", "sample"), by.y = c("ProteinID", "sample"), all.x = TRUE)

	results_tandem_NX_annot_disc_missing_uniqueness[is.na(results_tandem_NX_annot_disc_missing_uniqueness$NofDiscPeptides), "NofDiscPeptides"] <- 0
	results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep <- results_tandem_NX_annot_disc_missing_uniqueness[results_tandem_NX_annot_disc_missing_uniqueness$NofDiscPeptides > 1,]
	save(results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep_percolator.rda"))
	results_tandem_NX_annot_disc_missing_uniqueness$search_engine<-"TANDEM"
	#results_tandem_NX_annot_disc_missing_uniqueness_f <- unique(results_tandem_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "NextprotID", "search_engine", "Chr", "sample", "NofDiscPeptides")])
	results_tandem_NX_annot_disc_missing_uniqueness_f <- unique(results_tandem_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "ProteinID", "search_engine", "Chr", "sample", "NofDiscPeptides")])
	#results_tandem_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + NextprotID + search_engine + Chr + NofDiscPeptides, data = results_tandem_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")
	results_tandem_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + ProteinID + search_engine + Chr + NofDiscPeptides, data = results_tandem_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")

} else  {
	cat("Tandem has not identified any missing proteins!:")
	results_tandem_NX_annot_disc_missing_uniqueness <- results_tandem_NX_annot_disc_missing[0,]
	save(results_tandem_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_uniqueness_percolator.rda"))
	results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep <- results_tandem_NX_annot_disc_missing_uniqueness
	save(results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep_percolator.rda"))
	results_tandem_NX_annot_disc_missing_uniqueness_f <- results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep
	results_tandem_NX_annot_disc_missing_uniqueness_ff <- results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep
}

venn diagrams

pdf(file = paste0(projectDir, "Plots/results_Venn_MissingNXProteinsPerSearchEngine.pdf"), width = 10, height = 10)
#compare4List(paste(results_comet_NX_annot_disc_missing_uniqueness_2UniqPep$NextprotID), paste(results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep$NextprotID), paste(results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep$NextprotID), paste(results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep$NextprotID), "Comet results", "Mascot results", "Omssa results", "Tandem results",  "ProtFDRGlobal - Missing proteins (NX) for each search engine")
compare4List(paste(results_comet_NX_annot_disc_missing_uniqueness_2UniqPep$ProteinID), paste(results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep$ProteinID), paste(results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep$ProteinID), paste(results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep$ProteinID), "Comet results", "Mascot results", "Omssa results", "Tandem results",  "ProtFDRGlobal - Missing proteins (NX) for each search engine")

dev.off()

pdf(file = paste0(projectDir, "Plots/results_Venn_PeptidesOfMissingNXProteinsPerSearchEngine.pdf"), width = 10, height = 10)
compare4List(paste(results_comet_NX_annot_disc_missing_uniqueness_2UniqPep$PeptideSeq), paste(results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep$PeptideSeq), paste(results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep$PeptideSeq), paste(results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep$PeptideSeq), "Comet results", "Mascot results", "Omssa results", "Tandem results",  "ProtFDRGlobal - Peptides of missing proteins (NX) for each search engine")
dev.off()

# print data to table
#results_NX<-results_tandem_NX_annot_disc_missing_uniqueness_f
# results_NX_2<-results_tandem_NX_annot_disc_missing_uniqueness_ff
  results_NX <- rbind(results_comet_NX_annot_disc_missing_uniqueness_f, results_mascot_NX_annot_disc_missing_uniqueness_f, results_omssa_NX_annot_disc_missing_uniqueness_f, results_tandem_NX_annot_disc_missing_uniqueness_f)
  results_NX_2<-rbind(results_comet_NX_annot_disc_missing_uniqueness_ff, results_mascot_NX_annot_disc_missing_uniqueness_ff, results_omssa_NX_annot_disc_missing_uniqueness_ff, results_tandem_NX_annot_disc_missing_uniqueness_ff)
write.table(results_NX, file = paste0(projectDir, "Results/results_NX_percolator.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
write.table(results_NX_2, file = paste0(projectDir, "Results/results_NX_summaryby_percolator.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)

#####A PARTIR DE AQUI LO Q YO HE HECHO

# output_results<-results_NX_2[(results_NX_2$NofDiscPeptides>0),]
# write.table(output_results, file = paste0(projectDir, "Results/results_NX_filtered_NofDisc.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
#
# output_results_for_agg<-output_results2[,c("PeptideSeq","NexprotID","Chr")]
#
# output_results_for_agg_unique<-unique(output_results_for_agg)
#
# load("/home/nostromo/data/pepe/02_neXtprot_20180117_Feb18/nextProtXChrXENSP.RData")
#
# db_selected<-nextProtXChrXENSP[,c("NextprotID","GeneName")]
#
# output_final<-merge(output_results_for_agg, db_selected, by.x = c("NexprotID"), by.y = c("NextprotID"), all.x = TRUE)
# output_final<-unique(output_final)
# write.table(output_final,file=paste0(projectDir, "Results/peptides_NX_all_SE.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE))
