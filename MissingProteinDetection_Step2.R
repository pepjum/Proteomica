#!/usr/bin/Rscript

args=(commandArgs(TRUE))
projectDir <- args[1]		# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/28_PFortes_Shotgun_lncRNA_Feb18/
fastafileroot <- args[2]	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/28_PFortes_Shotgun_lncRNA_Feb18/nextProtDB20170801
nextprotfolder <- args[3]	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/23_neXtprot_20170801_Nov17/

library(Biostrings)
library(doBy)
library(reshape2)
library(ggplot2)
library(stringr)
source("/mnt/beegfs/agarin/dato-activo/01_Rscripts/A_Funciones/funcionesShotgun.R")
source("/mnt/beegfs/agarin/dato-activo/01_Rscripts/A_Funciones/funcionesVikv2.R")

load(paste0(projectDir, "Results/results_comet_NX_annot_disc_missing.rda")) 
load(paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing.rda"))
load(paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing.rda")) 
load(paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing.rda")) 

results_comet_NX_annot_disc_missing_peptides_uniqueness <- read.csv2(paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_peptides_uniqueness.csv"), header = TRUE, sep = ",") 
if(dim(results_comet_NX_annot_disc_missing_peptides_uniqueness)[1] != 0) {
	results_comet_NX_annot_disc_missing_uniqueness <- results_comet_NX_annot_disc_missing[paste(results_comet_NX_annot_disc_missing$PeptideSeq) %in% paste(results_comet_NX_annot_disc_missing_peptides_uniqueness$peptide), ]
	save(results_comet_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_uniqueness.rda"))
	tmp <- unique(results_comet_NX_annot_disc_missing_uniqueness[, c("sample", "NextprotID", "PeptideSeq")])
	tmp2 <- summaryBy(PeptideSeq~NextprotID + sample, data=tmp, FUN=length, keep.names=TRUE)
	colnames(tmp2)[3] <- "NofDiscPeptides"
	results_comet_NX_annot_disc_missing_uniqueness <- merge(results_comet_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("NextprotID", "sample"), by.y = c("NextprotID", "sample"), all.x = TRUE)
	results_comet_NX_annot_disc_missing_uniqueness[is.na(results_comet_NX_annot_disc_missing_uniqueness$NofDiscPeptides), "NofDiscPeptides"] <- 0
	results_comet_NX_annot_disc_missing_uniqueness_2UniqPep <- results_comet_NX_annot_disc_missing_uniqueness[results_comet_NX_annot_disc_missing_uniqueness$NofDiscPeptides > 1,]
	save(results_comet_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
	results_comet_NX_annot_disc_missing_uniqueness_f <- unique(results_comet_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "NextprotID", "search_engine", "Chr", "sample", "NofDiscPeptides")])
	results_comet_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + NextprotID + search_engine + Chr + NofDiscPeptides, data = results_comet_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")

} else {
	cat("Comet has not identified any missing proteins!:")
	results_comet_NX_annot_disc_missing_uniqueness <- results_comet_NX_annot_disc_missing[0,]
	save(results_comet_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_uniqueness.rda"))
	results_comet_NX_annot_disc_missing_uniqueness_2UniqPep <- results_comet_NX_annot_disc_missing_uniqueness
	save(results_comet_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
	results_comet_NX_annot_disc_missing_uniqueness_2UniqPep_f <- results_comet_NX_annot_disc_missing_uniqueness_2UniqPep
	results_comet_NX_annot_disc_missing_uniqueness_2UniqPep_ff <- results_comet_NX_annot_disc_missing_uniqueness_2UniqPep	
}

results_mascot_NX_annot_disc_missing_peptides_uniqueness <- read.csv2(paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_peptides_uniqueness.csv"), header = TRUE, sep = ",") 
if(dim(results_mascot_NX_annot_disc_missing_peptides_uniqueness)[1] != 0) {
	results_mascot_NX_annot_disc_missing_uniqueness <- results_mascot_NX_annot_disc_missing[paste(results_mascot_NX_annot_disc_missing$PeptideSeq) %in% paste(results_mascot_NX_annot_disc_missing_peptides_uniqueness$peptide), ]
	save(results_mascot_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_uniqueness.rda"))
	tmp <- unique(results_mascot_NX_annot_disc_missing_uniqueness[, c("sample", "NextprotID", "PeptideSeq")])
	tmp2 <- summaryBy(PeptideSeq~NextprotID + sample, data=tmp, FUN=length, keep.names=TRUE)
	colnames(tmp2)[3] <- "NofDiscPeptides"
	results_mascot_NX_annot_disc_missing_uniqueness <- merge(results_mascot_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("NextprotID", "sample"), by.y = c("NextprotID", "sample"), all.x = TRUE)
	results_mascot_NX_annot_disc_missing_uniqueness[is.na(results_mascot_NX_annot_disc_missing_uniqueness$NofDiscPeptides), "NofDiscPeptides"] <- 0
	results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep <- results_mascot_NX_annot_disc_missing_uniqueness[results_mascot_NX_annot_disc_missing_uniqueness$NofDiscPeptides > 1,]
	save(results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
	results_mascot_NX_annot_disc_missing_uniqueness_f <- unique(results_mascot_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "NextprotID", "search_engine", "Chr", "sample", "NofDiscPeptides")])
	results_mascot_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + NextprotID + search_engine + Chr + NofDiscPeptides, data = results_mascot_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")
} else  {
	cat("Mascot has not identified any missing proteins!:")
	results_mascot_NX_annot_disc_missing_uniqueness <- results_mascot_NX_annot_disc_missing[0,]
	save(results_mascot_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_uniqueness.rda"))
	results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep <- results_mascot_NX_annot_disc_missing_uniqueness
	save(results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
	results_mascot_NX_annot_disc_missing_uniqueness_f <- results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep
	results_mascot_NX_annot_disc_missing_uniqueness_ff <- results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep	
}			

results_omssa_NX_annot_disc_missing_peptides_uniqueness <- read.csv2(paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_peptides_uniqueness.csv"), header = TRUE, sep = ",") 
if(dim(results_omssa_NX_annot_disc_missing_peptides_uniqueness)[1] != 0) {
	results_omssa_NX_annot_disc_missing_uniqueness <- results_omssa_NX_annot_disc_missing[paste(results_omssa_NX_annot_disc_missing$PeptideSeq) %in% paste(results_omssa_NX_annot_disc_missing_peptides_uniqueness$peptide), ]
	save(results_omssa_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_uniqueness.rda"))
	tmp <- unique(results_omssa_NX_annot_disc_missing_uniqueness[, c("sample", "NextprotID", "PeptideSeq")])
	tmp2 <- summaryBy(PeptideSeq~NextprotID + sample, data=tmp, FUN=length, keep.names=TRUE)
	colnames(tmp2)[3] <- "NofDiscPeptides"
	results_omssa_NX_annot_disc_missing_uniqueness <- merge(results_omssa_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("NextprotID", "sample"), by.y = c("NextprotID", "sample"), all.x = TRUE)
	results_omssa_NX_annot_disc_missing_uniqueness[is.na(results_omssa_NX_annot_disc_missing_uniqueness$NofDiscPeptides), "NofDiscPeptides"] <- 0
	results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep <- results_omssa_NX_annot_disc_missing_uniqueness[results_omssa_NX_annot_disc_missing_uniqueness$NofDiscPeptides > 1,]
	save(results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
	results_omssa_NX_annot_disc_missing_uniqueness_f <- unique(results_omssa_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "NextprotID", "search_engine", "Chr", "sample", "NofDiscPeptides")])
	results_omssa_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + NextprotID + search_engine + Chr + NofDiscPeptides, data = results_omssa_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")
} else  {
	cat("Omssa has not identified any missing proteins!:")
	results_omssa_NX_annot_disc_missing_uniqueness <- results_omssa_NX_annot_disc_missing[0,]
	save(results_omssa_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_uniqueness.rda"))
	results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep <- results_omssa_NX_annot_disc_missing_uniqueness
	save(results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
	results_omssa_NX_annot_disc_missing_uniqueness_f <- results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep
	results_omssa_NX_annot_disc_missing_uniqueness_ff <- results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep	
}			

results_tandem_NX_annot_disc_missing_peptides_uniqueness <- read.csv2(paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_peptides_uniqueness.csv"), header = TRUE, sep = ",") 
if(dim(results_tandem_NX_annot_disc_missing_peptides_uniqueness)[1] != 0) {
	results_tandem_NX_annot_disc_missing_uniqueness <- results_tandem_NX_annot_disc_missing[paste(results_tandem_NX_annot_disc_missing$PeptideSeq) %in% paste(results_tandem_NX_annot_disc_missing_peptides_uniqueness$peptide), ]
	save(results_tandem_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_uniqueness.rda"))
	tmp <- unique(results_tandem_NX_annot_disc_missing_uniqueness[, c("sample", "NextprotID", "PeptideSeq")])
	tmp2 <- summaryBy(PeptideSeq~NextprotID + sample, data=tmp, FUN=length, keep.names=TRUE)
	colnames(tmp2)[3] <- "NofDiscPeptides"
	results_tandem_NX_annot_disc_missing_uniqueness <- merge(results_tandem_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("NextprotID", "sample"), by.y = c("NextprotID", "sample"), all.x = TRUE)
	results_tandem_NX_annot_disc_missing_uniqueness[is.na(results_tandem_NX_annot_disc_missing_uniqueness$NofDiscPeptides), "NofDiscPeptides"] <- 0
	results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep <- results_tandem_NX_annot_disc_missing_uniqueness[results_tandem_NX_annot_disc_missing_uniqueness$NofDiscPeptides > 1,]
	save(results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
	results_tandem_NX_annot_disc_missing_uniqueness_f <- unique(results_tandem_NX_annot_disc_missing_uniqueness[, c("PeptideSeq", "NextprotID", "search_engine", "Chr", "sample", "NofDiscPeptides")])
	results_tandem_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ PeptideSeq + NextprotID + search_engine + Chr + NofDiscPeptides, data = results_tandem_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")

} else  {
	cat("Tandem has not identified any missing proteins!:")
	results_tandem_NX_annot_disc_missing_uniqueness <- results_tandem_NX_annot_disc_missing[0,]
	save(results_tandem_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_uniqueness.rda"))
	results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep <- results_tandem_NX_annot_disc_missing_uniqueness
	save(results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
	results_tandem_NX_annot_disc_missing_uniqueness_f <- results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep
	results_tandem_NX_annot_disc_missing_uniqueness_ff <- results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep	
}

# venn diagrams 

pdf(file = paste0(projectDir, "Plots/results_Venn_MissingNXProteinsPerSearchEngine.pdf"), width = 10, height = 10)
compare4List(paste(results_comet_NX_annot_disc_missing_uniqueness_2UniqPep$NextprotID), paste(results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep$NextprotID), paste(results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep$NextprotID), paste(results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep$NextprotID), "Comet results", "Mascot results", "Omssa results", "Tandem results",  "ProtFDRGlobal - Missing proteins (NX) for each search engine")
dev.off()

pdf(file = paste0(projectDir, "Plots/results_Venn_PeptidesOfMissingNXProteinsPerSearchEngine.pdf"), width = 10, height = 10)
compare4List(paste(results_comet_NX_annot_disc_missing_uniqueness_2UniqPep$PeptideSeq), paste(results_mascot_NX_annot_disc_missing_uniqueness_2UniqPep$PeptideSeq), paste(results_omssa_NX_annot_disc_missing_uniqueness_2UniqPep$PeptideSeq), paste(results_tandem_NX_annot_disc_missing_uniqueness_2UniqPep$PeptideSeq), "Comet results", "Mascot results", "Omssa results", "Tandem results",  "ProtFDRGlobal - Peptides of missing proteins (NX) for each search engine")
dev.off()

# print data to table

results_NX <- rbind(results_comet_NX_annot_disc_missing_uniqueness_f, results_mascot_NX_annot_disc_missing_uniqueness_f, results_omssa_NX_annot_disc_missing_uniqueness_f, results_tandem_NX_annot_disc_missing_uniqueness_f)

write.table(results_NX, file = paste0(projectDir, "Results/results_NX.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)




