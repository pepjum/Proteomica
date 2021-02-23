#!/usr/bin/Rscript

#####EJECUTAR MANUALMENTE. LiNEA 24 del codigo modificar a mano o parametrizar

args=(commandArgs(TRUE))
projectDir <- args[1]		# /mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez.69/28_PFortes_Shotgun_lncRNA_Feb18/
fastafileroot <- args[2]	# /mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez.69/28_PFortes_Shotgun_lncRNA_Feb18/nextProtDB20170801
nextprotfolder <- args[3]	# /mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez.69/23_neXtprot_20170801_Nov17/
sample<-args[4] # NAME OF THE SAMPLE dataset

library(Biostrings)
library(doBy)
library(reshape2)
library(ggplot2)
library(stringr)
source("/mnt/beegfs/jgonzalez.69/dato-activo/01_Rscripts/A_Funciones/funcionesShotgun.R")
source("/mnt/beegfs/jgonzalez.69/dato-activo/01_Rscripts/A_Funciones/funcionesVikv2.R")

load(paste0(projectDir, "results/results_MSFragger_NX_annot_disc_missing.rda"))

results_MSFragger_NX_annot_disc_missing_peptides_uniqueness <- read.csv2(paste0(projectDir, "results/results_MSFragger_annot_disc_missing_peptides_uniqueness.csv"), header = TRUE, sep = ",")
if(dim(results_MSFragger_NX_annot_disc_missing_peptides_uniqueness)[1] != 0) {
	results_MSFragger_NX_annot_disc_missing_uniqueness <- results_MSFragger_NX_annot_disc_missing[paste(results_MSFragger_NX_annot_disc_missing$Peptide) %in% paste(results_MSFragger_NX_annot_disc_missing_peptides_uniqueness$peptide), ]
	save(results_MSFragger_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "results/results_MSFragger_NX_annot_disc_missing_uniqueness.rda"))
    results_MSFragger_NX_annot_disc_missing_uniqueness$sample<-sample
    ###ojo con la linea anterior.Cambiar para cada dataset
	tmp <- unique(results_MSFragger_NX_annot_disc_missing_uniqueness[, c("sample", "NextprotID", "Peptide")])
	tmp2 <- summaryBy(Peptide~NextprotID + sample, data=tmp, FUN=length, keep.names=TRUE)
	colnames(tmp2)[3] <- "NofDiscPeptides"
	results_MSFragger_NX_annot_disc_missing_uniqueness <- merge(results_MSFragger_NX_annot_disc_missing_uniqueness, tmp2, by.x = c("NextprotID", "sample"), by.y = c("NextprotID", "sample"), all.x = TRUE)
	results_MSFragger_NX_annot_disc_missing_uniqueness[is.na(results_MSFragger_NX_annot_disc_missing_uniqueness$NofDiscPeptides), "NofDiscPeptides"] <- 0
	results_MSFragger_NX_annot_disc_missing_uniqueness_2UniqPep <- results_MSFragger_NX_annot_disc_missing_uniqueness[results_MSFragger_NX_annot_disc_missing_uniqueness$NofDiscPeptides > 1,]
	save(results_MSFragger_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "results/results_MSFragger_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
	results_MSFragger_NX_annot_disc_missing_uniqueness_f <- unique(results_MSFragger_NX_annot_disc_missing_uniqueness[, c("Peptide", "NextprotID", "Chr", "sample", "NofDiscPeptides","Protein")])
	results_MSFragger_NX_annot_disc_missing_uniqueness_ff <- aggregate(sample ~ Peptide + NextprotID + Chr + NofDiscPeptides + Protein, data = results_MSFragger_NX_annot_disc_missing_uniqueness_f, FUN = paste, collapse = ",")

} else {
	cat("MSFragger has not identified any missing proteins!:")
	results_MSFragger_NX_annot_disc_missing_uniqueness <- results_MSFragger_NX_annot_disc_missing[0,]
	save(results_MSFragger_NX_annot_disc_missing_uniqueness, file = paste0(projectDir, "results/results_MSFragger_NX_annot_disc_missing_uniqueness.rda"))
	results_MSFragger_NX_annot_disc_missing_uniqueness_2UniqPep <- results_MSFragger_NX_annot_disc_missing_uniqueness
	save(results_MSFragger_NX_annot_disc_missing_uniqueness_2UniqPep, file = paste0(projectDir, "results/results_MSFragger_NX_annot_disc_missing_uniqueness_2UniqPep.rda"))
	results_MSFragger_NX_annot_disc_missing_uniqueness_2UniqPep_f <- results_MSFragger_NX_annot_disc_missing_uniqueness_2UniqPep
	results_MSFragger_NX_annot_disc_missing_uniqueness_2UniqPep_ff <- results_MSFragger_NX_annot_disc_missing_uniqueness_2UniqPep
}


results_NX<-results_MSFragger_NX_annot_disc_missing_uniqueness_f
results_NX_ff<-	results_MSFragger_NX_annot_disc_missing_uniqueness_ff

spec<-results_MSFragger_NX_annot_disc_missing_uniqueness[(results_MSFragger_NX_annot_disc_missing_uniqueness$Peptide=="KPRPMGIIAANVEK"),]
spec2<-results_MSFragger_NX_annot_disc_missing_uniqueness[(results_MSFragger_NX_annot_disc_missing_uniqueness$Peptide=="GYGLEVDMWAAGVILYILLCGFPPFR"),]
write.table(results_NX, file = paste0(projectDir, "results/results_NX.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
write.table(results_NX_ff, file = paste0(projectDir, "results/results_NX_ff.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
write.table(spec, file = paste0(projectDir, "results/peptide1_DCLK3_HUMAN.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
write.table(spec2, file = paste0(projectDir, "results/peptide2_DCLK3_HUMAN.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
