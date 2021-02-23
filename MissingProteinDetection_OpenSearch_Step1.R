#!/usr/bin/Rscript
###preparado para hacer dataset 1 a 1

args=(commandArgs(TRUE))
projectDir <- args[1]		# /mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez.69/28_PFortes_Shotgun_lncRNA_Feb18/
fastafileroot <- args[2]	# /mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez.69/28_PFortes_Shotgun_lncRNA_Feb18/db
nextprotfolder <- args[3]	# /mnt/beegfs/jgonzalez.69/dato-activo/03_Analysis/jgonzalez.69/23_neXtprot_20170801_Nov2017/

library(Biostrings)
library(doBy)
library(reshape2)
library(ggplot2)
library(stringr)
source("/mnt/beegfs/jgonzalez.69/dato-activo/01_Rscripts/A_Funciones/funcionesShotgun.R")
source("/mnt/beegfs/jgonzalez.69/dato-activo/01_Rscripts/A_Funciones/funcionesVikv2.R")

results_MSFragger<-read.csv2(paste0(projectDir, "psm.tsv"),sep="\t",header=T, fill=T, na.strings='')


db_peptidesXProtAll_filterAA <- get(load(paste0(fastafileroot, "_peptidesXProtAll_filterAA.rda")))

tmp_nextprot <- unique(db_peptidesXProtAll_filterAA[, c("Peptide", "NextprotID")])

results_MSFragger$Peptide <- paste(results_MSFragger$Peptide)
results_MSFragger<-unique(results_MSFragger)
results_MSFragger_NX <- merge(results_MSFragger,tmp_nextprot, by = "Peptide", all.x = TRUE)
results_MSFragger_NX <-unique(results_MSFragger_NX)
results_MSFragger_NX <- results_MSFragger_NX[!(is.na(results_MSFragger_NX$NextprotID)),]

dir.create(paste0(projectDir,"results"))

save(results_MSFragger_NX, file = paste0(projectDir, "results/results_MSFragger_NX.rda"))

# discriminant peptides

db_peptidesXProtAll_filterAA_disc <- get(load(paste0(fastafileroot,"_peptidesXProtAll_filterAA_disc.rda")))

load(paste0(nextprotfolder, "nextProtXChrXENSP.RData"))

results_MSFragger_NX_annot <- merge(results_MSFragger_NX, unique(nextProtXChrXENSP[, c("NextprotID", "Missing", "Chr")]), by.x = "NextprotID", by.y = "NextprotID")
results_MSFragger_NX_annot$discriminant <- 0
results_MSFragger_NX_annot[paste(results_MSFragger_NX_annot$Peptide) %in% paste(db_peptidesXProtAll_filterAA_disc$Peptide), "discriminant"] <- 1
results_MSFragger_NX_annot_disc <- results_MSFragger_NX_annot[results_MSFragger_NX_annot$discriminant == 1, ]

save(results_MSFragger_NX_annot_disc, file = paste0(projectDir, "results/results_MSFragger_NX_annot_disc.rda"))

results_MSFragger_NX_annot_disc_missing <- results_MSFragger_NX_annot_disc[results_MSFragger_NX_annot_disc$Missing == 1, ]

save(results_MSFragger_NX_annot_disc_missing, file = paste0(projectDir, "results/results_MSFragger_NX_annot_disc_missing.rda"))

write.table(unique(results_MSFragger_NX_annot_disc_missing$Peptide), file = paste0(projectDir, "results/results_MSFragger_NX_annot_disc_missing_peptides.txt"), col.names=FALSE, row.names = FALSE, sep="\t", quote=FALSE)


#### txt to peptide uniqueness checker
