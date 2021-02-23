#!/usr/bin/Rscript

args=(commandArgs(TRUE))
projectDir <- args[1]		# /home/nostromo/data/03_Analysis/agarin/28_PFortes_Shotgun_lncRNA_Feb18/
fastafileroot <- args[2]	# /home/nostromo/data/03_Analysis/agarin/28_PFortes_Shotgun_lncRNA_Feb18/nextProtDB20170801
nextprotfolder <- args[3]	# /home/nostromo/data/03_Analysis/agarin/23_neXtprot_20170801_Nov2017/

library(Biostrings)
library(doBy)
library(reshape2)
library(ggplot2)
library(stringr)
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesShotgun.R")
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesVikv2.R")


mascot_txtFiles <- list.files(path = paste0(mascot_rda_folder,dataset,"/percolator/"), pattern = "*pinresults_peptides.tsv.out")
comet_txtFiles <- list.files(path = paste0(comet_rda_folder,dataset,"/percolator/"), pattern = "*pinresults_peptides.tsv.out")
omssa_txtFiles <- list.files(path = paste0(omssa_rda_folder,dataset,"/percolator/"), pattern = "*pinresults_peptides.tsv.out")
tandem_txtFiles <- list.files(path = paste0(tandem_rda_folder,dataset,"/percolator/") = "*pinresults_peptides.tsv.out")


mascot_txtFiles<-paste0("/home/nostromo/data/pepe/EMBRIO_9_03/Dat_Files/PXD001381/percolator/",mascot_txtFiles)
comet_txtFiles<-paste0("/home/nostromo/data/pepe/EMBRIO_9_03/Omssa_files/PXD001381/percolator/",comet_txtFiles)
omssa_txtFiles<-paste0("/home/nostromo/data/pepe/EMBRIO_9_03/Omssa_files/PXD001381/percolator/",omssa_txtFiles)
tandem_txtFiles<-paste0("/home/nostromo/data/pepe/EMBRIO_9_03/Tandem_Files/PXD001381/percolator/",tandem_txtFiles)


mascot_file<-read.table(paste(mascot_txtFiles, row.names=NULL)
mascot_file$number<-seq(1:nrow(mascot_file))
names(mascot_file)<-c("PSMId","score","q.value","posterior_error_prob","peptide","proteinIds","Query","number")
PATTERN <- '\\[[0-9.-]+\\]'
PATTERN2 <- '[A-Z]{2,}'

mascot_file$PeptideSeq<-str_extract(str_remove_all(mascot_file$peptide, PATTERN), PATTERN2)

comet_file<-read.table(paste(comet_txtFiles, header=T)
comet_file$PeptideSeq<-str_extract(str_remove_all(comet_file$peptide, PATTERN), PATTERN2)

tandem_file<-read.table(paste(tandem_txtFiles, header=T)
tandem_file$PeptideSeq<-str_extract(str_remove_all(tandem_file$peptide, PATTERN), PATTERN2)
tmp_id<-paste(unlist(tandem_file$proteinIds))
tmp_id<-sapply(strsplit(tmp_id,"-"), "[", 1)
tandem_file$proteinIds<-tmp_id


results_comet<-comet_file
results_mascot<-mascot_file
results_tandem<-tandem_file

results_mascot<-results_mascot[(results_mascot$q.value<0.01),]
results_comet<-results_comet[(results_comet$q.value<0.01),]
results_omssa<-results_omssa[(results_omssa$q.value<0.01),]
results_tandem<-results_tandem[(results_tandem$q.value<0.01),]


db_peptidesXProtAll_filterAA <- get(load(paste0(fastafileroot, "_peptidesXProtAll_filterAA.rda")))

tmp_nextprot <- unique(db_peptidesXProtAll_filterAA[, c("Peptide", "ProteinID")])
tmp_id<-paste(unlist(tmp_nextprot$ProteinID))
tmp_id<-sapply(strsplit(tmp_id,"-"), "[", 1)
tmp_nextprot$ProteinID<-tmp_id

results_comet$PeptideSeq <- paste(results_comet$PeptideSeq)
results_comet <- unique(results_comet)
# quito cualquier referencia a Swissprot
results_comet$number <- NULL
results_comet$peptide <- NULL
results_comet$proteinIds <- NULL
results_comet <- unique(results_comet)
results_comet_NX <- merge(results_comet,tmp_nextprot, by.x = "PeptideSeq", by.y = "Peptide", all.x = TRUE)
results_comet_NX <- unique(results_comet_NX)

results_mascot$PeptideSeq <- paste(results_mascot$PeptideSeq)
results_mascot <- unique(results_mascot)
# quito cualquier referencia a Swissprot

results_mascot$number <- NULL
results_mascot$peptide <- NULL
results_mascot$proteinIds <- NULL

results_mascot <- unique(results_mascot)
results_mascot_NX <- merge(results_mascot,tmp_nextprot, by.x = "PeptideSeq", by.y = "Peptide", all.x = TRUE)
results_mascot_NX <- unique(results_mascot_NX)

#projectDir<-"/home/nostromo/data/pepe/EMBRIO_5_04/"
results_omssa$PeptideSeq <- paste(results_omssa$PeptideSeq)
results_omssa <- unique(results_omssa)
# quito cualquier referencia a Swissprot
results_omssa$number <- NULL
results_omssa$peptide <- NULL
results_omssa$proteinIds <- NULL

results_omssa <- unique(results_omssa)
results_omssa_NX <- merge(results_omssa,tmp_nextprot, by.x = "PeptideSeq", by.y = "Peptide", all.x = TRUE)
results_omssa_NX <- unique(results_omssa_NX)

results_tandem$PeptideSeq <- paste(results_tandem$PeptideSeq)
results_tandem <- unique(results_tandem)
# quito cualquier referencia a Swissprot
results_tandem$number <- NULL
results_tandem$peptide <- NULL
results_tandem$proteinIds <- NULL
results_tandem <- unique(results_tandem)
results_tandem_NX <- merge(results_tandem,tmp_nextprot, by.x = "PeptideSeq", by.y = "Peptide", all.x = TRUE)
results_tandem_NX <- unique(results_tandem_NX)

results_comet_NX <- results_comet_NX[!(is.na(results_comet_NX$ProteinID)),]
results_mascot_NX <- results_mascot_NX[!(is.na(results_mascot_NX$ProteinID)),]
results_omssa_NX <- results_omssa_NX[!(is.na(results_omssa_NX$ProteinID)),]
results_tandem_NX <- results_tandem_NX[!(is.na(results_tandem_NX$ProteinID)),]

results_comet_NX$PSM <- paste(results_comet_NX$datfile, results_comet_NX$Query, sep = "-")
results_mascot_NX$PSM <- paste(results_mascot_NX$datfile, results_mascot_NX$Query, sep = "-")
results_omssa_NX$PSM <- paste(results_omssa_NX$datfile, results_omssa_NX$Query, sep = "-")
results_tandem_NX$PSM <- paste(results_tandem_NX$datfile, results_tandem_NX$Query, sep = "-")

####cosecha del pepe

tmp_id<-paste(unlist(results_comet_NX$ProteinID))
tmp_id<-sapply(strsplit(tmp_id,"-"), "[", 1)
results_comet_NX$NextprotID<-tmp_id

tmp_id<-paste(unlist(results_mascot_NX$ProteinID))
tmp_id<-sapply(strsplit(tmp_id,"-"), "[", 1)
results_mascot_NX$NextprotID<-tmp_id

tmp_id<-paste(unlist(results_omssa_NX$ProteinID))
tmp_id<-sapply(strsplit(tmp_id,"-"), "[", 1)
results_omssa_NX$NextprotID<-tmp_id

tmp_id<-paste(unlist(results_tandem_NX$ProteinID))
tmp_id<-sapply(strsplit(tmp_id,"-"), "[", 1)
results_tandem_NX$NextprotID<-tmp_id


save(results_comet_NX, file = paste0(projectDir, "Results/results_comet_NX_percolator.rda"))
save(results_mascot_NX, file = paste0(projectDir, "Results/results_mascot_NX_percolator.rda"))
save(results_omssa_NX, file = paste0(projectDir, "Results/results_omssa_NX_percolator.rda"))
save(results_tandem_NX, file = paste0(projectDir, "Results/results_tandem_NX_percolator.rda"))

# Venn diagrams

pdf(file = paste0(projectDir, "Plots/results_Venn_NXProteinsPerSearchEngine.pdf"), width = 10, height = 10)
compare4List(paste(results_comet_NX$ProteinID), paste(results_mascot_NX$ProteinID), paste(results_omssa_NX$ProteinID), paste(results_tandem_NX$ProteinID), "Comet results", "Mascot results", "Omssa results", "Tandem results",  "Proteins (Nextprot) for each search engine")
dev.off()

# discriminant peptides

db_peptidesXProtAll_filterAA_disc <- get(load(paste0(fastafileroot,"_peptidesXProtAll_filterAA_disc.rda")))
load(paste0(nextprotfolder, "nextProtXChrXENSP.RData"))

results_comet_NX_annot <- merge(results_comet_NX, unique(nextProtXChrXENSP[, c("NextprotID", "Missing", "Chr")]), by.x = "ProteinID", by.y = "NextprotID")
results_comet_NX_annot$discriminant <- 0
results_comet_NX_annot[paste(results_comet_NX_annot$PeptideSeq) %in% paste(db_peptidesXProtAll_filterAA_disc$Peptide), "discriminant"] <- 1
results_comet_NX_annot_disc <- results_comet_NX_annot[results_comet_NX_annot$discriminant == 1, ]

results_mascot_NX_annot <- merge(results_mascot_NX, unique(nextProtXChrXENSP[, c("NextprotID", "Missing", "Chr")]), by.x = "ProteinID", by.y = "NextprotID")
results_mascot_NX_annot$discriminant <- 0
results_mascot_NX_annot[paste(results_mascot_NX_annot$PeptideSeq) %in% paste(db_peptidesXProtAll_filterAA_disc$Peptide), "discriminant"] <- 1
results_mascot_NX_annot_disc <- results_mascot_NX_annot[results_mascot_NX_annot$discriminant == 1, ]

results_omssa_NX_annot <- merge(results_omssa_NX, unique(nextProtXChrXENSP[, c("NextprotID", "Missing", "Chr")]), by.x = "ProteinID", by.y = "NextprotID")
results_omssa_NX_annot$discriminant <- 0
results_omssa_NX_annot[paste(results_omssa_NX_annot$PeptideSeq) %in% paste(db_peptidesXProtAll_filterAA_disc$Peptide), "discriminant"] <- 1
results_omssa_NX_annot_disc <- results_omssa_NX_annot[results_omssa_NX_annot$discriminant == 1, ]

results_tandem_NX_annot <- merge(results_tandem_NX, unique(nextProtXChrXENSP[, c("NextprotID", "Missing", "Chr")]), by.x = "ProteinID", by.y = "NextprotID")
results_tandem_NX_annot$discriminant <- 0
results_tandem_NX_annot[paste(results_tandem_NX_annot$PeptideSeq) %in% paste(db_peptidesXProtAll_filterAA_disc$Peptide), "discriminant"] <- 1
results_tandem_NX_annot_disc <- results_tandem_NX_annot[results_tandem_NX_annot$discriminant == 1, ]

save(results_comet_NX_annot, file = paste0(projectDir, "Results/results_comet_NX_annot_percolator.rda"))
save(results_comet_NX_annot_disc, file = paste0(projectDir, "Results/results_comet_NX_annot_disc_percolator.rda"))
save(results_mascot_NX_annot, file = paste0(projectDir, "Results/results_mascot_NX_annot_percolator.rda"))
save(results_mascot_NX_annot_disc, file = paste0(projectDir, "Results/results_mascot_NX_annot_disc_percolator.rda"))
save(results_omssa_NX_annot, file = paste0(projectDir, "Results/results_omssa_NX_annot_percolator.rda"))
save(results_omssa_NX_annot_disc, file = paste0(projectDir, "Results/results_omssa_NX_annot_disc_percolator.rda"))
save(results_tandem_NX_annot, file = paste0(projectDir, "Results/results_tandem_NX_annot_percolator.rda"))
save(results_tandem_NX_annot_disc, file = paste0(projectDir, "Results/results_tandem_NX_annot_disc_percolator.rda"))

results_comet_NX_annot_disc_missing <- results_comet_NX_annot_disc[results_comet_NX_annot_disc$Missing == 1, ]
results_mascot_NX_annot_disc_missing <- results_mascot_NX_annot_disc[results_mascot_NX_annot_disc$Missing == 1, ]
results_omssa_NX_annot_disc_missing <- results_omssa_NX_annot_disc[results_omssa_NX_annot_disc$Missing == 1, ]
results_tandem_NX_annot_disc_missing <- results_tandem_NX_annot_disc[results_tandem_NX_annot_disc$Missing == 1, ]

save(results_comet_NX_annot_disc_missing, file = paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_percolator.rda"))
save(results_mascot_NX_annot_disc_missing, file = paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_percolator.rda"))
save(results_omssa_NX_annot_disc_missing, file = paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_percolator.rda"))
save(results_tandem_NX_annot_disc_missing, file = paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_percolator.rda"))

write.table(unique(results_comet_NX_annot_disc_missing$PeptideSeq), file = paste0(projectDir, "Results/results_comet_NX_annot_disc_missing_peptides_percolator.txt"), col.names=FALSE, row.names = FALSE, sep="\t", quote=FALSE)

write.table(unique(results_mascot_NX_annot_disc_missing$PeptideSeq), file = paste0(projectDir, "Results/results_mascot_NX_annot_disc_missing_peptides_percolator.txt"), col.names=FALSE, row.names = FALSE, sep="\t", quote=FALSE)

write.table(unique(results_omssa_NX_annot_disc_missing$PeptideSeq), file = paste0(projectDir, "Results/results_omssa_NX_annot_disc_missing_peptides_percolator.txt"), col.names=FALSE, row.names = FALSE, sep="\t", quote=FALSE)

write.table(unique(results_tandem_NX_annot_disc_missing$PeptideSeq), file = paste0(projectDir, "Results/results_tandem_NX_annot_disc_missing_peptides_percolator.txt"), col.names=FALSE, row.names = FALSE, sep="\t", quote=FALSE)
