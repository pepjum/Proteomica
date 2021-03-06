#!/usr/bin/Rscript

args=(commandArgs(TRUE))
projectDir <- args[1]
databasename <- args[2]

library(Biostrings)
library(doBy)
library(reshape2)
library(ggplot2)
library(stringr)
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesShotgun.R")
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesVikv2.R")

db_peptidesXProtAll_filterAA_disc <- get(load(paste0(projectDir, "/", databasename, "_peptidesXProtAll_filterAA_disc.rda")))

mascot_rda_folder <- paste0(projectDir, "/Dat_Files/")
comet_rda_folder <- paste0(projectDir, "/Comet_files/")
omssa_rda_folder <- paste0(projectDir, "/Omssa_files/")
tandem_rda_folder <- paste0(projectDir, "/Tandem_Files/")

mascot_rdaFiles <- list.files(path = mascot_rda_folder, pattern = "*protFDR_Filter.rda")
comet_rdaFiles <- list.files(path = comet_rda_folder, pattern = "*protFDR_Filter.rda")
omssa_rdaFiles <- list.files(path = omssa_rda_folder, pattern = "*protFDR_Filter.rda")
tandem_rdaFiles <- list.files(path = tandem_rda_folder, pattern = "*protFDR_Filter.rda")

mascot_rdaFiles <- paste0(mascot_rda_folder,mascot_rdaFiles)
comet_rdaFiles <- paste0(comet_rda_folder,comet_rdaFiles)
omssa_rdaFiles <- paste0(omssa_rda_folder,omssa_rdaFiles)
tandem_rdaFiles <- paste0(tandem_rda_folder,tandem_rdaFiles)

resultDir <- paste0(projectDir, "/Results/")
y <- paste0("mkdir -p ", resultDir)
system(y)
y <- NULL

plotDir <- paste0(projectDir, "/Plots/")
y <- paste0("mkdir -p ", plotDir)
system(y)
y <- NULL

for (i in 1:length(mascot_rdaFiles)) {
	load(paste(mascot_rdaFiles[i]))
	cat(paste(mascot_rdaFiles[i]))
	cat("\n")
	cat(dim(dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter))
	cat("\n")
	if(i == 1) {
		results_mascot_TyD <- dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter
	} else {
		results_mascot_TyD <- rbind(results_mascot_TyD, dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter)
	}
}
save(results_mascot_TyD, file = paste0(resultDir, "results_mascot_TyD.rda"))

for (i in 1:length(comet_rdaFiles)) {
	load(paste(comet_rdaFiles[i]))
	cat(paste(comet_rdaFiles[i]))
	cat("\n")
	cat(dim(dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter))
	cat("\n")

	if(i == 1) {
		results_comet_TyD <- dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter
	} else {
		results_comet_TyD <- rbind(results_comet_TyD, dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter)
	}
}
save(results_comet_TyD, file = paste0(resultDir, "results_comet_TyD.rda"))

for (i in 1:length(omssa_rdaFiles)) {
	load(paste(omssa_rdaFiles[i]))
	cat(paste(omssa_rdaFiles[i]))
	cat("\n")
	cat(dim(dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter))
	cat("\n")

	if(i == 1) {
		results_omssa_TyD <- dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter
	} else {
		results_omssa_TyD <- rbind(results_omssa_TyD, dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter)
	}
}
save(results_omssa_TyD, file = paste0(resultDir, "results_omssa_TyD.rda"))

for (i in 1:length(tandem_rdaFiles)) {
	load(paste(tandem_rdaFiles[i]))
	cat(paste(tandem_rdaFiles[i]))
	cat("\n")
	cat(dim(dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter))

	cat("\n")
	if(i == 1) {
		results_tandem_TyD <- dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter
	} else {
		results_tandem_TyD <- rbind(results_tandem_TyD, dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter)
	}
}
save(results_tandem_TyD, file = paste0(resultDir, "results_tandem_TyD.rda"))

# check values
cat("Results Mascot:")
cat("\n")
cat("Range PSMFDR:", range(as.numeric(paste(results_mascot_TyD$psmFDR))))
cat("\n")
cat("Range pepFDR:", range(as.numeric(paste(results_mascot_TyD$pepFDR))))
cat("\n")
cat("Range protFDR:", range(as.numeric(paste(results_mascot_TyD$protFDR))))
cat("\n")
cat("Number of samples:", length(table(results_mascot_TyD$sample)))
cat("\n")
cat("Samples:", length(unique(results_mascot_TyD$sample)))
cat("\n")
cat("Number of target and decoy values:", table(results_mascot_TyD$database))
cat("\n")
cat("Example of target protein names:", paste(results_mascot_TyD[results_mascot_TyD$database == "T", "ProteinAccession"][1:5]))
cat("\n")
cat("Example of decoy protein names:", paste(results_mascot_TyD[results_mascot_TyD$database == "D", "ProteinAccession"][1:5]))
cat("\n")


cat("Results Comet:")
cat("\n")
cat("Range PSMFDR:", range(as.numeric(paste(results_comet_TyD$psmFDR))))
cat("\n")
cat("Range pepFDR:", range(as.numeric(paste(results_comet_TyD$pepFDR))))
cat("\n")
cat("Range protFDR:", range(as.numeric(paste(results_comet_TyD$protFDR))))
cat("\n")
cat("Number of samples:", length(table(results_comet_TyD$sample)))
cat("\n")
cat("Samples:", length(unique(results_comet_TyD$sample)))
cat("\n")
cat("Number of target and decoy values:", table(results_comet_TyD$database))
cat("\n")
cat("Example of target protein names:", paste(results_comet_TyD[results_comet_TyD$database == "T", "ProteinAccession"][1:5]))
cat("\n")
cat("Example of decoy protein names:", paste(results_comet_TyD[results_comet_TyD$database == "D", "ProteinAccession"][1:5]))
cat("\n")

cat("Results Omssa:")
cat("\n")
cat("Range PSMFDR:", range(as.numeric(paste(results_omssa_TyD$psmFDR))))
cat("\n")
cat("Range pepFDR:", range(as.numeric(paste(results_omssa_TyD$pepFDR))))
cat("\n")
cat("Range protFDR:", range(as.numeric(paste(results_omssa_TyD$protFDR))))
cat("\n")
cat("Number of samples:", length(table(results_omssa_TyD$sample)))
cat("\n")
cat("Samples:", length(unique(results_omssa_TyD$sample)))
cat("\n")
cat("Number of target and decoy values:", table(results_omssa_TyD$database))
cat("\n")
cat("Example of target protein names:", paste(results_omssa_TyD[results_omssa_TyD$database == "T", "ProteinAccession"][1:5]))
cat("\n")
cat("Example of decoy protein names:", paste(results_omssa_TyD[results_omssa_TyD$database == "D", "ProteinAccession"][1:5]))
cat("\n")

cat("Results Tandem:")
cat("\n")
cat("Range PSMFDR:", range(as.numeric(paste(results_tandem_TyD$psmFDR))))
cat("\n")
cat("Range pepFDR:", range(as.numeric(paste(results_tandem_TyD$pepFDR))))
cat("\n")
cat("Range protFDR:", range(as.numeric(paste(results_tandem_TyD$protFDR))))
cat("\n")
cat("Number of samples:", length(table(results_tandem_TyD$sample)))
cat("\n")
cat("Samples:", length(unique(results_tandem_TyD$sample)))
cat("\n")
cat("Number of target and decoy values:", table(results_tandem_TyD$database))
cat("\n")
cat("Example of target protein names:", paste(results_tandem_TyD[results_tandem_TyD$database == "T", "ProteinAccession"][1:5]))
cat("\n")
cat("Example of decoy protein names:", paste(results_tandem_TyD[results_tandem_TyD$database == "D", "ProteinAccession"][1:5]))
cat("\n")

# remove contaminants
results_mascot_TyD <- results_mascot_TyD[!(grepl("CRAP", results_mascot_TyD$ProteinAccession)),]
results_comet_TyD <- results_comet_TyD[!(grepl("CRAP", results_comet_TyD$ProteinAccession)),]
results_omssa_TyD <- results_omssa_TyD[!(grepl("CRAP", results_omssa_TyD$ProteinAccession)),]
results_tandem_TyD <- results_tandem_TyD[!(grepl("CRAP", results_tandem_TyD$ProteinAccession)),]

save(results_mascot_TyD, file = paste0(resultDir, "results_mascot_TyD.rda"))
save(results_comet_TyD, file = paste0(resultDir, "results_comet_TyD.rda"))
save(results_omssa_TyD, file = paste0(resultDir, "results_omssa_TyD.rda"))
save(results_tandem_TyD, file = paste0(resultDir, "results_tandem_TyD.rda"))

# remove decoy identifications

results_comet <- results_comet_TyD[results_comet_TyD$database == "T",]
results_mascot <- results_mascot_TyD[results_mascot_TyD$database == "T",]
results_omssa <- results_omssa_TyD[results_omssa_TyD$database == "T",]
results_tandem <- results_tandem_TyD[results_tandem_TyD$database == "T",]

results_comet$PSM <- paste(results_comet$datfile, results_comet$Query, sep = "-")
results_mascot$PSM <- paste(results_mascot$datfile, results_mascot$Query, sep = "-")
results_omssa$PSM <- paste(results_omssa$datfile, results_omssa$Query, sep = "-")
results_tandem$PSM <- paste(results_tandem$datfile, results_tandem$Query, sep = "-")

save(results_comet, file = paste0(resultDir, "results_comet.rda"))
save(results_mascot, file = paste0(resultDir, "results_mascot.rda"))
save(results_omssa, file = paste0(resultDir, "results_omssa.rda"))
save(results_tandem, file = paste0(resultDir, "results_tandem.rda"))

# Discriminant peptides

results_comet_annot <- results_comet
results_comet_annot$discriminant <- 0
results_comet_annot[paste(results_comet_annot$PeptideSeq) %in% paste(db_peptidesXProtAll_filterAA_disc$Peptide), "discriminant"] <- 1
results_comet_annot_disc <- results_comet_annot[results_comet_annot$discriminant == 1, ]

results_mascot_annot <- results_mascot
results_mascot_annot$discriminant <- 0
results_mascot_annot[paste(results_mascot_annot$PeptideSeq) %in% paste(db_peptidesXProtAll_filterAA_disc$Peptide), "discriminant"] <- 1
results_mascot_annot_disc <- results_mascot_annot[results_mascot_annot$discriminant == 1, ]

results_omssa_annot <- results_omssa
results_omssa_annot$discriminant <- 0
results_omssa_annot[paste(results_omssa_annot$PeptideSeq) %in% paste(db_peptidesXProtAll_filterAA_disc$Peptide), "discriminant"] <- 1
results_omssa_annot_disc <- results_omssa_annot[results_omssa_annot$discriminant == 1, ]

results_tandem_annot <- results_tandem
results_tandem_annot$discriminant <- 0
results_tandem_annot[paste(results_tandem_annot$PeptideSeq) %in% paste(db_peptidesXProtAll_filterAA_disc$Peptide), "discriminant"] <- 1
results_tandem_annot_disc <- results_tandem_annot[results_tandem_annot$discriminant == 1, ]

save(results_comet_annot, file = paste0(resultDir, "results_comet_annot.rda"))
save(results_comet_annot_disc, file = paste0(resultDir, "results_comet_annot_disc.rda"))
save(results_mascot_annot, file = paste0(resultDir, "results_mascot_annot.rda"))
save(results_mascot_annot_disc, file = paste0(resultDir, "results_mascot_annot_disc.rda"))
save(results_omssa_annot, file = paste0(resultDir, "results_omssa_annot.rda"))
save(results_omssa_annot_disc, file = paste0(resultDir, "results_omssa_annot_disc.rda"))
save(results_tandem_annot, file = paste0(resultDir, "results_tandem_annot.rda"))
save(results_tandem_annot_disc, file = paste0(resultDir, "results_tandem_annot_disc.rda"))

# Select proteins with at least 2 unique peptides

tmp <- unique(results_comet_annot_disc[, c("sample", "ProteinAccession", "PeptideSeq")])
tmp2 <- summaryBy(PeptideSeq~ProteinAccession + sample, data=tmp, FUN=length, keep.names=TRUE)
colnames(tmp2)[3] <- "NofDiscPeptides"
results_comet_annot_disc <- merge(results_comet_annot_disc, tmp2, by.x = c("ProteinAccession", "sample"), by.y = c("ProteinAccession", "sample"), all.x = TRUE)
results_comet_annot_disc[is.na(results_comet_annot_disc$NofDiscPeptides), "NofDiscPeptides"] <- 0

tmp <- unique(results_mascot_annot_disc[, c("sample", "ProteinAccession", "PeptideSeq")])
tmp2 <- summaryBy(PeptideSeq~ProteinAccession + sample, data=tmp, FUN=length, keep.names=TRUE)
colnames(tmp2)[3] <- "NofDiscPeptides"
results_mascot_annot_disc <- merge(results_mascot_annot_disc, tmp2, by.x = c("ProteinAccession", "sample"), by.y = c("ProteinAccession", "sample"), all.x = TRUE)
results_mascot_annot_disc[is.na(results_mascot_annot_disc$NofDiscPeptides), "NofDiscPeptides"] <- 0

tmp <- unique(results_omssa_annot_disc[, c("sample", "ProteinAccession", "PeptideSeq")])
tmp2 <- summaryBy(PeptideSeq~ProteinAccession + sample, data=tmp, FUN=length, keep.names=TRUE)
colnames(tmp2)[3] <- "NofDiscPeptides"
results_omssa_annot_disc <- merge(results_omssa_annot_disc, tmp2, by.x = c("ProteinAccession", "sample"), by.y = c("ProteinAccession", "sample"), all.x = TRUE)
results_omssa_annot_disc[is.na(results_omssa_annot_disc$NofDiscPeptides), "NofDiscPeptides"] <- 0

tmp <- unique(results_tandem_annot_disc[, c("sample", "ProteinAccession", "PeptideSeq")])
tmp2 <- summaryBy(PeptideSeq~ProteinAccession + sample, data=tmp, FUN=length, keep.names=TRUE)
colnames(tmp2)[3] <- "NofDiscPeptides"
results_tandem_annot_disc <- merge(results_tandem_annot_disc, tmp2, by.x = c("ProteinAccession", "sample"), by.y = c("ProteinAccession", "sample"), all.x = TRUE)
results_tandem_annot_disc[is.na(results_tandem_annot_disc$NofDiscPeptides), "NofDiscPeptides"] <- 0

results_comet_annot_disc_2UniqPep <- results_comet_annot_disc[results_comet_annot_disc$NofDiscPeptides > 1,]
results_mascot_annot_disc_2UniqPep <- results_mascot_annot_disc[results_mascot_annot_disc$NofDiscPeptides > 1,]
results_omssa_annot_disc_2UniqPep <- results_omssa_annot_disc[results_omssa_annot_disc$NofDiscPeptides > 1,]
results_tandem_annot_disc_2UniqPep <- results_tandem_annot_disc[results_tandem_annot_disc$NofDiscPeptides > 1,]

cat("ya estoy aqui")

save(results_comet_annot_disc_2UniqPep, file = paste0(resultDir, "results_comet_annot_disc_2UniqPep.rda"))
save(results_mascot_annot_disc_2UniqPep, file = paste0(resultDir, "results_mascot_annot_disc_2UniqPep.rda"))
results_omssa_annot_disc_2UniqPep$ProteinAccession<-sapply(strsplit(paste(results_omssa_annot_disc_2UniqPep$ProteinAccession)," "),"[",1)
save(results_omssa_annot_disc_2UniqPep, file = paste0(resultDir, "results_omssa_annot_disc_2UniqPep.rda"))
save(results_tandem_annot_disc_2UniqPep, file = paste0(resultDir, "results_tandem_annot_disc_2UniqPep.rda"))

# Venn diagrams

plotDir<-"/home/nostromo/data/pepe/09_EMBRIO_Abril18/Plots/"
databasename<-"2018_02_uniprot"
pdf(file = paste0(plotDir, "results_Venn_ProteinsPerSearchEngine.pdf"), width = 10, height = 10)
compare4List(paste(results_comet_annot_disc_2UniqPep$ProteinAccession), paste(results_mascot_annot_disc_2UniqPep$ProteinAccession), paste(results_omssa_annot_disc_2UniqPep$ProteinAccession), paste(results_tandem_annot_disc_2UniqPep$ProteinAccession), "Comet results", "Mascot results", "Omssa results", "Tandem results",  paste0("Proteins (", databasename, ") identified with 2 unique peptides. Results for each search engine."))
dev.off()

pdf(file = paste0(plotDir, "results_Venn_PeptidesPerSearchEngine.pdf"), width = 10, height = 10)
compare4List(paste(results_comet_annot_disc_2UniqPep$PeptideSeq), paste(results_mascot_annot_disc_2UniqPep$PeptideSeq), paste(results_omssa_annot_disc_2UniqPep$PeptideSeq), paste(results_tandem_annot_disc_2UniqPep$PeptideSeq), "Comet results", "Mascot results", "Omssa results", "Tandem results",  paste0("Peptides of proteins (", databasename, ") identified with 2 unique peptides. Results for each search engine"))
dev.off()

# Plot - peptides per sample - matrix

results_ALLSE_annot_disc_uniqueness_2UniqPep_f <- rbind(unique(results_comet_annot_disc_2UniqPep[, c("ProteinAccession","PeptideSeq","sample")]), unique(results_mascot_annot_disc_2UniqPep[, c("ProteinAccession","PeptideSeq","sample")]), unique(results_omssa_annot_disc_2UniqPep[, c("ProteinAccession","PeptideSeq","sample")]), unique(results_tandem_annot_disc_2UniqPep[, c("ProteinAccession","PeptideSeq","sample")]))

samples <- paste(unique(results_ALLSE_annot_disc_uniqueness_2UniqPep_f$sample))
proteins <- paste(unique(results_ALLSE_annot_disc_uniqueness_2UniqPep_f$ProteinAccession))
proteinsPerSamplesMatrix <- matrix(, nrow = length(proteins), ncol = length(samples))
#each column is a variation
for(i in 1:length(proteins)) {
	cat(proteins[i],"\n")
	for(j in 1:length(samples)){
		tmp <- results_ALLSE_annot_disc_uniqueness_2UniqPep_f[results_ALLSE_annot_disc_uniqueness_2UniqPep_f$ProteinAccession == proteins[i], ]
		proteinsPerSamplesMatrix[i,j] <- length(unique(tmp[(paste(tmp$sample) == samples[j]), "PeptideSeq"]))
	}
}
rownames(proteinsPerSamplesMatrix) <- paste(proteins)
colnames(proteinsPerSamplesMatrix) <- paste(samples)

resultDir<-"/home/nostromo/data/pepe/09_EMBRIO_Abril18/Results/"
write.table(proteinsPerSamplesMatrix, file = paste0(resultDir, "results_ALLSE_annot_disc_2UniqPep_nPeptidesPerProteinAndPerSample.txt"), col.names=TRUE, row.names = TRUE, sep="\t", quote=FALSE)

#### COMET

results_comet_annot_disc_2UniqPep_f <- unique(results_comet_annot_disc_2UniqPep[, c("ProteinAccession","PeptideSeq","sample")])

samples <- paste(unique(results_comet_annot_disc_2UniqPep_f$sample))
proteins <- paste(unique(results_comet_annot_disc_2UniqPep_f$ProteinAccession))
proteinsPerSamplesMatrix <- matrix(, nrow = length(proteins), ncol = length(samples))
#each column is a variation
for(i in 1:length(proteins)) {
	cat("comet ",proteins[i],"\n")
	for(j in 1:length(samples)){
		tmp <- results_comet_annot_disc_2UniqPep_f[results_comet_annot_disc_2UniqPep_f$ProteinAccession == proteins[i], ]
		proteinsPerSamplesMatrix[i,j] <- length(unique(tmp[(paste(tmp$sample) == samples[j]), "PeptideSeq"]))
	}
}
rownames(proteinsPerSamplesMatrix) <- paste(proteins)
colnames(proteinsPerSamplesMatrix) <- paste(samples)

write.table(proteinsPerSamplesMatrix, file = paste0(resultDir, "results_comet_annot_disc_2UniqPep_nPeptidesPerProteinAndPerSample.txt"), col.names=TRUE, row.names = TRUE, sep="\t", quote=FALSE)


#### MASCOT

results_mascot_annot_disc_2UniqPep_f <- unique(results_mascot_annot_disc_2UniqPep[, c("ProteinAccession","PeptideSeq","sample")])


samples <- paste(unique(results_mascot_annot_disc_2UniqPep_f$sample))
proteins <- paste(unique(results_mascot_annot_disc_2UniqPep_f$ProteinAccession))
proteinsPerSamplesMatrix <- matrix(, nrow = length(proteins), ncol = length(samples))
#each column is a variation
for(i in 1:length(proteins)) {
	cat("mascot ",proteins[i],"\n")
	for(j in 1:length(samples)){
		tmp <- results_mascot_annot_disc_2UniqPep_f[results_mascot_annot_disc_2UniqPep_f$ProteinAccession == proteins[i], ]
		proteinsPerSamplesMatrix[i,j] <- length(unique(tmp[(paste(tmp$sample) == samples[j]), "PeptideSeq"]))
	}
}
rownames(proteinsPerSamplesMatrix) <- paste(proteins)
colnames(proteinsPerSamplesMatrix) <- paste(samples)

write.table(proteinsPerSamplesMatrix, file = paste0(resultDir, "results_mascot_annot_disc_2UniqPep_nPeptidesPerProteinAndPerSample.txt"), col.names=TRUE, row.names = TRUE, sep="\t", quote=FALSE)

#### OMSSA

results_omssa_annot_disc_2UniqPep_f <- unique(results_omssa_annot_disc_2UniqPep[, c("ProteinAccession","PeptideSeq","sample")])

samples <- paste(unique(results_omssa_annot_disc_2UniqPep_f$sample))
proteins <- paste(unique(results_omssa_annot_disc_2UniqPep_f$ProteinAccession))
proteinsPerSamplesMatrix <- matrix(, nrow = length(proteins), ncol = length(samples))
#each column is a variation
cat("length(proteins)", length(proteins),"\n")
for(i in 1:length(proteins)) {
	cat("omssa ",proteins[i],"\n")
	for(j in 1:length(samples)){
		tmp <- results_omssa_annot_disc_2UniqPep_f[results_omssa_annot_disc_2UniqPep_f$ProteinAccession == proteins[i], ]
		proteinsPerSamplesMatrix[i,j] <- length(unique(tmp[(paste(tmp$sample) == samples[j]), "PeptideSeq"]))
	}
}
rownames(proteinsPerSamplesMatrix) <- paste(proteins)
colnames(proteinsPerSamplesMatrix) <- paste(samples)

write.table(proteinsPerSamplesMatrix, file = paste0(resultDir, "results_omssa_annot_disc_2UniqPep_nPeptidesPerProteinAndPerSample.txt"), col.names=TRUE, row.names = TRUE, sep="\t", quote=FALSE)

#### TANDEM

results_tandem_annot_disc_2UniqPep_f <- unique(results_tandem_annot_disc_2UniqPep[, c("ProteinAccession","PeptideSeq","sample")])

samples <- paste(unique(results_tandem_annot_disc_2UniqPep_f$sample))
proteins <- paste(unique(results_tandem_annot_disc_2UniqPep_f$ProteinAccession))
proteinsPerSamplesMatrix <- matrix(, nrow = length(proteins), ncol = length(samples))
#each column is a variation
for(i in 1:length(proteins)) {
	cat("tandem ",proteins[i],"\n")
	for(j in 1:length(samples)){
		tmp <- results_tandem_annot_disc_2UniqPep_f[results_tandem_annot_disc_2UniqPep_f$ProteinAccession == proteins[i], ]
		proteinsPerSamplesMatrix[i,j] <- length(unique(tmp[(paste(tmp$sample) == samples[j]), "PeptideSeq"]))
	}
}
rownames(proteinsPerSamplesMatrix) <- paste(proteins)
colnames(proteinsPerSamplesMatrix) <- paste(samples)

write.table(proteinsPerSamplesMatrix, file = paste0(resultDir, "results_tandem_annot_disc_2UniqPep_nPeptidesPerProteinAndPerSample.txt"), col.names=TRUE, row.names = TRUE, sep="\t", quote=FALSE)

# Print results to file

results_comet_annot_disc_2UniqPep_f <- unique(results_comet_annot_disc_2UniqPep[, c("search_engine", "sample", "ProteinAccession", "NofDiscPeptides", "PeptideSeq", "protFDR")])
results_mascot_annot_disc_2UniqPep_f <- unique(results_mascot_annot_disc_2UniqPep[, c("search_engine", "sample", "ProteinAccession", "NofDiscPeptides", "PeptideSeq", "protFDR")])
results_omssa_annot_disc_2UniqPep_f <- unique(results_omssa_annot_disc_2UniqPep[, c("search_engine", "sample", "ProteinAccession", "NofDiscPeptides", "PeptideSeq", "protFDR")])
results_tandem_annot_disc_2UniqPep_f <- unique(results_tandem_annot_disc_2UniqPep[, c("search_engine", "sample", "ProteinAccession", "NofDiscPeptides", "PeptideSeq", "protFDR")])

results_comet_annot_disc_2UniqPep_ff <- summaryBy(protFDR~search_engine + sample + ProteinAccession + NofDiscPeptides + PeptideSeq, data=results_comet_annot_disc_2UniqPep_f, FUN=min, keep.names=TRUE)
results_mascot_annot_disc_2UniqPep_ff <- summaryBy(protFDR~search_engine + sample + ProteinAccession + NofDiscPeptides + PeptideSeq, data=results_mascot_annot_disc_2UniqPep_f, FUN=min, keep.names=TRUE)
results_omssa_annot_disc_2UniqPep_ff <- summaryBy(protFDR~search_engine + sample + ProteinAccession + NofDiscPeptides + PeptideSeq, data=results_omssa_annot_disc_2UniqPep_f, FUN=min, keep.names=TRUE)
results_tandem_annot_disc_2UniqPep_ff <- summaryBy(protFDR~search_engine + sample + ProteinAccession + NofDiscPeptides + PeptideSeq, data=results_tandem_annot_disc_2UniqPep_f, FUN=min, keep.names=TRUE)

write.table(results_comet_annot_disc_2UniqPep_ff, file = paste0(resultDir, "results_comet_annot_disc_2UniqPep.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
write.table(results_mascot_annot_disc_2UniqPep_ff, file = paste0(resultDir, "results_mascot_annot_disc_2UniqPep.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
write.table(results_omssa_annot_disc_2UniqPep_ff, file = paste0(resultDir, "results_omssa_annot_disc_2UniqPep.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
write.table(results_tandem_annot_disc_2UniqPep_ff, file = paste0(resultDir, "results_tandem_annot_disc_2UniqPep.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
