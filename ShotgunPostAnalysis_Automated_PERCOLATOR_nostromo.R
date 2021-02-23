#!/usr/bin/Rscript

args=(commandArgs(TRUE))
projectDir <- args[1]
databasename <- args[2]
dataset<- args[3]

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

mascot_txtFiles <- list.files(path = paste0(mascot_rda_folder,dataset,"/percolator/"), pattern = "*pinresults_peptides.tsv.out")
comet_txtFiles <- list.files(path = paste0(comet_rda_folder,dataset,"/percolator/"), pattern = "*pinresults_peptides.tsv.out")
omssa_txtFiles <- list.files(path = paste0(omssa_rda_folder,dataset,"/percolator/"), pattern = "*pinresults_peptides.tsv.out")
tandem_txtFiles <- list.files(path = paste0(tandem_rda_folder,dataset,"/percolator/") = "*pinresults_peptides.tsv.out")

# mascot_rdaFiles <- paste0(mascot_rda_folder,mascot_rdaFiles)
# comet_rdaFiles <- paste0(comet_rda_folder,comet_rdaFiles)
# omssa_rdaFiles <- paste0(omssa_rda_folder,omssa_rdaFiles)
# tandem_rdaFiles <- paste0(tandem_rda_folder,tandem_rdaFiles)

mascot_txtFiles<-paste0("/home/nostromo/data/pepe/EMBRIO_9_03/Dat_Files/PXD001381/percolator/",mascot_txtFiles)
comet_txtFiles<-paste0("/home/nostromo/data/pepe/EMBRIO_9_03/Omssa_files/PXD001381/percolator/",comet_txtFiles)
omssa_txtFiles<-paste0("/home/nostromo/data/pepe/EMBRIO_9_03/Omssa_files/PXD001381/percolator/",omssa_txtFiles)
tandem_txtFiles<-paste0("/home/nostromo/data/pepe/EMBRIO_9_03/Tandem_Files/PXD001381/percolator/",tandem_txtFiles)

resultDir <- paste0(projectDir, "/Results/")
y <- paste0("mkdir -p ", resultDir)
system(y)
y <- NULL



for (i in 1:length(mascot_txtFiles)) {
	mascot_file<-read.table(paste(mascot_txtFiles[i]), row.names=NULL)
    mascot_file$number<-seq(1:nrow(mascot_file))
    names(mascot_file)<-c("PSMId","score","q.value","posterior_error_prob","peptide","proteinIds","Query","number")
	cat(paste(mascot_txtFiles[i]))
	cat("\n")
	cat(dim(mascot_file))
	cat("\n")
	if(i == 1) {
		results_mascot_TyD <- mascot_file
	} else {
		results_mascot_TyD <- rbind(results_mascot_TyD, mascot_file)
	}
}

save(results_mascot_TyD, file = paste0(resultDir, "results_mascot_TyD_percolator.rda"))

for (i in 1:length(comet_txtFiles)) {
	comet_file<-read.table(paste(comet_txtFiles[i]), header=T)
    comet_file$number<-seq(1:nrow(comet_file))
    #names(comet_file)<-c("PSMId","score","q.value","posterior_error_prob","peptide","proteinIds","Query","number")
	cat(paste(comet_txtFiles[i]))
	cat("\n")
	cat(dim(comet_file))
	cat("\n")
	if(i == 1) {
		results_comet_TyD <- comet_file
	} else {
		results_comet_TyD <- rbind(results_comet_TyD, mascot_file)
	}
}

save(results_comet_TyD, file = paste0(resultDir, "results_comet_TyD_percolator.rda"))

for (i in 1:length(omssa_txtFiles)) {
	omssa_file<-read.table(paste(omssa_txtFiles[i]), header=T)
    omssa_file$number<-seq(1:nrow(omssa_file))
    #names(omssa_file)<-c("PSMId","score","q.value","posterior_error_prob","peptide","proteinIds","Query","number")
	cat(paste(omssa_txtFiles[i]))
	cat("\n")
	cat(dim(omssa_file))
	cat("\n")
	if(i == 1) {
		results_omssa_TyD <- omssa_file
	} else {
		results_omssa_TyD <- rbind(results_omssa_TyD, omssa_file)
	}
}

save(results_omssa_TyD, file = paste0(resultDir, "results_omssa_TyD_percolator.rda"))

for (i in 1:length(tandem_txtFiles)) {
	tandem_file<-read.table(paste(tandem_txtFiles[i]), header=T)
    tandem_file$number<-seq(1:nrow(tandem_file))
    #names(tandem_file)<-c("PSMId","score","q.value","posterior_error_prob","peptide","proteinIds","Query","number")
	cat(paste(tandem_txtFiles[i]))
	cat("\n")
	cat(dim(tandem_file))
	cat("\n")
	if(i == 1) {
		results_tandem_TyD <- tandem_file
	} else {
		results_tandem_TyD <- rbind(results_tandem_TyD, tandem_file)
	}
}

save(results_tandem_TyD, file = paste0(resultDir, "results_tandem_TyD_percolator.rda"))

#filtrar por q value para quedarmos con los que pasan el test
results_mascot_TyD<-results_mascot_TyD[(results_mascot_TyD$q.value<0.01),]
results_comet_TyD<-results_comet_TyD[(results_comet_TyD$q.value<0.01),]
results_omssa_TyD<-results_omssa_TyD[(results_omssa_TyD$q.value<0.01),]
results_tandem_TyD<-results_tandem_TyD[(results_tandem_TyD$q.value<0.01),]

#clean peptide sequences

PATTERN <- '\\[[0-9.-]+\\]'
PATTERN2 <- '[A-Z]{2,}'
results_mascot_TyD$PeptideSeq<-str_extract(str_remove_all(results_mascot_TyD$peptide, PATTERN), PATTERN2)
results_comet_TyD$PeptideSeq<-str_extract(str_remove_all(results_comet_TyD$peptide, PATTERN), PATTERN2)
results_omssa_TyD$PeptideSeq<-str_extract(str_remove_all(results_omssa_TyD$peptide, PATTERN), PATTERN2)
results_tandem_TyD$PeptideSeq<-str_extract(str_remove_all(results_tandem_TyD$peptide, PATTERN), PATTERN2)

save(results_comet, file = paste0(resultDir, "results_comet.rda"))
save(results_mascot, file = paste0(resultDir, "results_mascot.rda"))
save(results_omssa, file = paste0(resultDir, "results_omssa.rda"))
save(results_tandem, file = paste0(resultDir, "results_tandem.rda"))

results_comet_annot <- results_comet_TyD
results_comet_annot$discriminant <- 0
results_comet_annot[paste(results_comet_annot$PeptideSeq) %in% paste(db_peptidesXProtAll_filterAA_disc$Peptide), "discriminant"] <- 1
results_comet_annot_disc <- results_comet_annot[results_comet_annot$discriminant == 1, ]

results_mascot_annot <- results_mascot_TyD
results_mascot_annot$discriminant <- 0
results_mascot_annot[paste(results_mascot_annot$PeptideSeq) %in% paste(db_peptidesXProtAll_filterAA_disc$Peptide), "discriminant"] <- 1
results_mascot_annot_disc <- results_mascot_annot[results_mascot_annot$discriminant == 1, ]

results_omssa_annot <- results_omssa_TyD
results_omssa_annot$discriminant <- 0
results_omssa_annot[paste(results_omssa_annot$PeptideSeq) %in% paste(db_peptidesXProtAll_filterAA_disc$Peptide), "discriminant"] <- 1
results_omssa_annot_disc <- results_omssa_annot[results_omssa_annot$discriminant == 1, ]

results_tandem_annot <- results_tandem_TyD
results_tandem_annot$discriminant <- 0
results_tandem_annot[paste(results_tandem_annot$PeptideSeq) %in% paste(db_peptidesXProtAll_filterAA_disc$Peptide), "discriminant"] <- 1
results_tandem_annot_disc <- results_tandem_annot[results_tandem_annot$discriminant == 1, ]

save(results_comet_annot, file = paste0(resultDir, "results_comet_annot_percolator.rda"))
save(results_comet_annot_disc, file = paste0(resultDir, "results_comet_annot_disc_percolator.rda"))
save(results_mascot_annot, file = paste0(resultDir, "results_mascot_annot_percolator.rda"))
save(results_mascot_annot_disc, file = paste0(resultDir, "results_mascot_annot_disc_percolator.rda"))
save(results_omssa_annot, file = paste0(resultDir, "results_omssa_annot_percolator.rda"))
save(results_omssa_annot_disc, file = paste0(resultDir, "results_omssa_annot_disc_percolator.rda"))
save(results_tandem_annot, file = paste0(resultDir, "results_tandem_annot_percolator.rda"))
save(results_tandem_annot_disc, file = paste0(resultDir, "results_tandem_annot_disc_percolator.rda"))

results_mascot_annot_disc$sample<-"PXD001381"
results_omssa_annot_disc$sample<-"PXD001381"
results_comet_annot_disc$sample<-"PXD001381"
results_tandem_annot_disc$sample<-"PXD001381"

tmp <- unique(results_comet_annot_disc[, c("sample", "proteinIds", "PeptideSeq")])
tmp2 <- summaryBy(PeptideSeq~proteinIds + sample, data=tmp, FUN=length, keep.names=TRUE)
colnames(tmp2)[3] <- "NofDiscPeptides"

results_comet_annot_disc <- merge(results_comet_annot_disc, tmp2, by.x = c("proteinIds", "sample"), by.y = c("proteinIds", "sample"), all.x = TRUE)
results_comet_annot_disc[is.na(results_comet_annot_disc$NofDiscPeptides), "NofDiscPeptides"] <- 0

tmp <- unique(results_mascot_annot_disc[, c("sample", "proteinIds", "PeptideSeq")])
tmp2 <- summaryBy(PeptideSeq~proteinIds + sample, data=tmp, FUN=length, keep.names=TRUE)
colnames(tmp2)[3] <- "NofDiscPeptides"

results_mascot_annot_disc <- merge(results_mascot_annot_disc, tmp2, by.x = c("proteinIds", "sample"), by.y = c("proteinIds", "sample"), all.x = TRUE)
results_mascot_annot_disc[is.na(results_mascot_annot_disc$NofDiscPeptides), "NofDiscPeptides"] <- 0

tmp <- unique(results_omssa_annot_disc[, c("sample", "proteinIds", "PeptideSeq")])
tmp2 <- summaryBy(PeptideSeq~proteinIds + sample, data=tmp, FUN=length, keep.names=TRUE)
colnames(tmp2)[3] <- "NofDiscPeptides"

results_omssa_annot_disc <- merge(results_omssa_annot_disc, tmp2, by.x = c("proteinIds", "sample"), by.y = c("proteinIds", "sample"), all.x = TRUE)
results_omssa_annot_disc[is.na(results_omssa_annot_disc$NofDiscPeptides), "NofDiscPeptides"] <- 0

tmp <- unique(results_tandem_annot_disc[, c("sample", "proteinIds", "PeptideSeq")])
tmp2 <- summaryBy(PeptideSeq~proteinIds + sample, data=tmp, FUN=length, keep.names=TRUE)
colnames(tmp2)[3] <- "NofDiscPeptides"

results_tandem_annot_disc <- merge(results_tandem_annot_disc, tmp2, by.x = c("proteinIds", "sample"), by.y = c("proteinIds", "sample"), all.x = TRUE)
results_tandem_annot_disc[is.na(results_tandem_annot_disc$NofDiscPeptides), "NofDiscPeptides"] <- 0

results_comet_annot_disc_2UniqPep <- results_comet_annot_disc[results_comet_annot_disc$NofDiscPeptides > 1,]
results_mascot_annot_disc_2UniqPep <- results_mascot_annot_disc[results_mascot_annot_disc$NofDiscPeptides > 1,]
results_omssa_annot_disc_2UniqPep <- results_omssa_annot_disc[results_omssa_annot_disc$NofDiscPeptides > 1,]
results_tandem_annot_disc_2UniqPep <- results_tandem_annot_disc[results_tandem_annot_disc$NofDiscPeptides > 1,]

save(results_comet_annot_disc_2UniqPep, file = paste0(resultDir, "results_comet_annot_disc_2UniqPep_percolator.rda"))
save(results_mascot_annot_disc_2UniqPep, file = paste0(resultDir, "results_mascot_annot_disc_2UniqPep_percolator.rda"))
save(results_omssa_annot_disc_2UniqPep, file = paste0(resultDir, "results_omssa_annot_disc_2UniqPep_percolator.rda"))
save(results_tandem_annot_disc_2UniqPep, file = paste0(resultDir, "results_tandem_annot_disc_2UniqPep_percolator.rda"))

#eliminar descripcion de tandem protein accesion

tmp_id<-paste(unlist(results_tandem_annot_disc_2UniqPep$proteinIds))
tmp_id<-sapply(strsplit(tmp_id,"-"), "[", 1)
results_tandem_annot_disc_2UniqPep$proteinIds<-tmp_id


results_ALLSE_annot_disc_uniqueness_2UniqPep_f <- rbind(unique(results_comet_annot_disc_2UniqPep[, c("proteinIds","PeptideSeq","sample")]), unique(results_mascot_annot_disc_2UniqPep[, c("proteinIds","PeptideSeq","sample")]), unique(results_omssa_annot_disc_2UniqPep[, c("proteinIds","PeptideSeq","sample")]), unique(results_tandem_annot_disc_2UniqPep[, c("proteinIds","PeptideSeq","sample")]))

samples <- paste(unique(results_ALLSE_annot_disc_uniqueness_2UniqPep_f$sample))
proteins <- paste(unique(results_ALLSE_annot_disc_uniqueness_2UniqPep_f$proteinIds))
proteinsPerSamplesMatrix <- matrix(, nrow = length(proteins), ncol = length(samples))

#each column is a variation
for(i in 1:length(proteins)) {
	for(j in 1:length(samples)){
		tmp <- results_ALLSE_annot_disc_uniqueness_2UniqPep_f[results_ALLSE_annot_disc_uniqueness_2UniqPep_f$proteinIds == proteins[i], ]
		proteinsPerSamplesMatrix[i,j] <- length(unique(tmp[(paste(tmp$sample) == samples[j]), "PeptideSeq"]))
	}
}

#### COMET

results_comet_annot_disc_2UniqPep_f <- unique(results_comet_annot_disc_2UniqPep[, c("ProteinAccession","PeptideSeq","sample")])

samples <- paste(unique(results_comet_annot_disc_2UniqPep_f$sample))
proteins <- paste(unique(results_comet_annot_disc_2UniqPep_f$ProteinAccession))
proteinsPerSamplesMatrix <- matrix(, nrow = length(proteins), ncol = length(samples))
#each column is a variation
for(i in 1:length(proteins)) {
	for(j in 1:length(samples)){
		tmp <- results_comet_annot_disc_2UniqPep_f[results_comet_annot_disc_2UniqPep_f$ProteinAccession == proteins[i], ]
		proteinsPerSamplesMatrix[i,j] <- length(unique(tmp[(paste(tmp$sample) == samples[j]), "PeptideSeq"]))
	}
}
rownames(proteinsPerSamplesMatrix) <- paste(proteins)
colnames(proteinsPerSamplesMatrix) <- paste(samples)

write.table(proteinsPerSamplesMatrix, file = paste0(resultDir, "results_comet_annot_disc_2UniqPep_nPeptidesPerProteinAndPerSample_percolator.txt"), col.names=TRUE, row.names = TRUE, sep="\t", quote=FALSE)


#### MASCOT

results_mascot_annot_disc_2UniqPep_f <- unique(results_mascot_annot_disc_2UniqPep[, c("proteinIds","PeptideSeq","sample")])


samples <- paste(unique(results_mascot_annot_disc_2UniqPep_f$sample))
proteins <- paste(unique(results_mascot_annot_disc_2UniqPep_f$proteinIds))
proteinsPerSamplesMatrix <- matrix(, nrow = length(proteins), ncol = length(samples))
#each column is a variation
for(i in 1:length(proteins)) {
	for(j in 1:length(samples)){
		tmp <- results_mascot_annot_disc_2UniqPep_f[results_mascot_annot_disc_2UniqPep_f$proteinIds == proteins[i], ]
		proteinsPerSamplesMatrix[i,j] <- length(unique(tmp[(paste(tmp$sample) == samples[j]), "PeptideSeq"]))
	}
}
rownames(proteinsPerSamplesMatrix) <- paste(proteins)
colnames(proteinsPerSamplesMatrix) <- paste(samples)

write.table(proteinsPerSamplesMatrix, file = paste0(resultDir, "results_mascot_annot_disc_2UniqPep_nPeptidesPerProteinAndPerSample_percolator.txt"), col.names=TRUE, row.names = TRUE, sep="\t", quote=FALSE)

#### OMSSA

results_omssa_annot_disc_2UniqPep_f <- unique(results_omssa_annot_disc_2UniqPep[, c("proteinIds","PeptideSeq","sample")])

samples <- paste(unique(results_omssa_annot_disc_2UniqPep_f$sample))
proteins <- paste(unique(results_omssa_annot_disc_2UniqPep_f$proteinIds))
proteinsPerSamplesMatrix <- matrix(, nrow = length(proteins), ncol = length(samples))
#each column is a variation
for(i in 1:length(proteins)) {
	for(j in 1:length(samples)){
		tmp <- results_omssa_annot_disc_2UniqPep_f[results_omssa_annot_disc_2UniqPep_f$proteinIds == proteins[i], ]
		proteinsPerSamplesMatrix[i,j] <- length(unique(tmp[(paste(tmp$sample) == samples[j]), "PeptideSeq"]))
	}
}
rownames(proteinsPerSamplesMatrix) <- paste(proteins)
colnames(proteinsPerSamplesMatrix) <- paste(samples)

write.table(proteinsPerSamplesMatrix, file = paste0(resultDir, "results_omssa_annot_disc_2UniqPep_nPeptidesPerProteinAndPerSample_percolator.txt"), col.names=TRUE, row.names = TRUE, sep="\t", quote=FALSE)

#### TANDEM

results_tandem_annot_disc_2UniqPep_f <- unique(results_tandem_annot_disc_2UniqPep[, c("proteinIds","PeptideSeq","sample")])

samples <- paste(unique(results_tandem_annot_disc_2UniqPep_f$sample))
proteins <- paste(unique(results_tandem_annot_disc_2UniqPep_f$proteinIds))
proteinsPerSamplesMatrix <- matrix(, nrow = length(proteins), ncol = length(samples))
#each column is a variation
for(i in 1:length(proteins)) {
	for(j in 1:length(samples)){
		tmp <- results_tandem_annot_disc_2UniqPep_f[results_tandem_annot_disc_2UniqPep_f$proteinIds == proteins[i], ]
		proteinsPerSamplesMatrix[i,j] <- length(unique(tmp[(paste(tmp$sample) == samples[j]), "PeptideSeq"]))
	}
}
rownames(proteinsPerSamplesMatrix) <- paste(proteins)
colnames(proteinsPerSamplesMatrix) <- paste(samples)

write.table(proteinsPerSamplesMatrix, file = paste0(resultDir, "results_tandem_annot_disc_2UniqPep_nPeptidesPerProteinAndPerSample_percolator.txt"), col.names=TRUE, row.names = TRUE, sep="\t", quote=FALSE)

results_comet_annot_disc_2UniqPep$search_engine<-"COMET"
results_mascot_annot_disc_2UniqPep$search_engine<-"MASCOT"
results_omssa_annot_disc_2UniqPep$search_engine<-"OMSSA"
results_tandem_annot_disc_2UniqPep$search_engine<-"TANDEM"


results_comet_annot_disc_2UniqPep_f <- unique(results_comet_annot_disc_2UniqPep[, c("search_engine","sample", "proteinIds", "NofDiscPeptides", "PeptideSeq","q.value")])
results_mascot_annot_disc_2UniqPep_f <- unique(results_mascot_annot_disc_2UniqPep[, c("search_engine","sample", "proteinIds", "NofDiscPeptides", "PeptideSeq","q.value")])
results_omssa_annot_disc_2UniqPep_f <- unique(results_omssa_annot_disc_2UniqPep[, c("search_engine","sample", "proteinIds", "NofDiscPeptides", "PeptideSeq","q.value")])
results_tandem_annot_disc_2UniqPep_f <- unique(results_tandem_annot_disc_2UniqPep[, c("search_engine","sample", "proteinIds", "NofDiscPeptides", "PeptideSeq","q.value")])

results_comet_annot_disc_2UniqPep_ff <- summaryBy(q.value~search_engine + sample + proteinIds + NofDiscPeptides + PeptideSeq, data=results_comet_annot_disc_2UniqPep_f, FUN=min, keep.names=TRUE)
results_mascot_annot_disc_2UniqPep_ff <- summaryBy(q.value~search_engine + sample + proteinids + NofDiscPeptides + PeptideSeq, data=results_mascot_annot_disc_2UniqPep_f, FUN=min, keep.names=TRUE)
results_omssa_annot_disc_2UniqPep_ff <- summaryBy(q.value~search_engine + sample + proteinIds + NofDiscPeptides + PeptideSeq, data=results_omssa_annot_disc_2UniqPep_f, FUN=min, keep.names=TRUE)
results_tandem_annot_disc_2UniqPep_ff <- summaryBy(q.value~search_engine + sample + proteinIds + NofDiscPeptides + PeptideSeq, data=results_tandem_annot_disc_2UniqPep_f, FUN=min, keep.names=TRUE)

write.table(results_comet_annot_disc_2UniqPep_ff, file = paste0(resultDir, "results_comet_annot_disc_2UniqPep_percolator.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
write.table(results_mascot_annot_disc_2UniqPep_ff, file = paste0(resultDir, "results_mascot_annot_disc_2UniqPep_percolator.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
write.table(results_omssa_annot_disc_2UniqPep_ff, file = paste0(resultDir, "results_omssa_annot_disc_2UniqPep_percolator.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
write.table(results_tandem_annot_disc_2UniqPep_ff, file = paste0(resultDir, "results_tandem_annot_disc_2UniqPep_percolator.txt"), col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
