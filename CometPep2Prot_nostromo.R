####################################################################################
###
args=(commandArgs(TRUE))

# params
currentFolder <- args[1]
# currentFolder <- "/mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/"
dataset <- args[2]
# dataset <- "349-028"
 txtFileName <- args[3]
# txtFileName <- "/mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/Comet_files/Shotgun_comet_349-028.txt"
DBtargetFile <- args[4]
# DBtargetFile <- "uniprot_sprot_2017_12.fasta"
decoy_id <- args[5]
# decoy_id <- "DECOY"
todo <- args[6]
# TOD; 1 si solo hay que hacer target; 2 si solo hay que hacer decoy; 0 si hay que hacer los dos.

# cat(currentFolder)
# cat("\n")
# cat(dataset)
# cat("\n")
# cat(txtFileName)
# cat("\n")
# cat(DBtargetFile)
# cat("\n")
# cat(decoy_id)
# cat("\n")
# cat(todo)
# cat("\n")

library(data.table)
library(Biostrings)
library(doBy)
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesShotgun.R")

mgfFolder <- paste(currentFolder, "MGFFiles/", sep = "/")
cometFolder <- paste(currentFolder, "Comet_files/", sep = "/")
datasetMgfFolder <- paste(mgfFolder, dataset, "/", sep = "")
datasetCometTargetFolder <- paste(cometFolder, dataset, "/", sep = "")
datasetCometDecoyFolder <- paste(cometFolder, dataset, "-D/",sep = "")
mgfFiles <- dir(path = datasetMgfFolder, pattern = ".mgf")
DBtargetFile <- paste(currentFolder, DBtargetFile, sep = "/")
dbName <- strsplit(basename(DBtargetFile), "\\.")[[1]][1]

if (todo == 0 || todo == 1) {
	targetCometFiles <- dir(path = datasetCometTargetFolder, pattern = ".txt")
	targetCometFiles <- targetCometFiles[!grepl("peptides.txt", targetCometFiles)]
	targetCometFiles <- targetCometFiles[!grepl("_peptides_out.txt", targetCometFiles)]
	targetCometFiles <- targetCometFiles[!grepl("_corrected.txt", targetCometFiles)]
	targetCometFiles <- paste(datasetCometTargetFolder, targetCometFiles, sep = "")
	targetOutCometFiles <- paste(substr(targetCometFiles, 1, nchar(targetCometFiles) - 4), "_peptides.txt", sep="")
	targetPMOutCometFiles <- paste(substr(targetCometFiles, 1, nchar(targetCometFiles) - 4), "_peptides_out.txt", sep="")
	targetCorrectedCometFiles <- paste(substr(targetCometFiles, 1, nchar(targetCometFiles) - 4), "_corrected.txt", sep="")
	targetLogCometFiles <-  paste(substr(targetCometFiles, 1, nchar(targetCometFiles) - 4), "_peptides_log.log", sep="")
	for(i in 1:length(targetCometFiles)) {

		cosmicResult <- read.table(targetCometFiles[i], skip = 2, header = FALSE, sep = "\t")
		if (dim(cosmicResult)[2] == 19) {
			cosmicResult[,19] <- NULL
			colnames(cosmicResult) <- c("scan", "num", "charge", "exp_neutral_mass", "calc_neutral_mass", "e-value", "xcorr", "delta_cn", "sp_score", "ions_matched", "ions_total", "plain_peptide", "modified_peptide", "prev_aa", "next_aa", "protein", "protein_count", "modifications")
			write.table(cosmicResult$plain_peptide, file = targetOutCometFiles[i], quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")

			y <- paste('java -XX:ParallelGCThreads=4 -jar /opt/PeptideMatchCMD_1.0.jar -a query -i ', currentFolder, '/PeptideMatch_', dbName, '_index -Q ',targetOutCometFiles[i],' -l -e -o ', targetPMOutCometFiles[i], ' > ', targetLogCometFiles[i], sep = "")
			system(y)
			y <- NULL

			cat(targetPMOutCometFiles[i])
			peptidesMatched <- read.csv2(targetPMOutCometFiles[i], header = FALSE, sep = "\t", skip = 2)
			peptidesMatched <- peptidesMatched[,1:2]
			peptidesMatched <- unique(peptidesMatched)
			colnames(peptidesMatched) <- c("peptide", "protein_id")
			cosmicResultMatched <- merge(cosmicResult, peptidesMatched, by.x = "plain_peptide", by.y = "peptide", all.x = TRUE)
			cosmicResultMatched <- unique(cosmicResultMatched)

			cosmicResultMatched$protein <- NULL
			cosmicResultMatched$protein <- cosmicResultMatched$protein_id
			cosmicResultMatched$protein_id <- NULL

			write.table(cosmicResultMatched, file = targetCorrectedCometFiles[i], quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")

			cosmicResult <- NULL
			cosmicResultMatched <- NULL
			peptidesMatched <- NULL
		} else {
			cat("Error in filename", targetCometFiles[i])
			error=1
			break
		}

	}
}

if (todo == 0 || todo == 2) {
	decoyCometFiles <- dir(path = datasetCometDecoyFolder, pattern = ".txt")
	decoyCometFiles <- decoyCometFiles[!grepl("peptides.txt", decoyCometFiles)]
	decoyCometFiles <- decoyCometFiles[!grepl("_peptides_out.txt", decoyCometFiles)]
	decoyCometFiles <- decoyCometFiles[!grepl("_corrected.txt", decoyCometFiles)]
	decoyCometFiles <- paste(datasetCometDecoyFolder, decoyCometFiles, sep = "")
	decoyOutCometFiles <- paste(substr(decoyCometFiles, 1, nchar(decoyCometFiles) - 4), "_peptides.txt", sep="")
	decoyPMOutCometFiles <- paste(substr(decoyCometFiles, 1, nchar(decoyCometFiles) - 4), "_peptides_out.txt", sep="")
	decoyCorrectedCometFiles <- paste(substr(decoyCometFiles, 1, nchar(decoyCometFiles) - 4), "_corrected.txt", sep="")
	decoyLogCometFiles <-  paste(substr(decoyCometFiles, 1, nchar(decoyCometFiles) - 4), "_peptides_log.log", sep="")
	for(i in 1:length(decoyCometFiles)) {
		cosmicResult <- read.table(decoyCometFiles[i], skip = 2, header = FALSE, sep = "\t")
		if (dim(cosmicResult)[2] == 19) {
			cosmicResult[,19] <- NULL
			colnames(cosmicResult) <- c("scan", "num", "charge", "exp_neutral_mass", "calc_neutral_mass", "e-value", "xcorr", "delta_cn", "sp_score", "ions_matched", "ions_total", "plain_peptide", "modified_peptide", "prev_aa", "next_aa", "protein", "protein_count", "modifications")
			write.table(cosmicResult$plain_peptide, file = decoyOutCometFiles[i], quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")

			y <- paste('java -XX:ParallelGCThreads=4 -jar /opt/PeptideMatchCMD_1.0.jar -a query -i ', currentFolder, '/PeptideMatch_', dbName, '_', decoy_id,'_index -Q ',decoyOutCometFiles[i],' -l -e -o ',decoyPMOutCometFiles[i], ' > ', decoyLogCometFiles[i], sep = "")
			system(y)
			y <- NULL

			cat(decoyPMOutCometFiles[i])
			peptidesMatched <- read.csv2(decoyPMOutCometFiles[i], header = FALSE, sep = "\t")
			peptidesMatched <- peptidesMatched[,1:2]
			peptidesMatched <- unique(peptidesMatched)
			colnames(peptidesMatched) <- c("peptide", "protein_id")
			cosmicResultMatched <- merge(cosmicResult, peptidesMatched, by.x = "plain_peptide", by.y = "peptide", all.x = TRUE)
			cosmicResultMatched <- unique(cosmicResultMatched)

			cosmicResultMatched$protein <- NULL
			cosmicResultMatched$protein <- cosmicResultMatched$protein_id
			cosmicResultMatched$protein_id <- NULL

			write.table(cosmicResultMatched, file = decoyCorrectedCometFiles[i], quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")

			cosmicResult <- NULL
			cosmicResultMatched <- NULL
			peptidesMatched <- NULL
		} else {
			cat("Error in filename", decoyCometFiles[i])
			error = 1
			break
		}

	}

}

targetCometFiles <- dir(path = datasetCometTargetFolder, pattern = "_corrected.txt")
targetCometFiles <- targetCometFiles[order(targetCometFiles, decreasing=FALSE)]
targetCometFiles <- paste(datasetCometTargetFolder, targetCometFiles, sep = "")
decoyCometFiles <- dir(path = datasetCometDecoyFolder, pattern = "_corrected.txt")
decoyCometFiles <- decoyCometFiles[order(decoyCometFiles, decreasing=FALSE)]
decoyCometFiles <- paste(datasetCometDecoyFolder, decoyCometFiles, sep = "")

tmp1 <- data.frame(a = targetCometFiles, b = "T", c = dataset, d = 1:length(targetCometFiles))
tmp2 <- data.frame(a = decoyCometFiles, b = "D", c = dataset, d = 1:length(decoyCometFiles))

tmp <- rbind(tmp1, tmp2)
write.table(tmp, file=txtFileName, col.names = FALSE, quote = FALSE, sep = " ", row.names = FALSE)
