####################################################################################
###
args=(commandArgs(TRUE))

# params
currentFolder <- args[1]
# currentFolder <- "/mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/"
dataset <- args[2]
# dataset <- "349-028"

# cat(currentFolder)
# cat("\n")
# cat(dataset)
# cat("\n")

library(data.table)
library(Biostrings)
library(doBy)
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesShotgun.R") #change for cluster

OmssaFolder <- paste0(currentFolder, "Omssa_files/")
datasetOmssaTargetFolder <- paste0(OmssaFolder, dataset, "/")
datasetOmssaDecoyFolder <- paste0(OmssaFolder, dataset, "-D/")
txtFileName <- paste0(currentFolder, "Omssa_files/Shotgun_omssa_", dataset, ".txt")
targetOmssaFiles <- dir(path = datasetOmssaTargetFolder, pattern = ".txt")
targetOmssaFiles <- targetOmssaFiles[order(targetOmssaFiles, decreasing=FALSE)]
targetOmssaFiles <- paste0(datasetOmssaTargetFolder, targetOmssaFiles)
decoyOmssaFiles <- dir(path = datasetOmssaDecoyFolder, pattern = ".txt")
decoyOmssaFiles <- decoyOmssaFiles[order(decoyOmssaFiles, decreasing=FALSE)]
decoyOmssaFiles <- paste0(datasetOmssaDecoyFolder, decoyOmssaFiles)

tmp1 <- data.frame(a = targetOmssaFiles, b = "T", c = dataset, d = 1:length(targetOmssaFiles))
tmp2 <- data.frame(a = decoyOmssaFiles, b = "D", c = dataset, d = 1:length(decoyOmssaFiles))

tmp <- rbind(tmp1, tmp2)
cat("Writing the ", txtFileName, " file")
write.table(tmp, file=txtFileName, col.names = FALSE, quote = FALSE, sep = " ", row.names = FALSE)
