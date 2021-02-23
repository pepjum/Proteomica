#!/usr/bin/Rscript

args=(commandArgs(TRUE))
inputPath <- args[1]
outputMGFFilePath <- args[2]
mgfFilePath <- args[3]
dataset <- args[4]
rdaFilePath <- args[5]
fase <- args[6]
search_engine <- args[7]
mgfTextFile <- args[8]
source("/home/bioinformatica/datos/01_Rscripts/A_Funciones/funcionesShotgun.R")
cat("Starting the script \n")

# cat(inputPath)
# cat("\n")
# cat(outputMGFFilePath)
# cat("\n")
# cat(mgfFilePath)
# cat("\n")
# cat(dataset)
# cat("\n")
# cat(rdaFilePath)
# cat("\n")
# cat(fase)
# cat("\n")
# cat(search_engine)
# cat("\n")
# cat(mgfTextFile)
# cat("\n")


if(identical(fase,"1")) {
	exp_metadata <- read.table(mgfTextFile, header = FALSE, sep = " ")
  	colnames(exp_metadata) <- c("file", "database", "sample", "fraction", "mgfFile")
  	cat("\n")
	resultObj <- get(load(rdaFilePath))
	resultTable <- resultObj$resProtGroups
	resultTable <- resultTable[resultTable$database == 'T',]
	mgfList <- ReadMGFFile(mgfFilePath)
	cat(mgfFilePath)
	cat("\n")
	cat("Length of MGFFile: ")
	cat(length(mgfList))
	cat("\n")
	mgf <- paste(dataset, ".mgf", sep = "")
	mgfFraction <- exp_metadata[paste(exp_metadata$mgfFile) == paste(mgf),"fraction"][1]
	if(any(nchar(paste(resultTable$fraction)) == 1) & all(nchar(paste(exp_metadata$fraction)) == 2)) {
			if ((mgfFraction < 10) && (nchar(mgfFraction) == 2)) {
				mgfFraction <- substring(mgfFraction, 1,1)
			} else {
				mgfFraction <- paste(mgfFraction)
			}
	} else if(all(nchar(paste(resultTable$fraction)) == 2) & any(nchar(paste(exp_metadata$fraction)) == 1)) {
			if ((mgfFraction < 10) && (nchar(mgfFraction) == 1)) {
				mgfFraction <- paste("0", mgfFraction, sep = "")
			} else {
				mgfFraction <- paste(mgfFraction)
			}
	} else {
			mgfFraction <- paste(mgfFraction)
	}

	resultTablef <- resultTable[paste(resultTable$fraction) == paste(mgfFraction),]

	if (dim(resultTablef)[1] == 0) {
		noHitMGF <- mgfList

	} else if (length(resultTablef$mgfIndex) == 0) {
		 ## método antiguo
		dats <- unique(paste(resultTablef$queryTitleDat))
		if (is.na(strsplit(names(mgfList[1])," ")[[1]][2]) == TRUE)	{
			mgfListf <- names(mgfList)
		} else {
			mgfListf <- unlist(sapply(1:length(names(mgfList)), FUN = function(x) strsplit(names(mgfList)," ")[[x]][1]))
		}
		if (is.na(strsplit(paste(dats), "\\%")[[1]][2]) == TRUE ) {
			datList <- paste(dats)
		} else {
			datList <- 	unlist(sapply(1:length(dats), FUN = function(x) strsplit(paste(dats), "\\%")[[x]][1]))
		}
		noHitMGF <- mgfList[!(mgfListf %in% datList)]
		mgfHits <- length(mgfList) - length(noHitMGF)
		cat("Number of hits: ")
		cat(mgfHits)
		cat("\n")
	} else {
	## nuevo método
		mgfHits <- resultTablef$mgfIndex
		hitMgfIndex <- as.numeric(paste(unique(mgfHits)))
		if (length(hitMgfIndex) != 0) {
			noHitMGF <- mgfList[-hitMgfIndex]
		} else {
			noHitMGF <- mgfList
		}
		cat("Number of hits: ")
		cat(length(hitMgfIndex))
		cat("\n")
	}
	cat("Length of new MGFFile: ")
	cat(length(noHitMGF))
	cat("\n")
	outputMGFFilePath <- paste(substr(outputMGFFilePath, 1, nchar(outputMGFFilePath) - 4), "_noHitKnown.mgf", sep="")
	WriteMGFFile(noHitMGF, outputMGFFilePath)


} else if (identical(fase,"2")) {
	exp_metadata <- read.table(mgfTextFile, header = FALSE, sep = " ")
  	colnames(exp_metadata) <- c("file", "database", "sample", "fraction", "mgfFile")
	resultTable <- get(load(rdaFilePath))
	resultTable <- resultTable[resultTable$database == 'T',]
	mgfList <- ReadMGFFile(mgfFilePath)
	cat(rdaFilePath)
	cat("\n")
	cat(mgfFilePath)
	cat("\n")
	cat("Length of MGFFile: ")
	cat(length(mgfList))
	cat("\n")
	mgf <- paste(dataset, ".mgf", sep = "")
	mgfFraction <- exp_metadata[paste(exp_metadata$mgfFile) == paste(mgf),"fraction"][1]
	if(any(nchar(paste(resultTable$fraction)) == 1) & all(nchar(paste(exp_metadata$fraction)) == 2)) {
			if ((mgfFraction < 10) && (nchar(mgfFraction) == 2)) {
				mgfFraction <- substring(mgfFraction, 1,1)
			} else {
				mgfFraction <- paste(mgfFraction)
			}
	} else if(all(nchar(paste(resultTable$fraction)) == 2) & any(nchar(paste(exp_metadata$fraction)) == 1)) {
			if ((mgfFraction < 10) && (nchar(mgfFraction) == 1)) {
				mgfFraction <- paste("0", mgfFraction, sep = "")
			} else {
				mgfFraction <- paste(mgfFraction)
			}
	} else {
			mgfFraction <- paste(mgfFraction)
	}

	resultTablef <- resultTable[paste(resultTable$fraction) == paste(mgfFraction),]
	if (dim(resultTablef)[1] != 0)
	{
		if (length(resultTablef$mgfIndex) == 0) {
			dats <- unique(paste(resultTablef$queryTitleDat))
			if (is.na(strsplit(names(mgfList[1])," ")[[1]][2]) == TRUE)	{
				mgfListf <- names(mgfList)
			} else {
				mgfListf <- unlist(sapply(1:length(names(mgfList)), FUN = function(x) strsplit(names(mgfList)," ")[[x]][1]))
			}
			if (is.na(strsplit(paste(dats), "\\%")[[1]][2]) == TRUE ) {
				datList <- paste(dats)
			} else {
				datList <- 	unlist(sapply(1:length(dats), FUN = function(x) strsplit(paste(dats), "\\%")[[x]][1]))
			}
			noHitMGF <- mgfList[!(mgfListf %in% datList)]
			mgfHits <- length(mgfList) - length(noHitMGF)
			cat("Number of hits: ")
			cat(mgfHits)
			cat("\n")
		} else {
		## nuevo método
			mgfHits <- resultTablef$mgfIndex
			hitMgfIndex <- as.numeric(paste(unique(mgfHits)))
			if (length(hitMgfIndex) != 0) {
				noHitMGF <- mgfList[-hitMgfIndex]
			} else {
				noHitMGF <- mgfList
			}
			cat("Number of hits: ")
			cat(length(hitMgfIndex))
			cat("\n")
		}
	} else
	{
		noHitMGF <- mgfList
		cat("Number of hits: 0")
		cat("\n")
	}

	cat("Length of new MGFFile: ")
	cat(length(noHitMGF))
	cat("\n")
	outputMGFFilePath <- paste(substr(outputMGFFilePath, 1, nchar(outputMGFFilePath) - 15), "_noHitUnknown.mgf", sep="")
	WriteMGFFile(noHitMGF, outputMGFFilePath)

} else {
	exp_metadata <- read.table(mgfTextFile, header = FALSE, sep = " ")
  	colnames(exp_metadata) <- c("file", "database", "sample", "fraction", "mgfFile")
	resultTable <- get(load(rdaFilePath))
	mgfList <- ReadMGFFile(mgfFilePath)
	cat(rdaFilePath)
	cat("\n")
	cat(mgfFilePath)
	cat("\n")
	cat("Length of MGFFile: ")
	cat(length(mgfList))
	cat("\n")

	if (dim(resultTable)[1] != 0)
	{
		mgf <- paste(dataset, ".mgf", sep = "")
		mgfFraction <- exp_metadata[paste(exp_metadata$mgfFile) == paste(mgf),"fraction"][1]
		if(any(nchar(paste(resultTable$fraction)) == 1) & all(nchar(paste(exp_metadata$fraction)) == 2)) {
			if ((mgfFraction < 10) && (nchar(mgfFraction) == 2)) {
				mgfFraction <- substring(mgfFraction, 1,1)
			} else {
				mgfFraction <- paste(mgfFraction)
			}
		} else if(all(nchar(paste(resultTable$fraction)) == 2) & any(nchar(paste(exp_metadata$fraction)) == 1)) {
			if ((mgfFraction < 10) && (nchar(mgfFraction) == 1)) {
				mgfFraction <- paste("0", mgfFraction, sep = "")
			} else {
				mgfFraction <- paste(mgfFraction)
			}
		} else {
			mgfFraction <- paste(mgfFraction)
		}

		resultTablef <- resultTable[paste(resultTable$fraction) == paste(mgfFraction),]
		if (dim(resultTablef)[1] != 0)
		{
			if (length(resultTablef$mgfIndex) == 0) {
				dats <- unique(paste(resultTablef$queryTitleDat))
				if (is.na(strsplit(names(mgfList[1])," ")[[1]][2]) == TRUE)	{
					mgfListf <- names(mgfList)
				} else {
					mgfListf <- unlist(sapply(1:length(names(mgfList)), FUN = function(x) strsplit(names(mgfList)," ")[[x]][1]))
				}
				if (is.na(strsplit(paste(dats), "\\%")[[1]][2]) == TRUE ) {
					datList <- paste(dats)
				} else {
					datList <- 	unlist(sapply(1:length(dats), FUN = function(x) strsplit(paste(dats), "\\%")[[x]][1]))
				}
				noHitMGF <- mgfList[!(mgfListf %in% datList)]
				mgfHits <- length(mgfList) - length(noHitMGF)
				cat("Number of hits: ")
				cat(mgfHits)
				cat("\n")
			} else {
			## nuevo método
				mgfHits <- resultTablef$mgfIndex
				hitMgfIndex <- as.numeric(paste(unique(mgfHits)))
				if (length(hitMgfIndex) != 0) {
					noHitMGF <- mgfList[-hitMgfIndex]
				} else {
					noHitMGF <- mgfList
				}
				cat("Number of hits: ")
				cat(length(hitMgfIndex))
				cat("\n")
			}
		} else
		{
			noHitMGF <- mgfList
			cat("Number of hits: 0")
			cat("\n")
		}

	} else
	{
		noHitMGF <- mgfList
			cat("Number of hits: 0")
			cat("\n")
	}
	cat("Length of new MGFFile: ")
		cat(length(noHitMGF))
		cat("\n")
		outputMGFFilePath <- paste(substr(outputMGFFilePath, 1, nchar(outputMGFFilePath) - 17), "_unassigned.mgf", sep="")
		WriteMGFFile(noHitMGF, outputMGFFilePath)

}
