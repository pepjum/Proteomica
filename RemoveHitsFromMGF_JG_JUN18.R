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
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesShotgun.R")
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
mgfFilePath<-"/home/nostromo/data/pepe/04_DataFelix_scORF_MASCOT_MAY18/MGFFiles/SUPERFELIX_nohits_PROTEINDB_MASCOT.mgf"   #path hasta el fichero mgf
dataset<-"SUPERFELIX_nohits_PROTEINDB_MASCOTs"    #nombre del mgf sin la extension
rdaFilePath<-"/home/nostromo/data/pepe/04_DataFelix_scORF_MASCOT_MAY18/Dat_Files/results_Peptides_mascot_FELIX_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter.rda"
fase="1"
search_engine<-"comet"
mgfTextFile<-"/home/nostromo/data/pepe/04_DataFelix_scORF_MASCOT_MAY18/Dat_Files/mascot_files.txt"
outputMGFFilePath<-"/home/nostromo/data/pepe/04_DataFelix_scORF_MASCOT_MAY18/MGFFiles/SUPERFELIX_nohits_PROTEINDB_MASCOT.mgf"
if(identical(fase,"1")) {
	exp_metadata <- read.table(mgfTextFile, header = FALSE, sep = " ")
	exp_metadata$V5<-paste0(dataset,".mgf")
  	colnames(exp_metadata) <- c("file", "database", "sample", "fraction", "mgfFile")
  	cat("\n")
	resultObj <- get(load(rdaFilePath))
	#resultTable <- resultObj$resProtGroups

	#resultTable <- resultTable[resultTable$database == 'T',]
	resultTable <- resultObj[resultObj$database=="T",]
	#ReadMGFFile por cada fichero mgf
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
	if(search_engine=="comet"){
		index<-c()
		list_unique_Scannumber<-paste(unique(resultTablef$ScanNumber))
		cat("length of list_unique_Scan","\n")
		length(list_unique_Scannumber)
		cat("\n")
		for (i in 1:length(list_unique_Scannumber)){

			ScanNr<-paste(list_unique_Scannumber[i])
			cat(ScanNr,i,"\n")
			for (j in 1:length(mgfList)){
				if ((ScanNr)==str_extract(mgfList[[j]][6], '(?<=SCANS\\=)[0-9]+')){  #en datos de felix es asi. en publicos hay que modificar esta linea ( 6 por 2 y SCANS por scan)
				# if ((ScanNr)==paste(str_extract(sapply(strsplit(mgfList[[j]][2]," "),"[" ,5),"[0-9]+"))){
					index<-c(index,j)
				}
			}
		}


		mgfHits <- as.numeric(paste(index))   #comet
	    mgfHits <- as.numeric(paste(resultTablef$mgfIndex))   #mascot

		hitMgfIndex <- as.numeric(paste(unique(mgfHits))) #mascot
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
	outputMGFFilePath <- paste(substr(outputMGFFilePath, 1, nchar(outputMGFFilePath) - 4), "_COMET.mgf", sep="")
	WriteMGFFile(noHitMGF, outputMGFFilePath)
	}else{
		mgfHits <- as.numeric(paste(resultTable$mgfIndex))   #mascot

		hitMgfIndex <- as.numeric(paste(unique(mgfHits))) #mascot
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
}

}


index<-c()
for (i in 1:nrow(resultTablef)){

	ScanNr<-paste(resultTablef[i,]$ScanNumber)
	cat(ScanNr,"\n")
	for (j in 1:length(mgfList)){
		if ((ScanNr)==(str_extract(sapply(strsplit(mgfList[[i]][2]," "),"[" ,5),"[0-9]+"))){
			index<-c(index,j)
		}
	}
}


##############################hasta aqui :)
else if (identical(fase,"2")) {
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

#cometsearch
index<-c()
for (i in 1:nrow(resultTablef)){

	ScanNr<-paste(resultTablef[i,]$ScanNumber)
	cat(ScanNr,"\n")
	for (j in 1:length(mgfList)){
		if ((ScanNr)==(str_extract(sapply(strsplit(mgfList[[i]][2]," "),"[" ,5),"[0-9]+"))){
			index<-c(index,j)
		}
	}
}
