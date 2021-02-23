#!/usr/bin/Rscript

args=(commandArgs(TRUE))
txtFileName <- args[1] #file Shotgun___.txt
dataset <- args[2] # NAME of sample
DBrdaFileName <- args[3] #0
decoy_id <- args[4]   #DECOY
proteinAnalysis <- args[5] # 0: psmFDR, 1: protFDRMAYU, 2: pepFDR, 3: protFDR
search_engine <- args[6] # 1: mascot, 2: comet, 3:tandem, 4: omssa
minn_AA <- args[7]   # 9
maxn_AA <- args[8]  # 30
dummyCodeFileName <- args[9]

cat(txtFileName)
cat("\n")
cat(dataset)
cat("\n")
cat(DBrdaFileName)
cat("\n")
cat(decoy_id)
cat("\n")
cat(proteinAnalysis)
cat("\n")
cat(search_engine)
cat("\n")
cat(minn_AA)
cat("\n")
cat(maxn_AA)
cat("\n")
cat(dummyCodeFileName)


library(Biostrings)
library(doBy)
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesShotgun.R")
cat("Starting the script \n")

psmFDRvalue = 0.01
protFDRvalue = 0.01
pepScore="score"
#decoy_id="decoy"
databaseLoaded <- 0
if (!is.na(dummyCodeFileName))
{
	cat("dummyCodeFileName ")
	cat(dummyCodeFileName)
	cat("...\n")
}
if (DBrdaFileName == 0)
{
	if (proteinAnalysis != 1)
	{
		error <- 1
	}
	error <- 0
}else
{
	error <- 0
	databaseLoaded <- 0
	cat("database: ")
	cat(DBrdaFileName)
	cat(" \n")
	if (nchar(DBrdaFileName) != 0)
	{
		tryCatch(
		{
			load(DBrdaFileName)
			if (exists("db_massfilter") == FALSE || exists("db_target_length") == FALSE)
			{
				message("Error: DBrdaFileName must contain two objects: db_massfilter and db_target_length")
				databaseLoaded = 0
			} else {
			# db_massfilter and db_target_length are loaded
			databaseLoaded = 1
			}
		},
		error=function(e)
		{
			message("Error:")
			message(e)
			databaseLoaded = 0
			error = 1
			return(NA)
		})
	  }

}
if (error == 0)
{
	prot_id = "ProteinAccession"
	pep_id="Pep_scan"
	pep_col="PeptideSeq"
	cat("Reading the file \n")
	currentDir = dirname(txtFileName)

	exp_metadata <- read.table(txtFileName, header = FALSE, sep = " ")
	colnames(exp_metadata) <- c("file", "database", "sample", "fraction")
	exp_metadata$datfile <- lapply(strsplit(paste(exp_metadata$file), "/"), FUN = function(x) x[length(x)])

	nExperiments <- length(unique(exp_metadata[,3]))
	nExperimentNames <- unique(exp_metadata[,3])

	pos <- match(dataset,nExperimentNames)

	if(proteinAnalysis == 1 && databaseLoaded == 0)
	{
		stop('Error: The database could not being loaded. \n You can either load a rda file using the DBrdaFileName parameter or set the DBTargetFileName and DBDecoyFileName to load the DB in the function. DBrdaFileName must contain two objects: db_massfilter and db_target_length. \n Try again')
	}else
	{
	   		experimentData = subset(exp_metadata, sample==nExperimentNames[pos])

				if (experimentData[1,2] == "TD") {
					#databases are concatenated
					concat_decoy = 1
				}else{
					#databases are not concatenated
					concat_decoy = 0
				}
			dataPSMMat <- data.frame()
				for (j in 1:nrow(experimentData))
				{
					if (j == 1)
					{
						if (search_engine == 1)
						{
							dataPSMMat <-  data.frame("file" = paste(experimentData[1,1]), parseMascotStandardOutput(paste(experimentData[1,1])), "database" = paste(experimentData[1,2]), "sample" = paste(experimentData[1,3]), "fraction" = paste(experimentData[1,4]), "datfile" = paste(experimentData[1,5]))
						}
						else if (search_engine == 2)
						{
							dataPSMMat <-  data.frame("file" = paste(experimentData[1,1]), parseCometStandardOutput(paste(experimentData[1,1])), "database" = paste(experimentData[1,2]), "sample" = paste(experimentData[1,3]), "fraction" = paste(experimentData[1,4]), "datfile" = paste(experimentData[1,5]))
						}
						else if (search_engine == 3)
						{
							dataPSMMat <-  data.frame("file" = paste(experimentData[1,1]), parseXTandemStandardOutput(paste(experimentData[1,1])), "database" = paste(experimentData[1,2]), "sample" = paste(experimentData[1,3]), "fraction" = paste(experimentData[1,4]), "datfile" = paste(experimentData[1,5]))
						}
						else if (search_engine == 4)
						{
							dataPSMMat <-  data.frame("file" = paste(experimentData[1,1]), parseOmssaStandardOutput(paste(experimentData[1,1])), "database" = paste(experimentData[1,2]), "sample" = paste(experimentData[1,3]), "fraction" = paste(experimentData[1,4]), "datfile" = paste(experimentData[1,5]))
						}
						else
						{
							cat("Error. Specify search_engine parameter correctly. 1 for mascot. 2 for comet. 3 for tandem. 4 for omssa")
							break
						}
					}else{
						if (search_engine == 1)
						{
							dataPSMMat <- rbind(dataPSMMat, data.frame("file" = paste(experimentData[j,1]), parseMascotStandardOutput(paste(experimentData[j,1])), "database" = paste(experimentData[j,2]), "sample" = paste(experimentData[j,3]), "fraction" = paste(experimentData[j,4]), "datfile" = paste(experimentData[j,5])))
						}
						else if (search_engine == 2)
						{
							dataPSMMat <- rbind(dataPSMMat, data.frame("file" = paste(experimentData[j,1]), parseCometStandardOutput(paste(experimentData[j,1])), "database" = paste(experimentData[j,2]), "sample" = paste(experimentData[j,3]), "fraction" = paste(experimentData[j,4]), "datfile" = paste(experimentData[j,5])))
						}
						else if (search_engine == 3)
						{
							dataPSMMat <- rbind(dataPSMMat, data.frame("file" = paste(experimentData[j,1]), parseXTandemStandardOutput(paste(experimentData[j,1])), "database" = paste(experimentData[j,2]), "sample" = paste(experimentData[j,3]), "fraction" = paste(experimentData[j,4]), "datfile" = paste(experimentData[j,5])))
						}
						else if (search_engine == 4)
						{
							dataPSMMat <- rbind(dataPSMMat, data.frame("file" = paste(experimentData[j,1]), parseOmssaStandardOutput(paste(experimentData[j,1])), "database" = paste(experimentData[j,2]), "sample" = paste(experimentData[j,3]), "fraction" = paste(experimentData[j,4]), "datfile" = paste(experimentData[j,5])))
						}
					}
				}

			if  (search_engine == 1)
			{
				dataPSMMat$search_engine <- "mascot"
				search_engine_str = "mascot"
			} else if (search_engine == 2)
			{
				dataPSMMat$search_engine <- "comet"
				search_engine_str = "comet"
			} else if (search_engine == 3)
			{
				dataPSMMat$search_engine <- "tandem"
				search_engine_str = "tandem"
			} else if (search_engine == 4)
			{
				dataPSMMat$search_engine <- "omssa"
				search_engine_str = "omssa"
			}

			save(dataPSMMat, file=paste(currentDir, "/results_Peptides_", search_engine_str, "_", dataset, "_dataPSMMat.rda", sep = ""))

			cat("Filtering AA length ...\n")
			thr_AA <- as.numeric(paste(minn_AA)) - 1
			dataPSMMat_filterAA <- dataPSMMat[nchar(paste(dataPSMMat$PeptideSeq)) > thr_AA,]
			thr_AA <- as.numeric(paste(maxn_AA)) + 1
			dataPSMMat_filterAA <- dataPSMMat_filterAA[nchar(paste(dataPSMMat_filterAA$PeptideSeq)) < thr_AA,]
			save(dataPSMMat_filterAA, file=paste(currentDir, "/results_Peptides_", search_engine_str, "_", dataset, "_dataPSMMat_filterAA.rda", sep = ""))
			cat("Calculating psmFDR and filtering ...\n")
			dataPSMMat_filterAA_psmFDR <- psmFDR(dataPSMMat_filterAA, pepScore=pepScore, decoy_id=decoy_id, concat_decoy = concat_decoy, prot_id = prot_id)
			dataPSMMat_filterAA_psmFDR_Filter <- dataPSMMat_filterAA_psmFDR[dataPSMMat_filterAA_psmFDR$psmFDR < psmFDRvalue,]
			save(dataPSMMat_filterAA_psmFDR_Filter,file=paste(currentDir, "/results_Peptides_", search_engine_str, "_", dataset, "_dataPSMMat_filterAA_PSMFDR_filter.rda", sep = "") )
			if (proteinAnalysis == 1) # protFDRMAYU
			{
				tryCatch(
				{
					cat("Calculating MAYU protFDR and filtering ...\n")
					db_npep_massFilter_res <- db_massfilter[(names(db_massfilter) %in% paste(unique(dataPSMMat_filterAA_psmFDR_Filter[, prot_id])))]
					db_bySize_500 <- cut(as.numeric(db_npep_massFilter_res)[which(as.numeric(db_npep_massFilter_res)<500)], breaks=18)
					db_bySize_10000 <- cut(as.numeric(db_npep_massFilter_res)[which(as.numeric(db_npep_massFilter_res)>=500)], breaks=2)
					DB_partitions_res <- vector(mode = "numeric", length(db_npep_massFilter_res))
					DB_partitions_res[which(as.numeric(db_npep_massFilter_res)<500)] <- as.numeric(db_bySize_500)
					DB_partitions_res[which(as.numeric(db_npep_massFilter_res)>=500)] <- as.numeric(db_bySize_10000)+18
					names(DB_partitions_res) <- names(db_npep_massFilter_res)

					pepscan <- paste(dataPSMMat_filterAA_psmFDR_Filter$sample, dataPSMMat_filterAA_psmFDR_Filter$fraction, dataPSMMat_filterAA_psmFDR_Filter$Query,sep="_")
					dataPSMMat_filterAA_psmFDR_Filter$Pep_scan <- pepscan
					dataPSMMat_filterAA_psmFDR_Filter_protFDR <- protFDRMayu(dataPSMMat_filterAA_psmFDR_Filter, DB_partitions_res, length(db_massfilter), decoy_char=decoy_id, prot_id=prot_id, pep_id=pep_id, pep_score=pepScore)
					dataPSMMat_filterAA_psmFDR_Filter_protFDR_Res <- dataPSMMat_filterAA_psmFDR_Filter_protFDR[[2]]
					dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter <- dataPSMMat_filterAA_psmFDR_Filter_protFDR_Res[dataPSMMat_filterAA_psmFDR_Filter_protFDR_Res$prot_shFDR < protFDRvalue,]
					cat("dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter")
					cat(dim(dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter))
					# No hacer para comet
					#if (search_engine == 1)
					#{
						cat("Protein Grouping ...\n")
						dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping <- PAnalyzer(PID_res = dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter, pep_col=pep_col, prot_col=prot_id)
						cat("dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping")
						cat(dim(dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping))
						cat("...\n")
						cat(" conclusive \n ")
						conclusive <- unique(paste(dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping[dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping$PGClass == 1, prot_id]))
						if (length(grep(decoy_id, conclusive)) != 0) {
							conclusive_sindecoy <- conclusive[-grep(decoy_id, conclusive)]
						} else {
							conclusive_sindecoy <- conclusive
						}
						cat(" indistinguible \n ")
						indistinguible <- unique(paste(dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping[dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping$PGClass == 0.5, prot_id]))
						if (length(grep(decoy_id, indistinguible)) != 0) {
							indistinguible_sindecoy <- indistinguible[-grep(decoy_id, indistinguible)]
						} else {
							indistinguible_sindecoy <- indistinguible
						}
						cat(" ambiguous \n ")
						ambiguous <- unique(paste(dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping[dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping$PGClass == 0.25, prot_id]))
						if (length(grep(decoy_id, ambiguous)) != 0) {
							ambiguous_sindecoy <- ambiguous[-grep(decoy_id, ambiguous)]
						} else {
							ambiguous_sindecoy <- ambiguous
						}
						cat(" nonconclusive \n ")
						nonconclusive <- unique(paste(dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping[dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping$PGClass == 0, prot_id]))
						if (length(grep(decoy_id, nonconclusive)) != 0) {
							nonconclusive_sindecoy <- nonconclusive[-grep(decoy_id, nonconclusive)]
						} else {
							nonconclusive_sindecoy <- nonconclusive
						}
					#}
					if (!is.na(dummyCodeFileName))
					{
						tryCatch(
						{
							cat("loading dummy code table \n")
							load(dummyCodeFileName)
							prot_id_dummy = paste(prot_id, "Dummy", sep="")
							dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter[,prot_id_dummy] <- dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter[,prot_id]

							dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter[,prot_id] <- dummyCodeTable$name[match(dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter[,prot_id_dummy], dummyCodeTable$dummyCodes)]

							dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping[,prot_id_dummy] <- dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping[,prot_id]

							dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping[,prot_id] <- dummyCodeTable$name[match(dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping[,prot_id_dummy], dummyCodeTable$dummyCodes)]

							conclusive_sindecoy_dummy <- conclusive_sindecoy
							conclusive_sindecoy <- dummyCodeTable$name[match(conclusive_sindecoy_dummy, dummyCodeTable$dummyCodes)]
							indistinguible_sindecoy_dummy <- indistinguible_sindecoy
							indistinguible_sindecoy <- dummyCodeTable$name[match(indistinguible_sindecoy_dummy, dummyCodeTable$dummyCodes)]
							ambiguous_sindecoy_dummy <- ambiguous_sindecoy
							ambiguous_sindecoy <- dummyCodeTable$name[match(ambiguous_sindecoy_dummy, dummyCodeTable$dummyCodes)]
							nonconclusive_sindecoy_dummy <- nonconclusive_sindecoy
							nonconclusive_sindecoy	<- dummyCodeTable$name[match( nonconclusive_sindecoy_dummy, dummyCodeTable$dummyCodes)]
						},
						error=function(e)
						{
							message("Error with dummyCodeFileName:")
							message(e)
						})
					}

						dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping <- dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping[dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping$database == "T",]

						cat("writing ...\n")
						write.table(dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping[,-3], file=paste(currentDir, "/results_", search_engine_str, "_",  dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping$sample[1], ".txt", sep=''), quote = FALSE, sep = "\t", row.names = FALSE)

						results <- list(dataPSMMat_filterAA_psmFDR_Filter = dataPSMMat_filterAA_psmFDR_Filter, resProtGroups =dataPSMMat_filterAA_psmFDR_Filter_protFDR_Filter_Grouping, conclusive = conclusive_sindecoy, indistinguible = indistinguible_sindecoy, ambiguous = ambiguous_sindecoy, nonconclusive = nonconclusive_sindecoy, resProtNSAF = data_NSAF_QProtCounts)

						save(results, file=paste(currentDir, "/results_", search_engine_str, "_", dataPSMMat_filterAA_psmFDR_Filter$sample[1], ".rda", sep = ""))
						cat("Results save in: ")
						cat(paste(currentDir, "/results_", search_engine_str, "_", dataset, ".rda", "\n", sep = ""))
						cat("Object's name is: results \n")


					},
				error=function(e)
				{
					message("Error:")
					message(e)
					message("...\n")
					return(NA)
				})
			} else if (proteinAnalysis == 2) # pepFDR
			{
				cat("Calculating pepFDR and filtering ...\n")
				dataPSMMat_filterAA_psmFDR_Filter_pepFDR <- pepFDR(dataPSMMat_filterAA_psmFDR_Filter, pepScore=pepScore, concat_decoy = concat_decoy, pep_id = pep_col)
				dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter <- dataPSMMat_filterAA_psmFDR_Filter_pepFDR[dataPSMMat_filterAA_psmFDR_Filter_pepFDR$pepFDR < 0.01,]
				save(dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter, file = paste(currentDir, "/results_Peptides_", search_engine_str, "_", dataset, "_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter.rda", sep = ""))


			} else if (proteinAnalysis == 3) # pepFDR, protFDR besthit
			{
				cat("Calculating pepFDR and filtering ...\n")
				dataPSMMat_filterAA_psmFDR_Filter_pepFDR <- pepFDR(dataPSMMat_filterAA_psmFDR_Filter, pepScore=pepScore, concat_decoy = concat_decoy, pep_id = pep_col)
				dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter <- dataPSMMat_filterAA_psmFDR_Filter_pepFDR[dataPSMMat_filterAA_psmFDR_Filter_pepFDR$pepFDR < 0.01,]
				save(dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter, file = paste(currentDir, "/results_Peptides_", search_engine_str, "_", dataset, "_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter.rda", sep = ""))
				cat("Calculating protFDR and filtering ...\n")
				dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR <- protFDR(dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter, pepScore=pepScore, concat_decoy = concat_decoy, prot_id = prot_id)
				dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter <- dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR[dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR$protFDR < 0.01,]
				save(dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter, file=paste(currentDir, "/results_Peptides_", search_engine_str, "_", dataset, "_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter.rda", sep = ""))

			}
			else
			{
					dataPSMMat_filterAA_psmFDR_Filter <- dataPSMMat_filterAA_psmFDR_Filter[dataPSMMat_filterAA_psmFDR_Filter$database == "T",]

					cat("writing ...\n")
					write.table(dataPSMMat_filterAA_psmFDR_Filter, file=paste(currentDir, "/results_Peptides_", search_engine_str, "_", dataset, "_dataPSMMat_filterAA_psmFDR_Filter.txt", sep=''), quote = FALSE, sep = "\t", row.names = FALSE)
					cat("Saving ...\n")

					cat("Results save in: ")
					cat(paste(currentDir, "/results_Peptides_", search_engine_str, "_", dataset, "_dataPSMMat.rda", sep = ""))
					cat("Object's name is: dataPSMMat \n")
					cat(paste(currentDir, "/results_Peptides_", search_engine_str, "_", dataset, "_dataPSMMat_filterAA.rda", sep = ""))
					cat("Object's name is: dataPSMMat_filterAA \n")
					cat(paste(currentDir, "/results_Peptides_", search_engine_str, "_", dataset, "_dataPSMMat_filterAA_psmFDR_Filter.rda", "\n", sep = ""))
					cat("Object's name is: dataPSMMat_filterAA_psmFDR_Filter \n")
			}


	  }

	cat("Ending the script \n")
}else
{
	message("Error: Check the parameters")
	message("...\n")
	return(NA)
}
