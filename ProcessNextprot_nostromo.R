#!/usr/bin/Rscript

args=(commandArgs(TRUE))
nextprot_folder <- args[1]				# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/23_neXtprot_20170801_Nov17/

library(Biostrings)
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesEli.R")
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesShotgun.R")
library(gdata)

chrList <- c(1:22, "X", "Y", "MT")
NextProtXChr <- data.frame()

for (i in 1:length(chrList)) {
	cat(i, "\n")
	filetmp <- paste0(nextprot_folder, "chr_reports/nextprot_chromosome_", chrList[i], ".txt")
	con <- file(filetmp, "r", blocking = FALSE)
	geneListv2 <- readLines(con)
	geneListv2.lst <- strsplit(geneListv2[20:(length(geneListv2)-6)], " ")
	geneListv2.lst <- lapply(geneListv2.lst, FUN = function(x) x[x != ""])
	# take the gene name
	geneListv2.genename <- t(as.data.frame(lapply(geneListv2.lst, FUN = function(x) x[[1]])))
	# take the nextprotid
	geneListv2.NXID <- t(as.data.frame(lapply(geneListv2.lst, FUN = function(x) x[[2]])))
	# take the PE value
	geneListv2.PE <- t(as.data.frame(lapply(geneListv2.lst, FUN = function(x) x[[7]])))

	if (i == 1) {
		nextProtXChr <- data.frame(chrList[i], geneListv2.NXID, geneListv2.genename, geneListv2.PE)
	} else {
		nextProtXChr <- rbind(nextProtXChr, data.frame(chrList[i], geneListv2.NXID, geneListv2.genename, geneListv2.PE))
	}
}
rownames(nextProtXChr) <- NULL
colnames(nextProtXChr) <- c("Chr", "NextprotID", "GeneName", "ProteinEvidence")
tmp <- paste(nextProtXChr$ProteinEvidence)
tmp[tmp == "protein"] = "PE1"
tmp[tmp == "transcript"] = "PE2"
tmp[tmp == "homology"] = "PE3"
tmp[tmp == "predicted"] = "PE4"
tmp[tmp == "uncertain"] = "PE5"
nextProtXChr$PE <- tmp
nextProtXChr$Missing <- (nextProtXChr$PE %in% c("PE2", "PE3", "PE4"))*1

write.table(nextProtXChr, file= paste0(nextprot_folder, "NextProt.txt"), col.names=TRUE, row.names=FALSE, sep="\t")
save(nextProtXChr, file= paste0(nextprot_folder, "nextProt.RData"))

nextprotXensp <- read.csv2(file = paste0(nextprot_folder, "nextprot_ensp.txt"), header = FALSE, sep = "\t")

nextprotXensp$NextprotID  <- paste(unlist(lapply(strsplit(paste(nextprotXensp[,1]), "-"), FUN = function(x) x[1])))

nextProtXChrXENSP <- merge(nextProtXChr, nextprotXensp, by.x = 2, by.y = 3, all.x = TRUE)

colnames(nextProtXChrXENSP)[7] <- "IsoformID"
colnames(nextProtXChrXENSP)[8] <- "ENSP"

save(nextProtXChrXENSP, file = paste0(nextprot_folder, "nextProtXChrXENSP.RData"))

nextprotXensg <- read.csv2(file = paste0(nextprot_folder, "nextprot_ensg.txt"), header = FALSE, sep = "\t")

nextProtXChrXENSPXENSG <- merge(nextProtXChrXENSP, nextprotXensg, by.x = 1, by.y = 1, all.x = TRUE)

colnames(nextProtXChrXENSPXENSG)[9] <- "ENSG"

save(nextProtXChrXENSPXENSG, file = paste0(nextprot_folder, "nextProtXChrXENSPXENSG.RData"))
