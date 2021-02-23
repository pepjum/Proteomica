args=(commandArgs(TRUE))
fileroot <- args[1]

initialfile <- paste0(fileroot, "_simple.txt")
con <- file(initialfile, open="r")
db_raw <- readLines(con)
close(con)

db_protsAll <- db_raw[grep("Protein Name", db_raw)]
db_protsAll <- t(data.frame(strsplit(unlist(lapply(strsplit(db_protsAll, ": "), FUN = function(x) x[2])), "\\|")))
colnames(db_protsAll) <- NULL
rownames(db_protsAll) <- NULL

db_npeptidesAll <- db_raw[grep("Number", db_raw)]
db_NpeptidesAll <- as.numeric(paste(unlist(lapply(strsplit(db_npeptidesAll, " = "), FUN = function(x) x[2]))))

tmp <- db_raw[-grep("Protein Name", db_raw)]
tmp <- tmp[tmp != ""]
tmp <- tmp[-grep("Number", tmp)]
tmp <- tmp[-grep("No", tmp)]
tmp <- tmp[-c(1,2)]
tmp2 <- strsplit(tmp, " ")

tmp2 <- lapply(tmp2, FUN = function(x) x[x != ""])

db_peptidesAll.df <- data.frame("PepNo" = unlist(lapply(tmp2, FUN = function(x) x[1])), "Range" = unlist(lapply(tmp2, FUN = function(x) x[2])), "IsotopicMass" = unlist(lapply(tmp2, FUN = function(x) x[3])), "AverageMass" = unlist(lapply(tmp2, FUN = function(x) x[4])), "Peptide" = unlist(lapply(tmp2, FUN = function(x) x[5])))

save(db_peptidesAll.df, file = paste0(fileroot, "_peptidesAll.df.rda"))

save(db_NpeptidesAll, file = paste0(fileroot, "_NpeptidesAll.rda"))

npepind <- c()
for (i in 1:length(db_NpeptidesAll)) {
	npepind <- c(npepind, rep(i, db_NpeptidesAll[i]))
}

db_peptidesXProtAll.df <- data.frame(db_protsAll[npepind,], db_peptidesAll.df)
colnames(db_peptidesXProtAll.df)[1:3] <- c("ProteinID", "Name", "Description")
db_peptidesXProtAll.df$Range <- paste("'", db_peptidesXProtAll.df$Range, sep = "")
write.table(db_peptidesXProtAll.df, paste0(fileroot, "_simple_digested.txt"), col.names=TRUE, sep="\t", quote=FALSE)

save(db_peptidesXProtAll.df, file = paste0(fileroot, "_peptidesXProtAll.df.rda"))

n_AA = 8
max_n_AA = 51
db_peptidesXProtAll_filterAA <- db_peptidesXProtAll.df[nchar(paste(db_peptidesXProtAll.df$Peptide)) > n_AA,]
db_peptidesXProtAll_filterAA <- db_peptidesXProtAll_filterAA[nchar(paste(db_peptidesXProtAll_filterAA$Peptide)) < max_n_AA,]

save(db_peptidesXProtAll_filterAA, file = paste0(fileroot, "_peptidesXProtAll_filterAA.rda"))
# luego cogemos el df de la bd y quitamos la columna secuencia, y cogemos las filas que contengan transcritos únicos.
db_peptidesXProtAll_filterAA_unique <- unique(db_peptidesXProtAll_filterAA[,c("ProteinID","Peptide")])

write.table(db_peptidesXProtAll_filterAA_unique$Peptide, file = paste0(fileroot, "_peptidesXProtAll_filterAA_unique_peptides.txt"), col.names=FALSE, row.names = FALSE, sep="\t", quote=FALSE)

inputFile <- paste0(fileroot, "_peptidesXProtAll_filterAA_unique_peptides.txt")
outputFile <- paste0(fileroot, "_peptidesXProtAll_filterAA_unique_peptides_unique.txt")
y <- paste('sort -k1,1 ', inputFile,' | uniq -u > ', outputFile, sep = "")
system(y)
y <- NULL

db_peptidesXProtAll_filterAA_unique_peptides_unique <- read.csv2(paste0(fileroot, "_peptidesXProtAll_filterAA_unique_peptides_unique.txt"), header = FALSE, sep = "\t", fill = TRUE)

db_peptidesXProtAll_filterAA_disc <- merge(db_peptidesXProtAll_filterAA_unique_peptides_unique, db_peptidesXProtAll_filterAA, by.x = "V1", by.y = "Peptide", all.x = TRUE)
colnames(db_peptidesXProtAll_filterAA_disc)[1] <- "Peptide"

save(db_peptidesXProtAll_filterAA_disc, file = paste0(fileroot, "_peptidesXProtAll_filterAA_disc.rda"))
