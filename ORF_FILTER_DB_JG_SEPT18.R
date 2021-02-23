
library(Biostrings)
library(data.table)
library(doBy)
library(gdata)
library(plyr)
library(ggplot2)
library(reshape2)
library(scales)

source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesVikv2.R")
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesEli.R")
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesAlba.R")

#### FIRST we make digestion of PROTEINDB

perl /opt/proteogest/proteogest.pl -i PROTEINDB_ORFEOME.fasta -d -a -c trypsin -R M -W 15.99 -G 0 -L 1


######
#LOAD PROTEINDB DIGESTED

PROTEIN_DB_dig_930FileName<-"PROTEINDB_ORFEOME_simple.txt"

con <- file(PROTEIN_DB_dig_930FileName, open="r")
PROTEIN_DB_dig_930_raw <- readLines(con)
close(con)

PROTEIN_DB_dig_930_protsAll <- PROTEIN_DB_dig_930_raw[grep("Protein Name", PROTEIN_DB_dig_930_raw)]
PROTEIN_DB_dig_930_protsAll <- t(data.frame(strsplit(unlist(lapply(strsplit(PROTEIN_DB_dig_930_protsAll, ": "), FUN = function(x) x[2])), "\\|")))
colnames(PROTEIN_DB_dig_930_protsAll) <- NULL
rownames(PROTEIN_DB_dig_930_protsAll) <- NULL

library(stringr)
PROTEIN_DB_dig_930_npeptidesAll <- PROTEIN_DB_dig_930_raw[grep("Number", PROTEIN_DB_dig_930_raw)]

PROTEIN_DB_dig_930_NpeptidesAll <- as.integer(str_extract(PROTEIN_DB_dig_930_npeptidesAll, '[0-9]+$'))

tmp <- PROTEIN_DB_dig_930_raw[-grep("Protein Name", PROTEIN_DB_dig_930_raw)]
tmp <- tmp[tmp != ""]
tmp <- tmp[-grep("Number", tmp)]
tmp <- tmp[-grep("No", tmp)]
tmp <- tmp[-c(1,2)]
tmp2 <- strsplit(tmp, " ")

tmp2 <- lapply(tmp2, FUN = function(x) x[x != ""])

PROTEIN_DB_dig_930_peptidesAll.df <- data.frame("PepNo" = unlist(lapply(tmp2, FUN = function(x) x[1])), "Range" = unlist(lapply(tmp2, FUN = function(x) x[2])), "IsotopicMass" = unlist(lapply(tmp2, FUN = function(x) x[3])), "AverageMass" = unlist(lapply(tmp2, FUN = function(x) x[4])), "Peptide" = unlist(lapply(tmp2, FUN = function(x) x[5])))


#ESTO NO ME HACE FALTA PoRQUE NO NECESITO identificar a que proteina viene cada peptido, solo comparar las secuencias

# npepind <- c()
# for (i in 1:length(PROTEIN_DB_dig_930_NpeptidesAll)) {
#  	npepind <- c(npepind, rep(i, PROTEIN_DB_dig_930_NpeptidesAll[i]))
# }

#PROTEIN_DB_dig_930_peptidesXProtAll.df <- data.frame(PROTEIN_DB_dig_930_protsAll[npepind,], PROTEIN_DB_dig_930_peptidesAll.df)

PROTEIN_DB_dig_930_peptidesAll.df$Name<-c(1:(nrow(PROTEIN_DB_dig_930_peptidesAll.df)))

#colnames(PROTEIN_DB_dig_930_peptidesXProtAll.df)[1] <- c("transcript_ID")
#PROTEIN_DB_dig_930_peptidesXProtAll.df$Range <- paste("'", PROTEIN_DB_dig_930_peptidesXProtAll.df$Range, sep = "")
#PROTEIN_DB_dig_930_peptidesXProtAll.df$Transcript_id_dig <- paste0(PROTEIN_DB_dig_930_peptidesXProtAll.df$transcript_ID, "-",PROTEIN_DB_dig_930_peptidesXProtAll.df$PepNo)

n_AA = 8
max_n_AA = 31
PROTEIN_DB_dig_930_peptides_filterAA <- PROTEIN_DB_dig_930_peptidesAll.df[nchar(paste(PROTEIN_DB_dig_930_peptidesAll.df$Peptide)) > n_AA,]

PROTEIN_DB_dig_930_peptides_filterAA <- PROTEIN_DB_dig_930_peptides_filterAA[nchar(paste(PROTEIN_DB_dig_930_peptides_filterAA$Peptide)) < max_n_AA,]



PROTEIN_DB_dig_930_peptides_filterAA_unique <- PROTEIN_DB_dig_930_peptides_filterAA[,c(5,6)]
colnames(PROTEIN_DB_dig_930_peptides_filterAA_unique)[1:2]<-c("peptideSeq","translation_id2")

#collapse peptideSeq values and concatenate translation_id values
library(plyr)
PROTEIN_DB_translated_unique_df <- ddply(PROTEIN_DB_dig_930_peptides_filterAA_unique, .(peptideSeq), summarise, paste(translation_id2, collapse = "|"))
#26358
####
PROTEINDB<-PROTEIN_DB_translated_unique_df

colnames(PROTEIN_DB_translated_unique_df) <- c("seq", "name")

save(PROTEIN_DB_translated_unique_df, file="PROTEIN_DB_translated_unique_df_digested.rda")

#LOAD UNKkNOWN DB


UNKNOWN_DB_filename<-"/home/nostromo/data/pepe/04_DataFelix_scORF_MASCOT_MAY18/All_in_one/Known_lnc/ORF_NONCODING_DB.fasta"
#UNKNOWN_DB_filename<-"~/dato-activo/03_Analysis/jgonzalez69/11_ORFEOME_SEPT18/ORF_UNKNOWN_FILTERED.fasta"
#UNKNOWN_DB<-readAAStringSet(UNKNOWN_DB_filename, format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
UNKNOWN_DB<-readAAStringSet(UNKNOWN_DB_filename, format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)

UNKNOWN_DF<-data.frame("seq"=paste(UNKNOWN_DB), "name"=names(UNKNOWN_DB))
### COMPARAR LAS DOS BD Y QUEDARSE SOLO CON LOS PEPTIDOS UNKNOWN QUE NO ESTAN EN PROTEINDB
peptides_UNKNOWN_not_in_PROTEINDB<-UNKNOWN_DF[!(  UNKNOWN_DF$seq %in% PROTEIN_DB_translated_unique_df$seq),]

UNKOWN<-peptides_UNKNOWN_not_in_PROTEINDB
#####GENERAR FASTA TARGET DEFINITIVO Y SU DECOY

UNKNOWN_FINAL_db <- AAStringSet(peptides_UNKNOWN_not_in_PROTEINDB$seq)

source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesShotgun.R")

DB_unknown_TargetFileName = "ORF_UNKNOWN_FILTERED.fasta"

names(UNKNOWN_FINAL_db) <- peptides_UNKNOWN_not_in_PROTEINDB$name
writeXStringSet(UNKNOWN_FINAL_db,DB_unknown_TargetFileName)

db_target<-UNKNOWN_FINAL_db
db_decoy <- db_target
names(db_decoy) <- gsub("T", "D", names(db_decoy))
seq <- paste(as.character(db_decoy))

seq_decoy <- sapply(seq,  FUN=function(x) .strPseudoreverseDecoy(x))
seq_decoy_fa <- AAStringSet(seq_decoy)

DB_unknown_DecoyFileName = "ORF_UNKNOWN_FILTERED_D.fasta"

names(seq_decoy_fa) <-names(db_decoy)
writeXStringSet(seq_decoy_fa, DB_unknown_DecoyFileName)

###### CARGAR CODING novel

CODING_DB_FILENAME<-"/home/nostromo/data/pepe/04_DataFelix_scORF_MASCOT_MAY18/All_in_one/ORF_unknown_translations_Coding_novel_digested.fasta"
CODING_DB_FILENAME<-"~/dato-activo/03_Analysis/jgonzalez69/11_ORFEOME_SEPT18/ORF_CODING_NOVEL_FILTERED.fasta"
CODING_NOVEL_DB<-readAAStringSet(CODING_DB_FILENAME, format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)

#1)COMPARAR CONTRA PROTEIN_DB Y QUEDARSE CON LOS QUE NO ESTEN EN PROTEIN DB
CODING_NOVEL_DF<-data.frame("seq"=paste(CODING_NOVEL_DB), "name"=names(CODING_NOVEL_DB))
#113111
peptides_CODING_NOVEL_not_in_PROTEINDB<-CODING_NOVEL_DF[!(CODING_NOVEL_DF$seq %in% PROTEIN_DB_translated_unique_df$seq ),]
#113042
#2)COMPARAR CONTRA UNKNOWN DB Y QUEDARSE CON LO QUE TAMPOCO ESTEN EN UNKNOWN DB
peptides_CODING_not_in_PROTEINDB_not_in_UNKNOWN_DB<- peptides_CODING_NOVEL_not_in_PROTEINDB[!( peptides_CODING_NOVEL_not_in_PROTEINDB$seq %in% peptides_UNKNOWN_not_in_PROTEINDB$seq ),]
#110043
# GENERAR FASTA TARGET Y DECOY DEFINITVOS

CODING_NOVEL_FINAL_db <- AAStringSet(peptides_CODING_not_in_PROTEINDB_not_in_UNKNOWN_DB$seq)

CODING<-peptides_CODING_not_in_PROTEINDB_not_in_UNKNOWN_DB

DB_unknown_TargetFileName = "ORF_CODING_NOVEL_FILTERED.fasta"

names(CODING_NOVEL_FINAL_db) <- peptides_CODING_not_in_PROTEINDB_not_in_UNKNOWN_DB$name
writeXStringSet(CODING_NOVEL_FINAL_db,DB_unknown_TargetFileName)

db_target<-CODING_NOVEL_FINAL_db
db_decoy <- db_target
names(db_decoy) <- gsub("T", "D", names(db_decoy))
seq <- paste(as.character(db_decoy))

seq_decoy <- sapply(seq,  FUN=function(x) .strPseudoreverseDecoy(x))
seq_decoy_fa <- AAStringSet(seq_decoy)

DB_unknown_DecoyFileName = "ORF_CODING_NOVEL_FILTERED_D.fasta"

names(seq_decoy_fa) <- names(db_decoy)
writeXStringSet(seq_decoy_fa, DB_unknown_DecoyFileName)

###### CARGAR NONCODING_NOVEL_DB

non_CODING_DB_FILENAME<-"/home/nostromo/data/pepe/04_DataFelix_scORF_MASCOT_MAY18/All_in_one/ORF_unknown_translations_nonCoding_novel_digested.fasta"
non_CODING_DB_FILENAME<-"~/dato-activo/03_Analysis/jgonzalez69/11_ORFEOME_SEPT18/ORF_NON_CODING_NOVEL_FILTERED.fasta"

non_CODING_NOVEL_DB<-readAAStringSet(non_CODING_DB_FILENAME, format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)

#1)COMPARAR CONTRA PROTEIN_DB Y QUEDARSE CON LOS QUE NO ESTEN EN PROTEIN DB
non_CODING_NOVEL_DF<-data.frame("seq"=paste(non_CODING_NOVEL_DB), "name"=names(non_CODING_NOVEL_DB))

peptides_NON_CODING_NOVEL_not_in_PROTEINDB<-non_CODING_NOVEL_DF[!( non_CODING_NOVEL_DF$seq %in% PROTEINDB$seq ),]

#2)COMPARAR CONTRA UNKNOWN DB Y QUEDARSE CON LO QUE TAMPOCO ESTEN EN UNKNOWN DB
peptides_NON_CODING_not_in_PROTEINDB_not_in_UNKNOWN_DB<- peptides_NON_CODING_NOVEL_not_in_PROTEINDB[!(peptides_NON_CODING_NOVEL_not_in_PROTEINDB$seq %in%  UNKOWN$seq  ),]

#3)COMPARAR CONTRA CODING NOVEL DB Y QUEDARSE CON LO QUE TAMPOCO ESTÁ aqui

peptides_NON_CODING_not_in_PROTEINDB_not_in_UNKNOWN_DB_not_in_CODING_NOVEL<- peptides_NON_CODING_not_in_PROTEINDB_not_in_UNKNOWN_DB[!(peptides_NON_CODING_not_in_PROTEINDB_not_in_UNKNOWN_DB$seq %in% CODING$seq),]

#generar FASTA noncoding definitivo TARGET Y DECOY_

non_CODING_NOVEL_FINAL_db <- AAStringSet(peptides_NON_CODING_not_in_PROTEINDB_not_in_UNKNOWN_DB_not_in_CODING_NOVEL$seq)

DB_unknown_TargetFileName = "ORF_NON_CODING_NOVEL_FILTERED.fasta"

names(non_CODING_NOVEL_FINAL_db) <- peptides_NON_CODING_not_in_PROTEINDB_not_in_UNKNOWN_DB_not_in_CODING_NOVEL$name
writeXStringSet(non_CODING_NOVEL_FINAL_db,DB_unknown_TargetFileName)

db_target<-non_CODING_NOVEL_FINAL_db
db_decoy <- db_target
names(db_decoy) <- gsub("T", "D", names(db_decoy))
seq <- paste(as.character(db_decoy))

seq_decoy <- sapply(seq,  FUN=function(x) .strPseudoreverseDecoy(x))
seq_decoy_fa <- AAStringSet(seq_decoy)

DB_unknown_DecoyFileName = "ORF_NON_CODING_NOVEL_FILTERED_D.fasta"

names(seq_decoy_fa) <- names(db_decoy)
writeXStringSet(seq_decoy_fa, DB_unknown_DecoyFileName)
