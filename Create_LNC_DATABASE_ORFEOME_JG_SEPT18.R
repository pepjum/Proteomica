library(Biostrings)
library(data.table)
library(doBy)
library(gdata)
library(plyr)
library(ggplot2)
library(reshape2)
library(scales)

setwd("/home/nostromo/data/pepe/04_DataFelix_scORF_MASCOT_MAY18/All_in_one/Known_lnc/")

source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesVikv2.R")
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesEli.R")
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesAlba.R")

Gingeras_orig_fileName = "known_lnc_features.gtf"

Gingeras_orig <- read.table(Gingeras_orig_fileName, header = FALSE, sep = "\t")
#TAKE OUT THE FIRST ELEMENT OF SOME COLUMNS (starting from the 9th)
tmp <- t(apply(as.data.frame(Gingeras_orig[,9]), 1, parseMITRANSorig))
# Quito los TCONS contenidos en otro (contained_in)
Gingeras_original <- data.frame(Gingeras_orig[,1:8], tmp)

####known_only
colnames(Gingeras_original) <- c("chr", "database", "type", "start", "end", "c1", "strand", "c2","NA", "gene_id", "transcript_id", "exon_number", "gene_name", "oId", "ENSEMBL_transcript_id", "tstatus", "tss_id")

tmp <- NULL
rm(Gingeras_orig)
#save(Gingeras_original,file="nonCoding_novel_gtf.rda")
save(Gingeras_original,file="known_lnc_gtf.rda")

df_Gingeras<-data.frame("gene_id"=Gingeras_original$gene_id,"miTransc_chr"=Gingeras_original$chr,"G19_transcript_id"=Gingeras_original$transcript_id,"G19_oId"=Gingeras_original$oId,"miTransc_start"=Gingeras_original$start,"miTransc_end"=Gingeras_original$end,"miTransc_strand"=Gingeras_original$strand,"miTransc_transcript_id"=Gingeras_original$transcript_id,"miTransc_tstatus"=Gingeras_original$tstatus,"miTtrans_tss_id"=Gingeras_original$tss_id,"gene_name"=Gingeras_original$gene_name,"ENST"=Gingeras_original$ENSEMBL_transcript_id)

ORF_unknown_translations <- df_Gingeras
rm(df_Gingeras)

save(ORF_unknown_translations, file="ORF_unknown_translations_without_translated_known_lnc_features.rda")
#save(ORF_unknown_translations, file="ORF_unknown_translations_without_translated_Coding_novel.rda")


## THREE FRAME TRANSLATION (in the cluster)
library(GenomicRanges)
library(BSgenome)

ORF_unknown_translations_known_strand <- ORF_unknown_translations[-(grep("\\.", ORF_unknown_translations$miTransc_strand)),]
ORF_unknown_translations_unknown_strand <- ORF_unknown_translations[grep("\\.", ORF_unknown_translations$miTransc_strand),]
ORF_unknown_translations_known_strand$miTransc_strand <- factor(paste(ORF_unknown_translations_known_strand$miTransc_strand))

save(ORF_unknown_translations_known_strand, file= "ORF_unknown_translations_known_strand_lnc_features.rda")
#43671

save(ORF_unknown_translations_unknown_strand, file= "ORF_unknown_translations_unknown_strand_lnc_features.rda")
#7381


DBTargetFileName <-"known_lnc_features_reloaded.fa"


Gingeras_FASTA <- readAAStringSet(DBTargetFileName, format="fasta",nrec=-1L, skip=0L, use.names=TRUE)

nameFASTAF = names(Gingeras_FASTA)
sequence = paste(Gingeras_FASTA)
df_FASTA <- data.frame(nameFASTAF, sequence)

df_FASTA$TCONS<-sapply(strsplit(paste(df_FASTA$nameFASTAF)," "), "[",1)

df_TEMP_w_seq<-merge(ORF_unknown_translations,df_FASTA, by.x="miTransc_transcript_id", by.y="TCONS")
#51142 nrows
###separar el dataframe por strand

df_strand_plus<-df_TEMP_w_seq[which(df_TEMP_w_seq$miTransc_strand =="+"),]
#23032
df_strand_minus<-df_TEMP_w_seq[which(df_TEMP_w_seq$miTransc_strand =="-"),]
#20729
df_strand_unknown<-df_TEMP_w_seq[which(df_TEMP_w_seq$miTransc_strand =="."),]
#7381
#### HAY que hacer reversa complementaria de las secuencias con strand "-"
seqs_plus<-df_strand_plus$sequence
seqs_minus<-df_strand_minus$sequence
seqs_unknown<-df_strand_unknown$sequence

seq_plus_obj<-DNAStringSet(seqs_plus)
seq_minus_obj <-DNAStringSet(seqs_minus)
seq_unknown_obj<-DNAStringSet(seqs_unknown)

seq_minus_ord <- reverseComplement(seq_minus_obj)

#elimino las secuencias originales
df_strand_minus$sequence<-NULL

#añado secuencias complemento_reverso
df_strand_minus$sequence<-paste(seq_minus_ord)

seqs_known_strand <- c(seq_plus_obj, seq_minus_ord)

tmp_plus<-data.frame("names"=df_strand_plus$G19_transcript_id,"sequences"=seqs_plus,"strand"=df_strand_plus$miTransc_strand)
#tmp_plus<-unique(tmp_plus)
#7714
tmp_minus<-data.frame("names"=df_strand_minus$G19_transcript_id,"sequences"=seq_minus_ord,"strand"=df_strand_minus$miTransc_strand)
#tmp_minus<-unique(tmp_minus)
#7094
tmp_unknown<-data.frame("names"=df_strand_unknown$G19_transcript_id,"sequences"=seqs_unknown,"strand"=df_strand_unknown$miTransc_strand)
#tmp_unknown<-unique(tmp_unknown)

ORF_unknown_translations_known_strand_translated <-rbind(tmp_plus,tmp_minus)
##43671
ORF_unknown_translations_unknown_strand_translated <- tmp_unknown
#7381
##Remove sequences contains NNN
Nindx <- grep('N', seqs_known_strand)
if(length(Nindx) > 0)
{
    seqs_known_strand <- seqs_known_strand[-Nindx]
	ORF_unknown_translations_known_strand_translated <- ORF_unknown_translations_known_strand_translated[-Nindx,]
}

save(ORF_unknown_translations_known_strand_translated, file =  "ORF_unknown_translations_known_strand_Known_lnc.rda")

save(ORF_unknown_translations_unknown_strand_translated, file =  "ORF_unknown_translations_unknown_strand_Known_lnc.rda")

save(seqs_known_strand, file =  "ORF_unknown_translations_known_strand-seqs_known_lnc.rda")

save(seq_unknown_obj, file =  "ORF_unknown_translations_unknown_strand-seqs_known_lnc.rda")

###TRADUCCION A PROTEINA

#three frame translation, as we know the strand
peptides_r1 <- translate(seqs_known_strand)
peptides_r2 <- translate(subseq(seqs_known_strand, start=2))
peptides_r3 <- translate(subseq(seqs_known_strand, start=3))

all_pep <- rbind(cbind("1", paste(ORF_unknown_translations_known_strand_translated[, 'names']), paste(ORF_unknown_translations_known_strand_translated[, 'strand']), paste(as.data.frame(peptides_r1)[,1])), cbind("2", paste(ORF_unknown_translations_known_strand_translated[, 'names']), paste(ORF_unknown_translations_known_strand_translated[, 'strand']), paste(as.data.frame(peptides_r2)[,1])), cbind("3", paste(ORF_unknown_translations_known_strand_translated[, 'names']), paste(ORF_unknown_translations_known_strand_translated[, 'strand']), paste(as.data.frame(peptides_r3)[,1])))



all_pep[,3] <- gsub("+", "f53", all_pep[,3], fixed=TRUE)
all_pep[,3] <- gsub("-", "f35", all_pep[,3], fixed=TRUE)

all_pep <- cbind(all_pep, paste(paste(all_pep[,2]), "_", paste(all_pep[,3]), "_", paste(all_pep[,1]), sep=""))

colnames(all_pep) <- c("3frame_id", "transcript_id", "strand", "peptideSeq", "translation_id")


save(all_pep, file =  "ORF_unknown_translations_known_strand-all_pep_known_lnc.rda")


### remove the sequence after the peptide stop codon. TAMBIEN ELIMINA LO QUE HAYA DELANTE DE LA METIONINA
all_pep[,4] <- unlist(lapply(strsplit(paste(all_pep[,4]), "\\*"), FUN = function(y) paste(unlist(lapply(1:length(y), function(x) if(regexpr("(M[A-Z]+)", y[x]) != "-1") substr(y[x], regexpr("(M[A-Z]+)", y[x]), nchar(y[x])))), collapse="*") ))

#remove empty sequences
all_pep <- all_pep[-(which(all_pep[,4] == "")),]
#646140

save(all_pep, file= "ORF_unknown_translations_known_strand-all_pep-f_known_lnc.rda")

#hacer el 6 frame de los que no sabemos el strand

seq_plus <- seq_unknown_obj
seq_minus <- reverseComplement(seq_plus)

	##Remove sequences contains NNN
Nindx_p <- grep('N', seq_plus)
Nindx_m <- grep('N', seq_minus)
ORF_unknown_translations_unknown_strand_translated_plus <- tmp_unknown
ORF_unknown_translations_unknown_strand_translated_minus <- tmp_unknown

ORF_unknown_translations_unknown_strand_translated_minus$sequences<-NULL
ORF_unknown_translations_unknown_strand_translated_minus$sequences<-paste(seq_minus)

if(length(Nindx_p) > 0)
{
    seq_plus <- seq_plus[-Nindx_p]
	ORF_unknown_translations_unknown_strand_translated_plus <- ORF_unknown_translations_unknown_strand_translated_plus[-Nindx_p]
}

if(length(Nindx_m) > 0)
{
    seq_minus <- seq_minus[-Nindx_m]
	ORF_unknown_translations_unknown_strand_translated_minus <- ORF_unknown_translations_unknown_strand_translated_minus[-Nindx_m,]
}

save(seq_plus, file =  "ORF_unknown_translations_unknown_strand-seq_plus_known_lnc.rda")


save(seq_minus, file =  "ORF_unknown_translations_unknown_strand-seq_minus_known_lnc.rda")

save(ORF_unknown_translations_unknown_strand_translated_plus, file =  "ORF_unknown_translations_unknown_strand_translated_plus_known_lnc.rda")

save(ORF_unknown_translations_unknown_strand_translated_minus, file =  "ORF_unknown_translations_unknown_strand_translated_minus_known_lnc.rda")


#six frame translation, as we dont know the strand
peptides_r1 <- translate(seq_plus)
peptides_r2 <- translate(subseq(seq_plus, start=2))
peptides_r3 <- translate(subseq(seq_plus, start=3))
peptides_r4 <- translate(seq_minus)
peptides_r5 <- translate(subseq(seq_minus, start=2))
peptides_r6 <- translate(subseq(seq_minus, start=3))


peptides_r1_df <- cbind("1", paste(ORF_unknown_translations_unknown_strand_translated_plus[, 'names']), paste(ORF_unknown_translations_unknown_strand_translated_plus[, 'strand']), paste(as.data.frame(peptides_r1)[,1]), paste(paste(ORF_unknown_translations_unknown_strand_translated_plus[, 'names']), "_f53_1", sep=""))
peptides_r2_df <- cbind("2", paste(ORF_unknown_translations_unknown_strand_translated_plus[, 'names']), paste(ORF_unknown_translations_unknown_strand_translated_plus[, 'strand']), paste(as.data.frame(peptides_r2)[,1]), paste(paste(ORF_unknown_translations_unknown_strand_translated_plus[, 'names']), "_f53_2", sep=""))
peptides_r3_df <- cbind("3", paste(ORF_unknown_translations_unknown_strand_translated_plus[, 'names']), paste(ORF_unknown_translations_unknown_strand_translated_plus[, 'strand']), paste(as.data.frame(peptides_r3)[,1]), paste(paste(ORF_unknown_translations_unknown_strand_translated_plus[, 'names']), "_f53_3", sep=""))
peptides_r4_df <- cbind("1", paste(ORF_unknown_translations_unknown_strand_translated_minus[, 'names']), paste(ORF_unknown_translations_unknown_strand_translated_minus[, 'strand']), paste(as.data.frame(peptides_r4)[,1]), paste(paste(ORF_unknown_translations_unknown_strand_translated_minus[, 'names']), "_f35_1", sep=""))
peptides_r5_df <- cbind("2", paste(ORF_unknown_translations_unknown_strand_translated_minus[, 'names']), paste(ORF_unknown_translations_unknown_strand_translated_minus[, 'strand']), paste(as.data.frame(peptides_r5)[,1]), paste(paste(ORF_unknown_translations_unknown_strand_translated_minus[, 'names']), "_f35_2", sep=""))
peptides_r6_df <- cbind("3", paste(ORF_unknown_translations_unknown_strand_translated_minus[, 'names']), paste(ORF_unknown_translations_unknown_strand_translated_minus[, 'strand']), paste(as.data.frame(peptides_r6)[,1]), paste(paste(ORF_unknown_translations_unknown_strand_translated_minus[, 'names']), "_f35_3", sep=""))

all_pep <- rbind(peptides_r1_df, peptides_r2_df, peptides_r3_df, peptides_r4_df, peptides_r5_df, peptides_r6_df)

colnames(all_pep) <- c("3frame_id", "transcript_id", "strand", "peptideSeq", "translation_id")

all_pep[,4] <- unlist(lapply(strsplit(paste(all_pep[,4]), "\\*"), FUN = function(y) paste(unlist(lapply(1:length(y), function(x) if(regexpr("(M[A-Z]+)", y[x]) != "-1") substr(y[x], regexpr("(M[A-Z]+)", y[x]), nchar(y[x])))), collapse="*") ))


save(all_pep, file =  "ORF_unknown_translations_unknown_strand-all_pep_known_lnc.rda")


#remove empty sequences
all_pep <- all_pep[-(which(all_pep[,4] == "")),]
#215085

save(all_pep, file =  "ORF_unknown_translations_unknown_strand-all_pep-f_known_lnc.rda")


##COMBINE THE THREE FRAME TRANLASTIONS WITH THE SIX FRAME TRANSLATIONS

all_pep_6f <- all_pep
all_pep <- NULL

load("ORF_unknown_translations_known_strand-all_pep-f_known_lnc.rda")


all_pep_3f <- all_pep
all_pep <- NULL
all_pep <- rbind(all_pep_3f, all_pep_6f)

#861225

ORF_unknown_translations_translated <- merge(ORF_unknown_translations, all_pep, by.x = 'miTransc_transcript_id', by.y = 'transcript_id', all.x=TRUE)
#980754

#luego se sobreescribe, por eso guardo la data.table
library(data.table)
ORF_unknown_translations_translated_dt <- setDT(ORF_unknown_translations_translated, keep.rownames=TRUE)

save(ORF_unknown_translations_translated_dt, file =  "ORF_unknown_translations_translated_dt_known_lnc.rda")


library(plyr)
ORF_unknown_translations_translated_dt_sorted <- arrange(ORF_unknown_translations_translated_dt, translation_id)

ORF_unknown_translations_non_translated <- ORF_unknown_translations_translated_dt_sorted[which(is.na(ORF_unknown_translations_translated_dt_sorted$translation_id)),]
#37
ORF_unknown_translations_translated_dt_sorted <- ORF_unknown_translations_translated_dt_sorted[which(!(is.na(ORF_unknown_translations_translated_dt_sorted$translation_id))),]
#980717

ORF_unknown_translations_translated_dt_sorted<-unique(ORF_unknown_translations_translated_dt_sorted)

tmp1 <- strsplit(paste(ORF_unknown_translations_translated_dt_sorted$peptideSeq), "\\*")
#980717
ORF_unknown_translations_translated_dt_sorted$freq <- unlist(lapply(strsplit(paste(ORF_unknown_translations_translated_dt_sorted$peptideSeq), "\\*"), FUN = function(x) length(x)))

tmp2 <- data.table(translation_id = rep(ORF_unknown_translations_translated_dt_sorted$translation_id, sapply(tmp1, length)), peptideSeq = unlist(tmp1))
#28599622

ORF_unknown_translations_translated_dt_sorted$peptideSeq <- NULL
ORF_unknown_translations_translated_dt_sorted$'3frame_id' <- NULL
ORF_unknown_translations_translated_dt_sorted$miTransc_start<-NULL
ORF_unknown_translations_translated_dt_sorted$miTransc_end<-NULL
ORF_unknown_translations_translated_dt_sorted$miTransc_tstatus<-NULL
ORF_unknown_translations_translated_dt_sorted$miTtrans_tss_id<-NULL
ORF_unknown_translations_translated_dt_sorted$G19_transcript_id<-NULL
ORF_unknown_translations_translated_dt_sorted$G19_oId<-NULL
tmp1 <- NULL
ORF_unknown_translations_translated_dt_sorted$rn<-NULL


tmp2<-tmp2[which(nchar(tmp2$peptideSeq) >="8"),]
tmp2<-unique(tmp2)
#2533350
#120981
setkey(ORF_unknown_translations_translated_dt_sorted, miTransc_transcript_id)
ORF_unknown_translations_translated_dt_sorted<-unique(ORF_unknown_translations_translated_dt_sorted)
#101227
ORF_unknown_translations_translated_dt_sorted_mod <- merge(ORF_unknown_translations_translated_dt_sorted, tmp2, by="translation_id", allow.cartesian=TRUE)


save(ORF_unknown_translations_translated_dt_sorted_mod, file = "ORF_unknown_translations_translated_dt_sorted_mod_known_lnc.rda")



save(ORF_unknown_translations_translated_dt_sorted, file = "ORF_unknown_translations_translated_dt_sorted_known_lnc.rda")
rm(ORF_unknown_translations_translated_dt_sorted)

save(tmp2,  file = "ORF_unknown_translations_translated_dt_sorted-tmp2_known_lnc.rda")
rm(tmp2)

ORF_unknown_translations_non_translated <- ORF_unknown_translations_translated[!(ORF_unknown_translations_translated$miTransc_transcript_id %in% ORF_unknown_translations_translated_dt_sorted_mod$miTransc_transcript_id),]


save(ORF_unknown_translations_non_translated, file="ORF_unknown_translations_non_translated_known_cod_features_not_in_GENCODE.rda")


ORF_unknown_translations_translated <- ORF_unknown_translations_translated_dt_sorted_mod[!(ORF_unknown_translations_translated_dt_sorted_mod$miTransc_transcript_id %in% ORF_unknown_translations_non_translated$miTrans_transcript_id),]

rm(ORF_unknown_translations_translated_dt_sorted_mod)
rm(ORF_unknown_translations_non_translated)

ORF_unknown_translations_translated<-ORF_unknown_translations_translated[which(!is.na(ORF_unknown_translations_translated$peptideSeq)),]

save(ORF_unknown_translations_translated, file="ORF_unknown_translations_translated_known_lnc.rda")


##with more cluster memory borrar objetos
rm(list(ls))

ORF_unknown_translations_translated$ind <- c(1:(nrow(ORF_unknown_translations_translated)))


ORF_unknown_translations_translated$translation_id2 <- paste(ORF_unknown_translations_translated$translation_id, "_", ORF_unknown_translations_translated$ind, sep = "")

#save(ORF_unknown_translations_translated, file="ORF_unknown_translations_translated_1.rda")
save(ORF_unknown_translations_translated, file="ORF_unknown_translations_translated_1_known_lnc.rda")
#175751

## CREATE THE DATABASE
library(plyr)
#remove the sequences with less than 8

ORF_unknown_translations_translated_f <- ORF_unknown_translations_translated[which(nchar(ORF_unknown_translations_translated$peptideSeq) >= 8),]


aa<-ORF_unknown_translations_translated_f

load("ORF_unknown_translationsKNOWNCOD_translated_f_NOT_GENCODE.rda")   #los que vienen de KNOWN CODING


ORF_unknown_translations_translated_unique_all <- rbind(aa[,c(10,12)],ORF_unknown_translations_translated_f_NOT_GENCODE[,c(11,13)])


ORF_unknown_translations_translated_unique_df <- ddply(ORF_unknown_translations_translated_unique_all, .(peptideSeq), summarise, paste(translation_id2, collapse = "|"))
#195028
colnames(ORF_unknown_translations_translated_unique_df) <- c("seq", "name")



ORF_unknown_translations_translated_unique_df$dummyNamesTarget <- paste("T_", row.names(ORF_unknown_translations_translated_unique_df), sep = "")


ORF_unknown_translations_translated_unique_df$dummyNamesDecoy <- paste("D_", row.names(ORF_unknown_translations_translated_unique_df), sep = "")


ORF_unknown_translations_dummyCodeTable <- data.frame(dummyCodes = c(ORF_unknown_translations_translated_unique_df$dummyNamesTarget, ORF_unknown_translations_translated_unique_df$dummyNamesDecoy), name = c(ORF_unknown_translations_translated_unique_df$name, ORF_unknown_translations_translated_unique_df$name))


DB_unknown_translations_DummyCodeFileName = "ORF_unknown_translations_dummyCodeTable_known_cod_features_ALL.rda"

save(ORF_unknown_translations_dummyCodeTable,file = DB_unknown_translations_DummyCodeFileName)


source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesShotgun.R")

ORF_unknown_translations_db <- AAStringSet(ORF_unknown_translations_translated_unique_df$seq)
DB_unknown_TargetFileName = "ORF_unknown_translations_known_lnc_features_ALL.fasta"


names(ORF_unknown_translations_db) <- ORF_unknown_translations_translated_unique_df$dummyNamesTarget
writeXStringSet(ORF_unknown_translations_db,DB_unknown_TargetFileName)


names(ORF_unknown_translations_db) <- ORF_unknown_translations_translated_unique_df$name
writeXStringSet(ORF_unknown_translations_db,paste(substr(DB_unknown_TargetFileName, 1, nchar(DB_unknown_TargetFileName) - 6), "_old.fasta", sep=""))

####hacer PROTEOGEST

perl /opt/proteogest/proteogest.pl -i ORF_unknown_translations_known_lnc_features_ALL.fasta -d -a -c trypsin -R M -W 15.99 -G 1 

#####################

ORF_unknown_translations_dig_930FileName<-"ORF_unknown_translations_known_lnc_features_ALL_simple.txt"
#ORF_unknown_translations_dig_930FileName<-"ORF_unknown_translations_known_cod_features_IN_GENCODE_simple.txt"
con <- file(ORF_unknown_translations_dig_930FileName, open="r")
ORF_unknown_translations_dig_930_raw <- readLines(con)
close(con)

ORF_unknown_translations_dig_930_protsAll <- ORF_unknown_translations_dig_930_raw[grep("Protein Name", ORF_unknown_translations_dig_930_raw)]
ORF_unknown_translations_dig_930_protsAll <- t(data.frame(strsplit(unlist(lapply(strsplit(ORF_unknown_translations_dig_930_protsAll, ": "), FUN = function(x) x[2])), "\\|")))
colnames(ORF_unknown_translations_dig_930_protsAll) <- NULL
rownames(ORF_unknown_translations_dig_930_protsAll) <- NULL

library(stringr)
ORF_unknown_translations_dig_930_npeptidesAll <- ORF_unknown_translations_dig_930_raw[grep("Number", ORF_unknown_translations_dig_930_raw)]
#ORF_unknown_translations_dig_930_NpeptidesAll <- as.numeric(paste(unlist(lapply(strsplit(ORF_unknown_translations_dig_930_npeptidesAll, " = "), FUN = function(x) x[2]))))
ORF_unknown_translations_dig_930_NpeptidesAll <- as.integer(str_extract(ORF_unknown_translations_dig_930_npeptidesAll, '[0-9]+$'))

tmp <- ORF_unknown_translations_dig_930_raw[-grep("Protein Name", ORF_unknown_translations_dig_930_raw)]
tmp <- tmp[tmp != ""]
tmp <- tmp[-grep("Number", tmp)]
tmp <- tmp[-grep("No", tmp)]
tmp <- tmp[-c(1,2)]
tmp2 <- strsplit(tmp, " ")

tmp2 <- lapply(tmp2, FUN = function(x) x[x != ""])

ORF_unknown_translations_dig_930_peptidesAll.df <- data.frame("PepNo" = unlist(lapply(tmp2, FUN = function(x) x[1])), "Range" = unlist(lapply(tmp2, FUN = function(x) x[2])), "IsotopicMass" = unlist(lapply(tmp2, FUN = function(x) x[3])), "AverageMass" = unlist(lapply(tmp2, FUN = function(x) x[4])), "Peptide" = unlist(lapply(tmp2, FUN = function(x) x[5])))

save(ORF_unknown_translations_dig_930_peptidesAll.df, file = "ORF_unknown_translations_dig_930_peptidesAll.df_known_lnc_features_ALL.rda")
#1823425
save(ORF_unknown_translations_dig_930_NpeptidesAll, file = "ORF_unknown_translations_dig_930_NpeptidesAll_known_lnc_features_ALL.rda")
#165813

npepind <- c()
for (i in 1:length(ORF_unknown_translations_dig_930_NpeptidesAll)) {
 	npepind <- c(npepind, rep(i, ORF_unknown_translations_dig_930_NpeptidesAll[i]))
}


save(npepind, file = "ORF_known_lnc_features_npepind_ALL.rda")

ORF_unknown_translations_dig_930_peptidesXProtAll.df <- data.frame(ORF_unknown_translations_dig_930_protsAll[npepind,], ORF_unknown_translations_dig_930_peptidesAll.df)
colnames(ORF_unknown_translations_dig_930_peptidesXProtAll.df)[1] <- c("transcript_ID")
ORF_unknown_translations_dig_930_peptidesXProtAll.df$Range <- paste("'", ORF_unknown_translations_dig_930_peptidesXProtAll.df$Range, sep = "")
ORF_unknown_translations_dig_930_peptidesXProtAll.df$Transcript_id_dig <- paste0(ORF_unknown_translations_dig_930_peptidesXProtAll.df$transcript_ID, "-",ORF_unknown_translations_dig_930_peptidesXProtAll.df$PepNo)
write.table(ORF_unknown_translations_dig_930_peptidesXProtAll.df, file="ORF_unknown_translations_dig_930_simple_digested_known_lnc_features.txt", col.names=TRUE, sep="\t", quote=FALSE)

save(ORF_unknown_translations_dig_930_peptidesXProtAll.df, file = "ORF_unknown_translations_dig_930_peptidesXProtAll_known_lnc_features.rda")



n_AA = 8
max_n_AA = 31
ORF_unknown_translations_dig_930_peptides_filterAA <- ORF_unknown_translations_dig_930_peptidesXProtAll.df[nchar(paste(ORF_unknown_translations_dig_930_peptidesXProtAll.df$Peptide)) > n_AA,]
#31274

ORF_unknown_translations_dig_930_peptides_filterAA <- ORF_unknown_translations_dig_930_peptides_filterAA[nchar(paste(ORF_unknown_translations_dig_930_peptides_filterAA$Peptide)) < max_n_AA,]
#28182


save(ORF_unknown_translations_dig_930_peptides_filterAA, file = "ORF_unknown_translations_dig_930_peptides_filterAA_known_lnc_features.rda")


ORF_unknown_translations_dig_930_peptides_filterAA_unique <- ORF_unknown_translations_dig_930_peptides_filterAA[,c(6,7)]
colnames(ORF_unknown_translations_dig_930_peptides_filterAA_unique)[1:2]<-c("peptideSeq","translation_id2")

#collapse peptideSeq values and concatenate translation_id values
library(plyr)
ORF_unknown_translations_translated_unique_df <- ddply(ORF_unknown_translations_dig_930_peptides_filterAA_unique, .(peptideSeq), summarise, paste(translation_id2, collapse = "|"))
#26358
####
#####create FASTA Db

colnames(ORF_unknown_translations_translated_unique_df) <- c("seq", "name")
ORF_unknown_translations_translated_unique_df$dummyNamesTarget <- paste("T_", row.names(ORF_unknown_translations_translated_unique_df), sep = "")
ORF_unknown_translations_translated_unique_df$dummyNamesDecoy <- paste("D_", row.names(ORF_unknown_translations_translated_unique_df), sep = "")

ORF_unknown_translations_dummyCodeTable <- data.frame(dummyCodes = c(ORF_unknown_translations_translated_unique_df$dummyNamesTarget, ORF_unknown_translations_translated_unique_df$dummyNamesDecoy), name = c(ORF_unknown_translations_translated_unique_df$name, ORF_unknown_translations_translated_unique_df$name))

#DB_unknown_translations_DummyCodeFileName = "ORF_unknown_translations_dummyCodeTable.rda"
DB_unknown_translations_DummyCodeFileName = "ORF_unknown_translations_dummyCodeTable_known_lnc_features_digested.rda"

save(ORF_unknown_translations_dummyCodeTable,file = DB_unknown_translations_DummyCodeFileName)

source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesShotgun.R")
ORF_unknown_translations_db <- AAStringSet(ORF_unknown_translations_translated_unique_df$seq)
#DB_unknown_TargetFileName = "ORF_unknown_translations_GINGERAS.fasta"
DB_unknown_TargetFileName = "ORF_NONCODING_DB.fasta"

names(ORF_unknown_translations_db) <- ORF_unknown_translations_translated_unique_df$dummyNamesTarget
writeXStringSet(ORF_unknown_translations_db,DB_unknown_TargetFileName)
names(ORF_unknown_translations_db) <- ORF_unknown_translations_translated_unique_df$name
writeXStringSet(ORF_unknown_translations_db,paste(substr(DB_unknown_TargetFileName, 1, nchar(DB_unknown_TargetFileName) - 6), "_old.fasta", sep=""))

db_target<-ORF_unknown_translations_db
db_decoy <- db_target
names(db_decoy) <- gsub("T", "D", names(db_decoy))
seq <- paste(as.character(db_decoy))

seq_decoy <- sapply(seq,  FUN=function(x) .strPseudoreverseDecoy(x))
seq_decoy_fa <- AAStringSet(seq_decoy)
#DB_unknown_DecoyFileName = "ORF_unknown_translations_GINGERAS_D.fasta"
DB_unknown_DecoyFileName = "ORF_NONCODING_D.fasta"


names(seq_decoy_fa) <- ORF_unknown_translations_translated_unique_df$dummyNamesDecoy
writeXStringSet(seq_decoy_fa, DB_unknown_DecoyFileName)
names(seq_decoy_fa) <- ORF_unknown_translations_translated_unique_df$name
writeXStringSet(seq_decoy_fa,paste(substr(DB_unknown_DecoyFileName, 1, nchar(DB_unknown_DecoyFileName) - 6), "_old.fasta", sep=""))
