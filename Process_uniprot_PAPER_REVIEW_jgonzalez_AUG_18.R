#descargo uniprot_swissprot_2017_11. HAY que filtrar solo HUMAN
library(GenomicRanges)
library(Biostrings)

DBTargetFileName<-"/home/nostromo/data/pepe/uniprot_2017_11/uniprot_sprot.fasta"
uniprot_FASTA <- readAAStringSet(DBTargetFileName, format="fasta",nrec=-1L, skip=0L, use.names=TRUE)

names_uniprot<-names(uniprot_FASTA)
sequence = paste(uniprot_FASTA)

uniprot_df_all<-data.frame("Protein_ID"=names_uniprot,"sequence"=sequence)


uniprot_df_all$organism<-str_extract(uniprot_df_all$Protein_ID,"(?<=_)[A-Z]+")
uniprot_df_human3<-uniprot_df_all[(uniprot_df_all$organism=="HUMAN"),]

uniprot_human_db <- AAStringSet(uniprot_df_human$sequence)
names(uniprot_human_db) <- uniprot_df_human$Protein_ID
#DB_unknown_TargetFileName = "ORF_unknown_translations_GINGERAS.fasta"
DB_unknown_TargetFileName = "uniprot_swissprot_2017_11_human.fasta"

writeXStringSet(uniprot_human_db,DB_unknown_TargetFileName)


#### PROTEOGEST
perl /opt/proteogest/proteogest.pl -i uniprot_swissprot_2017_11_human.fasta -d -a -c trypsin -R M -W 15.99 -G 1

####### after proteogest

ORF_unknown_translations_dig_930FileName <- "uniprot_swissprot_2017_11_human_simple.txt"

con <- file(ORF_unknown_translations_dig_930FileName, open="r")
ORF_unknown_translations_dig_930_raw <- readLines(con)
close(con)

ORF_unknown_translations_dig_930_protsAll <- ORF_unknown_translations_dig_930_raw[grep("Protein Name", ORF_unknown_translations_dig_930_raw)]
ORF_unknown_translations_dig_930_protsAll <- t(data.frame(strsplit(unlist(lapply(strsplit(ORF_unknown_translations_dig_930_protsAll, ": "), FUN = function(x) x[2])), "\\|")))

colnames(ORF_unknown_translations_dig_930_protsAll) <- NULL
rownames(ORF_unknown_translations_dig_930_protsAll) <- NULL

ORF_unknown_translations_dig_930_npeptidesAll <- ORF_unknown_translations_dig_930_raw[grep("Number", ORF_unknown_translations_dig_930_raw)]
ORF_unknown_translations_dig_930_NpeptidesAll <- as.numeric(paste(unlist(lapply(strsplit(ORF_unknown_translations_dig_930_npeptidesAll, " = "), FUN = function(x) x[2]))))

tmp <- ORF_unknown_translations_dig_930_raw[-grep("Protein Name", ORF_unknown_translations_dig_930_raw)]
tmp <- tmp[tmp != ""]
tmp <- tmp[-grep("Number", tmp)]
tmp <- tmp[-grep("No", tmp)]
tmp <- tmp[-c(1,2)]
tmp2 <- strsplit(tmp, " ")

tmp2 <- lapply(tmp2, FUN = function(x) x[x != ""])

ORF_unknown_translations_dig_930_peptidesAll.df <- data.frame("PepNo" = unlist(lapply(tmp2, FUN = function(x) x[1])), "Range" = unlist(lapply(tmp2, FUN = function(x) x[2])), "IsotopicMass" = unlist(lapply(tmp2, FUN = function(x) x[3])), "AverageMass" = unlist(lapply(tmp2, FUN = function(x) x[4])), "Peptide" = unlist(lapply(tmp2, FUN = function(x) x[5])))


npepind <- c()
for (i in 1:length(ORF_unknown_translations_dig_930_NpeptidesAll)) {
 	npepind <- c(npepind, rep(i, ORF_unknown_translations_dig_930_NpeptidesAll[i]))
}

ORF_unknown_translations_dig_930_peptidesXProtAll.df <- data.frame(ORF_unknown_translations_dig_930_protsAll[npepind,], ORF_unknown_translations_dig_930_peptidesAll.df)

colnames(ORF_unknown_translations_dig_930_peptidesXProtAll.df)[1:3] <- c("Protein_ID","Name", "Description")

ORF_unknown_translations_dig_930_peptidesXProtAll.df$Protein_ID<-ORF_unknown_translations_dig_930_peptidesXProtAll.df$Name

library(stringr)
ORF_unknown_translations_dig_930_peptidesXProtAll.df$Name <-str_extract(paste(ORF_unknown_translations_dig_930_peptidesXProtAll.df$Description), '^[A-Z0-9_]+(?= )')

write.table(ORF_unknown_translations_dig_930_peptidesXProtAll.df, file="uniprot_swissprot_dig_930_simple_digested.txt", col.names=TRUE, sep="\t", quote=FALSE)

save(ORF_unknown_translations_dig_930_peptidesXProtAll.df, file = "uniprot_swissprot_dig_930_peptidesXProtAll.df.rda")


ORF_unknown_translations_dig_930_peptidesXProtAll.df<-unique(ORF_unknown_translations_dig_930_peptidesXProtAll.df)
#2089716

n_AA = 8
max_n_AA = 50
ORF_unknown_translations_dig_930_peptides_filterAA <- ORF_unknown_translations_dig_930_peptidesXProtAll.df[nchar(paste(ORF_unknown_translations_dig_930_peptidesXProtAll.df$Peptide)) > n_AA,]

 ORF_unknown_translations_dig_930_peptides_filterAA <- ORF_unknown_translations_dig_930_peptides_filterAA[nchar(paste(ORF_unknown_translations_dig_930_peptides_filterAA$Peptide)) <= max_n_AA,]
# 2210285(filtrados)

fileroot<-"uniprot_swissprot"
save(ORF_unknown_translations_dig_930_peptides_filterAA, file = paste0(fileroot, "_peptidesXProtAll_filterAA_9-50.rda"))

db_peptidesXProtAll_filterAA_unique <- unique(ORF_unknown_translations_dig_930_peptides_filterAA[,c("Protein_ID","Peptide")])

write.table(db_peptidesXProtAll_filterAA_unique$Peptide, file = paste0(fileroot, "_peptidesXProtAll_filterAA_unique_peptides_9-50.txt"), col.names=FALSE, row.names = FALSE, sep="\t", quote=FALSE)
#1264099
####ALBA FUERA DE R

inputFile <- paste0(fileroot, "_peptidesXProtAll_filterAA_unique_peptides_9-50.txt")
outputFile <- paste0(fileroot, "_peptidesXProtAll_filterAA_unique_peptides_unique_9-50.txt")
y <- paste('sort -k1,1 ', inputFile,' | uniq -u > ', outputFile, sep = "")
system(y)
y <- NULL



db_peptidesXProtAll_filterAA_unique_peptides_unique <- read.csv2(paste0(fileroot, "_peptidesXProtAll_filterAA_unique_peptides_unique_9-50.txt"), header = FALSE, sep = "\t", fill = TRUE)
#1182517 (fuera de R)

db_peptidesXProtAll_filterAA_disc <- merge(db_peptidesXProtAll_filterAA_unique_peptides_unique, ORF_unknown_translations_dig_930_peptides_filterAA, by.x = "V1", by.y = "Peptide", all.x = TRUE)
colnames(db_peptidesXProtAll_filterAA_disc)[1] <- "Peptide"  ####este tiene los peptidos totales que son unicos
#2075756
save(db_peptidesXProtAll_filterAA_disc, file = paste0(fileroot, "_peptidesXProtAll_filterAA_disc_9-50.rda"))

db_peptidesXProtAll_filterAA_disc$nextProt_ID<-paste0("NX_",db_peptidesXProtAll_filterAA_disc$Protein_ID)
db_peptidesXProtAll_filterAA_disc<-unique(db_peptidesXProtAll_filterAA_disc)
load("/home/nostromo/data/pepe/02_neXtprot_20180117_Feb18/nextProt.RData")
#nextProtXChrXENSPXENSG

tmp_nextProt<-unique(nextProtXChr[,c(2,5,6)])
unique_nexProt<-unique(nextProtXChr)

#db_peptidesXProtAll_filterAA_disc_annot<-merge(db_peptidesXProtAll_filterAA_disc,tmp_nextProt,by.x="nextProt_ID",by.y="NextprotID", all.x=TRUE)
db_peptidesXProtAll_filterAA_disc_annot<-merge(db_peptidesXProtAll_filterAA_disc,unique_nexProt,by.x="nextProt_ID",by.y="NextprotID", all.x=TRUE)

#1182517  Se repiten lineas por algunos gene_names

db_peptidesXProtAll_filterAA_disc_annot_Unannotated_NX<-db_peptidesXProtAll_filterAA_disc_annot[which(paste(db_peptidesXProtAll_filterAA_disc_annot$Missing)=="NA"),]


db_peptidesXProtAll_filterAA_disc_annot_missing<-db_peptidesXProtAll_filterAA_disc_annot[(db_peptidesXProtAll_filterAA_disc_annot$Missing==1),]


#numero de proteinas de Nextprot que no tienen peptidos unicos
nextProtID_names<-unique(paste(tmp_nextProt$NextprotID))
triptics_names<-unique(paste((db_peptidesXProtAll_filterAA_disc$nextProt_ID)))
proteins_without_triptics_peptides<-setdiff(nextProtID_names, triptics_names) # = 229
proteins_without_triptics_peptides<-as.data.frame(proteins_without_triptics_peptides)
names(proteins_without_triptics_peptides)<-"NextprotID"
proteins_without_triptics_peptides_annot<-merge(proteins_without_triptics_peptides,nextProtXChr, by="NextprotID")

#comprobar names fasta uniprot y nextprot
library(GenomicRanges)
library(Biostrings)

DBTargetFileName<-"/home/nostromo/data/pepe/uniprot_2017_11/uniprot_swissprot_2017_11_human.fasta"
uniprot_FASTA <- readAAStringSet(DBTargetFileName, format="fasta",nrec=-1L, skip=0L, use.names=TRUE)

names_uniprot<-names(uniprot_FASTA)

names_uniprot_cleaned<-str_extract( names_uniprot,"(?<=\\|)[A-Z0-9]+(?=\\|)")
names_uniprot_NX<-paste0("NX_",names_uniprot_cleaned)

###check data
proteins_uniprot_not_in_nextprot<-setdiff(names_uniprot_NX, nextProtID_names) #13
proteins_nextprot_not_in_uniprot<-setdiff(nextProtID_names,names_uniprot_NX)  #0
#one hit wonders totales
tmp<-data.frame(Peptide= db_peptidesXProtAll_filterAA_disc_annot$Peptide,Nextprot_ID= db_peptidesXProtAll_filterAA_disc_annot$nextProt_ID)
library(plyr)
one_hit_wonders_df<-aggregate(Peptide ~ Nextprot_ID, data=tmp, FUN=paste0)
one_hit_wonders_df$number_of_unique_peptides<-str_count(paste(one_hit_wonders_df$Peptide),",")+1

one_hit_wonders<-one_hit_wonders_df[(one_hit_wonders_df$number_of_unique_peptides)==1,]
#57

####values by PE level one hit one_hit_wonders_df

one_hit_wonders_annot<-merge(one_hit_wonders,tmp_nextProt,by.x="Nextprot_ID", by.y="NextprotID", all.x=T)

table(one_hit_wonders_annot["PE"])
