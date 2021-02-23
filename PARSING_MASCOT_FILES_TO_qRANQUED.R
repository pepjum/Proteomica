args=(commandArgs(TRUE))

targetFILE <- args[1]
#targetFILE<-"/home/nostromo/data/pepe/08_PERCOLATOR_JUN18/Dat_files/PXD001737/F293351.dat"

#source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesShotgun.R")


#fileMascot<-parseMascotDATFileComplete(targetFILE)
#fileMascot<-parseMascotDATFileSimple(targetFILE)

#targetFILE<-"/home/nostromo/data/pepe/08_PERCOLATOR_JUN18/Dat_files/PXD001737/F293351.dat"
nameFile<-basename(targetFILE)
dirFile<-dirname(targetFILE)

parseMascotDATFileSimple_G<-function(fileName) {

  cat("Reading file ", fileName, " ... \n");

  conn=file(fileName,open="r")
  linn=readLines(conn)

  indQuery <- grep("^q[0-9]", linn)

  queries <- linn[indQuery]

  cat("Processing queries ... \n");

  queries.lst <- strsplit(queries, "=")

  queries.lst_1 <- lapply(queries.lst, FUN = function(x) x[1])
  queries.lst_2 <- lapply(queries.lst, FUN = function(x) x[2])

  queries.lst_1.df <- strsplit(unlist(queries.lst_1), "_")

  queries.lst_1.df2 <- data.frame("query" = unlist(lapply(queries.lst_1.df, FUN = function(x) x[1])), "rank" = unlist(lapply(queries.lst_1.df, FUN = function(x) x[2])), "ob1" = unlist(lapply(queries.lst_1.df, FUN = function(x) x[3])), "ob2" = unlist(lapply(queries.lst_1.df, FUN = function(x) x[4])))

  queries.lst_2_str <- unlist(queries.lst_2)
  queries.lst_2_strlst <- strsplit(queries.lst_2_str, ";")

  queries.lst_2_df <- data.frame("assignPep" = unlist(lapply(queries.lst_2_strlst, FUN = function(x) x[1])),  "assignProt" = unlist(lapply(queries.lst_2_strlst, FUN = function(x) x[2])))

  queries.df <- data.frame(queries.lst_1.df2, queries.lst_2_df)
  queries.assign.df <- queries.df[!is.na(queries.df$assignProt),]


  queries.assign.df.p1 <- as.data.frame(strsplit(paste(queries.assign.df$assignPep), ","))
  queries.assign.df.p1 <- t(queries.assign.df.p1)
  queries.assign.df.p1 <- data.frame(queries.assign.df[, c(1,2)], queries.assign.df.p1)
  queries.assign.df.p2 <- strsplit(paste(queries.assign.df$assignProt), ",")

  queries.assign.prot.df <- unlist(strsplit(unlist(queries.assign.df.p2), ";"))
  pepProt <- unlist(lapply(queries.assign.df.p2, length))
  queries.assign.prot.df2 <- strsplit(queries.assign.prot.df, ":")
  queries.assign.prot.df3 <- data.frame(unlist(lapply(queries.assign.prot.df2, FUN = function(x) x[1])), unlist(lapply(queries.assign.prot.df2, FUN = function(x) x[2])), unlist(lapply(queries.assign.prot.df2, FUN = function(x) x[3])), unlist(lapply(queries.assign.prot.df2, FUN = function(x) x[4])), unlist(lapply(queries.assign.prot.df2, FUN = function(x) x[5])))
  colnames(queries.assign.prot.df3) <- NULL
  rownames(queries.assign.prot.df3) <- NULL

  indexPep <- list()

  for (i in 1:length(pepProt)) {
	  indexPep[i] <- list(rep(i, pepProt[i]))
  }

  indexPep.vec <- unlist(indexPep)

  queriesPepProt <- data.frame(queries.assign.df.p1[indexPep.vec,], queries.assign.prot.df3)

  colnames(queriesPepProt) <- c("Query", "Rank", "MissedCleavages", "PeptideMr", "Delta", "NumberOfIonsMatched", "PeptideSeq", "PeaksFromIons1", "VarMod", "ion_score", "IonSeries", "PeaksFromIons2", "PeaksFromIons3", "ProteinAccession", "FrameNumber", "Start", "End", "Multiplicity")
  rownames(queriesPepProt) <- NULL
  queriesPepProt[, "ProteinAccession"] <- sub("\\\"", "", queriesPepProt[, "ProteinAccession"])
  queriesPepProt[, "ProteinAccession"] <- sub("\\\"", "", queriesPepProt[, "ProteinAccession"])

  # Merged to add the Enzime cuts
  enzimeCuts_df <- setNames(queries.df[,grep(paste(c("query", "rank", "ob1", "assignPep"), collapse="|"), colnames(queries.df))], c("Query", "Rank", "ob1", "enzCut"))
  enzimeCuts_df <- enzimeCuts_df[c(which(enzimeCuts_df$ob1 == "terms")),]
  enzimeCuts_df <- enzimeCuts_df[,grep(paste(c("Query", "Rank", "enzCut"), collapse="|"), colnames(enzimeCuts_df))]

  queriesPepProt <- merge(queriesPepProt, enzimeCuts_df, by=c("Query", "Rank"), all=T)

  return(queriesPepProt)

}

fileMascot<-parseMascotDATFileSimple_G(targetFILE)

terms<-paste(fileMascot$enzCut)

N_term<-substr(terms,1,1)
C_term<-substr(terms,3,3)

fileMascot$Nterm<-N_term
fileMascot$C_term<-C_term

###add enzyme cuts columns to output
fileMascot$enzN <- (fileMascot$Nterm =="R" | fileMascot$Nterm =="K")*1
fileMascot$enzC <- (fileMascot$C_term =="K"| fileMascot$C_term =="R")*1
fileMascot$enzInt<-fileMascot$MissedCleavages
fileMascot$CalcMass<-fileMascot$PeptideMr
fileMascot$dM <- fileMascot$Delta


tmp<-fileMascot[,c("CalcMass","dM")]

tmp$CalcMass <- as.numeric(paste(tmp$CalcMass))
tmp$dM <- as.numeric(paste(tmp$dM))
fileMascot$ExpMass<-apply(tmp, 1, sum)
fileMascot$absdM<-abs(tmp$dM)
fileMascot$PepLen<-nchar(as.character(fileMascot$PeptideSeq))
fileMascot$Label=1
fileMascot$PeptideSeq<-paste0(fileMascot$Nterm,".",fileMascot$PeptideSeq,".",fileMascot$C_term)

output_df<-fileMascot[,c("Query","Label","ExpMass","CalcMass","ion_score","dM","absdM","PepLen","enzN","enzC","enzInt","PeptideSeq","ProteinAccession")]
names(output_df)<-c("SpecId","Label","ExpMass","CalcMass","ion_score","dM","absdM","PepLen","enzN","enzC","enzInt","Peptide","Proteins")

output_df$ScanNr<-output_df$SpecId
output_df$SpecId<-paste0(targetFILE,"_",output_df$SpecId)
#dir.create(paste0(dirFile,"/percolator/",showWarnings=FALSE))

write.table(output_df, file=paste0(dirFile,"/",nameFile,".pin"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
