args=(commandArgs(TRUE))

targetFILE <- args[1]
path<-dirname(targetFILE)

percolator_file<-read.table(targetFILE, header=T,fill=T)

percolator_df<-percolator_file[c("SpecId","Label","ScanNr","ExpMass","CalcMass","hyperscore","deltaScore","frac_ion_b","frac_ion_y","Mass","dM","absdM","PepLen","Charge2","Charge3","Charge4","enzN","enzC","enzInt","Peptide"),]

write.table(percolator_df,file=paste0(path,"/","ALL_TANDEM_rev.pin"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
