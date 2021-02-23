args=(commandArgs(TRUE))

targetFileName <- args[1]
decoyFileName  <- args[2]
# targetFileName <- '/home/nostromo/data/pepe/08_PERCOLATOR_JUN18/Omssa_files/PXD001737/qe1_2014feb18_19_uc_at12_1_1rev.pin'
# decoyFileName  <- '/home/nostromo/data/pepe/08_PERCOLATOR_JUN18/Omssa_files/PXD001737-D//qe1_2014feb18_19_uc_at12_1_1rev.pin'




print(paste0('Target file: ' , targetFileName))
print(paste0('Decoy file: '  , decoyFileName))
#reading files
# omssatarget<-read.table(targetFileName,fill = T, header = F, stringsAsFactors = F, na.strings ="")
# omssadecoy<-read.table(decoyFileName, fill = T, header = F, stringsAsFactors = F, na.strings ="")
omssatarget <- read.table(targetFileName,fill = T, header = F, col.names = paste0('V', seq(1:max(count.fields(targetFileName)))), stringsAsFactors = F, na.strings ="")
omssadecoy  <- read.table(decoyFileName, fill = T, header = F, col.names = paste0('V', seq(1:max(count.fields(decoyFileName)))),stringsAsFactors = F, na.strings ="")



omssatarget<-omssatarget[which(!duplicated(omssatarget)),]
omssadecoy<-omssadecoy[which(!duplicated(omssadecoy)),]
#mycolumns<-paste("V",1:17,sep="")

omssatarget<-omssatarget[,1:17]
omssadecoy<-omssadecoy[,1:17]

header <- as.character(unlist(omssatarget[1,]))


omssatarget = omssatarget[-1,]
omssadecoy = omssadecoy[-1,]

names(omssatarget)<-header
names(omssadecoy)<-header

#omssatarget<-omssatarget[which(complete.cases(omssatarget[,c(1:16)])),]
#omssadecoy<-omssadecoy[which(complete.cases(omssadecoy[,c(1:16)])),]

omssadecoy$SpecId<-sub('target', 'decoy',omssadecoy$SpecId)
#change label column to decoyFileName
omssadecoy$Label="-1"
omssadecoy$Proteins<-paste0("DECOY_",omssadecoy$Proteins)
#merging Files

output_to_percolator<-rbind(omssatarget,omssadecoy)

#create outputdir

directory<-dirname(targetFileName)
namefile<-basename(targetFileName)

dir.create(paste0(directory,"/","percolator/"), showWarnings = FALSE)

write.table(output_to_percolator, file=paste0(directory,"/","percolator/",namefile), col.names=TRUE, row.names=FALSE,quote=FALSE, sep="\t")
