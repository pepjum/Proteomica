args=(commandArgs(TRUE))

library(plyr)
targetDirectory <- args[1]

print(targetDirectory)

aa<-list.files(targetDirectory,pattern="*.pin")
#print(aa)

path<-targetDirectory
#print(path)
bb<-paste0(path,aa)
#print(bb)
counter <- 0
for (fl in bb){
    print(fl)
    tmp <- read.table(fl, fill = T, header = F, stringsAsFactors = F)
    counter <- counter + (nrow(tmp)-1)
    if((sum(paste(tmp[1,]) == '') != 0)){
        newNames <- c(paste(tmp[1,])[1:ncol(tmp)-1], paste0(paste(tmp[1,])[ncol(tmp)-1], seq(1:sum(paste(tmp[1,]) == ''))))
        colnames(tmp) <- newNames
        tmp <- tmp[-c(1),]
    }else{
        colnames(tmp) <- tmp[1,]
        tmp <- tmp[-c(1),]
    }
    if (fl==bb[1]){
        df<-tmp
    }else{
        df<-rbind.fill(df,tmp)
    }
}
ifelse(counter == nrow(df), 'Ok!', ':c')

# tables <- lapply(bb, read.table, header = TRUE, fill = T, row.names = NULL)
#
# combined.df <- do.call(rbind , tables)
df$lnNumSP<-NULL

write.table(df,file=paste0(path,"ALL_OMSSA.pin"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
