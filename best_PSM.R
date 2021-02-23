args=(commandArgs(TRUE))

library(stringr)

tandem_file<-args[1]
#

#batch_number<-str_extract(tandem_file, "[0-9]+")


tandem_file_load<-get(load(tandem_file))

tandem_best_psm<-data.frame()
for (query in unique(paste(tandem_file_load$Query))){
    #print(query)
    tmp<-tandem_file_load[which(tandem_file_load$Query==query),]
    if((nrow(tmp) >=2 )& length(unique(paste(tmp$score) >1))){
        tmp_selected<-tmp[which.max(as.numeric(paste(tmp$score))),]
        tandem_best_psm<-rbind(tandem_best_psm,tmp_selected)
    }else if(nrow(tmp)==1){
        tandem_best_psm<-rbind(tandem_best_psm,tmp)
    }

}

output_file<-paste0(lapply(strsplit(paste(tandem_file),"\\."),"[",1),"_bestPSM.Rdata")
save(tandem_best_psm, file=output_file)
