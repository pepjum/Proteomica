
database<-args[1]   #sptxt
spread_value<-args[2] #0
Dalton<-args[3]   #1

source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesShotgun.R")

#mgf_dir="."
options(stringsAsFactors = FALSE)

library(stringr)
library(stringi)



#spread_value<-0


readSpectraSTdb<-function(x){
    con <- file(db_spectra, open="r")
    db_raw <- readLines(con)
    close(con)

#db_dig_protsAll <- db_dig_raw[grep("Protein Name", db_dig_raw)]
    Names_spectra_ALL<-db_raw[grep("^Name",db_raw)]
    lib_id_spectra_ALL<-db_raw[grep("^LibID",db_raw)]
    PrecursorMZ_ALL<-db_raw[grep("^PrecursorMZ",db_raw)]
    Status_ALL<-db_raw[grep("^Status",db_raw)]
    Full_name_ALL<-db_raw[grep("^FullName",db_raw)]
    Comments_ALL<-db_raw[grep("^Comment",db_raw)]
    NumPeaks<-db_raw[grep("^NumPeaks",db_raw)]
    Names_spectra_ALL <- t(data.frame(strsplit(unlist(lapply(strsplit(Names_spectra_ALL, ": "), FUN = function(x) x[2])), "\\|")))

    colnames(Names_spectra_ALL) <- NULL
    rownames(Names_spectra_ALL) <- NULL
    colnames(Names_spectra_ALL)<-"Name"

    lib_id_spectra_ALL<-as.integer(str_extract(lib_id_spectra_ALL,'[0-9]+$'))
    PrecursorMZ_ALL<-as.numeric(str_extract(PrecursorMZ_ALL,'[0-9\\.]+'))
    Status_ALL<-sapply(strsplit(paste(Status_ALL)," "),'[',2)
    Full_name_ALL<-sapply(strsplit(paste(Full_name_ALL)," "),'[',2)
    NumPeaks<-as.numeric(str_extract(NumPeaks,'[0-9\\.]+'))


    tmp<-db_raw[-c(1:6)]  #elimino cabecera
    tmp<-lapply(tmp, FUN=function(x) x=paste0(x,"####")) #meto un fin de linea reconocible
    tmp<-lapply(tmp, FUN=function(x) if(grepl("^NumPeaks",x)){x=paste0(x,"@")}else{x=x})
    tmp<-stri_paste(tmp,collapse='') #uno todo en un string
    tmp<-strsplit(paste(tmp),"#Name")  #separo por name y asi separamos los espectros
    tmp <- unlist(lapply(tmp, FUN = function(x) x[x != ""]))
    tmp<-sapply(strsplit(paste(tmp),"@"),'[',3)

    SpectraST_df <- cbind(Names_spectra_ALL,data.frame("Lib_ID" = paste(lib_id_spectra_ALL), "PrecursorMZ" = paste(PrecursorMZ_ALL), "Status" = paste(Status_ALL), "FullName" = paste(Full_name_ALL), "NumPeaks" = paste(NumPeaks), "Comments"=paste(Comments_ALL), "Peaks"=paste(tmp)))

    return(SpectraST_df)
}


#####CALL TO DB
db_spectra<-readSpectraSTdb(database)

list_all_peaks<-paste(db_spectra$Peaks)

bining<-function(dataframe,Dalton){
    dataframe<-as.data.frame(dataframe)
	new = 0
	old = nrow(dataframe)

	dataframe$M<-as.numeric(str_extract(paste(dataframe$M),'^[0-9]*'))
	while( old != new){

		tmp = data.frame('M' = NULL, 'Z' = NULL)
		SKIP_LINE = F
		iteration = nrow(dataframe)- 1
		for (i in 1:(iteration)){
			if(SKIP_LINE){
				SKIP_LINE = F
				next
			}
			dataframe<-dataframe[order(as.numeric(paste(dataframe$M))),]
			if((dataframe[i+1, 'M'] - dataframe[i, 'M'] ) <= Dalton){
				#print('shrinking')
				maxM = max(dataframe[i+1, 'M'], dataframe[i, 'M'] )
				maxZ = max(dataframe[i+1, 'Z'], dataframe[i, 'Z'] )
				tmp = rbind(tmp, data.frame('M' = maxM, 'Z' = maxZ))
				SKIP_LINE = T
			}else{
				tmp = rbind(tmp, dataframe[i,])
				SKIP_LINE = F
			}
		}
		tmp <- rbind(tmp, dataframe[i+1,])
		new <- nrow(tmp)
		old <- nrow(dataframe)
		dataframe <- tmp
	}
	if((dataframe[nrow(dataframe), 'M'] - dataframe[nrow(dataframe)-1, 'M'] ) <= 1){
		maxM = max(dataframe[i, 'M'], dataframe[i-1 , 'M'] )
		maxZ = max(dataframe[i, 'Z'], dataframe[i-1, 'Z'] )
		dataframe2 <- dataframe[c(seq(1,nrow(dataframe)-2)), ]
		dataframe2 <- rbind(dataframe2, data.frame('M' = maxM, 'Z' = maxZ))
	}else{
        dataframe2 <- rbind(dataframe, data.frame('M' = maxM, 'Z' = maxZ))

    }
    dataframe2<-unique(dataframe2)
	return(dataframe2)

}

####transform peaks info in dataframes to work in it

list_df_peaks_not_binned<-list()
for(j in 1:length(list_all_peaks)){
    cat(j, "of",length(list_all_peaks),"\n")
    tmp<-strsplit(paste(list_all_peaks[j]), "####")
    tmp<-tmp[!grepl('^###',tmp)]
    tmp<-lapply(tmp, FUN= function(x) x[!grepl('^###',x)])
    m<-paste(lapply(strsplit(unlist(tmp),"\\t"),'[',1))
    z<-paste(lapply(strsplit(unlist(tmp),"\\t"),'[',2))
    df_peaks<-data.frame("M"=m,"Z"=z)
    print(j)
    print(dim(df_peaks))
    list_df_peaks_not_binned[[j]]<-df_peaks

}

list_df_peaks_binned<-list()
for (k in 1: length(list_df_peaks_not_binned)){
    cat(k,"\n")
    df_peaks_binned<-bining(list_df_peaks_not_binned[k],1)
    list_df_peaks_binned[[k]]<-df_peaks_binned
}

list_df_peaks_binned_and_spread_out<-list()
for (i in 1:length(list_df_peaks_binned)){
    cat(i, "of ", length(list_df_peaks_binned), "\n" )
    df<-as.data.frame(list_df_peaks_binned[i][[1]])
    df$percentage<-0

    df$percentage<-(as.numeric(paste(df$Z)))*as.numeric(paste(spread_value))/2    #el porcentaje del pico a transferir hay que dividirlo entre sus dos vecinos
    ###ahora quito el porcentaje de los picos originales
    df$peak_less<-as.numeric(paste(df$Z))-as.numeric(paste(df$percentage))    #el porcentaje del pico a transferir hay que dividirlo entre sus dos vecinos

    #sumar los pedazos de los adyacentes
    df$peaks_final<-0
    df$peaks_final[1]<-as.numeric(paste(df$peak_less[1]))+as.numeric(paste(df$percentage[2]))
    df$peaks_final[nrow(df)]<-as.numeric(paste(df$peak_less[nrow(df)]))+as.numeric(paste(df$percentage[nrow(df)-1]))

    for (j in 2:(nrow(df)-1)){
        df$peaks_final[j]<-as.numeric(paste(df$peak_less[j]))+as.numeric(paste(df$percentage[j-1]))+as.numeric(paste(df$percentage[j+1]))    #el porcentaje del pico a transferir hay que dividirlo entre sus dos vecinos
    }
    df_final<-df[,-c(2,3,4)]
    #df_final<-df
    print(i)
    print(dim(df_final))
    list_df_peaks_binned_and_spread_out[[i]]<-df_final
}

save(list_df_peaks_binned_and_spread_out, file="/home/margaret/data/pepe/HUMAN_NIST_25_7_2018/HUMAN_best.Rdata")

list_df_normalized<-list()
for (i in 1:length(list_df_peaks_binned_and_spread_out)){
    cat(i, "of ", length(list_df_peaks_binned_and_spread_out), "\n" )
    df<-as.data.frame(list_df_peaks_binned_and_spread_out[i][[1]])
    df$square<-(df$peaks_final)^2
    sqrt_sum_tot_peaks_squared<-sqrt(sum(as.numeric(df$square)))
    df$normalized<-df$peaks_final/sqrt_sum_tot_peaks_squared
    df<-df[,c(1,4)]
    list_df_normalized[[i]]<-df
}

save(list_df_normalized, file="/home/margaret/data/pepe/HUMAN_NIST_25_7_2018/HUMAN_best_spread_out_normalized.Rdata")


cat("end of processing database")
