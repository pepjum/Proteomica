args=(commandArgs(TRUE))


mgf<-args[1]      #mgf file
database<-args[2]   #database spectaST
range<-args[3]     #for selecting spectra by mz difference of .... (0.1 by default)
spread_value<-args[4]    #spread out
Dalton<-args[5]       # Value of daltons for bining function purposes

source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesShotgun.R")

#mgf_dir="."
options(stringsAsFactors = FALSE)

library(stringr)
library(stringi)
#setwd(paste0(mgf_dir,"/"))

cat("reading MGF file","\n")
mgf_file<-ReadMGFFile(mgf)

name_spectrums<-c()
mass_values<-c()
peaks_list<-c()

#READ SPECTRUMS AND KEEP NAMES AND M/Z CHARGE

for (i in 1:length(mgf_file)){
    name_file<-mgf_file[[i]][2]
    name_spectrums<-c(name_spectrums,name_file)
    mass<-mgf_file[[i]][4]
    mass_values<-c(mass_values,mass)
    peaks<-mgf_file[[i]][6:length(mgf_file[[i]])-1]
    peaks_cod<-paste(peaks,collapse="###")
    peaks_list<-c(peaks_list, peaks_cod)

}
mass_values<-sapply(strsplit(paste(mass_values),"="),`[`,2)
mass_values<-sapply(strsplit(paste(mass_values)," "),`[`,1)

name_files<-sapply(strsplit(paste(name_spectrums),"="),`[`,2)
name_spectrums<-sapply(strsplit(paste(name_files)," "),`[`,1)

#creating df from mgf file

df_spectrums_mgf<-data.frame("spectrum"=name_spectrums,"mz_value"=mass_values, "peaks"=paste(peaks_list))
rm(mgf_file)

cat("loading database","\n")
db_spectra<-get(load(database))   #_500-1500.Rdata

selectSpectraByMZrange<-function(SpectraST_df,PrecursorMGF,range, ind){

    threshold_up <- (as.numeric(PrecursorMGF)+ as.numeric(range))
    threshold_down <- (as.numeric(PrecursorMGF)- as.numeric(range))
    l_up <- which(as.numeric(SpectraST_df$PrecursorMZ) < threshold_up)
    l_down <- which(as.numeric(SpectraST_df$PrecursorMZ) > threshold_down)
    tmp <- SpectraST_df[intersect(l_up,l_down ),]
    if(nrow(tmp)!=0){
        tmp$idx <- ind
        return(tmp)
    }else{
        return(tmp)
    }
}
#range=0.5

cat("selecting spectrums for MZ range","\n")
df_spectrums_selected_for_each_spectrum<-data.frame()
for (i in 1:nrow(df_spectrums_mgf)){
    #cat(i,"\n")
    PrecursorMGF<-as.numeric(paste(df_spectrums_mgf$mz_value[i]))
    last <- selectSpectraByMZrange(db_spectra,PrecursorMGF,range, i)
    if((nrow(last) != 0 && nrow(df_spectrums_selected_for_each_spectrum) == 0)){
        #print('First?')
        df_spectrums_selected_for_each_spectrum<-last
    }else if (nrow(last) !=0){
        #print('Adding one!')
        df_spectrums_selected_for_each_spectrum <- rbind(df_spectrums_selected_for_each_spectrum, selectSpectraByMZrange(db_spectra,PrecursorMGF,range, i))
    }
}

#####binado de señal del mgf

list_all_peaks<-paste(df_spectrums_mgf$peaks)

cat("bining peaks","\n")
bining<-function(dataframe,Dalton){
    dataframe<-as.data.frame(dataframe)
	new = 0
	old = nrow(dataframe)

	dataframe$M<-as.numeric(str_extract(paste(dataframe$M),'^[0-9]*'))
    dataframe$Z<-as.numeric(paste(dataframe$Z))
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
			if(as.numeric(paste(dataframe[i+1, 'M'])) - as.numeric(paste(dataframe[i, 'M']))  <= as.numeric(Dalton)){
				#print('shrinking')
				maxM = max(as.numeric(paste(dataframe[i+1, 'M'])), as.numeric(paste(dataframe[i, 'M'] )))
				maxZ = max(as.numeric(paste(dataframe[i+1, 'Z'])), as.numeric(paste(dataframe[i, 'Z'] )))
				tmp = rbind(tmp, data.frame('M' = as.numeric(maxM), 'Z' = as.numeric(maxZ)))
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
	if(as.numeric(paste(dataframe[nrow(dataframe), 'M'])) - as.numeric(paste(dataframe[nrow(dataframe)-1, 'M'])) <= as.numeric(Dalton)){
		maxM = max(as.numeric(paste(dataframe[nrow(dataframe), 'M'])), as.numeric(paste(dataframe[nrow(dataframe)-1 , 'M'] )))
		maxZ = max(as.numeric(paste(dataframe[nrow(dataframe), 'Z'])), as.numeric(paste(dataframe[nrow(dataframe)-1, 'Z'] )))
		dataframe2 <- dataframe[c(seq(1,nrow(dataframe)-2)), ]
		dataframe2 <- rbind(dataframe2, data.frame('M' = as.numeric(maxM), 'Z' = as.numeric(maxZ)))
	}else{
        dataframe2 <- dataframe

    }
    dataframe2<-unique(dataframe2)
	return(dataframe2)

}
####transform peaks info in dataframes to work in it
cat("transform peaks info to a dataframe to work in it","\n")

list_df_peaks_not_binned<-list()
for(j in 1:length(list_all_peaks)){
    tmp<-strsplit(paste(list_all_peaks[j]), "###")
    m<-paste(lapply(strsplit(unlist(tmp)," "),'[',1))
    z<-paste(lapply(strsplit(unlist(tmp)," "),'[',2))
    df_peaks<-data.frame("M"=m,"Z"=z)
    #print(j)
    list_df_peaks_not_binned[[j]]<-df_peaks

}
##### resampling before bining. Despues hacer el binado con alguna funcion de algun paquete, a ver si es mas rapida

### bining peaks
cat("bining peaks","\n")
list_df_peaks_binned<-list()
for (k in 1: length(list_df_peaks_not_binned)){
    #cat(k,"\n")
    df_peaks_binned<-bining(list_df_peaks_not_binned[k],1)
    list_df_peaks_binned[[k]]<-df_peaks_binned
}

###### enabling matching (trasladar porcentaje de picos a sus vecinos ######
#spread_value<-0
cat("spread_out","\n")
list_df_peaks_binned_and_spread_out<-list()
for (i in 1:length(list_df_peaks_binned)){
    if (i %% 100 ==0){
        #cat(i, "of ", length(list_df_peaks_binned), "\n" )
    }
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
    #print(dim(df_final))
    list_df_peaks_binned_and_spread_out[[i]]<-df_final
}

cat("selecting range of peak list","\n")
list_selected<-list()
for (i in 1:length(list_df_peaks_binned_and_spread_out)){
    #cat(i, "of ", length(list_df_peaks_binned_and_spread_out), "\n" )
    df<-as.data.frame(list_df_peaks_binned_and_spread_out[i][[1]])
    df_ranges_tmp<-df[as.numeric(paste(df$M))>=500,]
    df_ranges_final<-df_ranges_tmp[as.numeric(paste(df_ranges_tmp$M))<=1500,]
    list_selected[[i]]<-df_ranges_final
}



cat("normalization","\n")

#####normalization
list_df_normalized<-list()
for (i in 1:length(list_selected)){
    df<-as.data.frame(list_selected[i][[1]])
    df$normalized<-as.numeric(df$peaks_final)/sqrt(sum(as.numeric(df$peaks_final)^2))
    df<-df[,c(1,3)]
    list_df_normalized[[i]]<-df
}

#creating a vector adding zeros to positions where we have not mz values
cat("creating a vector adding zeros to positions where we have not mz values","\n")
list_df_matrix<-list()
for (i in 1:length(list_df_normalized)){
    #cat(i, "of ", length(list_df_normalized), "\n" )
    df<-as.data.frame(list_df_normalized[i][[1]])
    matrix<-data.frame(M=as.numeric(seq(1:1500)),normalized=0)
    if (nrow(df)==0){
        df<- matrix
    }
    matrix[which(as.numeric(matrix$M) %in% as.numeric(df$M)),]<-df$normalized
    matrix$M<-seq(1:nrow(matrix))
    df_ranges_tmp<-matrix[as.numeric(paste(matrix$M))>=500,]
    df_ranges_final<-df_ranges_tmp[as.numeric(paste(df_ranges_tmp$M))<=1500,]

    list_df_matrix[[i]]<-df_ranges_final
}


#convert dataframe into a string to indexing in the info dataframe
cat("convert peaks dataframe into string to indexing in the info dataframe","\n")
list_linecoded<-c()
for (i in 1:length(list_df_matrix)){
    #cat(i,"/",length(list_df_matrix),"\n")
    df<-as.data.frame(list_df_matrix[i][[1]])
    df$line<-paste0(df$M,"\t",df$normalized,"###")
    linecod<-paste0(df$line)
    linecod<-paste(linecod,collapse="")
    list_linecoded<-c(list_linecoded, linecod)
}

df_spectrums_mgf$peaks<-NULL
df_spectrums_mgf$Peaks_mgf<-paste(list_linecoded)

#indexando identificaciones a los espectros

df_spectrums_mgf$idx<-seq(1:(nrow(df_spectrums_mgf)))
#names(df_spectrums_mgf)[3]<-"Peaks_mgf"
df_results_tmp<-merge(df_spectrums_mgf,df_spectrums_selected_for_each_spectrum, by="idx", all.x=TRUE )

#eliminar no hits
df_results_tmp2<-df_results_tmp[which(!is.na(df_results_tmp$Name)),]
rm(df_results_tmp)

#### calcular dot product y dot bias
cat("calculating dot product and dot bias","\n")
dot_product<-c()
Dot_bias<-c()
for (i in 1:nrow(df_results_tmp2)){
    #cat(i,"/",nrow(df_results_tmp2),"\n")
    tmp_mgf<-strsplit(paste(df_results_tmp2$Peaks_mgf[i]), "###")
    m_mgf<-paste(lapply(strsplit(unlist(tmp_mgf),"\t"),'[',1))
    z_mgf<-paste(lapply(strsplit(unlist(tmp_mgf),"\t"),'[',2))
    df_peaks_mgf<-data.frame("M"=as.numeric(m_mgf),"Z"=as.numeric(z_mgf))
    tmp_database<-strsplit(paste(df_results_tmp2$Peaks[i]), "###")
    m_database<-paste(lapply(strsplit(unlist(tmp_database),"\t"),'[',1))
    z_database<-paste(lapply(strsplit(unlist(tmp_database),"\t"),'[',2))
    df_peaks_database<-data.frame("M"=as.numeric(m_database),"Z"=as.numeric(z_database))
    df_peaks_mgf_tmp1<-df_peaks_mgf[(df_peaks_mgf$M >500),]
    df_peaks_mgf_tmp2<-df_peaks_mgf_tmp1[(df_peaks_mgf_tmp1$M <1500),]
    df_peaks_database_tmp1<-df_peaks_database[(df_peaks_database$M >500),]
    df_peaks_database_tmp2<-df_peaks_database_tmp1[(df_peaks_database_tmp1$M <1500),]
    names(df_peaks_database_tmp2)[2]<-"Z_database"
    df_tmp_res<-merge(df_peaks_mgf_tmp2,df_peaks_database_tmp2, by="M")
    df_tmp_res$dot<-(as.numeric(df_tmp_res$Z))*(as.numeric(df_tmp_res$Z_database))
    D<-sum(as.numeric(df_peaks_mgf_tmp2$Z)*(as.numeric(df_peaks_database_tmp2$Z_database)))
    #df_tmp_res$pot<-(as.numeric(df_tmp_res$Z)^2)*(as.numeric(df_tmp_res$Z_database)^2)
    DB<-sqrt(sum(as.numeric(as.numeric(df_tmp_res$Z)^2)*(as.numeric(df_tmp_res$Z_database)^2)))/as.numeric(D)
    dot_product<-c(dot_product,D)
    Dot_bias<-c(Dot_bias,DB)
}

df_results_tmp2$D<-as.numeric(dot_product)
df_results_tmp2$DB<-as.numeric(Dot_bias)
df_results_tmp2$Peaks_mgf<-NULL
df_results_tmp2$Peaks<- NULL

df_results_tmp2[with(df_results_tmp2, order(idx,-D)),]

####calculating Delta
cat("calculating DeltaD","\n")
df_results_tmp2$Delta<-NA
for (i in unique(df_results_tmp2$idx)){
    tmp<-df_results_tmp2[which(df_results_tmp2$idx == i),]
    tmp_max1<-tmp[(which((tmp$D) ==max(tmp$D))),]
    tmp_max2<-tmp[-c(which((tmp$D) ==max(tmp$D))),]
    tmp_max3<-tmp_max2[(which((tmp_max2$D) ==max(tmp_max2$D))),]
    D1<-as.numeric(tmp_max1$D)
    D2<-as.numeric(tmp_max3$D)
    D3<-(D1-D2)/D1
    df_results_tmp2[df_results_tmp2$idx == i, 'Delta']<-as.numeric(D3)
}

df_results_tmp2$Fscore<- NA
df_results_tmp2[which(is.nan(df_results_tmp2$DB)), "DB"]<-0
#df_results_tmp2[which(is.na(df_results_tmp2$Delta)), "Delta"]<-1

cat("calculating Fscore","\n")
for (i in 1:nrow(df_results_tmp2)){
    if (as.numeric(df_results_tmp2$DB[i]) < 0.1){
        df_results_tmp2$Fscore[i]<- (0.6*(as.numeric(df_results_tmp2$D[i])))+(0.4*(as.numeric(df_results_tmp2$Delta[i]))) - 0.12
    }else if (as.numeric(df_results_tmp2$DB[i]) > 0.35 && as.numeric(df_results_tmp2$DB[i]) <= 0.4 ){
        df_results_tmp2$Fscore[i]<- (0.6*(as.numeric(df_results_tmp2$D[i])))+(0.4*(as.numeric(df_results_tmp2$Delta[i]))) - 0.12
    }else if(as.numeric(df_results_tmp2$DB[i]) > 0.4 && as.numeric(df_results_tmp2$DB[i]) <= 0.45){
        df_results_tmp2$Fscore[i]<- (0.6*(as.numeric(df_results_tmp2$D[i])))+(0.4*(as.numeric(df_results_tmp2$Delta[i]))) - 0.18
    }else if(as.numeric(df_results_tmp2$DB[i]) > 0.45){
        df_results_tmp2$Fscore[i]<- (0.6*(as.numeric(df_results_tmp2$D[i])))+(0.4*(as.numeric(df_results_tmp2$Delta[i]))) - 0.18
    }else{
        df_results_tmp2$Fscore[i]<- (0.6*(as.numeric(df_results_tmp2$D[i])))+(0.4*(as.numeric(df_results_tmp2$Delta[i]))) - 0
    }
}


df_results_tmp3<-df_results_tmp2[which(df_results_tmp2$Fscore > 0),]
save(df_results_tmp3, file =paste0(mgf,"results.Rdata"))
write.table(df_results_tmp3, file=paste0(mgf,"results.txt"), row.names=F, col.names=F, quote = FALSE, sep="\t")
