args=(commandArgs(TRUE))


mgf<-args[1]      #mgf file
database<<-args[2]   #database spectaST
rango<<-args[3]     #for selecting spectra by mz difference of .... (0.1 by default)
spread_value<<-args[4]    #spread out
Dalton<<-args[5]       # Value of daltons for bining function purposes
output<-args[6]   # 1 for Rdata, 2 for txt, 3 for feather

source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesShotgun.R")
source("/home/margaret/data/pepe/scripts/functions_SE.R")
#mgf_dir="."
cat("started at","\n")
old<-Sys.time()
print(old)

options(stringsAsFactors = FALSE)

library(stringr)
library(stringi)
library(plyr)
library(future.apply)
library(compiler)
library(feather)
#setwd(paste0(mgf_dir,"/"))


bining_c<-cmpfun(bining)

enableJIT(3)
cat("reading MGF file","\n")

ReadMGFFile_c<-cmpfun(ReadMGFFile)

plan(multicore, workers=40)

options(future.globals.maxSize= +Inf)

mgf_file<-ReadMGFFile_c(mgf)

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
df_spectrums_mgf$idx<-seq(1:nrow(df_spectrums_mgf))

rm(mgf_file)

cat("loading database","\n")
db_spectra<<-get(load(database))   #_500-1500.Rdata


#for (i in 1:nrow(df_spectrums_mgf)){
searching<-function(row_dataframe){
    spectrum<-row_dataframe

    PrecursorMZ<-as.numeric(paste(spectrum["mz_value"]))
    threshold_up <- (as.numeric(PrecursorMZ)+ as.numeric(rango))
    threshold_down <- (as.numeric(PrecursorMZ)- as.numeric(rango))
    l_up <- which(as.numeric(db_spectra$PrecursorMZ) < threshold_up)
    l_down <- which(as.numeric(db_spectra$PrecursorMZ) > threshold_down)
    tmp <- db_spectra[intersect(l_up,l_down ),]
    if (nrow(tmp) != 0){
        list_all_peaks<-paste(spectrum["peaks"])
        tmp_mgf<-strsplit(paste(list_all_peaks), "###")
        m<-paste(lapply(strsplit(unlist(tmp_mgf)," "),'[',1))
        z<-paste(lapply(strsplit(unlist(tmp_mgf)," "),'[',2))
        df_peaks<-data.frame("M"=m,"Z"=z)
        df_peaks_binned<-bining_c(df_peaks,1)
        if (spread_value !=0){
        cat("spread_out","\n")
        df_peaks_binned_and_spread_out<-spread_out(df_peaks_binned,spread_value)
        }else{
        df_peaks_binned_and_spread_out<-df_peaks_binned
        }
        df_ranges_tmp<-df_peaks_binned_and_spread_out[as.numeric(paste(df_peaks_binned_and_spread_out$M))>=500,]
        df_ranges_final<-df_ranges_tmp[as.numeric(paste(df_ranges_tmp$M))<=1500,]
        df_ranges_final$normalized<-as.numeric(df_ranges_final$Z)/sqrt(sum(as.numeric(df_ranges_final$Z)^2))
        matrix<-data.frame(M=as.numeric(seq(1:1500)),normalized=0)
        if (nrow(df_ranges_final)==0){
            df_ranges_final<- matrix
        }
        matrix[which(as.numeric(matrix$M) %in% as.numeric(df_ranges_final$M)),]<-df_ranges_final$normalized
        names(matrix)[2]<-"Z"
        df_ranges_tmp<-matrix[as.numeric(paste(matrix$M))>=500,]
        df_ranges_final<-df_ranges_tmp[as.numeric(paste(df_ranges_tmp$M))<=1500,]
        df_ranges_final$line<-paste0(df_ranges_final$M,"\t",df_ranges_final$Z,"###")
        linecod<-paste0(df_ranges_final$line)
        linecod<-paste(linecod,collapse="")
        spectrum<-spectrum[,-3]
        spectrum["Peaks_mgf"]<-paste(linecod)
        number_reps<-nrow(tmp)
        list_name<-rep(spectrum["spectrum"], number_reps)
        list_mz<-rep(spectrum["mz_value"], number_reps)
        list_peaks<-rep(spectrum["Peaks_mgf"], number_reps)
        df_mgfs<-data.frame("spectrum" = paste(list_name), "PrecursorMGF" = paste(list_mz), "Peaks_mgf" = paste(list_peaks))
        df_results_tmp<-cbind(df_mgfs,tmp)
        res<-ldply(apply(df_results_tmp,1, FUN=function(x) Dot_and_DB(x)),rbind)
        res<-res[,c(2,3)]
        df_results_tmp<-cbind(df_results_tmp,res)
        df_results_tmp$Peaks_mgf<-NULL
        df_results_tmp$Peaks<- NULL
        df_results_tmp<-df_results_tmp[which(df_results_tmp$D > 0),]
        df_results_tmp<-df_results_tmp[with(df_results_tmp, order(-D)),]
        df_results_tmp[which(is.nan(df_results_tmp$DB)), "DB"]<-0
        tmp_max1<-df_results_tmp[(as.numeric(paste(df_results_tmp$D)) ==max(as.numeric(paste(df_results_tmp$D)))),]
        tmp_max2<-df_results_tmp[-c(as.numeric(paste(df_results_tmp$D)) ==max(as.numeric(paste(df_results_tmp$D)))),]
        tmp_max3<-tmp_max2[(as.numeric(paste(tmp_max2$D)) ==max(as.numeric(paste(tmp_max2$D)))),]
        D1<-as.numeric(tmp_max1$D)
        D2<-as.numeric(tmp_max3$D)
        D3<-(D1-D2)
        df_results_tmp$Delta<-as.numeric(D3)
        res_2<-ldply(apply(df_results_tmp,1, FUN = function(x) F_score(x)),rbind)
        names(res_2)[2]<-"Fscore"
        df_results_tmp<-cbind(df_results_tmp,res_2)

    }else{
        return(NA)
    }
}

df_results_final<-ldply(future_apply(df_spectrums_mgf,1, FUN=function(x) searching(x),future.seed=TRUE),rbind)



new_name<-paste(substr(mgf, 1, nchar(mgf) - 3), "results", sep="")

if(output == 1){
    save(df_results_final, file = paste0(new_name,".Rdata"))
}else if(output == 2){
    write.table(df_results_final, file= paste0(new_name,".txt"), row.names=F, col.names=F, quote = FALSE, sep="\t")
}else if(output == 3){
    write_feather(df_results_final, paste0(new_name,".feather"))
}

new<-Sys.time() - old
print(new)
cat("end of search")
