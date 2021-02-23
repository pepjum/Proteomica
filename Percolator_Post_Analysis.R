library(stringr)
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesVikv2.R")


fileinput_comet_percolator<-read.table("/home/nostromo/data/pepe/08_PERCOLATOR_JUN18/Comet_files/PXD001381/percolator/ALL_COMET_REVISED.pinresults_peptides.tsv.out", header=T)
fileinput_tandem_percolator<-read.table("/home/nostromo/data/pepe/EMBRIO_9_03/Tandem_Files/PXD001381/percolator/ALL_TANDEM_rev.pinresults_peptides.tsv.out", header=T)
fileinput_mascot_percolator<-read.table("/home/nostromo/data/pepe/EMBRIO_9_03/Dat_Files/PXD001381/percolator/ALL_MASCOT.pinresults_peptides.tsv.out", row.names=NULL)

fileinput_mascot_percolator$number<-seq(1:nrow(fileinput_mascot_percolator))
names(fileinput_mascot_percolator)<-c("PSMId","score","q.value","posterior_error_prob","peptide","ProteinAccession","query","number")




df_results_comet<-fileinput_comet_percolator[(fileinput_comet_percolator$q.value <0.01),]
df_results_tandem<-fileinput_tandem_percolator[(fileinput_tandem_percolator$q.value <0.01),]
df_results_mascot<-fileinput_mascot_percolator[(fileinput_mascot_percolator$q.value <0.01),]


list_peptides_comet<-paste(df_results_comet$peptide)
list_peptides_tandem<-paste(df_results_tandem$peptide)
list_peptides_mascot<-paste(df_results_mascot$peptide)


PATTERN <- '\\[[0-9.-]+\\]'
PATTERN2 <- '[A-Z]{2,}'
comet_good_list<-str_extract(str_remove_all(list_peptides_comet, PATTERN), PATTERN2)
tandem_good_list<-str_extract(str_remove_all(list_peptides_tandem, PATTERN), PATTERN2)
mascot_good_list<-str_extract(str_remove_all(list_peptides_mascot, PATTERN), PATTERN2)

unique_peptides_comet<-unique(comet_good_list)
unique_peptides_tandem<-unique(tandem_good_list)
unique_peptides_mascot<-unique(mascot_good_list)
#plots
pdf("COMPARE3SE_PERCOLATOR_PXD001383.pdf")
compare3List(unique_peptides_comet,unique_peptides_tandem,unique_peptides_mascot,"comet_percolator","tandem_percolator","mascot_percolator","PERCOLATOR PXD001383")
dev.off()

intersec_list_tmp<-intersect(unique_peptides_comet,unique_peptides_tandem)
intersect_list_def<-intersect(intersec_list_tmp,unique_peptides_mascot)


write.table(intersect_list_def, file="INTERSECT_3SE_PXD001383_PERCOLATOR.txt", col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\n")

####ALBA FDR

load("/home/nostromo/data/pepe/08_PERCOLATOR_JUN18/Comet_files/results_Peptides_comet_PXD001381_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter.rda")
comet<-dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter

load("/home/nostromo/data/pepe/08_PERCOLATOR_JUN18/Tandem_Files/results_Peptides_tandem_PXD001381_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter.rda")
tandem<-dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter

load("/home/nostromo/data/pepe/EMBRIO_9_03/Dat_Files/results_Peptides_mascot_PXD001381_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter.rda")
mascot<-dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter

mascot_ok<-mascot[-(grep("^DECOY",mascot$ProteinAccession)),]
tandem_ok<-tandem[-(grep("^DECOY",tandem$ProteinAccession)),]
comet_ok<-comet[-(grep("^DECOY",comet$ProteinAccession)),]


pep_comet<-unique(paste(comet_ok$PeptideSeq))
pep_tandem<-unique(paste(tandem_ok$PeptideSeq))
pep_mascot<-unique(paste(mascot_ok$PeptideSeq))

pdf("COMPARE3SE_HPP_PXD001383.pdf")
compare3List(pep_comet,pep_tandem,pep_mascot,"comet","tandem","mascot","FDR PXD001383")
dev.off()

pep_ALL_fdr<-unique(c(pep_comet,pep_tandem,pep_mascot))
pep_ALL_PERCOLATOR<-unique(c(unique_peptides_comet,unique_peptides_tandem,unique_peptides_mascot))

pdf("COMPARE_HPP_PERCOLATOR_PXD001383.pdf")
compare2List(pep_ALL_fdr,pep_ALL_PERCOLATOR,"FDR","PERCOLATOR","PXD001383")
dev.off()
intersect_PERC_HPP<-intersect(pep_ALL_fdr,pep_ALL_PERCOLATOR)
write.table(intersect_PERC_HPP, file="INTERSECT_PXD001383__HPP_PERCOLATOR.txt", col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\n")
diff_HPP<-setdiff(pep_ALL_fdr,pep_ALL_PERCOLATOR)
write.table(diff_HPP, file="PXD001383_HITS_HPP_NOHITS_PERCOLATOR.txt", col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\n")
diff_PERCOLATOR<-setdiff(pep_ALL_PERCOLATOR,pep_ALL_fdr)
write.table(diff_PERCOLATOR, file="PXD001383_HITS_PERCOLATOR_NOHITS_HPP.txt", col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\n")
