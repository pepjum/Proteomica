### step 3
args=(commandArgs(TRUE))

folder_output<-args[1]
datos<-args[2]


files<-list.files(folder_output, pattern=".Rdata")
files<-paste0(folder_output,files)

TCGA_uPE1_enrichedGO_pval<-data.frame()
TCGA_uPE1_enrichedGO_geneXenrichedGO<-data.frame()
for(i in 1:length(files)){
    cat(i,"\n")
    loaded<-get(load(files[[i]]))

    TCGA_uPE1_enrichedGO_pval<-rbind(TCGA_uPE1_enrichedGO_pval,loaded$pval)
    TCGA_uPE1_enrichedGO_geneXenrichedGO<-rbind(TCGA_uPE1_enrichedGO_geneXenrichedGO, loaded$geneXenrichedGO)
    cat(dim(TCGA_uPE1_enrichedGO_pval),"\n")
}

TCGA_uPE1_enrichedGO_pval_sel <- (apply(TCGA_uPE1_enrichedGO_pval,2,FUN=function(x) (x<0.01)*1))
TCGA_uPE1_enrichedGO_pval_sel[is.na(TCGA_uPE1_enrichedGO_pval_sel)] <- 0

write.table(cbind(GO_id=rownames(TCGA_uPE1_enrichedGO_pval_sel),TCGA_uPE1_enrichedGO_pval_sel), file=paste0(folder_output,"TCGA_enrichedGOmat.txt"), row.names=F, sep="\t", quote=F)
