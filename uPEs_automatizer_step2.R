args=(commandArgs(TRUE))
library(icesTAF)

source("~/data/01_Rscripts/A_Funciones/rlain.R")

directory<-args[1]
database<-args[2]
datos<-args[3]
outname_chunk<-args[4]

if(datos=="TCGA"){

        file_loaded<-get(load(paste0(directory,"TCGA_c_Ene20.rda")))
        file_loaded<-as.matrix(file_loaded)

        dir_querys<-paste0(directory,"chunks_",database)
        namecleaned<-unlist(lapply(strsplit(paste(outname_chunk),"/"),"[",2))
        chunk_query<-paste0(dir_querys,"/",namecleaned)
        chunk<-get(load(chunk_query))
        name_chunk<-basename(chunk_query)
        name_chunk<-lapply(strsplit(paste(name_chunk),"\\."),"[",1)

        uPE1_enrichedGO <- enrichGeneDB(as.matrix(file_loaded), as.matrix(chunk), fdr = FALSE, p_mat=0.01)
        output_dir<-paste0(directory,"OUTPUT_",datos,"_",database)
        dir.create(output_dir, showWarnings = FALSE)
        save(uPE1_enrichedGO, file=paste0(output_dir,"/",name_chunk,"_output.Rdata"))

}else if(datos=="GTEX"){

        file_loaded<-get(load(paste0(directory,"GTEX_c_Ene20.rda")))

        dir_querys<-paste0(directory,"chunks_",database)
        namecleaned<-unlist(lapply(strsplit(paste(outname_chunk),"/"),"[",2))
        chunk_query<-paste0(dir_querys,"/",namecleaned)
        chunk<-get(load(chunk_query))
        name_chunk<-basename(chunk_query)
        name_chunk<-lapply(strsplit(paste(name_chunk),"\\."),"[",1)

        uPE1_enrichedGO <- enrichGeneDB(as.matrix(file_loaded), chunk, fdr = FALSE, p_mat=0.01)
        output_dir<-paste0(directory,"OUTPUT_",datos,"_",database)
        dir.create(output_dir, showWarnings = FALSE)
        save(uPE1_enrichedGO, file=paste0(output_dir,"/",name_chunk,"_output.Rdata"))

}else if(datos=="CCLE"){

        file_loaded<-get(load(paste0(directory,"CCLE_c_Ene20.rda")))

        dir_querys<-paste0(directory,"chunks_",database)
        namecleaned<-unlist(lapply(strsplit(paste(outname_chunk),"/"),"[",2))
        chunk_query<-paste0(dir_querys,"/",namecleaned)
        chunk<-get(load(chunk_query))
        name_chunk<-basename(chunk_query)
        name_chunk<-lapply(strsplit(paste(name_chunk),"\\."),"[",1)

        uPE1_enrichedGO <- enrichGeneDB(as.matrix(file_loaded), as.matrix(chunk), fdr = FALSE, p_mat=0.01)
        output_dir<-paste0(directory,"OUTPUT_",datos,"_",database)
        dir.create(output_dir, showWarnings = FALSE)
        save(uPE1_enrichedGO, file=paste0(output_dir,"/",name_chunk,"_output.Rdata"))

}
