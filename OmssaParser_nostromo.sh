#!/bin/bash

CURRENTFOLDER=$1	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/
DATASETNAME=$2		# 349-028

omssaTargetFolder=$CURRENTFOLDER'Omssa_files/'$DATASETNAME'/'

xmlFilesT=$(find "${omssaTargetFolder}" -maxdepth 1 -iname "*.xml")
xmlFilesT=( $xmlFilesT )

for ((i=0; i < ${#xmlFilesT[@]}; i ++))
do
		xmlFilesTfile=${xmlFilesT[$i]}
		xmlFilesTfilenoext="${xmlFilesTfile%.*}"
		Rscript /home/nostromo/data/pepe/jgonzalez.69/OmssaParser.R $xmlFilesTfilenoext
done

omssaDecoyFolder=$CURRENTFOLDER'Omssa_files/'$DATASETNAME'-D/'

xmlFilesD=$(find "${omssaDecoyFolder}" -maxdepth 1 -iname "*.xml")
xmlFilesD=( $xmlFilesD )

for ((i=0; i < ${#xmlFilesD[@]}; i ++))
do
		xmlFilesDfile=${xmlFilesD[$i]}
		xmlFilesDfilenoext="${xmlFilesDfile%.*}"
		Rscript /home/nostromo/data/pepe/jgonzalez.69/OmssaParser.R $xmlFilesDfilenoext 
done

wait
