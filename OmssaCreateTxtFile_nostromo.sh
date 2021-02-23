#! /bin/bash


CURRENTFOLDER=$1	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/
DATASETNAME=$2		# 349-028

Rscript /home/nostromo/data/pepe/jgonzalez.69/OmssaCreateTxtFile_nostromo.R $CURRENTFOLDER $DATASETNAME
