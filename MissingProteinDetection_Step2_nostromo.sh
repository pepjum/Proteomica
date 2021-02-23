#!/bin/bash


CURRENTFOLDER=$1	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/28_PFortes_Shotgun_lncRNA_Feb18/
FASTAFILE=$2		# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/23_neXtprot_20170801_Nov2017/nextProtDB20170801.fasta

FASTAFILEFOLDER=$(dirname "${FASTAFILE}")
FASTAFILENAME=$(basename "${FASTAFILE}")
FASTAFILENAMENOEXT="${FASTAFILENAME%.*}"
FASTAFILENAMEEXT="${FASTAFILENAME##*.}"

fileroot=$FASTAFILEFOLDER'/'$FASTAFILENAMENOEXT

Rscript /home/nostromo/data/pepe/scripts/MissingProteinDetection_Step2_nostromo.R $CURRENTFOLDER $fileroot
