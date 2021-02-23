#!/bin/bash


CURRENTFOLDER=$1	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/28_PFortes_Shotgun_lncRNA_Feb18/
FASTAFILE=$2		# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/23_neXtprot_20170801_Nov2017/nextProtDB20170801.fasta
NEXTPROTFOLDER=$3	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/23_neXtprot_20170801_Nov2017/

FASTAFILEFOLDER=$(dirname "${FASTAFILE}")
FASTAFILENAME=$(basename "${FASTAFILE}")
FASTAFILENAMENOEXT="${FASTAFILENAME%.*}"
FASTAFILENAMEEXT="${FASTAFILENAME##*.}"

perl /opt/proteogest/proteogest.pl -i $FASTAFILE -d -a -c trypsin -R M -W 15.99 -G 1 -L 2

fileroot=$FASTAFILEFOLDER'/'$FASTAFILENAMENOEXT

echo "starting the digestion of the database"

Rscript /home/nostromo/data/pepe/scripts/DigestDatabase.R $fileroot

echo "starting to process the database"

Rscript /home/nostromo/data/pepe/scripts/ProcessNextprot_nostromo.R $NEXTPROTFOLDER

echo "starting to detection of missing proteins"

Rscript /home/nostromo/data/pepe/scripts/MissingProteinDetection_Step1_nostromo.R $CURRENTFOLDER $fileroot $NEXTPROTFOLDER

echo "The detection of missing proteins has finished."
