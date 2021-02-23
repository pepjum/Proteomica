#! /bin/bash


FASTAFILE=$1	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/28_PFortes_Shotgun_lncRNA_Feb18/uniprot_sprot_2017_12_CRAP.fasta

FASTAFILEFOLDER=$(dirname "${FASTAFILE}")
FASTAFILENAME=$(basename "${FASTAFILE}")
FASTAFILENAMENOEXT="${FASTAFILENAME%.*}"
FASTAFILENAMEEXT="${FASTAFILENAME##*.}"

perl /opt/proteogest/proteogest.pl -i $FASTAFILE -d -a -c trypsin -R M -W 15.99 -G 1

fileroot=$FASTAFILEFOLDER'/'$FASTAFILENAMENOEXT

Rscript /home/nostromo/data/pepe/scripts/DigestDatabase.R $fileroot
