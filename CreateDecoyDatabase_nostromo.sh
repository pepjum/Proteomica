#! /bin/bash


FASTAFILE=$1	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/28_PFortes_Shotgun_lncRNA_Feb18/uniprot_sprot_2017_12_CRAP.fasta
DECOYID=$2  	# DECOY

FASTAFILEFOLDER=$(dirname "${FASTAFILE}")
FASTAFILENAME=$(basename "${FASTAFILE}")
FASTAFILENAMENOEXT="${FASTAFILENAME%.*}"
FASTAFILENAMEEXT="${FASTAFILENAME##*.}"
DECOYFASTAFILE=$FASTAFILEFOLDER'/'$FASTAFILENAMENOEXT'_'$DECOYID'.'$FASTAFILENAMEEXT

echo $DECOYFASTAFILE
Rscript /home/nostromo/data/pepe/jgonzalez.69/CreateDecoyDatabase_nostromo.R $FASTAFILE $DECOYFASTAFILE $DECOYID
