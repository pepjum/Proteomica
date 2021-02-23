#!/bin/bash


RFILE=$1 			# /mnt/beegfs/agarin/dato-activo/02_ClusterScripts/agarin/CometPep2Prot.R
CURRENTFOLDER=$2	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/
DATASET=$3			# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/MGFFiles/349-028
FILENAME=$4 		# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/Comet_files/Shotgun_comet_349-028.txt
DBTARGETFILE=$5 	# uniprot_sprot_2017_12.fasta
DECOYID=$6			# DECOY

datasetname=$(basename "${DATASET}")

mgfFiles=$(find "${DATASET}" -maxdepth 1 -iname "*.mgf")
mgfFiles=( $mgfFiles )
nMGFFiles=${#mgfFiles[@]}

datasetCometTargetFolder=$( echo $DATASET | sed -e  "s/MGFFiles/Comet_files/g")
datasetCometDecoyFolder=$datasetCometTargetFolder'-D/'
datasetCometTargetFolder=$datasetCometTargetFolder'/'
CometTargetFiles=$(find "${datasetCometTargetFolder}" -iname "*.txt" -not -name "*_peptides.txt" -not -name "*_peptides_out.txt" -not -name "*_corrected.txt" -not -name "*_peptides_log.log")
CometTargetFiles=( $CometTargetFiles )
nCometTargetFiles=${#CometTargetFiles[@]}
CorrectedTargetFiles=$(find "${datasetCometTargetFolder}" -iname "*_corrected.txt")
CorrectedTargetFiles=( $CorrectedTargetFiles )
nCorrectedTargetFiles=${#CorrectedTargetFiles[@]}
PeptideTargetFiles=$(find "${datasetCometTargetFolder}" -iname "*_peptides.txt")
PeptideTargetFiles=( $PeptideTargetFiles )
nPeptideTargetFiles=${#PeptideTargetFiles[@]}
PeptideOutTargetFiles=$(find "${datasetCometTargetFolder}" -iname "*_peptides_out.txt")
PeptideOutTargetFiles=( $PeptideOutTargetFiles )
nPeptideOutTargetFiles=${#PeptideOutTargetFiles[@]}
PeptideLogTargetFiles=$(find "${datasetCometTargetFolder}" -iname "*_peptides_log.log")
PeptideLogTargetFiles=( $PeptideLogTargetFiles )
nPeptideLogTargetFiles=${#PeptideLogTargetFiles[@]}


CometDecoyFiles=$(find "${datasetCometDecoyFolder}" -iname "*.txt" -not -name "*_peptides.txt" -not -name "*_peptides_out.txt" -not -name "*_corrected.txt" -not -name "*_peptides_log.log")
CometDecoyFiles=( $CometDecoyFiles )
nCometDecoyFiles=${#CometDecoyFiles[@]}
CorrectedDecoyFiles=$(find "${datasetCometDecoyFolder}" -iname "*_corrected.txt")
CorrectedDecoyFiles=( $CorrectedDecoyFiles )
nCorrectedDecoyFiles=${#CorrectedDecoyFiles[@]}
PeptideDecoyFiles=$(find "${datasetCometDecoyFolder}" -iname "*_peptides.txt")
PeptideDecoyFiles=( $PeptideDecoyFiles )
nPeptideDecoyFiles=${#PeptideDecoyFiles[@]}
PeptideOutDecoyFiles=$(find "${datasetCometDecoyFolder}" -iname "*_peptides_out.txt")
PeptideOutDecoyFiles=( $PeptideOutDecoyFiles )
nPeptideOutDecoyFiles=${#PeptideOutDecoyFiles[@]}
PeptideLogDecoyFiles=$(find "${datasetCometDecoyFolder}" -iname "*_peptides_log.log")
PeptideLogDecoyFiles=( $PeptideLogDecoyFiles )
nPeptideLogDecoyFiles=${#PeptideLogDecoyFiles[@]}

target=0
if [ $nMGFFiles != $nCometTargetFiles ]; then
	echo "There are no search files in the folder. Try agarin"
	exit
elif [ $nMGFFiles != $nCorrectedTargetFiles ]; then
	target=1
elif [ $nMGFFiles != $nPeptideTargetFiles ]; then
	target=1
elif [ $nMGFFiles != $nPeptideOutTargetFiles ]; then
	target=1
elif [ $nMGFFiles != $nPeptideLogTargetFiles ]; then
	target=1
else
	target=0
fi

decoy=0
if [ $nMGFFiles != $nCometDecoyFiles ]; then
	echo "There are no search files in the folder. Try agarin DECOY"
	exit
elif [ $nMGFFiles != $nCorrectedDecoyFiles ]; then
	decoy=1
elif [ $nMGFFiles != $nPeptideDecoyFiles ]; then
	decoy=1
elif [ $nMGFFiles != $nPeptideOutDecoyFiles ]; then
	decoy=1
elif [ $nMGFFiles != $nPeptideLogDecoyFiles ]; then
	decoy=1
else
	decoy=0
fi

echo $datasetname
if [ $target -eq 1 ]; then
	if [ $decoy -eq 1 ]; then
		TOD=0
		echo "doing PeptideMatch for target and decoy files."
	else
		TOD=1
		echo "doing PeptideMatch for target files. Decoy files are already processed."
	fi
else
	if [ $decoy -eq 1 ]; then
		TOD=2
		echo "doing PeptideMatch for decoy files. Target files are already processed."
	else
		TOD=3
	fi
fi

if [ $TOD -ne 3 ]; then
	# TOD; 1 si solo hay que hacer target; 2 si solo hay que hacer decoy; 0 si hay que hacer los dos.
  Rscript $RFILE $CURRENTFOLDER $datasetname $FILENAME $DBTARGETFILE $DECOYID $TOD
else
	echo "PeptideMatch is already done. Continue."
fi

# echo $datasetname
# echo $datasetCometTargetFolder
# echo $datasetCometDecoyFolder
# echo $target
# echo $decoy
# echo $TOD

# echo $nMGFFiles
# echo $nCometTargetFiles
# echo $nCorrectedTargetFiles
# echo $nPeptideTargetFiles
# echo $nPeptideOutTargetFiles
# echo $nPeptideLogTargetFiles
# echo $nCometDecoyFiles
# echo $nCorrectedDecoyFiles
# echo $nPeptideDecoyFiles
# echo $nPeptideOutDecoyFiles
# echo $nPeptideLogDecoyFiles
