#!/bin/bash



CURRENTFOLDER=$1	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/
DATASET=$2			# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/MGFFiles/349-028
COMETTARGETFILE=$3 	# comet_Swissprot201712.txt
COMETDECOYFILE=$4 	# comet_Swissprot201712_D.txt
DBTARGETFILE=$5 	# uniprot_sprot_2017_12.fasta
DBDECOYFILE=$6 		# uniprot_sprot_2017_12_D.fasta
DBRDAFILE=$7 		# 3
DECOYID=$8			# DECOY

#module load comet
#spack load jdk
#spack load java

datasetname=$(basename "${DATASET}")
cometfolder=$CURRENTFOLDER'Comet_files'
mgfFolder=$CURRENTFOLDER'/MGFFiles/'
mkdir -p $cometfolder

mgfFiles=$(find "${DATASET}" -maxdepth 1 -iname "*.mgf")
mgfFiles=( $mgfFiles )

datasetCometTargetFolder=$cometfolder'/'$datasetname'/'
mkdir -p $datasetCometTargetFolder
datasetCometDecoyFolder=$cometfolder'/'$datasetname'-D/'
mkdir -p $datasetCometDecoyFolder

nMGFFiles=${#mgfFiles[@]}
CometFilesTarget=$(find "${datasetCometTargetFolder}" -iname "*.txt" -not -name "*_peptides.txt" -not -name "*_peptides_out.txt" -not -name "*_corrected.txt" -not -name "*_peptides_log.log")
CometFilesTarget=( $CometFilesTarget )
nCometFilesTarget=${#CometFilesTarget[@]}

CometFilesDecoy=$(find "${datasetCometDecoyFolder}" -iname "*.txt" -not -name "*_peptides.txt" -not -name "*_peptides_out.txt" -not -name "*_corrected.txt" -not -name "*_peptides_log.log")
CometFilesDecoy=( $CometFilesDecoy )
nCometFilesDecoy=${#CometFilesDecoy[@]}

if [ $nMGFFiles != $nCometFilesTarget ]; then
	echo 'Comet target DB search...'
	param='-P'$CURRENTFOLDER$COMETTARGETFILE
	for ((i=0; i < ${#mgfFiles[@]}; i ++))
	do
		#echo $i;
		/usr/bin/comet.exe $param ${mgfFiles[$i]}     #en secuencial. En cluster poner &
	done
	wait
	cometFiles=$(find "${DATASET}" -maxdepth 1 -iname "*.txt")
	#cometFiles2=$(find "${DATASET}" -maxdepth 1 -iname "*.pep.xml")
	cometFiles3=$(find "${DATASET}" -maxdepth 1 -iname "*.pin")

	mv -t $datasetCometTargetFolder $cometFiles
	#mv -t $datasetCometTargetFolder $cometFiles2
	mv -t $datasetCometTargetFolder $cometFiles3





else
	echo "Target searches were already done. Continue"
fi

if [ $nMGFFiles != $nCometFilesDecoy ]; then
	echo 'Comet decoy DB search...'
	param='-P'$CURRENTFOLDER$COMETDECOYFILE
	for ((i=0; i < ${#mgfFiles[@]}; i ++))
	do
		/usr/bin/comet.exe $param ${mgfFiles[$i]}
	done
	wait
	cometFiles=$(find "${DATASET}" -maxdepth 1 -iname "*.txt")
	cometFiles2=$(find "${DATASET}" -maxdepth 1 -iname "*.pep.xml")
	cometfiles3=$(find "${DATASET}" -maxdepth 1 -iname "*.pin")
	mv -t $datasetCometDecoyFolder $cometFiles
	#mv -t $datasetCometDecoyFolder $cometFiles2
	mv -t $datasetCometDecoyFolder $cometFiles3

else
	echo "Decoy searches were already done. Continue"
fi

wait
