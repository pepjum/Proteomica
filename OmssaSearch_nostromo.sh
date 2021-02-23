#!/bin/bash

CURRENTFOLDER=$1	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/
DATASET=$2			# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/MGFFiles/349-028
DBINDEXROOT=$3	 	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/26_Navarrabiomed_Missing_Ene18/uniprot_sprot_2017_12
QUALITY=$4 			# HIGH; MEDIUM OR LOW.

echo $DATASET

datasetname=$(basename "${DATASET}")
omssafolder=$CURRENTFOLDER'Omssa_files'
mgfFolder=$CURRENTFOLDER'/MGFFiles/'
mkdir -p $omssafolder

mgfFiles=$(find "${DATASET}" -maxdepth 1 -iname "*.mgf")
mgfFiles=( $mgfFiles )

datasetOmssaTargetFolder=$omssafolder'/'$datasetname'/'
mkdir -p $datasetOmssaTargetFolder
datasetOmssaDecoyFolder=$omssafolder'/'$datasetname'-D/'
mkdir -p $datasetOmssaDecoyFolder

nMGFFiles=${#mgfFiles[@]}
OmssaFilesTarget=$(find "${datasetOmssaTargetFolder}" -iname "*.txt" -not -name "*_peptides.txt" -not -name "*_peptides_out.txt" -not -name "*_corrected.txt" -not -name "*_peptides_log.log")
OmssaFilesTarget=( $OmssaFilesTarget )
nOmssaFilesTarget=${#OmssaFilesTarget[@]}

OmssaFilesDecoy=$(find "${datasetOmssaDecoyFolder}" -iname "*.txt" -not -name "*_peptides.txt" -not -name "*_peptides_out.txt" -not -name "*_corrected.txt" -not -name "*_peptides_log.log")
OmssaFilesDecoy=( $OmssaFilesDecoy )
nOmssaFilesDecoy=${#OmssaFilesDecoy[@]}

if [ $nMGFFiles != $nOmssaFilesTarget ]; then
	echo 'Omssa target DB search...'
	for ((i=0; i < ${#mgfFiles[@]}; i ++))
	do
		mgfile=${mgfFiles[$i]}
		mgffileroot="${mgfile%.*}"
		echo $mgfile
		DBINDEXFILE=$DBINDEXROOT'_index'
		# e trypsin; -i b,y ions; -mf carbamidomethyl C static mod; -mv oxidized M variable mod;
		if  [ $QUALITY == "HIGH" ]; then
			echo "high quality" # precursor tolerance 10 ppm; fragment tolerance 0.05; -hl retain top 10 hits; -v 1 cleavage
			/opt/omssa-2.1.9.linux/omssacl -w -il -e 0 -i 1,4 -mf 3 -mv 1 -te 10 -teppm -to 0.05 -hl 10 -v 1 -fm $mgffileroot".mgf" -d $DBINDEXFILE -o $mgffileroot".txt" -fomx $mgffileroot".omx"  -op $mgffileroot".pepXML" -ox $mgffileroot".xml" -oc $mgffileroot"_summary.csv" -logfile$mgffileroot".log"
		elif [ $QUALITY == "MEDIUM" ]; then
			echo "medium quality" # precursor tolerance 20 ppm; fragment tolerance 0.05; -hl retain top 10 hits; -v 1 cleavage
			/opt/omssa-2.1.9.linux/omssacl -w -il -e 0 -i 1,4 -mf 3 -mv 1 -te 20 -teppm -to 0.05 -hl 10 -v 1 -fm $mgffileroot".mgf" -d $DBINDEXFILE -o $mgffileroot".txt" -fomx $mgffileroot".omx" -op $mgffileroot".pepXML" -ox $mgffileroot".xml" -oc $mgffileroot"_summary.csv" -logfile$mgffileroot".log"
		else
			echo "low quality" # precursor tolerance 20 ppm; fragment tolerance 0.5; -hl retain top 10 hits; -v 1 cleavage
			/opt/omssa-2.1.9.linux/omssacl -w -il -e 0 -i 1,4 -mf 3 -mv 1 -te 20 -teppm -to 0.5 -hl 10 -v 1 -fm $mgffileroot".mgf" -d $DBINDEXFILE -o $mgffileroot".txt" -fomx $mgffileroot".omx" -op $mgffileroot".pepXML" -ox $mgffileroot".xml" -oc $mgffileroot"_summary.csv" 	-logfile$mgffileroot".log"
		fi
		mv $mgffileroot".omx" $datasetOmssaTargetFolder
		mv $mgffileroot".txt" $datasetOmssaTargetFolder
		mv $mgffileroot".pepXML" $datasetOmssaTargetFolder
		mv $mgffileroot".xml" $datasetOmssaTargetFolder
		mv $mgffileroot"_summary.csv" $datasetOmssaTargetFolder
		mv $mgffileroot".log" $datasetOmssaTargetFolder
	done
else
	echo "Target searches were already done. Continue"
fi

if [ $nMGFFiles != $nOmssaFilesDecoy ]; then
	echo 'Omssa decoy DB search...'
	for ((i=0; i < ${#mgfFiles[@]}; i ++))
	do
		mgfile=${mgfFiles[$i]}
		mgffileroot="${mgfile%.*}"
		echo $mgfile
		DBINDEXFILE=$DBINDEXROOT'_DECOY_index'
		# e trypsin; -i b,y ions; -mf carbamidomethyl C static mod; -mv oxidized M variable mod;
		if  [ $QUALITY == "HIGH" ]; then
			echo "high quality" # precursor tolerance 10 ppm; fragment tolerance 0.05; -hl retain top 10 hits; -v 1 cleavage
			/opt/omssa-2.1.9.linux/omssacl -w -il -e 0 -i 1,4 -mf 3 -mv 1 -te 10 -teppm -to 0.05 -hl 10 -v 1 -fm $mgffileroot".mgf" -o $mgffileroot".txt" -fomx $mgffileroot".omx" -d $DBINDEXFILE -op $mgffileroot".pepXML" -ox $mgffileroot".xml" -oc $mgffileroot"_summary.csv" -logfile$mgffileroot".log"
		elif [ $QUALITY == "MEDIUM" ]; then
			echo "medium quality" # precursor tolerance 20 ppm; fragment tolerance 0.05; -hl retain top 10 hits; -v 1 cleavage
			/opt/omssa-2.1.9.linux/omssacl -w -il -e 0 -i 1,4 -mf 3 -mv 1 -te 20 -teppm -to 0.05 -hl 10 -v 1 -fm $mgffileroot".mgf" -o $mgffileroot".txt" -fomx $mgffileroot".omx" -d $DBINDEXFILE -op $mgffileroot".pepXML" -ox $mgffileroot".xml" -oc $mgffileroot"_summary.csv" -logfile$mgffileroot".log"
		else
			echo "low quality" # precursor tolerance 20 ppm; fragment tolerance 0.5; -hl retain top 10 hits; -v 1 cleavage
			/opt/omssa-2.1.9.linux/omssacl -w -il -e 0 -i 1,4 -mf 3 -mv 1 -te 20 -teppm -to 0.5 -hl 10 -v 1 -fm $mgffileroot".mgf" -o $mgffileroot".txt" -fomx $mgffileroot".omx" -d $DBINDEXFILE -op $mgffileroot".pepXML" -ox $mgffileroot".xml" -oc $mgffileroot"_summary.csv" -logfile$mgffileroot".log"
		fi
		mv $mgffileroot".omx" $datasetOmssaTargetFolder
		mv $mgffileroot".txt" $datasetOmssaTargetFolder
		mv $mgffileroot".pepXML" $datasetOmssaDecoyFolder
		mv $mgffileroot".xml" $datasetOmssaDecoyFolder
		mv $mgffileroot"_summary.csv" $datasetOmssaDecoyFolder
		mv $mgffileroot".log" $datasetOmssaTargetFolder
	done
else
	echo "Decoy searches were already done. Continue"
fi

wait
