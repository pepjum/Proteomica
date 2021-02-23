#!/bin/bash


currentfolder=$1	# folder where all the data is stored. There must be one folder: MGFFiles
comettargetfile=$2	# comet parameter file for target searches
cometdecoyfile=$3	# comet parameter file for decoy searches
dbtargetfile=$4		# database file for target searches
dbdecoyfile=$5		# database file for decoy searches
dbrdafile=$6		# several options (rda file, 0, 1, 2, 3) => rda file: database rda file for MAYU protFDR analyses, 0: psmFDR, 1: protFDRMAYU, 2: pepFDR, 3: protFDR
decoyid=$7			# decoy sequence identifier.
MGFFOLDER=$8 		# mgfFolder with all the experiments to be done

EXP=$MGFFOLDER'mgfFolderList.txt'
#find $MGFFOLDER -mindepth 1 -maxdepth 1 -type d -print > $EXP   #descomentar luego

cd /home/nostromo/data/pepe/scripts/
nParameters=$#
if [[ $nParameters -eq 8 ]]
then
	./PeptideMatch_nostromo.sh "${currentfolder}" "${dbtargetfile}" "${dbdecoyfile}" "${decoyid}"
	echo 'antes'
	while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line
		#line="${line%"${line##*[![:space:]]}"}"
		echo 'entro en Cometsearch'            # remove spaces from line
		dataset=$line
		./CometSearch_nostromo.sh "${currentfolder}" "${dataset}" "${comettargetfile}" "${cometdecoyfile}" "${dbtargetfile}" "${dbdecoyfile}" "${dbrdafile}" "${decoyid}"
	done < "${EXP}"
	echo "termina cometsearch"
	while IFS='' read -r line || [[ -n "$line" ]]; do
		echo "entro en CometPep2Prot"
		dataset=$line
		datasetname=$(basename "${dataset}")
		filename=$currentfolder"Comet_files/Shotgun_comet_"$datasetname".txt"
		./CometPep2Prot_nostromo.sh /home/nostromo/data/pepe/scripts/CometPep2Prot_nostromo.R "${currentfolder}" "${dataset}" "${filename}" "${dbtargetfile}" "${decoyid}"
	done < "${EXP}"
	# ##### voy por aqui
	echo "termina cometpep2prot"
	while IFS='' read -r line || [[ -n "$line" ]]; do
		dataset=$line
		datasetname=$(basename "${dataset}")
		filename=$currentfolder"Comet_files/Shotgun_comet_"$datasetname".txt"
	  Rscript /home/nostromo/data/pepe/scripts/ShotgunAnalysisAutomated_nostromo.R $filename $datasetname 0 $decoyid 3 2 9 50

	done < "${EXP}"

else
	echo "Parameters are missing."
	echo "USAGE: ./CometLauncher_nostromo.sh currentfolder comettargetfile cometdecoyfile dbtargetfile dbdecoyfile dbrdafile decoyid mgfFolderPath "
	echo "EXAMPLE: ./CometLauncher_nostromo.sh /home/nostromo/data/pepe/26_Navarrabiomed_Missing_Ene18/ comet_Swissprot201712.txt comet_Swissprot201712_D.txt uniprot_sprot_2017_12.fasta uniprot_sprot_2017_12_D.fasta 3 DECOY /home/nostromo/data/pepe/26_Navarrabiomed_Missing_Ene18/Prueba2/MGFFiles/"
fi
