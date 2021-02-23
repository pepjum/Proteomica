#! /bin/bash


currentfolder=$1	# folder where all the data is stored. There must be one folder: MGFFiles
dbtargetfile=$2		# database file for target searches
dbdecoyfile=$3		# database file for decoy searches
dbrdafile=$4		# several options (rda file, 0, 1, 2, 3) => rda file: database rda file for MAYU protFDR analyses, 0: psmFDR, 1: protFDRMAYU, 2: pepFDR, 3: protFDR
decoyid=$5			# decoy sequence identifier.
MGFFOLDER=$6 		# mgfFolder with all the experiments to be done
QUALITY=$7			# HIGH; MEDIUM OR LOW.

# echo $currentfolder
# echo $dbtargetfile
# echo $dbdecoyfile
# echo $dbrdafile
# echo $decoyid
# echo $MGFFOLDER
# echo $QUALITY

if [[ $nParameters -eq 0 ]]
then
	EXP=$MGFFOLDER'mgfFolderList.txt'
	find $MGFFOLDER -mindepth 1 -maxdepth 1 -type d -print > $EXP

	dbtarget=$currentfolder$dbtargetfile
	dbtargetname="${dbtargetfile%.*}"
	dbdecoy=$currentfolder$dbdecoyfile
	dbdecoyname="${dbdecoyfile%.*}"

	# PREPARE THE DATABASE
	blastdbtargetout=$currentfolder$dbtargetname'_index'
	/home/nostromo/lotus_pipeline/bin/ncbi-blast-2.2.29+/bin/makeblastdb -in $dbtarget -dbtype prot
	/home/nostromo/lotus_pipeline/bin/ncbi-blast-2.2.29+/bin/blastdb_aliastool -dblist $dbtarget -dbtype prot -title $dbtargetname -out $blastdbtargetout
	echo $dbtarget
	echo $dbtargetname
	echo $blastdbtargetout
	echo $dbdecoy
	echo $dbdecoyname
	blastdbdecoyout=$currentfolder$dbdecoyname'_index'
	echo $blastdbdecoyout


	/home/nostromo/lotus_pipeline/bin/ncbi-blast-2.2.29+/bin/makeblastdb -in $dbdecoy -dbtype prot
	/home/nostromo/lotus_pipeline/bin/ncbi-blast-2.2.29+/bin/blastdb_aliastool -dblist $dbdecoy -dbtype prot -title $dbdecoyname -out $blastdbdecoyout

	exit   #OJO!!!
	while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line
		dataset=$line
		dbrootpath=$currentfolder$dbtargetname
		datasetname=$(basename "${dataset}")
		txtFileName=$currentfolder"Omssa_files/Shotgun_omssa_"$datasetname".txt"
		# SEARCH
		sbatchout=$(sbatch --job-name OSEARCH /mnt/beegfs/agarin/dato-activo/02_ClusterScripts/agarin/OmssaSearch.sbs "${currentfolder}" "${dataset}" "${dbrootpath}" "${QUALITY}")
		echo $sbatchout
		jobid1=$(echo ${sbatchout} | grep -oh "[1-9][0-9]*$")
		echo 'waiting for:'$jobid1
		# PARSE FROM XML TO TXT
		sbatchout2=$(sbatch --depend=afterok:${jobid1} --job-name OPARSE /mnt/beegfs/agarin/dato-activo/02_ClusterScripts/agarin/OmssaParser.sbs "${currentfolder}" "${datasetname}")
		echo $sbatchout2
		jobid2=$(echo ${sbatchout2} | grep -oh "[1-9][0-9]*$")
		echo 'waiting for:'$jobid2
		# CREATE THE TXT FILE
		sbatchout3=$(sbatch --depend=afterok:${jobid2} --job-name OTXT /mnt/beegfs/agarin/dato-activo/02_ClusterScripts/agarin/OmssaCreateTxtFile.sbs "${currentfolder}" "${datasetname}")
		echo $sbatchout3
		jobid3=$(echo ${sbatchout3} | grep -oh "[1-9][0-9]*$")
		echo 'waiting for:'$jobid3
		# SHOTGUN ANALYSIS
		sbatch --depend=afterok:${jobid3} --job-name OSHOT /mnt/beegfs/agarin/dato-activo/02_ClusterScripts/agarin/ShotgunAnalysisAutomated.sh /mnt/beegfs/agarin/dato-activo/02_ClusterScripts/agarin/ShotgunAnalysisAutomated.R $txtFileName $datasetname 0 $decoyid 3 4 9 30

	done < "${EXP}"



else
	echo "Parameters are missing."
	echo "USAGE: sbatch /mnt/beegfs/agarin/dato-activo/02_ClusterScripts/agarin/OmssaLauncher.sbs currentfolder dbtargetfile dbdecoyfile analysistype decoyid mgfFolder quality"
	echo "USAGE: sbatch /mnt/beegfs/agarin/dato-activo/02_ClusterScripts/agarin/OmssaLauncher.sbs"
fi
