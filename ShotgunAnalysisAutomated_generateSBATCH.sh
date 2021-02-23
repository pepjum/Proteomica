#!/bin/bash 

#### USAGE
## sh ShotgunAnalysisAutomated_generateSBATCH.sh txtFullPath databaseFullPath decoy_id proteinAnalysis search_engine minnAA maxnAA [dummyCodeTableFullPath]

pathToFile=$1
database=$2
decoy_id=$3
proteinAnalysis=$4 # 0: psmFDR, 1: protFDRMAYU, 2: pepFDR, 3: protFDR
search_engine=$5 # 1: mascot, 2: comet, 3:tandem, 4: omssa
minnAA=$6
maxnAA=$7
dummyCodeTable=$8
nParameters=$#

if [ $nParameters -gt 1 ] 
then
	datasets=($(awk '{print $3}' $pathToFile))
	unique_datasets=($(echo "${datasets[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
	echo 'path to file '$pathToFile
	echo 'database '$database
	echo 'decoy_id '$decoy_id
	if [ $search_engine -eq 1 ] 
	then
		echo 'search_engine mascot'
	elif [ $search_engine -eq 2 ] 
	then
		echo 'search_engine comet'
	elif [ $search_engine -eq 3 ] 
	then
		echo 'search_engine xtandem'
	elif [ $search_engine -eq 4 ] 
	then
		echo 'search_engine omssa'
	fi
	if [ $nParameters -eq 8 ] 
	then
		echo 'dummycodetable '$dummyCodeTable
	fi

	if [ $proteinAnalysis -eq 0 ]	
	then
		database=0		
	fi
	for i in ${!unique_datasets[@]};
	  do 	
			if [ $nParameters -eq 7 ] 
			then			
				echo "sbatch /mnt/beegfs/agarin/dato-activo/02_ClusterScripts/agarin/ShotgunAnalysisAutomated.sh /mnt/beegfs/agarin/dato-activo/02_ClusterScripts/agarin/ShotgunAnalysisAutomated.R $pathToFile ${unique_datasets[i]} $database $decoy_id $proteinAnalysis $search_engine $minnAA $maxnAA "
				

			elif [ $nParameters -eq 8 ] 
			then
				echo "sbatch /mnt/beegfs/agarin/dato-activo/02_ClusterScripts/agarin/ShotgunAnalysisAutomated.sh /mnt/beegfs/agarin/dato-activo/02_ClusterScripts/agarin/ShotgunAnalysisAutomated.R $pathToFile ${unique_datasets[i]} $database $decoy_id $proteinAnalysis $search_engine $minnAA $maxnAA $dummyCodeTable"
			 else
				echo "Check the parameters"
				echo "USAGE sh /mnt/beegfs/agarin/dato-activo/02_ClusterScripts/agarin/ShotgunAnalysisAutomated_launcher.sh txtFullPath databaseFullPath decoy_id proteinAnalysis search_engine minAA maxAA [dummyCodeTableFullPath]"
			 fi
	done  > /mnt/beegfs/agarin/dato-activo/02_ClusterScripts/agarin/ShotgunAnalysisAutomated_launcher.cmd
	
else
	echo "Parameters are missing. Check them"
	echo "USAGE sh /mnt/beegfs/agarin/dato-activo/02_ClusterScripts/agarin/ShotgunAnalysisAutomated_generateSBATCH.sh txtFullPath databaseFullPath decoy_id proteinAnalysis minAA maxAA search_engine [dummyCodeTableFullPath]"
fi

