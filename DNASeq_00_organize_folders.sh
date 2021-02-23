#!/bin/bash

# usage 
# sh 00_organize_folders.sh rootFolder outputFolder extension 
# example executed in /mnt/beegfs/agarin/dato-activo/03_Analysis/apatino/13_Exomes_SmokersExtreme_Ene18/CUN_01/20180108/
# sh /mnt/beegfs/agarin/dato-activo/02_ClusterScripts/agarin/DNASeq_00_organize_folders.sh FASTQ FASTQ fastq.gz

nParameters=$#
if [ $nParameters -eq 3 ] 
then
	folderName=${1%/}
	outputFolderName=${2%/}
	extension=${3%/}
	extPattern="*_1."$extension
	find $folderName -name $extPattern -type f | while read file;
	do
	    f1=$(basename "$file")
	   	# remove the characters after "_1"
	   	FASTQfolderName=$(sed 's/\(.*\)_1.*/\1/' <<< $f1)
	   	# remove the trimming info. remove all the character after the last "_" 
	   	#FASTQfolderName=$(sed 's/\(.*\)_.*/\1/' <<< $FASTQfolderName)
	   	# get the also _2 filename
	   	f2=$(sed 's/_1/_2/' <<< $f1)
	   	mkdir -p $outputFolderName'/'$FASTQfolderName
	    # mv $folderName'/'$f1 $outputFolderName'/'$FASTQfolderName'/'
	    # mv $folderName'/'$f2 $outputFolderName'/'$FASTQfolderName'/'
	    # copy also other files such as fastqcs.
	    input=$folderName'/'
	    pattern=$FASTQfolderName'*'
	    files=$(find "${input}" -maxdepth 1 -iname "${pattern}")
	    echo $pattern
	    echo $outputFolderName'/'$FASTQfolderName'/'
	    echo $files
		mv -t $outputFolderName'/'$FASTQfolderName'/' $files
	    ## begin check code
	    # echo $outputFolderName'/'$FASTQfolderName
	    # echo $folderName'/'$f1 $outputFolderName'/'$FASTQfolderName'/'
	    # echo $folderName'/'$f2 $outputFolderName'/'$FASTQfolderName'/'
	   	## end check code
	   	
	done
	
else
	echo "Parameters are missing. Check the usage inside the script"
fi