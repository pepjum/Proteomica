#!/bin/bash

while getopts "f:" OPTION
do
	case $OPTION in
	f) sdir=$OPTARG
    ;;
	?)
	usage
	exit
	esac
done
cd $sdir
EXP=$sdir'files_to_best_psm_corrected.txt'

 #find . -type f -name "*batch*" -print > $EXP
#sed "s/[./]//" $sdir'files_to_best_psm.txt' > $sdir'files_to_best_psm_corrected.txt'

#EXP=$sdir'files_to_best_psm_corrected.txt'
while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line

	var2=$(echo $line | cut -f2 -d /)
	printf "$var2\n"
    Rscript /home/nostromo/data/pepe/scripts/best_PSM.R $var2 &

done < "${EXP}"
