dataset=$1 		#/Tandem_files/datasetquesea

echo 'changing target .dat files to .pin format'

EXP=$dataset'/inputFilesTargetDat.log'
find $dataset -mindepth 1 -maxdepth 1 -type f -iname "*.dat"  -not -name "*_corrected.tsv" -not -name "*_peptides.tsv" -not -name "*_out.tsv" -not -name "*_peptides_log.log" > $EXP  #descomentar luego

cd /home/nostromo/data/pepe/scripts/MASCOT_PERCOLATOR/

while IFS='' read -r line || [[ -n "$line" ]]; do
  file=$line
  #echo $file "files"
  DIR=$(dirname "${file}")
  #echo $DIR "DIR"
  namefile=$(basename "${file}")
  namefilefilenoext="${namefile%.*}"
  Rscript PARSING_MASCOT_FILES_TO_PERCOLATOR.R $line
  #python /home/nostromo/data/pepe/scripts/Reoplace.py $DIR

done < "${EXP}"
wait

echo 'Merge target files in one .pin file'

Rscript MERGE_FILES_MASCOT_PERCOLATOR_INPUT_JUL18.R $dataset'/'


echo 'changing decoy .dat files to .pin format'


EXP2=$dataset'/inputFilesDecoyDat.log'
find $dataset'-D/' -mindepth 1 -maxdepth 1 -type f -iname "*.dat"  -not -name "*_corrected.tsv" -not -name "*_peptides.tsv" -not -name "*_out.tsv" -not -name "*_peptides_log.log" > $EXP2  #descomentar luego


while IFS='' read -r line || [[ -n "$line" ]]; do
  file=$line
  #echo $file "files"
  DIR=$(dirname "${file}")
  #echo $DIR "DIR"
  namefile=$(basename "${file}")
  namefilefilenoext="${namefile%.*}"
  Rscript PARSING_MASCOT_FILES_TO_PERCOLATOR.R $line
  #python /home/nostromo/data/pepe/scripts/Reoplace.py $DIR

done < "${EXP2}"
wait

echo 'Merge decoy files in one .pin file'

Rscript MERGE_FILES_MASCOT_PERCOLATOR_DECOYINPUT_JUL18.R $dataset'-D/'

 percolator_dir=$dataset'/percolator/'

echo 'Creating a one file .pin input'

MascotTargetFile=$percolator_dir'MASCOT_TARGET.pin'
MascotDecoyFile=$percolator_dir'MASCOT_DECOY.pin'

Rscript BIND_TYD_MASCOT_PERCOLATOR_JUL18.R $MascotTargetFile $MascotDecoyFile

 inputfilePerc=$percolator_dir'ALL_MASCOT.pin'

echo 'launching percolator'

percolator -A -Y -P DECOY_ -r $inputfilePerc'results_peptides.tsv.out' -l $inputfilePerc'results_proteins.tsv.out' -J $inputfilePerc'.out' -m $inputfilePerc'results_psm.out' --override -c -q -d 2 -T 0.2 -i 30 --spectral-counting-fdr 0.01 $inputfilePerc
