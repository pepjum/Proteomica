#!/bin/bash

CURRENTFOLDER=$1	# /mnt/beegfs/agarin/dato-activo/03_Analysis/agarin/28_PFortes_Shotgun_lncRNA_Feb18/
DATABASENAME=$2		# nextProtDB20170801

Rscript /home/nostromo/data/pepe/scripts/ShotgunPostAnalysisAutomated_nostromo.R $CURRENTFOLDER $DATABASENAME
