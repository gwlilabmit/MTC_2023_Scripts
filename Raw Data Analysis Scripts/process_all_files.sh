#!/bin/sh
#$ -S/bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ _M mirae@mit.edu
#$ -pe whole_nodes 1
#########################################

#we source many variables from another shell script with needed variables.  Pass it when you call this script.
#do this by running sbatch process_seq.sh ./config.sh (or sub path to the config file if not in same folder).
source $1

for f in ${Raw_Data_Location}*$Original_End; do
	echo "Starting Script for $f ..."
	sbatch fastq_to_dataframe.sh $1 $f
done
