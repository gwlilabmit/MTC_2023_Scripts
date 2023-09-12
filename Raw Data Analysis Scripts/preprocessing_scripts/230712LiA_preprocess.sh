#!/bin/sh
#########################################

module load seqtk
#Pass the name of the bowtie file we are aligning to this when running.
#Run this script using: sbatch process_bowtie_files.sh path_to_bowtie_file
file=$1
new_file=$2
echo "Trimming fastq files based on quality scores. "
seqtk trimfq $1 > $2
echo "Done trimming!"

