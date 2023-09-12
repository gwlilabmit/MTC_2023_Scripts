#!/bin/sh
#########################################

module load python
#Pass the name of the bowtie file we are aligning to this when running.
#Run this script using: sbatch process_bowtie_files.sh path_to_bowtie_file
echo "Trimming PolyA tail from reads before aligning. "
python trim_polyA_and_adapter.py $1 $2
echo "Done trimming!"

