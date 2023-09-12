## Publication Raw Data Processing. 

This folder contains all of the scripts used to process fastq sequencing files used in this publication.  We did all of our analysis on a SLURM-run server.  Each experiment has its own "config.sh" file which contains all of the exerpiment specific commands.  To run the analysis, submit a job by running `sbatch process_all_files.sh __config.sh`.  To do this yourself, please modify the config files to point to your local version of the raw sequencing data. Several modules are required by the analysis pipeline.  These will be automatically loaded by the jobs, but please be sure your environment has them.  These are:
- python 
- samtools
- bowtie (not bowtie2)
- cutadapt
- bedtools
- seqtk

