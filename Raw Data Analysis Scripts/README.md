# Raw data analysis
### Going from fastq files to wigs and gene-count dataframes.

The raw data processing for this project is standard for the analysis done in the Li lab and follows the following steps:

- Trimming of each read as needed to remove any adapter sequence or nucleotides that may have been added as part of the library preparation. This is done via a custom python script `trim_reads.py`.
- Alignment to the bowtie index of choice. This is done using [bowtie](https://bowtie-bio.sourceforge.net/index.shtml) (not bowtie2) and provided indices.
- Sorting of aligned output (sam file), and compression to BAM file using `samtools`.
- Generate a "depth file" using the 5' end of each read as a read count at that location. (here we also seperate reads from the + and - strands). This is down using the `bedtools genomecov` command.
- Conversion of density files to wigs (the default format used in our lab to view sequencing results). This is done through a custom python script `density_to_wig.py`.
- Conversion of wigs to read/gene dataframe files.  This is done through a custom python file that requires a CDS file for the genome annotation that the reads were aligned to: `wig_to_df.py`.

* The `analysis_helper_functions.pyc` file just contains some helper files for some of the other python files. 

If you wish to do this yourself all you need to do is run the bash script `process_seq.sh` and provide a relevant `config.sh` file as the argument to set the relevant parameters.  You will likely need to modify the `config.sh` file to reflect the locations of files to analyze on your machine (and where you would like intermediate files to be stored).

Also provided here are the bowtie indices and the cds files that were used in our analysis.
