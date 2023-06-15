# Analysis Scripts for Molecular Time Capsule Project

### What you will find here:

- Copies of the scripts/exact parameters I used to convert my raw sequencing reads to processed read/gene dataframe files used in downstream analysis.
- The bowtie indices I aligned my raw data with and their corresponding CDS files.
- The dataframe read/gene files for each of the experiments referenced in my paper.
- The wig files for each of the experiments referenced in the paper.
- Jupyter notebook files for each figure in my paper which transform the data in the read/gene files into the exact figures seen in the paper.

### More Details:

<details>
<summary>Raw data Processing:</summary>
<br>

The raw data processing consists of the following steps:
- Trimming of each read as needed to remove any adapter sequence or nucleotides that may have been added as part of the library preparation. This is done via a custom python script `trim_reads.py`.
- Alignment to the bowtie index of choice. This is done using bowtie (not bowtie2) and provided indices.
- Sorting of aligned reads, and compression to BAM file using samtools.
- Generate a "depth file" using the 5' end of each read as a read count at that location. (here we also seperate reads from the + and - strands). This is down using bedtools.
- Conversion of density files to wigs (the default format used in our lab to view sequencing results). This is done through a custom python script `density_to_wig.py`.
- Conversion of wigs to read/gene dataframe files.  This is done through a custom python file that requires a CDS file for the genome annotation that the reads were aligned to: `wig_to_df.py`.

All of these steps are collected in a single bash shell program called `process_seq.sh`.  This program takes a single argument from the command line: - another shell file (denoted as a `config.sh` file) which contains all the experiment specific parameters.  For each experiment in this paper I have a seperate `config.sh` file available with the exact parameters used for the pulished analysis. If you wish to redo this analysis yourself you need only modify the relevant `config.sh` file to update the parameters for the relevant location of the raw reads and (optionally) where you want the processed reads and intermediates to be stored.

In order to run this analysis you will require the following:

- python
  - pandas
- samtools
- bedtools
- bowtie

</details>

<details>
<summary>Raw data Processing:</summary>
<br>

The raw data processing consists of the following steps:
- Trimming of each read as needed to remove any adapter sequence or nucleotides that may have been added as part of the library preparation. This is done via a custom python script `trim_reads.py`.
- Alignment to the bowtie index of choice. This is done using bowtie (not bowtie2) and provided indices.
- Sorting of aligned reads, and compression to BAM file using samtools.
- Generate a "depth file" using the 5' end of each read as a read count at that location. (here we also seperate reads from the + and - strands). This is down using bedtools.
- Conversion of density files to wigs (the default format used in our lab to view sequencing results). This is done through a custom python script `density_to_wig.py`.
- Conversion of wigs to read/gene dataframe files.  This is done through a custom python file that requires a CDS file for the genome annotation that the reads were aligned to: `wig_to_df.py`.

All of these steps are collected in a single bash shell program called `process_seq.sh`.  This program takes a single argument from the command line: - another shell file (denoted as a `config.sh` file) which contains all the experiment specific parameters.  For each experiment in this paper I have a seperate `config.sh` file available with the exact parameters used for the pulished analysis. If you wish to redo this analysis yourself you need only modify the relevant `config.sh` file to update the parameters for the relevant location of the raw reads and (optionally) where you want the processed reads and intermediates to be stored.

In order to run this analysis you will require the following:

- python
  - pandas
- samtools
- bedtools
- bowtie

</details>
<br>

<details>
<summary> Figure Analysis </summary>
<br>

Each figure or subfigure has its own folder which contains:
- The final version of each figure that was included in the paper.  Where possible there will be both the png of the file that was included in the paper, and an interactive version of the figure as well.
- A jupyter notebook which converts the processed data into the figure with annotations.

</details>