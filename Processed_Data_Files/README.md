# Processed Data Files

### What you will find here:

- wig files for each experiment.  After aligning to the correct genome of choice we generate wig files for each experiment.  These are tab seperated files that consist of a header indicating what genome was used for bowtie alignment of the reads, and then two columns of data, the first indicating the genome position and the second indicating the number of reads whose 5' end aligned to that position. Each experiment will have 2 wig files for each of the genomes it was aligned to, one for each strand.
- dataframe files.  These files take the data in each wig files and compare it to a CDS file to generate the number of reads found in each coding body sequence for each gene. They are written as simple tab deliminated text files with the gene name, read count, and several pieces of information about each gene.  We like to read in this file and work with it using pandas.
