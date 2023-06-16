#!/bin/sh
#$ -S/bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ _M mirae@mit.edu
#$ -pe whole_nodes 1
#########################################

module load python
module load bowtie
module load samtools
module load bedtools

#we source many variables from another shell script with needed variables.  Pass it when you call this script.
#do this by running sbatch process_seq.sh ./config.sh (or sub path to the config file if not in same folder).
source $1

#First we trim the files if that has been requested:
if $Trimming; then
	echo -e "\n\n--------------------------------------------\nTrimming Files"
	for f in ${Raw_Data_Location}*.fastq; do
		echo "Trimming $f ..."
		python /home/mirae/data/Sequencing/Sequencing_Analysis/trim_reads.py $f $Trim_Seq $Trim_Save $Trim_End
	done
else
	Trim_Save=$Raw_Data_Location
	Trim_End=".fastq"
fi

#Next we will align all files to the provided genome
echo -e "\n--------------------------------------------\nAligning Files using ${Genome##*/} Bowtie index."
for f in ${Trim_Save}*$Trim_End; do
	f_no_start=${f##*/}
	Sorted_File=${Sorted_Save}${f_no_start/$Trim_End/$Sorted_End}
	if [ ! -f "$Sorted_File" ]; then
		echo "Aligning $f"
		bowtie --threads 2 --sam $Bowtie_Args $Genome $f | samtools view -u | samtools sort -T temp_sort_file.bam -o $Sorted_File
	fi
done

#We take the aligned and sorted genome and generate a "depth coverage" txt file
echo -e "\n--------------------------------------------\nGenerating Genome Coverage txt files using bedtools."
for f in ${Sorted_Save}*$Sorted_End; do
	f_no_start=${f##*/}
	Plus_Depth_File=${Depth_Save}${f_no_start/$Sorted_End/_plus_depth.txt}
	Minus_Depth_File=${Depth_Save}${f_no_start/$Sorted_End/_minus_depth.txt}
	if [ ! -f "$Plus_Depth_File" ]; then
		echo "Generating depth files for $f"
		bedtools genomecov -d -5 -strand + -ibam $f > $Plus_Depth_File
		bedtools genomecov -d -5 -strand - -ibam $f > $Minus_Depth_File
	fi
done

#Convert those depth files to wigs
echo -e "\n--------------------------------------------\nUsing Depth Files to Make wig files."
for f in ${Depth_Save}*_depth.txt; do
	f_no_start=${f##*/}
	Wig_File=${Wig_Save}${f_no_start/_depth.txt/.wig}
	if [ ! -f "$Wig_File" ]; then
		echo "Converting ${f_no_start/_depth.txt/} depth file into a wig file..."
		python /home/mirae/data/Sequencing/Sequencing_Analysis/density_to_wig.py $f ${Wig_Save}${f_no_start/_depth.txt/.wig}
	fi
done


#If a valid CDS_Dir is provided then make a dataframe file with counts for each gene.
if [ -d "$CDS_Dir" ]; then
	echo -e "\n--------------------------------------------\nUsing CDS Files to Generate Gene Specifc Counts."
	for f in ${Wig_Save}*_minus.wig; do
		f_no_start=${f##*/}
		echo "Making ${DataFrame_Save}${f_no_start/_minus.wig/} dataframe.."
		python /home/mirae/data/Sequencing/Sequencing_Analysis/wig_to_df.py ${f/_minus.wig/} $CDS_Dir $CDS_Files $CDS_Genomes $CDS_Names ${DataFrame_Save}${f_no_start/_minus.wig/}
	done
else
	echo "No valid CDS directory was provided so not generating gene count CDS"
fi

