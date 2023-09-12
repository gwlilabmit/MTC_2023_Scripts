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
raw_file=$2
f_name_start=${raw_file##*/}
bowtie_name=${Sorted_Save}${f_name_start/$Original_End/$Sorted_End}
plus_depth_file_name=${Depth_Save}${f_name_start/$Original_End/_plus_depth.txt}
minus_depth_file_name=${Depth_Save}${f_name_start/$Original_End/_minus_depth.txt}
plus_wig_name=${Wig_Save}${f_name_start/$Original_End/_plus.wig}
minus_wig_name=${Wig_Save}${f_name_start/$Original_End/_minus.wig}
dataframe_file=${DataFrame_Save}${f_name_start/$Original_End/}

#Next we will align all files to the provided genome
empty_str=""

if [ "$Preprocess_File" != "$empty_str" ]; then
	trim_name=${Trim_Save}${f_name_start/$Original_End/$Trim_End}
	if [ ! -f "$trim_name" ]; then
		filename=$(readlink -f $Preprocess_File)
		srun $filename $raw_file $trim_name
	else
		echo "$trim_name already exists"
	fi
else
	trim_name=$raw_file
fi


if [ ! -f "$bowtie_name" ]; then
	echo "Aligning $trim_name to $Genome"
	bowtie --threads 2 --sam $Bowtie_Args $Genome $trim_name | samtools view -u | samtools sort -T ${f_name_start/$Original_End/_sort_file.bam} -o $bowtie_name
fi

echo "Generating depth files for ${f_name_start/$Original_End/}"
if [ ! -f "$plus_depth_file_name" ]; then
	bedtools genomecov -d -5 -strand + -ibam $bowtie_name > $plus_depth_file_name
fi
if [ ! -f "$minus_depth_file_name" ]; then
	bedtools genomecov -d -5 -strand - -ibam $bowtie_name > $minus_depth_file_name
fi

echo "Converting depth files to wigs for ${f_name_start/$Original_End/}"
if [ ! -f "$plus_wig_name" ]; then
	python /home/mirae/data/Publication_Raw_Data_Processing/density_to_wig.py $plus_depth_file_name $plus_wig_name
fi
if [ ! -f "$minus_wig_name" ]; then
	python /home/mirae/data/Publication_Raw_Data_Processing/density_to_wig.py $minus_depth_file_name $minus_wig_name
fi


#If a valid CDS_Dir is provided then make a dataframe file with counts for each gene.
if [ -d "$CDS_Dir" ]; then
	echo -e "\n--------------------------------------------\nUsing CDS Files to Generate Gene Specifc Counts."
	echo "Making ${DataFrame_Save}${f_no_start/_minus.wig/} dataframe.."
	if [ ! -f "$dataframe_file" ]; then
		python /home/mirae/data/Publication_Raw_Data_Processing/wig_to_df.py ${minus_wig_name/_minus.wig/} $CDS_Dir $CDS_Files $CDS_Genomes $CDS_Names $dataframe_file
	fi
else
	echo "No valid CDS directory was provided so not generating gene count CDS"
fi



