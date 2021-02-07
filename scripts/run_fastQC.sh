#!/bin/bash
# run_fastQC.sh
# Ashley Mae Conard
# Last Mod: Nov 5, 2018
# Runs FastQC in /research/compbio/users/aconard/FastQC for all .fastq.gz files in a given directory

# Check to make sure input is correct
if [ $# -ne 3 ]; then
	echo $0: Usage: ./run_fastQC.sh /PATH/TO/DIR/ \(containing fastq.gz files\) /PATH/TO/OUTPUT_DIR/ \(suggest placing in /results/fastQC\) NUM_THREADS
	exit 1
fi

mkdir -p $2

# Move into directory and create fastQC script
cd $1
find $PWD > file.txt # grab all .fastq.gz files in DIR
sed -i "1d" file.txt # remove dir itself from file.txt
sed -i "$ d" file.txt # remove "file.txt" from file.txt
sed -i "1 i\--outdir=$2 :::" file.txt
sed -i '1 i\/research/compbio/users/aconard/FastQC/fastqc {1}\n' file.txt  
sed -i "1 i\parallel -j $3" file.txt

tr '\n' ' ' < file.txt > commands.txt # replace new line with space
rm -rf file.txt
bash commands.txt
rm -rf commands.txt
