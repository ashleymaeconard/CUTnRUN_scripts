#!/bin/bash
# get_SRA.sh
# Ashley Mae Conard
# Last Mod: June 26, 2019

# Check to make sure input is correct
if [ $# -ne 6 ]; then
	echo $0: "Usage: ./get_SRA.sh /PATH/TO/SRA_ACCESSION_FILE/ (the .csv file option from NCBI's SRA) PROCESSORS PAIRED_OR_NOT (YES is 1, NO is 0)" 
	exit 1
fi

# Assign input to variable names
INPUT_DIR=$1
RESULTS_DIR=$2
PROCESSORS=$3
PAIRED_OR_NOT=$4
 
# Check if RNA-seq data is paired or not
if [ "$PAIRED_OR_NOT" -ne "1" ]; then
   echo "This script only works for paired-end RNA-seq data processing at the moment.";
   exit;
fi

OUTPUT_DIR=$INPUT_DIR

parallel -j 10 --bar fastq-dump {1} --split-files -O /ltmp/aconard/SingleCell/Wang_2014/ -gzip ::: 
    
START=0
i=$START

echo "GET SRA files from .xls file from NCBI's SRA on $PROCESSORS processors."
# Add commands to a script.txt file
for dir in $INPUT_DIR/*/
    do
    echo ${dir}
    # iterate through input directory
    if [ -d ${dir} ] ; then
    
        # iterate through only the R1 replicates
        for R1 in ${dir}/*R1*
            do
                echo "Getting paired .fastq for $R1"
                # iterate through each R1 to determine when to add 'wait' to script
                ((i = i + 1))
                
                # if number of processors reached, add wait to form a batch that will finish, and then process the next batch
                if (( ${i}%${PROCESSORS}==0 )); then
                    echo "wait" >> $COMMAND_SCRIPT
                fi 
                
                # generate the bowtie2 command for the next R1 and R2
                fileName=$(echo `basename $R1`) # get filename from .fq.gz
                fileParent=$(echo `basename $(dirname $R1)`)

                # create folder for all Bowtie2 outputs per file by removing R1 and R2 (e.g. MTb8-8)
                fName=$(echo ${fileName::$NUM_CHAR_KEEP_NAME}) 
                folderName=${OUTPUT_DIR}${fileParent}"/"${fName}

                # get 2nd read pair
                R2=${R1//R1/R2}

                # output directory for BOWTIE2/sample_name
                mkdir -p ${folderName}

                # write Bowtie2 command to a file, followed by .bam creation and sorting
                echo "Adding ${fileName} to run_Bowtie.txt script"
                echo " "
                echo "(bowtie2 -p ${THREADS} -x /data/compbio/aconard/rnaseq_larschan/BDGP6/Drosophila_melanogaster/Ensembl/BDGP6/Sequence/Bowtie2Index/genome -1 $R1 -2 $R2 --un $folderName/out_un.sam --al $folderName/out_al.sam --un-conc $folderName/out_unconc.sam --met-file $folderName/out_met-file.tsv -S $folderName/out.sam 2> $folderName/alignment_info.txt; samtools view -bS $folderName/out.sam > $folderName/out.bam; rm -rf $folderName/out.sam; samtools sort $folderName/out.bam -o $folderName/out.sorted.bam) &" >> $COMMAND_SCRIPT  
        done
    fi
done

# run command_script (.txt file saved in /results/bowtie2/
bash $COMMAND_SCRIPT
