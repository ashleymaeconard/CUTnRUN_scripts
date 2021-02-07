#!/bin/bash
# create_avg_prof.sh
# Ashley Mae Conard
# Last Mod: Dec 29, 2019
# Resources: https://www.biostars.org/p/42844/
#            http://seqanswers.com/forums/archive/index.php/t-45817.html (remove chr from file and remove header)

# Check to make sure input is correct
if [ $# -ne 4 ]; then
	echo $0: "Usage: ./create_avg_prof.sh /PATH/TO/SORTED.BAMS/FOLDER(s)/ (containing folders of .sorted.bam files, assuming at most 2 subdirectories (e.g. only write up to /bowtie2/ here: .../results/bowtie2/clamp_sg_fm/MR_1-1/out.sorted.bam ) PROCESSORS MFOLD_LOWER MFOLD_UPPER (default is 5 50) \n" # NOTE: MAKE SURE source ~/venv/bin/activate FOR PROPER PYTHON ENVIRONMENT.
	exit 1
fi

# Assign input to variable names
INPUT_DIR=$1
TF=$2
# MFOLD_LOWER=$3
# MFOLD_UPPER=$4

samtools sort ${TF}_control.bam -o ${TF}_control.sorted.bam &
samtools sort ${TF}_1.bam -o ${TF}_1.sorted.bam & 
samtools sort ${TF}_2.bam -o ${TF}_2.sorted.bam & 
samtools sort ${TF}_3.bam -o ${TF}_3.sorted.bam & 

samtools index ${TF}_control.sorted.bam & 
samtools index ${TF}_1.sorted.bam &
samtools index ${TF}_2.sorted.bam & 
samtools index ${TF}_3.sorted.bam & 

bamCompare -b1 ${TF}_1.sorted.bam -b2 ${TF}_control.sorted.bam -o ${TF}_1_log2ratio.bigWig; bamCompare -b1 ${TF}_2.sorted.bam -b2 ${TF}_control.sorted.bam -o ${TF}_2_log2ratio.bigWig; bamCompare -b1 ${TF}_3.sorted.bam -b2 ${TF}_control.sorted.bam -o ${TF}_3_log2ratio.bigWig

wiggletools mean ${TF}_1_log2ratio.bigWig ${TF}_2_log2ratio.bigWig ${TF}_3_log2ratio.bigWig | wigToBigWig stdin /data/compbio/aconard/BDGP6/dm6.chrom.sizes ${TF}_avg_allReps_log2ratio.bigWig
# computeMatrix - avg_allReps_log2ratio
computeMatrix scale-regions -S ${TF}_avg_allReps_log2ratio.bigWig -R /data/compbio/aconard/BDGP6/cluster.*.bed --binSize 250 --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --regionBodyLength 5000 -o matrix.genes.clusters.avg_allReps_log2ratio.mat.gz --skipZeros --smartLabels --sortRegions descend
# plotHeatmap
plotHeatmap -m matrix.genes.clusters.avg_allReps_log2ratio.mat.gz -out heatmap_genes.clusters.avg_allReps_log2ratio_se_${TF}.png --heatmapHeight 25 --heatmapWidth 15 --labelRotation 45 --missingDataColor red --plotType=se

plotProfile -m matrix.genes.clusters.avg_allReps.mat.gz -out profile_genes.clusters.avg_allReps_se_${TF}.png --labelRotation 45 --plotType=se --colors salmon olivedrab turquoise dodgerblue goldenrod --plotHeight 10

#######################
OR WITHOUT LOG2RATIO

wiggletools mean ${TF}_1.bigWig ${TF}_2.bigWig ${TF}_3.bigWig | wigToBigWig stdin /data/compbio/aconard/BDGP6/dm6.chrom.sizes ${TF}_avg_allReps.bigWig
# computeMatrix - avg_allReps
computeMatrix scale-regions -S ${TF}_avg_allReps.bigWig -R /data/compbio/aconard/BDGP6/cluster.*.bed --binSize 250 --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --regionBodyLength 5000 -o matrix.genes.clusters.avg_allReps.mat.gz --skipZeros --smartLabels --sortRegions descend
# plotHeatmap
plotHeatmap -m matrix.genes.clusters.avg_allReps.mat.gz -out heatmap_genes.clusters.avg_allReps_se_${TF}.png --heatmapHeight 25 --heatmapWidth 15 --labelRotation 45 --missingDataColor red --plotType=se


# create a bdg without header 
cp $bdg $bdgSorted

# Sort the bedgraph file
# bedSort out.sorted.bam_treat_pileup_noheader.bdg out.sorted.bam_treat_pileup.sorted.bdg

# bedGraphToBigWig
# bedGraphToBigWig out.sorted.bam_treat_pileup.sorted.bdg ../../../BDGP6/bdgp6.chrom.size out.sorted.bam_summits.bigWig

# computeMatrix - WHOLE GENE 1,2,3
#computeMatrix scale-regions -S cwo_1.bigWig cwo_2.bigWig cwo_3.bigWig -R ../../../BDGP6/cluster.*.bed --binSize 250 --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --regionBodyLength 5000 -o matrix.genes.clusters.reps123.mat.gz --skipZeros --smartLabels --sortRegions descend
# plotHeatmap
#plotHeatmap -m matrix.genes.clusters.reps123.mat.gz -out heatmap_genes.clusters.reps123_se_cwo.png --heatmapHeight 25 --heatmapWidth 15 --labelRotation 45 --missingDataColor red --plotType=se

# computeMatrix - TSS 1,2,3 and merged
# (base) micmacs /data/compbio/aconard/test_avg_prof/perturbedTFs/cwo $ computeMatrix reference-point -S cwo*.bigWig -R ../../../BDGP6/cluster.*.bed --binSize 250 -o matrix.ref-point.clusters.reps123nMerged.mat.gz --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros --smartLabels --sortRegions descend
# plotHeatmap TSS 1,2,3 and merged
# (base) micmacs /data/compbio/aconard/test_avg_prof/perturbedTFs/cwo $ plotHeatmap -m matrix.ref-point.clusters.reps123nMerged.mat.gz -out heatmap_ref-point.clusters.reps123nMerged_se_cwo.png --heatmapHeight 25 --heatmapWidth 15 --labelRotation 45 --missingDataColor red --plotType=se

# computeMatrix - TSS 1,2,3
# computeMatrix reference-point -S cwo_1.bigWig cwo_2.bigWig cwo_3.bigWig -R ../../../BDGP6/cluster.*.bed --binSize 250 -o matrix.ref-point.clusters.reps123.mat.gz --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros --smartLabels --sortRegions descend
# plotHeatmap TSS 1,2,3
# (base) micmacs /data/compbio/aconard/test_avg_prof/perturbedTFs/cwo $ plotHeatmap -m matrix.ref-point.clusters.reps123.mat.gz -out heatmap_ref-point.clusters.reps123_se_cwo.png --heatmapHeight 25 --heatmapWidth 15 --labelRotation 45 --missingDataColor red --plotType=se

###########################
ls *.sam | parallel -j4 -k bash convert2bam.sh {}

sample=$1
describer=$(echo ${sample} | sed 's/.sam//')  
   
# Convert file from SAM to BAM format  
samtools view -b $sample > ${describer}.uns.bam  
   
# Sort BAM file  
samtools sort ${describer}.uns.bam ${describer}   
   
# index the bam file  
samtools index ${describer}.bam  
   
# Revove intermediate files  
rm ${describer}.uns.bam  
 
###########################

# plotProfile - with confidence intervals
# plotProfile -m matrix.clusters.mat.gz -out exampleprofilecwo_bin200.png --plotType=se --labelRotation 45 --plotHeight 15 --plotWidth 15


bedGraph=$(echo) 
# RESULTS_DIR=$(echo `dirname $INPUT_DIR`)

# # Check to see if results/bowtie2/ directory exists
# if [ -d "$RESULTS_DIR/macs2/" ]; then
# 	echo "${RESULTS_DIR}/macs2/ already exists."
#     OUTPUT_DIR=${RESULTS_DIR}"/macs2/"
# else
#     OUTPUT_DIR=${RESULTS_DIR}"/macs2/"
#     mkdir $OUTPUT_DIR
# fi

# # Create command script for Bowtie2
# COMMAND_SCRIPT=${OUTPUT_DIR}"/run_macs2.txt" 
# echo $COMMAND_SCRIPT

# # Remove command script if already made
# if [ -f ${COMMAND_SCRIPT} ]; then
# 	echo -e "REPLACING ${COMMAND_SCRIPT}\n"
# 	rm -rf ${COMMAND_SCRIPT}
# fi

# START=0
# i=$START

# echo "Processing .sorted.bam files with MACS2 on $PROCESSORS processors."

# # Add commands to a script.txt file
# for dir in $INPUT_DIR/*/
#     do
#     echo ${dir}
#     # iterate through input directory (with sample names)
#     if [ -d ${dir} ] ; then
#         declare -a my_array
#         unset my_array
        
#         # iterate through only replicates to add them to list to be merged
#         for file in ${dir}/*/*.sorted.bam
#             do       
   
#             # iterate through each sample to determine when to add 'wait' to script
#             ((i = i + 1))

#             # if number of processors reached, add wait to form a batch that will finish, and then process the next batch
#             if (( ${i}%${PROCESSORS}==0 )); then
#                 echo "wait" >> $COMMAND_SCRIPT
#             fi 

#             # generate folder structure 
#             fileName=$(echo `basename $file`) # get filename from .fq.gz
#             fileReplicate=$(echo `basename $(dirname $file)`)
#             fileSample=$(echo `basename $(dirname $(dirname $file))`)

#             # create folder for all outputs per file by removing replicate identifer
#             folderName=${OUTPUT_DIR}${fileSample}/${fileReplicate}/

#             # output directory for BOWTIE2/sample_name
#             mkdir -p ${folderName}

#             # write Bowtie2 command to a file, followed by .bam creation and sorting
#             #echo "Adding ${my_array[@]} to run_macs2.txt script."
#             echo " "
#             # ${my_array[@]} instead of $file
#             echo "(macs2 callpeak -t $file -B -f AUTO --nomodel --SPMR --keep-dup all -g dm --trackline -n $fileName --cutoff-analysis --call-summits -p 0.01 -m ${MFOLD_LOWER} ${MFOLD_UPPER} --outdir $folderName 2> $folderName/$fileName.log) &" >> $COMMAND_SCRIPT  
#         done
#     fi
# done

# # run command_script (.txt file saved in /results/bowtie2/
# echo "RUNNING ${COMMAND_SCRIPT}"
# bash ${COMMAND_SCRIPT}

