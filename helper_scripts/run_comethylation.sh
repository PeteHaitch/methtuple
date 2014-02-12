#!/bin/bash

#### DESCRIPTION ####
# Peter Hickey (peter.hickey@gmail.com)
# 04/02/2014
# A example bash script for processing a BAM file with comethylation in a pseudo-parallel chromosome-by-chromsome fashion.

#### WARNINGS #####
# WARNING: This is not a robust, general purpose script. Rather it works for my particular setup and suggests how comethylation might be sped up by processing chromosome-level (or some other sub-genome level) units
# WARNING: Script will produce output as subdirectories of OUTDIR
# WARNING: This script runs 1 job per chromosome and each job can use 10-20Gb of memory. Roughly, the larger chromosomes require more memory and more computational time.
# WARNING: The SAMPLE_NAME variable must not contain the "." (period) character
# REQUIRES: comethylation, samtools and Rscript must be in your $PATH
# REQUIRES: tabulate_hist.r to be in the same folder as run_comethylation.sh

#### ALGORITHM ####
# The idea is as follows:
# (1) Split the BAM by chromosome and create chromosome-level BAMs. The BAM must be sorted in coordinate order and indexed
#		(a) if (Single-end data): Continue because there is no need to sort single-end data
# 		(b) if (Paired-end data): Sort the chromosome-level BAMs by queryname order
# (2) For each m-tuple value:
#		(a) Process each (sorted) chromosome-level BAM with comethylation with that m-tuple value. Each chromosome is processed as a background job. 
#		(b) Do not move onto the next m-tuple value until all chromosome jobs finish.
#		(c) Concetanate all chromosome-level m-tuples files

#### Parameters ####
BAM= # The path of the BAM file, e.g. "ex.bam". The BAM must be sorted in coordinate order and indexed. 
OUTDIR= # The path to the directory used for output, e.g. "my_outdir". Multiple subdirectories and files will be created in this folder.
SAMPLE_NAME= # The sample name, which will be the prefix of all output files, e.g. "exciting_sample". WARNING: Must not contain the "." (period) character
M= # A bash array of the "m" in m-tuples, e.g. (1 2 3).
CHROMS= # A bash array of the chromosomes to be processed. It is recommended, but not required, that this array be in karyotypic order e.g. (chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM chrL).
COMETH_PARAMS= # Additional parameters to be passed to comethylation e.g. '--methylationType CG --ignoreDuplicates --minQual 0 --ignoreStart_r1 5 --ignoreStart_r2 $3 --ignoreEnd_r1 10 --ignoreEnd_r2 3 --methylationType CG --overlappingPairedEndFilter XM'.
SEQ_TYPE= # "SE" (single-end data) or "PE" (paired-end data).

#### Step 0 ####
# Check parameters
if [[ -n "${BAM}" && -n "${OUTDIR}" && -n "${SAMPLE_NAME}" && -n "${M}" && -n "${CHROMS}" && -n "${COMETH_PARAMS}" && -n "${SEQ_TYPE}" ]]
then
	if [[ "$SAMPLE_NAME" != "${SAMPLE_NAME//./}" ]]
	then
		echo "Sorry, the SAMPLE_NAME variable must not contain a '.' character"
		exit	
	else
		echo "BAM file = ${BAM}"
		echo "Output directory = ${OUTDIR}"
		echo "Sample name = ${SAMPLE_NAME}"
		echo "'m' in m-tuples = ${M[@]}"
		echo "Chromosomes = ${CHROMS[@]}"
		echo "Additional comethylation parameters = ${COMETH_PARAMS}"
		echo "Sequencing type = ${SEQ_TYPE}"
	fi

else
	echo "Please check all parameters are valid"
	exit
fi

#### Step 1 ####
# Create QS BAM file for each chromosome
if [ ${SEQ_TYPE} == 'PE' ]
then
	echo "Paired-end sequencing"
	echo "Creating chromosome-level BAMs and sorting by queryname"
	for CHROM in ${CHROMS[@]}
	do
		CMD="samtools view -u ${BAM} ${CHROM} | samtools sort -n -f - ${OUTDIR}/${SAMPLE_NAME}_${CHROM}.bam &"
		eval ${CMD}
	done
	wait
elif [ ${SEQ_TYPE} == 'SE'  ] 
then
	echo "Single-end sequencing"
	echo "Creating chromosome-level BAMs"
	for CHROM in ${CHROMS[@]}
	do
		CMD="samtools view -b ${BAM} ${CHROM} > ${OUTDIR}/${SAMPLE_NAME}_${CHROM}.bam &"
		eval ${CMD}
	done
	wait
else
	echo "The sequencing type (SEQ_TYPE parameter) must be one of 'SE' (for single-end sequencing) or 'PE' (paired-end for paired-end sequencing)"
	exit
fi

#### Step 2 ####
# Loop over "m" in m-tuples
for m in ${M[@]}
do
	## mkdir for output
	CMD="mkdir -p ${OUTDIR}/${m}_tuples/log ${OUTDIR}/${m}_tuples/hist"
	echo ${CMD}
	eval ${CMD}

	## Extract methylation loci m-tuples for each chromosome
	for CHROM in ${CHROMS[@]}
	do      
		CMD="comethylation ${COMETH_PARAMS} --mTuple ${m} ${OUTDIR}/${SAMPLE_NAME}_${CHROM}.bam ${OUTDIR}/${m}_tuples/${SAMPLE_NAME}_${CHROM} &> ${OUTDIR}/${m}_tuples/log/${SAMPLE_NAME}_${CHROM}.${m}.log &"
		eval ${CMD}
	done
	# Wait until all chromosomes have finished before moving to next m-tuple
	wait

	## Concatenate all chromosome-level data at the genome level. The only chromosome-level files that are retained are the .hist files.
	# Move all .hist files to a 'hist' directory
	mv ${OUTDIR}/${m}_tuples/*.hist ${OUTDIR}/${m}_tuples/hist
	# Create the genome-level files, i.e. the concatenated files
	EXTENSION=$(basename ${OUTDIR}/${m}_tuples/${SAMPLE_NAME}_${CHROMS[0]}.*.tsv | cut -d "." -f 2-)
	head -n 1 ${OUTDIR}/${m}_tuples/${SAMPLE_NAME}_${CHROMS[@]}.${EXTENSION} > ${OUTDIR}/${m}_tuples/${SAMPLE_NAME}.${EXTENSION} # Add the header
	touch ${OUTDIR}/${m}_tuples/${SAMPLE_NAME}.reads_that_failed_QC.txt
	Rscript tabulate_hist.r ${OUTDIR}/${m}_tuples/hist
	for CHROM in ${CHROMS[@]}
	do      
		cat ${OUTDIR}/${m}_tuples/${SAMPLE_NAME}_${CHROM}.${EXTENSION} | grep -v -P  "^chr\tpos1" \ 
		>> ${OUTDIR}/${m}_tuples/${SAMPLE_NAME}.${EXTENSION} && rm ${OUTDIR}/${m}_tuples/${SAMPLE_NAME}_${CHROM}.${EXTENSION}
		cat ${OUTDIR}/${m}_tuples/${SAMPLE_NAME}_${CHROM}.reads_that_failed_QC.txt \ 
		>> ${OUTDIR}/${m}_tuples/${SAMPLE_NAME}.reads_that_failed_QC.txt && rm ${OUTDIR}/${m}_tuples/${SAMPLE_NAME}_${CHROM}.reads_that_failed_QC.txt
	done
	# The final .tsv and "reads_that_failed_QC.txt" files are gzipped
	gzip ${OUTDIR}/${m}_tuples/${SAMPLE_NAME}.${EXTENSION}
	gzip ${OUTDIR}/${m}_tuples/${SAMPLE_NAME}.reads_that_failed_QC.txt
done

#### Tidy up by removing chromosome-level BAM files ####
for CHROM in ${CHROMS[@]}
	do
		rm ${OUTDIR}/${SAMPLE_NAME}_${CHROM}.bam
	done