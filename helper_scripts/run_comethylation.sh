#!/usr/bin/env bash

#### DESCRIPTION ####
# Peter Hickey (peter.hickey@gmail.com)
# 04/02/2014
# A example bash script for processing a BAM file with comethylation in a pseudo-parallel chromosome-by-chromsome fashion.
# This script makes extensive use of GNU parallel [O. Tange (2011): GNU Parallel - The Command-Line Power Tool, ;login: The USENIX Magazine, February 2011:42-47.; https://www.gnu.org/software/parallel/]. 

#### TODOs ####
# TODO: Allow for "." characters in $SAMPLE_NAME (instead, don't allow "/" characters)

#### LICENSE ####
## Copyright (C) 2012 - 2014 Peter Hickey (peter.hickey@gmail.com)

## This file is part of comethylation.

## comethylation is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
## comethylation is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with comethylation.  If not, see <http://www.gnu.org/licenses/>.
#### REQUIREMENTS ####
# REQUIRES: comethylation, SAMtools, GNU parallel and Rscript must be in your $PATH
# REQUIRES: tabulate_hist.r to be in the same directory as run_comethylation.sh

#### WARNINGS #####
# WARNING: Script will produce output as subdirectories of $OUTDIR
# WARNING: Roughly, this script runs 1 job per chromosome and each job can use 1-10Gb of memory. Larger chromosomes and larger values of "m" require more memory and more computational time.
# WARNING: The $SAMPLE_NAME variable must not contain the "." (period) character
# WARNING: You may wish to increase the memory available to SAMtools sort via the -m flag, which will decrease the number of intermediate files (requires a recent-ish version of SAMtools).

#### Required variables ####
BAM= # The path of the BAM file, e.g. data/cs_pe_directional.bam. The BAM must be sorted in coordinate order and indexed.
OUTDIR= # The path to the directory used for output, e.g. my_outdir. Multiple subdirectories and files will be created in this folder.
SAMPLE_NAME= # The sample name, which will be the prefix of all output files, e.g. exciting_sample. WARNING: Must not contain the "." (period) character
M= # A bash array of the "m" in m-tuples, e.g. (1 2 3).
CHROMS= # A bash array of the chromosomes to be processed. It is recommended, but not required, that this array be in karyotypic order e.g. (chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM chrL).
COMETH_PARAMS= # Additional parameters to be passed to comethylation e.g. "--methylationType CG --ignoreDuplicates --minQual 0 --ignoreStart_r1 5 --ignoreStart_r2 $3 --ignoreEnd_r1 10 --ignoreEnd_r2 3 --methylationType CG --overlappingPairedEndFilter XM".
SEQ_TYPE= # SE (single-end data) or PE (paired-end data).
NUM_CORES= # The number of cores to be used by GNU parallel, e.g. 8
TABULATE_HIST= # The path to the helper script "tabulate_hist.R"

#### Algorithm ####
# The idea is as follows:
# (1) Split the BAM by chromosome and create chromosome-level BAM files. The original BAM must be sorted in coordinate order and indexed.
#               (a) if (Single-end data): Continue because there is no need to sort single-end data.
#               (b) if (Paired-end data): Sort the chromosome-level BAMs by queryname order.
# (2) For each m-tuple value:
#               (a) Process each chromosome-level BAM with comethylation.
#               (c) Concetanate all chromosome-level m-tuples files.
# (3) Clean up by removing chromosome-level BAM files

#### Step 0 ####
# Check variables
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
# Create chromosome-level BAM files (queryname sorted if paired-end sequencing)
# NB: {.} takes the value of {} without the file extension in parallel commands
if [ ${SEQ_TYPE} == 'PE' ]
then
        echo "Paired-end sequencing"
        echo "Creating chromosome-level BAMs and sorting by queryname"
        parallel --joblog create_chromosome_level_BAMs.log -j ${NUM_CORES} "samtools view -u ${BAM} {} | samtools sort -m 10G -n - ${OUTDIR}/${SAMPLE_NAME}_{.}" ::: ${CHROMS[@]} 
elif [ ${SEQ_TYPE} == 'SE'  ] 
then
        echo "Single-end sequencing"
        echo "Creating chromosome-level BAMs"
        parallel --joblog create_chromosome_level_BAMs.log -j ${NUM_CORES} "samtools view -u ${BAM} {} > ${OUTDIR}/${SAMPLE_NAME}_{.}.bam" ::: ${CHROMS[@]}
else
        echo "The sequencing type (SEQ_TYPE parameter) must be one of 'SE' (for single-end sequencing) or 'PE' (paired-end for paired-end sequencing)"
        exit
fi

#### Step 2 ####
## mkdir for output
parallel -j ${NUM_CORES} "mkdir -p ${OUTDIR}/{}_tuples/log ${OUTDIR}/{}_tuples/hist" ::: ${M[@]}

## Extract methylation loci m-tuples for each chromosome
echo "Extracting methylation loci m-tuples for each chromosome..."
parallel --joblog comethylation.log -j ${NUM_CORES} "comethylation ${COMETH_PARAMS} --mTuple {1} ${OUTDIR}/${SAMPLE_NAME}_{2}.bam ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}_{2} &> ${OUTDIR}/{1}_tuples/log/${SAMPLE_NAME}_{1}.{2}.log" ::: ${M[@]} ::: ${CHROMS[@]}

## Concatenate all chromosome-level data at the genome level. The only chromosome-level files that are retained are the .hist files.
echo "Concatenating chromosome-level files into genome-level files..."
# Move all .hist files to a 'hist' directory
parallel -j ${NUM_CORES} "mv ${OUTDIR}/{}_tuples/*.hist ${OUTDIR}/{}_tuples/hist" ::: ${M[@]}
# Create the header for the genome-wide .tsv file (for given m)
parallel -j ${NUM_CORES} "EXTENSION=\$(basename ${OUTDIR}/{}_tuples/${SAMPLE_NAME}_${CHROMS[0]}.*.tsv | cut -d '.' -f 2-);
head -n 1 ${OUTDIR}/{}_tuples/${SAMPLE_NAME}_${CHROMS[0]}.\${EXTENSION} > ${OUTDIR}/{}_tuples/${SAMPLE_NAME}.\${EXTENSION}" ::: ${M[@]} # Note the need to escape the $ character for variables defined within the command executed by command to avoid the shell interpreting them 
# Touch the file storing all reads_that_failed_QC (for a given m)
parallel -j ${NUM_CORES} "touch ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}.reads_that_failed_QC.txt;
Rscript ${TABULATE_HIST} ${OUTDIR}/{1}_tuples/hist" ::: ${M[@]}
# Concatenate the .tsv files and gzip (for given m)
parallel --joblog cat_tsvs.log -j ${NUM_CORES} "EXTENSION=$(basename ~/methylation_m-tuples/${SAMPLE_NAME}/{1}/${SAMPLE_NAME}_${CHROMS[0]}.*.tsv | cut -d '.' -f 2-);
	for CHROM in ${CHROMS[@]}
	do
		cat ~/methylation_m-tuples/${SAMPLE_NAME}/{1}/${SAMPLE_NAME}_${CHROM}.${EXTENSION} | grep -v -P  '^chr\tpos1'
	done | gzip -c - >> ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}.\${EXTENSION}" ::: ${M[@]}
# Concatenate the reads_that_failed_QC files and gzip (for given m)
parallel --joblog cat_reads_that_failed_QC_files.log -j ${NUM_CORES} "
	for CHROM in ${CHROMS[@]}
	do
		cat ~/methylation_m-tuples/${SAMPLE_NAME}/{1}/${SAMPLE_NAME}_${CHROM}.reads_that_failed_QC.txt
	done | gzip -c - >> ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}.reads_that_failed_QC.txt" ::: ${M[@]}

#### Step 3 ####
## Tidy up by removing chromosome-level BAM files
echo "Removing chromosome-level BAMs..."
parallel -j ${NUM_CORES} "rm ${OUTDIR}/${SAMPLE_NAME}_{1}.bam" ::: ${CHROMS[@]}