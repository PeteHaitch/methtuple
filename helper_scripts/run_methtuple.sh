#!/usr/bin/env bash

# DESCRIPTION
#-----------------------------------------------------------------
# Peter Hickey (peter.hickey@gmail.com)
# 28/08/2014
# A example bash script for processing a BAM file with methtuple
# in parallelon a chromosome-by-chromsome basis.
# This script makes extensive use of GNU parallel
# [O. Tange (2011): GNU Parallel - The Command-Line Power Tool, ;login:
# The USENIX Magazine, February 2011:42-47.
# https://www.gnu.org/software/parallel/].
#-----------------------------------------------------------------

# The MIT License (MIT)
#
# Copyright (c) [2012-2015] [Peter Hickey (peter.hickey@gmail.com)]
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# REQUIREMENTS
#-------------------------------------------------------------------------------
# REQUIRES: methtuple, SAMtools (>= v0.1.19), GNU parallel and
# Rscript must be in your $PATH.
#-------------------------------------------------------------------------------

# WARNINGS
#-------------------------------------------------------------------------------
# WARNING: Script will produce output as sub-directories of $OUTDIR
# WARNING: Roughly, this script runs 1 job per chromosome and each job can use
# 1-10Gb of memory. Larger chromosomes and larger values of "m" require more
# memory and more computational time.
#-------------------------------------------------------------------------------

# NOTES
#-------------------------------------------------------------------------------
# NOTE: You may wish to increase the maximum memory per thread when running
# samtools sort via the -m flag. I have not explored the effects of increasing
# the number of threads available for sorting and compression (via the -@
# flag); in particular, I do not know whether this plays nicely with GNU
# parallel.

# Required variables
#-------------------------------------------------------------------------------
# The path of the BAM file, e.g. data/cs_pe_directional.bam.
# The BAM must be sorted in coordinate order and indexed.
BAM=
# The path to the directory used for output, e.g. my_outdir.
# Multiple subdirectories and files will be created in this folder.
OUTDIR=
# A bash array of the "m" in m-tuples, e.g. (1 2 3).
M=
# A bash array of the chromosomes to be processed.
# It is recommended, but not required, that this array be in karyotypic order
# e.g. (chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13
#         chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)
CHROMS=
# Additional parameters to be passed to methtuple e.g.
# "--methylation-type CG --ignore-duplicates --min-mapq 0 --ir1p 1-5,98-100
#  --ir2p 1-10,98-100 --overlap-filter XM"
METHTUPLE_OPTIONS=
# SE (single-end data) or PE (paired-end data).
SEQ_TYPE=
# The number of cores to be used by GNU parallel, e.g. 8
NUM_CORES=
# The path to the helper script "tabulate_hist.R"
TABULATE_HIST=
#-------------------------------------------------------------------------------

# Description of algorithm
#-------------------------------------------------------------------------------
# The idea is as follows:
# (1) Split the BAM by chromosome and create chromosome-level BAM files.
# The original BAM must be sorted in coordinate order and indexed.
#       (a) if (Single-end data): Continue because there is no need to sort
#       single-end data.
#       (b) if (Paired-end data): Sort the chromosome-level BAMs by queryname
#       order.
# (2) For each m-tuple value:
#       (a) Process each chromosome-level BAM with methtuple.
#       (b) Concetanate all chromosome-level m-tuples files.
# (3) Clean up by removing chromosome-level BAM files
#-------------------------------------------------------------------------------

# Step 0
#-------------------------------------------------------------------------------
# Check variables
if [[ -n "${BAM}" && -n "${OUTDIR}" && -n "${M}" && -n "${CHROMS}" \
  && -n "${METHTUPLE_OPTIONS}" && -n "${SEQ_TYPE}" ]]
then
  echo "BAM file = ${BAM}"
  echo "Output directory = ${OUTDIR}"
  echo "'m' in m-tuples = ${M[@]}"
  echo "Chromosomes = ${CHROMS[@]}"
  echo "Additional methtuple options = ${METHTUPLE_OPTIONS}"
  echo "Sequencing type = ${SEQ_TYPE}"
else
  echo "Please check all parameters are valid"
  exit
fi
SAMPLE_NAME=$(basename "$BAM")
SAMPLE_NAME="${SAMPLE_NAME%.*}"
#-------------------------------------------------------------------------------

# Step 1
#-------------------------------------------------------------------------------
# Create chromosome-level BAM files (queryname-sorted if paired-end sequencing)
if [ ${SEQ_TYPE} == 'PE' ]
then
  echo "Paired-end sequencing"
  echo "Creating chromosome-level BAMs and sorting by queryname"
  parallel --joblog create_chromosome_level_BAMs.log -j ${NUM_CORES} "
  samtools view -u ${BAM} {1} | samtools sort -n - \
  ${OUTDIR}/${SAMPLE_NAME}_{1}" ::: ${CHROMS[@]}
elif [ ${SEQ_TYPE} == 'SE'  ]
then
  echo "Single-end sequencing"
  echo "Creating chromosome-level BAMs"
        parallel --joblog create_chromosome_level_BAMs.log -j ${NUM_CORES} "
        samtools view -u ${BAM} {1} > ${OUTDIR}/${SAMPLE_NAME}_{1}.bam" \
        ::: ${CHROMS[@]}
else
  echo "The sequencing type (SEQ_TYPE parameter) must be one of 'SE' (for \
  single-end sequencing) or 'PE' (paired-end for paired-end sequencing)"
  exit
fi
#-------------------------------------------------------------------------------

# Step 2
#-------------------------------------------------------------------------------
# mkdir for output
parallel -j ${NUM_CORES} "mkdir -p ${OUTDIR}/{1}_tuples/log \
${OUTDIR}/{1}_tuples/hist" ::: ${M[@]}

# Extract methylation loci m-tuples for each chromosome
echo "Extracting methylation loci m-tuples for each chromosome..."
parallel --joblog methtuple.log -j ${NUM_CORES} "methtuple \
${METHTUPLE_OPTIONS} -m {1} -o {1}_tuples/${SAMPLE_NAME}_{2} \
${OUTDIR}/${SAMPLE_NAME}_{2}.bam &> \
${OUTDIR}/{1}_tuples/log/${SAMPLE_NAME}_{1}.{2}.log" ::: ${M[@]} ::: \
${CHROMS[@]}

# The only chromosome-level files that are retained are the .hist files.
echo "Concatenating chromosome-level files into genome-level files..."
# Move all .hist files to a 'hist' directory and create the genome-wide
# .hist file.
parallel -j ${NUM_CORES} "mv ${OUTDIR}/{1}_tuples/*.hist \
${OUTDIR}/{1}_tuples/hist;
Rscript ${TABULATE_HIST} ${SAMPLE_NAME} ${OUTDIR}/{1}_tuples/hist" ::: ${M[@]}

# Create the header for the genome-wide .tsv file (for given m).
# Note the need to escape the $ character for variables defined within the
# command to avoid the shell interpreting them (i.e. EXTENSION)
parallel -j ${NUM_CORES} "
EXTENSION=\$(basename ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}_${CHROMS[0]}.*.tsv \
| rev | cut -d '.' -f -3 | rev);
head -n 1 ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}_${CHROMS[0]}.\${EXTENSION} > \
${OUTDIR}/{1}_tuples/${SAMPLE_NAME}.\${EXTENSION}" ::: ${M[@]}

# Concatenate the .tsv files and gzip (for given m).
# Note the need to escape the $ character for variables defined within the
# command executed to avoid the shell interpreting them (i.e. EXTENSION and
# CHROM)
parallel --joblog cat_tsvs.log -j ${NUM_CORES} "
EXTENSION=\$(basename ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}_${CHROMS[0]}.*.tsv \
| rev | cut -d '.' -f -3 | rev);
for CHROM in ${CHROMS[@]}
  do
		cat ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}_\${CHROM}.\${EXTENSION}  \
    | grep -v -P  '^chr\tstrand' >> \
    ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}.\${EXTENSION}
  done
gzip ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}.\${EXTENSION}" ::: ${M[@]}

# Concatenate the reads_that_failed_QC files and gzip (for given m).
# Note the need to escape the $ character for variables defined within the
# command executed to avoid the shell interpreting them (i.e. CHROM)
parallel --joblog cat_reads_that_failed_QC_files.log -j ${NUM_CORES} "
for CHROM in ${CHROMS[@]}
  do
		cat ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}_\${CHROM}.reads_that_failed_QC.txt \
    >> ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}.reads_that_failed_QC.txt
	done
gzip ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}.reads_that_failed_QC.txt" ::: ${M[@]}
#-------------------------------------------------------------------------------

# Step 3
#-------------------------------------------------------------------------------
# Tidy up by removing chromosome-level BAM files
echo "Removing chromosome-level files..."
parallel -j ${NUM_CORES} "rm ${OUTDIR}/${SAMPLE_NAME}_{1}.bam" ::: ${CHROMS[@]}
parallel -j ${NUM_CORES} "
rm ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}_*.tsv;
rm ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}_*.reads_that_failed_QC.txt" ::: ${M[@]}
#-------------------------------------------------------------------------------
