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
# Copyright (c) [2012-2014] [Peter Hickey (peter.hickey@gmail.com)]
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

# Variables that are specified per-samples
#-------------------------------------------------------------------------------
# The path of the BAM file, e.g. data/cs_pe_directional.bam.
# The BAM must be sorted in coordinate order and indexed.
BAM=
# The path to the directory used for output, e.g. my_outdir.
# Multiple subdirectories and files will be created in this folder.
OUTDIR=
# A bash array of the chromosomes to be processed.
# It is recommended, but not required, that this array be in karyotypic order
# e.g. (chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13
#         chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)
# Check chromosome names, e.g. chrM vs. chrMt, chrL, etc.
CHROMS=
# Read positions to ignore based on analysis of M-bias
# Set IR1P="--ir1p 0" if not ignoring any positions
IR1P=
# Not used if data are SE
IR2P=
# Single-end (SE) or paired-end (PE)
SEQ_TYPE=
# The number of cores to be used by GNU parallel, e.g. 8
NUM_CORES=

# Variables fixed for all samples analysed in thesis
#-------------------------------------------------------------------------------
# A bash array of the "m" in m-tuples, e.g. (1 2 3).
M=(1 2 3 4 5 6 7 8)
# --methylation-type
METHYLATION_TYPE="--methylation-type CG"
# All of the samples I'm analysing have mapQ set to 255, i.e. unavailable.
# This is because until recently Bismark didn't supply mapQ, and even now it
# only does it with Bowtie2.
MIN_MAPQ="--min-mapq 0"
OVERLAP_FILTER="--overlap-filter XM_ol"
# All samples use Phred33 qualities.
# Not all samples have real base qualities, e.g. most of the Lister data
MIN_BASE_QUAL="--min-base-qual 3"
# Additional parameters to be passed to methtuple e.g.
# "--methylation-type CG --ignore-duplicates --min-mapq 0 --ir1p 1-5,98-100
#  --ir2p 1-10,98-100 --overlap-filter XM"
OTHER_METHTUPLE_OPTIONS="--ignore-duplicates"
# The path to the helper script "tabulate_hist.R"
TABULATE_HIST=~/methtuple/helper_scripts/tabulate_hist.R
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
if [[ -n "${BAM}" && -n "${OUTDIR}" && -n "${CHROMS}" && -n "${IR1P}" && \
   -n "${IR2P}" && -n "${SEQ_TYPE}" && -n "${NUM_CORES}" && -n "${M}" && \
   -n "${METHYLATION_TYPE}" && -n "${MIN_MAPQ}" && -n "${OVERLAP_FILTER}" && \
   -n "${OVERLAP_FILTER}" && -n "${MIN_BASE_QUAL}" && \
   -n "${OTHER_METHTUPLE_OPTIONS}" && -n "${TABULATE_HIST}" ]]
then
  echo "BAM file = ${BAM}"
  echo "Sequencing type = ${SEQ_TYPE}"
  echo "Output directory = ${OUTDIR}"
  echo "Chromosomes = ${CHROMS[@]}"
  if [ ${SEQ_TYPE} == 'SE' ]
  then
    echo  "methtuple options = ${IR1P} ${METHYLATION_TYPE} ${MIN_MAPQ}\
    ${OVERLAP_FILTER} ${MIN_BASE_QUAL} ${OTHER_METHTUPLE_OPTIONS}"
  else [ ${SEQ_TYPE} == 'PE' ]
    echo  "methtuple options = ${IR1P} ${IR2P} ${METHYLATION_TYPE} ${MIN_MAPQ}\
    ${OVERLAP_FILTER} ${MIN_BASE_QUAL} ${OTHER_METHTUPLE_OPTIONS}"
  fi
  echo "'m' in m-tuples = ${M[@]}"
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
  # TODO: Don't exit cause that quits a ssh session, but do something else.
  exit
fi
#-------------------------------------------------------------------------------

# Step 2
#-------------------------------------------------------------------------------
# mkdir for output
parallel -j ${NUM_CORES} "mkdir -p ${OUTDIR}/{1}_tuples/log \
${OUTDIR}/{1}_tuples/hist" ::: ${M[@]}
mkdir -p ${OUTDIR}/2ac_tuples/log ${OUTDIR}/2ac_tuples/hist

# Extract methylation loci m-tuples for each chromosome
echo "Extracting methylation loci m-tuples for each chromosome..."
if [ ${SEQ_TYPE} == 'SE' ]
then
  parallel --joblog methtuple.log -j ${NUM_CORES} "/usr/bin/time -v methtuple \
  ${METHYLATION_TYPE} -m {1} ${OTHER_METHTUPLE_OPTIONS} ${MIN_MAPQ} \
  ${OVERLAP_FILTER} ${IR1P} ${MIN_BASE_QUAL} -o {1}_tuples/${SAMPLE_NAME}_{2} \
  ${OUTDIR}/${SAMPLE_NAME}_{2}.bam &> \
  ${OUTDIR}/{1}_tuples/log/${SAMPLE_NAME}_{1}.{2}.log" ::: ${M[@]} ::: \
  ${CHROMS[@]}
  # "2ac"
  parallel --joblog methtuple_2ac.log -j ${NUM_CORES} "/usr/bin/time -v \
  methtuple ${METHYLATION_TYPE} -m 2 --all-combinations \
  ${OTHER_METHTUPLE_OPTIONS} ${MIN_MAPQ} ${OVERLAP_FILTER} ${IR1P} \
  ${MIN_BASE_QUAL} -o 2ac_tuples/${SAMPLE_NAME}_{1} \
  ${OUTDIR}/${SAMPLE_NAME}_{1}.bam &> \
  ${OUTDIR}/2ac_tuples/log/${SAMPLE_NAME}_2ac.{1}.log" ::: ${CHROMS[@]}
else
  parallel --joblog methtuple.log -j ${NUM_CORES} "/usr/bin/time -v methtuple \
  ${METHYLATION_TYPE} -m {1} ${OTHER_METHTUPLE_OPTIONS} ${MIN_MAPQ} \
  ${OVERLAP_FILTER} ${IR1P} ${IR2P} ${MIN_BASE_QUAL} -o \
  {1}_tuples/${SAMPLE_NAME}_{2} ${OUTDIR}/${SAMPLE_NAME}_{2}.bam &> \
  ${OUTDIR}/{1}_tuples/log/${SAMPLE_NAME}_{1}.{2}.log" ::: ${M[@]} ::: \
  ${CHROMS[@]}
  # "2ac"
  parallel --joblog methtuple_2ac.log -j ${NUM_CORES} "/usr/bin/time -v \
  methtuple ${METHYLATION_TYPE} -m 2 --all-combinations \
  ${OTHER_METHTUPLE_OPTIONS} ${MIN_MAPQ} ${OVERLAP_FILTER} ${IR1P} ${IR2P} \
  ${MIN_BASE_QUAL} -o 2ac_tuples/${SAMPLE_NAME}_{1} \
  ${OUTDIR}/${SAMPLE_NAME}_{1}.bam &> \
  ${OUTDIR}/2ac_tuples/log/${SAMPLE_NAME}_2ac.{1}.log" ::: ${CHROMS[@]}
fi

# The only chromosome-level files that are retained are the .hist files.
echo "Concatenating chromosome-level files into genome-level files..."
# Move all .hist files to a 'hist' directory and create the genome-wide
# .hist file.
parallel -j ${NUM_CORES} "mv ${OUTDIR}/{1}_tuples/*.hist \
${OUTDIR}/{1}_tuples/hist;
Rscript ${TABULATE_HIST} ${SAMPLE_NAME} ${OUTDIR}/{1}_tuples/hist" ::: ${M[@]}
# "2ac"
mv ${OUTDIR}/2ac_tuples/*.hist ${OUTDIR}/2ac_tuples/hist
Rscript ${TABULATE_HIST} ${SAMPLE_NAME} ${OUTDIR}/2ac_tuples/hist

# Create the header for the genome-wide .tsv file (for given m).
# Note the need to escape the $ character for variables defined within the
# command to avoid the shell interpreting them (i.e. EXTENSION)
parallel -j ${NUM_CORES} "
EXTENSION=\$(basename ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}_${CHROMS[0]}.*.tsv \
| rev | cut -d '.' -f -3 | rev);
head -n 1 ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}_${CHROMS[0]}.\${EXTENSION} > \
${OUTDIR}/{1}_tuples/${SAMPLE_NAME}.\${EXTENSION}" ::: ${M[@]}
# "2ac"
EXTENSION=$(basename ${OUTDIR}/2ac_tuples/${SAMPLE_NAME}_${CHROMS[0]}.*.tsv \
| rev | cut -d '.' -f -3 | rev)
head -n 1 ${OUTDIR}/2ac_tuples/${SAMPLE_NAME}_${CHROMS[0]}.${EXTENSION} > \
${OUTDIR}/2ac_tuples/${SAMPLE_NAME}.${EXTENSION}

# Concatenate the .tsv files and gzip (for given m).
# Note the need to escape the $ character for variables defined within the
# command executed to avoid the shell interpreting them (i.e. EXTENSION and
# CHROM)
parallel --joblog cat_tsvs.log -j ${NUM_CORES} "
EXTENSION=\$(basename ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}_${CHROMS[0]}.*.tsv \
| rev | cut -d '.' -f -3 | rev);
for CHROM in ${CHROMS[@]}
  do
		cat ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}_\${CHROM}.\${EXTENSION} \
    | grep -v -P  '^chr\tstrand' >> \
    ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}.\${EXTENSION}
  done
gzip ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}.\${EXTENSION}" ::: ${M[@]}
# "2ac"
EXTENSION=$(basename ${OUTDIR}/2ac_tuples/${SAMPLE_NAME}_${CHROMS[0]}.*.tsv \
| rev | cut -d '.' -f -3 | rev)
for CHROM in ${CHROMS[@]}
  do
    cat ${OUTDIR}/2ac_tuples/${SAMPLE_NAME}_${CHROM}.${EXTENSION} \
    | grep -v -P  '^chr\tstrand' >> \
    ${OUTDIR}/2ac_tuples/${SAMPLE_NAME}.${EXTENSION}
  done
gzip ${OUTDIR}/2ac_tuples/${SAMPLE_NAME}.${EXTENSION}

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
# "2ac"
for CHROM in ${CHROMS[@]}
  do
    cat ${OUTDIR}/2ac_tuples/${SAMPLE_NAME}_${CHROM}.reads_that_failed_QC.txt \
    >> ${OUTDIR}/2ac_tuples/${SAMPLE_NAME}.reads_that_failed_QC.txt
  done
gzip ${OUTDIR}/2ac_tuples/${SAMPLE_NAME}.reads_that_failed_QC.txt
#-------------------------------------------------------------------------------

# Step 3
#-------------------------------------------------------------------------------
# Tidy up by removing chromosome-level BAM files
echo "Removing chromosome-level files..."
parallel -j ${NUM_CORES} "rm ${OUTDIR}/${SAMPLE_NAME}_{1}.bam" ::: ${CHROMS[@]}
parallel -j ${NUM_CORES} "
rm ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}_*.tsv;
rm ${OUTDIR}/{1}_tuples/${SAMPLE_NAME}_*.reads_that_failed_QC.txt" ::: ${M[@]}
rm ${OUTDIR}/2ac_tuples/${SAMPLE_NAME}_*.tsv
rm ${OUTDIR}/2ac_tuples/${SAMPLE_NAME}_*.reads_that_failed_QC.txt
#-------------------------------------------------------------------------------
