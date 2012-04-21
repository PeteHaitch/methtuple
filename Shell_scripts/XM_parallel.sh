# Peter Hickey 
# 22/04/2012

# A shell script to process a BAM files with PE_SAM2MS.py on a chromosome-by-chromosome basis
# Produces one ".ms" file for each chromosome/contig
# These .ms files should be concatenated at the completion of this script

# Project-specific variables
BAM=DM_ADS-adipose_Bismark_Lister_parameters.bam
N_PARALLEL=22 # The maximum number of subprocesses to run at any one time

# Get the list of chromosomes
samtools view -H ${BAM} | grep "@SQ" | awk '{ print $2 }' | awk -F "SN:" '{ print $2}' > chrnames.tmp
CHRS=(` cat "chrnames.tmp" `)

# Process the files in batches of size ${N_PARALLEL}
N_BATCHES=$(echo "(${#CHRS[@]} + ${N_PARALLEL} - 1)/${N_PARALLEL}" | bc)  # Number of batches required to process the FASTQ-splits. While this expression looks ugly, it simply computes CEILING(number of chromosomes / ${N_PARALLEL})

for (( i=0; i < ${N_BATCHES}; i++ ))
do
    echo "Processing batch $[$i + 1] of ${N_BATCHES} with PE_SAM2MS_v2.py..."

    # Create a bam file for each chromosome
    for (( j=0; j < ${N_PARALLEL}; j++ ))
    do
	idx=$[$i * ${N_PARALLEL} + $j]
	samtools view -b ${BAM} ${CHRS[${idx}]} > ${CHRS[${idx}]}.bam & 
    done
    wait

    # Index each bam file
    for (( j=0; j < ${N_PARALLEL}; j++ ))
    do
	idx=$[$i * ${N_PARALLEL} + $j]
	samtools index ${CHRS[${idx}]}.bam & 
    done
    wait

    # Processes BAM file with PE_SAM2MS_v2.py
    for (( j=0; j < ${N_PARALLEL}; j++ ))
    do
	idx=$[$i * ${N_PARALLEL} + $j]
	python PE_SAM2MS_v2.py ${CHRS[${idx}]}.bam ${CHRS[${idx}]}.ms & 
    done
    wait    
    
    # Remove temporary BAM files and indices
    for (( j=0; j < ${N_PARALLEL}; j++ ))
    do
	idx=$[$i * ${N_PARALLEL} + $j]
	rm ${CHRS[${idx}]}.bam ${CHRS[${idx}]}.bam.bai 
    done
    wait

done

