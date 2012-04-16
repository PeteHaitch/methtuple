#!/usr/bin/python
import argparse
import sys
import pysam

## Fix the FLAG values in a Bismark paired-end SAM file. 
## Specifically, correct the strand information in the FLAG and add a tag (XS:Z:<tag> to encode which DNA-strand the read is informative for, where <tag> = OT, CTOT, OB or CTOB.
## See accompanying Word document "Paired_end_read_orientation.docx" for details.

###### WARNING - CURRENTLY ONLY HANDLE DIRECTIONAL-LIBRARIES ########
###### WARNING - RENAME OLD AND NEW FILE HANDLES PRIOR TO USE #######

## TODO: Add check that TLEN is positive for read1 (read2) if read is informative for OT (OB) strand. This is necessary to ensure read1 (read2) is leftmost compared to read2 (read1) for the OT (OB) strand.
## TODO: Add program name and settings to @PG tag in header, rather than as a comment (@CO)
## TODO: Add argparse options

# Open the old SAM/BAM file
OLD = pysam.Samfile('trimmed_paired_SRR097428_1_HWI-BRUNOP20X_0637:1.fastq.gz.val_paired_1.fq.gz_bismark_pe.bam', 'rb')

#### Standard Bismark FLAG values
#  67 - read paired, read mapped in proper pair, first in pair
# 131 - read paired, read mapped in proper pair, second in pair
# 115 - read paired, read mapped in proper pair, read reverse strand, mate reverse strand, first in pair
# 179 - read paired, read mapped in proper pair, read reverse strand, mate reverse strand, second in pair
####

#### SAM-spec-compliant FLAG values
#  99 - read paired, read mapped in proper pair, mate reverse strand, first in pair
# 147 - read paired, read mapped in proper pair, read reverse strand, second in pair
#  83 - read paired, read mapped in proper pair, read reverse strand, first in pair
# 163 - read paired, read mapped in proper pair, mate reverse strand, second in pair
####

# Create new SAM/BAM file
header = OLD.header
header['CO'] = ['Bismark SAM/BAM file corrected by correct_Bismark_PE_SAM.py']
NEW = pysam.Samfile("fixed_trimmed_paired_SRR097428_1_HWI-BRUNOP20X_0637:1.fastq.gz.val_paired_1.fq.gz_bismark_pe.bam", "wb", header = header)

for read in OLD:
    read.qname = read.qname.rstrip('/[1-2]') # Remove the /1 or /2 suffix attached by Bowtie/Bismark to paired-end reads
    if read.flag == 67L and read.opt('XR') == 'CT' and read.opt('XG') == 'CT':
        read.flag = 99L
        read.tags = read.tags + [('XS', 'OT')]
    elif read.flag == 131L and read.opt('XR') == 'GA' and read.opt('XG') == 'CT':
        read.flag = 147L
        read.tags = read.tags + [('XS', 'OT')]
    elif read.flag == 115L and read.opt('XR') == 'CT' and read.opt('XG') == 'GA':
        read.flag = 83L
        read.tags = read.tags + [('XS', 'OB')]
    elif read.flag == 179L and read.opt('XR') == 'GA' and read.opt('XG') == 'GA':
        read.flag = 163L
        read.tags = read.tags + [('XS', 'OB')]
    else:
        sys.exit('Error: Unexpected FLAG/XR/XG combination. Is your sample non-directional? This script only works with for directional-librarires.')
    NEW.write(read)

NEW.close()
OLD.close()
