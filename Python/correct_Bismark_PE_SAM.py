#!/usr/bin/python
import argparse
import sys
import pysam

## Fix the FLAG values in a Bismark paired-end SAM file. 
## Specifically, correct the strand information in the FLAG and add a tag (XS:Z:<tag> to encode which DNA-strand the read is informative for, where <tag> = OT, CTOT, OB or CTOB.
## See accompanying Word document "Paired_end_read_orientation.docx" for details.

## TODO: Add funtion to correct single-end data
## TODO: Change XS-tag to XI-tag so it doesn't clash with Bowtie2's XS-tag

###### WARNING - CURRENTLY ONLY HANDLE DIRECTIONAL-LIBRARIES ########

# Command line passer
parser = argparse.ArgumentParser(description='Fix the FLAG values in a Bismark paired-end SAM file. Specifically, correct the strand information in the FLAG and add a tag (XS:Z:<tag> to encode which DNA-strand the read is informative for, where <tag> = OT, CTOT, OB or CTOB. See accompanying Word document "Paired_end_read_orientation.docx" for details.\nWARNING: Currently only supports directional-libraries.')
parser.add_argument('infile', metavar = 'in.bam',
                  help='The path to the original Bismark SAM/BAM file')
parser.add_argument('outfile', metavar = 'out.bam',
                  help='The path to the new SAM/BAM file')
args = parser.parse_args()

# Open the old SAM/BAM file
OLD = pysam.Samfile(args.infile)

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
id = 'correct_Bismark_PE_SAM.py'
vn = '0.1'
cl = ' '.join(['python correct_Bismark_PE_SAM.py', args.infile, args.outfile])
XM_tag_PG = {'ID': id, 'VN': vn, 'CL': cl}
header['PG'].append(XM_tag_PG)
header['CO'] = ['Bismark SAM/BAM file strand-FLAGs corrected by correct_Bismark_PE_SAM.py']
NEW = pysam.Samfile(args.outfile, "wb", header = header)

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
