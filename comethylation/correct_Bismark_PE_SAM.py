#!/usr/bin/env python
import argparse
import sys
import pysam

#### LICENSE ####
## Copyright (C) 2012 - 2014 Peter Hickey (peter.hickey@gmail.com)

## This file is part of Comethylation.

## Comethylation is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
## Comethylation is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with Comethylation.  If not, see <http://www.gnu.org/licenses/>.

#### IMPORTANT ####
# The typical use case for this script is to fix the FLAG for each read prior to running comethylation.py for SAM/BAM files created with Bismark version < 0.8.3. However, this fix is now done on-the-fly by comethylation.py and so this script is no longer required for this use case.
# Someone may still want to run this script in order to correct the FLAG for each read in a SAM/BAM files created with Bismark version < 0.8.3, e.g. so that paired-end reads are properly displayed in IGV.
#### IMPORTANT ####


## Fix the FLAG values in a Bismark paired-end SAM file. 
## Specifically, correct the strand information in the FLAG and add a tag (XS:Z:<tag> to encode which DNA-strand the read is informative for, where <tag> = OT, CTOT, OB or CTOB.
## See accompanying Word document "Paired_end_read_orientation.docx" for details.

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
