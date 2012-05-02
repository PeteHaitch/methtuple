#!/usr/bin/python
import re
import random
from numpy import array
import argparse
import sys
import csv
import pysam

## Extract comethylation signal for reads overlapping multiple CpGs
## SAM input file is the output of Bismark (no pre-processing required)
## CpG overlap is as determined by Bismark output, i.e. only positions with Z or z in the methylation string are considered CpGs
## This program is a re-write of SAM2MS_deprecated.py in order to utilise the pysam library. This version contains much additional functionality.


## This program is Copyright (C) 2012, Peter Hickey (hickey@wehi.edu.au)

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.


# Order of fields in Bismark SAM file ###
  #####################################################################################
  ### (1) QNAME                                                                     ###
  ### (2) FLAG                                                                      ### 
  ### (3) RNAME                                                                     ###
  ### (4) POS                                                                       ###
  ### (5) MAPQ                                                                      ###
  ### (6) CIGAR                                                                     ###
  ### (7) RNEXT                                                                     ###
  ### (8) PNEXT                                                                     ###
  ### (9) TLEN                                                                      ###
  ### (10) SEQ                                                                      ###
  ### (11) QUAL                                                                     ###
  ### (12) NM-tag (edit distance to reference)                                      ###
  ### (13) XX-tag (base-by-base mismatches to the reference, not including indels)  ###
  ### (14) XM-tag (methylation call string)                                         ###
  ### (15) XR-tag (read conversion state for the alignment)                         ###
  ### (16) XG-tag (genome conversion state for the alignment)                       ###
  #####################################################################################

# Explanation of methylation string (XM) 
  #################################################################
  ### . for bases not involving cytosines                       ###
  ### X for methylated C in CHG context (was protected)         ###
  ### x for not methylated C in CHG context (was converted)     ###
  ### H for methylated C in CHH context (was protected)         ###
  ### h for not methylated C in CHH context (was converted)     ###
  ### Z for methylated C in CpG context (was protected)         ###
  ### z for not methylated C in CpG context (was converted)     ###
  #################################################################

### Stranded-ness in Bismark SAM files
# A read that aligns to the reverse strand has its read sequence reverse-complemented prior to being stored in the SEQ field of the SAM file (consistent with Bowtie)
# This means that the 5' end of the original read (aka the 'good quality end') is at the leftmost position in the SEQ field regardless of which strand the read aligned against.

### A note on duplicates
# True PCR duplicates may be missed by MarkDuplicates if the reads have undergone different amounts of trimming since the start and end co-ordinates of the duplicate reads no longer need coincide even if both reads have the same POS field in the SAM file.

## TODO: Extract strand information from XS tag, rather than relying on read orientation (which is the incorrect way to infer strand information).
## TODO: Add a count of the number of CpGs between the CpG-pair to the output for each CpG pair
## TODO: Add normalising factor to "position in read of (un)methylated CpGs" based on number of reads of length n
## TODO: Use "read.positions" to extract position of CpG in reference genome rather than my method of read.start + index. This method may(?) generalise better to reads containing indels, i.e. using the Bowtie2 pipeline
## TODO: Clarifiy how "read.positions" works if there is an insertion or deletion in the read; example read HWI-ST567_0245:2:65:16410:91305#GGCTAC in 46210.bam
## TODO: Discuss "random choice of CpGs in a read" with Terry - naive approach results in pairs with intra-pair distance > 40 (< 30%). Prove this analytically. Might choose outermost pair instead.
## TODO: Paired-end version. Be careful of orientation of read pairs and stranded-ness for ignore5, ignore3 options and parsing of XM tag
## TODO: Add flags to skip 'non-primary alignment' reads, 'failed reads', etc.

# Command line passer
parser = argparse.ArgumentParser(description='Extract the methylation calls for a CpG pair from reads that overlap multiple CpGs from a Bismark SAM/BAM file. The output file contains the positions of the CpG pair using 1-based co-ordinates on the forward strand of the cytosine in each CpG. If a read overlaps more than two CpGs there are several ways to construct the pairs (see the --pairChoice argument). The output of this file can be used for analysing comethylation along a read.')
parser.add_argument('SAM', metavar = 'SAM',
                  help='The path to the Bismark SAM/BAM file')
parser.add_argument('read1_output', type=argparse.FileType('w'), nargs=1,
                   help='The output filename for read1')
parser.add_argument('read2_output', type=argparse.FileType('w'), nargs=1,
                   help='The output filename for read2')
parser.add_argument('-5', '--ignore5', metavar = '<int>', type = int,
                  default=0,
                  help='Ignore <int> bases from 5\' (left) end of reads (default: 0). WARNING: Parameter value not sanity checked by program.')
parser.add_argument('-3', '--ignore3', metavar = '<int>', type = int,
                  default=0,
                  help='Ignore <int> bases from 3\' (right) end of reads (default: 0). WARNING: Parameter value not sanity checked by program.')
parser.add_argument('--maxReadLength', metavar = '<int>', type = int,
                  default=150,
                  help='Maximum read length in bp (default: 150)')
parser.add_argument('--pairChoice', metavar = '<string>',
                  default="random",
                  help='Method for constructing CpG pairs in reads that contain more than two CpGs: random, leftmost or outermost (default: random)')
parser.add_argument('--ignoreDuplicates', action='store_true',help='Ignore reads that have been flagged as PCR duplicates by Picard\'s MarkDuplicates function')

args = parser.parse_args()

# Create (possibly redundant) pointers to files from command line arguments
OUT1 = args.read1_output[0]
OUT2 = args.read2_output[0]
maxReadLength = args.maxReadLength
pairChoice = args.pairChoice
ignore5 = args.ignore5
ignore3 = args.ignore3

# Variable initialisations
CpG_pattern = re.compile(r"[Zz]")
methylated_CpG_pattern = re.compile(r"[Z]")
unmethylated_CpG_pattern = re.compile(r"[z]")
read_counter = 0 # The number of reads processed by SAM2MS.py
duplicate_read_counter = 0 # The number of duplicate reads ignore by SAM2MS.py
CpGs_per_read = [0] * (maxReadLength/2 + 2) # The number of reads containing x CpGs (x = 0, 1, ...; the minimum number of CpGs is 0 and the maximum number of CpGs is half the read length + 1)
CpG_positions = [0] * maxReadLength # Count how often each position on the read was a CpG
methylated_CpG_positions = [0] * maxReadLength # Count how many times each position on the read was a methylated CpG
unmethylated_CpG_positions = [0] * maxReadLength # Count how many times each position on the read was an unmethylated CpG

# Function to write tab-separated output files
tabWriter1 = csv.writer(OUT1, delimiter='\t', quotechar=' ', quoting=csv.QUOTE_MINIMAL)
tabWriter2 = csv.writer(OUT2, delimiter='\t', quotechar=' ', quoting=csv.QUOTE_MINIMAL)

# Opens SAM/BAM file
sam = pysam.Samfile(args.SAM) 

# Process each line of SAM file
for read in sam:
    # Skip duplicate reads if requested
    if args.ignoreDuplicates and read.is_duplicate:
        read_counter += 1
        duplicate_read_counter += 1
        continue
    else:
        read_counter += 1 # Increment the counter
        chrom = sam.getrname(read.tid) # 
        start = read.pos + 1 # read.pos is 0-based in accordance with BAM file specificiations and Python standards, regardless of whether the file is in BAM or SAM format. I convert this to 1-based positions.
        XM = [tag[1] for tag in read.tags if tag[0]=='XM'][0] # Tags are stored as a Python list of 2-tuples [("TAG_NAME", "TAG_VALUE"), ...]. This assumes there is one, and one only, XM tag.
        CpG_index = [m.start() for m in re.finditer(CpG_pattern, XM)] # Stores index of CpGs in read (C for Watson strand, G for Crick strand)
        methylated_CpG_index = [m.start() for m in re.finditer(methylated_CpG_pattern, XM)] # Stores index of methylated CpGs in read
        unmethylated_CpG_index = [m.start() for m in re.finditer(unmethylated_CpG_pattern, XM)] # Stores index of unmethylated CpGs in read
        n_CpGs = len(CpG_index) # The number of CpGs the read overlaps
        CpGs_per_read[n_CpGs] += 1 # Record how many CpGs were covered by the read
        # Increment the counters
        for i in CpG_index:
            CpG_positions[i] += 1
        for i in methylated_CpG_index:
            methylated_CpG_positions[i] += 1
        for i in unmethylated_CpG_index:
            unmethylated_CpG_positions[i] += 1
        # Ignore --ignore5 and --ignore3 bases from read
        CpG_index = [pos for pos in CpG_index if pos > ignore5 and pos <= (len(XM) - ignore3)]
        n_CpGs = len(CpG_index) # Update the number of CpGs the read overlaps post-"ignore"
        # If there is more than one CpG in the read then we want to process that read
        if n_CpGs > 1:
            # Create a CpG pair - the choice is random, leftmost or outermost.
            # The outermost pair ensures the greatest intra-pair distance.
            if pairChoice=="random":
                CpG_pair = random.sample(CpG_index, 2)
                CpG_pair.sort() # Important to sort the result to ensure CpG_1 < CpG_2 by position
            elif pairChoice=="outermost":
                CpG_pair = [CpG_index[0], CpG_index[-1]]
            elif pairChoice=="leftmost":
                CpG_pair = [CpG_index[0], CpG_index[1]]
            else:
                sys.exit("Error: pairChoice must be one of random, leftmost or outermost. Please retry.")
            index = array(CpG_pair)
            # If a read maps to the reverse strand the Z/z characters in the XM string point to the G in the CpG - I want to point to the C in the CpG so I move the position coordinates 1bp to the left
            if read.is_reverse :
                positions = start + index - 1 # positions are 1-based positions of the cytosines on the forward strand in each CpG
                strand = "-"
            else:
                positions = start + index # positions are 1-based positions of the cytosines on the forward strand in each CpG
                strand = "+"
            output=[chrom, positions[0], positions[1], XM[index[0]], XM[index[1]], strand]
            if read.is_read1:
                tabWriter1.writerow(output)
            elif read.is_read2:
                tabWriter2.writerow(output)
                

# Close the file connections
sam.close()
OUT1.close()
OUT2.close()

# Print summary statistics to standard output
print '%d reads processed' % read_counter
if args.ignoreDuplicates:
    print '%d duplicate reads ignored' % duplicate_read_counter
print '\nHistogram of CpGs per read'
sys.stdout.write('x')
sys.stdout.write('\t')
sys.stdout.write('Count\n')
print '------------------'
for index in range(len(CpGs_per_read)):
    print index, '\t', CpGs_per_read[index]
print '%d (%d%%) reads contained > 1 CpG (this includes ignored CpGs!)' % (sum(CpGs_per_read) - CpGs_per_read[0], 100 * (sum(CpGs_per_read) - CpGs_per_read[0])/read_counter)
print 'Position in read of methylated CpGs'
print methylated_CpG_positions
print 'Position in read of unmethylated CpGs'
print unmethylated_CpG_positions
print 'Position in read of CpGs'
print CpG_positions