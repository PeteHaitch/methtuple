#!/usr/bin/python
import re # Python's regular expression module
import random
from numpy import array
import argparse
import sys
import csv

## Extract comethylation signal for reads overlapping multiple CpGs
## SAM input file is the output of Bismark (no pre-processing required)
## CpG overlap is as determined by Bismark output, i.e. only positions with Z or z in the methylation string are considered CpGs

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
######################################################################################

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

## TODO: [options] might include ignoring positions at end of read, ignoring low quality reads/positions
## TODO: Pretty up the file and add error checks & warnings
## TODO: Write a basic version of this program implemented in C. Compare run times of python vs. C implementation
## TODO: Discuss "random choice of CpGs in a read" with Terry - naive approach results in pairs with intra-pair distance > 40 (< 30%). Prove this analytically. Might choose outermost pair.
## TODO: Implement command line option to choose how CpG pairs are selected, e.g. random, outermost, innermost, other(?)
## TODO: Write Version 2 that looks at read content, rather than XM tag, to determine methylation status. NB: Will need to be very careful with reads aligning to Crick-strand (NB: unmethylated reverse strand reads are A at the G in the CpG and methylated reads are G at the G in the CpG.)

# Command line passer
parser = argparse.ArgumentParser(description='Extract the methylation call at two CpGs for reads that overlap multiple CpG from a Bismark SAM file. If a read overlaps more than two CpGs, select two CpGs at random. The output of this file can be used for analysing comethylation along a read.')
parser.add_argument('samfile', type=argparse.FileType('r'), nargs=1,
                   help='The (sorted) Bismark SAM file to be processed. Can be a Bismark BAM files piped via SAMtools')
parser.add_argument('output', type=argparse.FileType('w'), nargs=1,
                   help='The output filename')
parser.add_argument('--ignore5', metavar = '<int>',
                  default=0,
                  help='NOT YET IMPLEMENTED ignore <int> bases from 5\' (left) end of reads (default: 0)') ## TODO: Add 'type' argument
parser.add_argument('--ignore3', metavar = '<int>',
                  default=0,
                  help='NOT YET IMPLEMENTED ignore <int> bases from 3\' (right) end of reads (default: 0)') ## TODO: Add 'type' argument
parser.add_argument('--maxReadLength', metavar = '<int>', type = int,
                  default=150,
                  help='NOT YET IMPLEMENTED Maximum read length in bp (default: 150)')
parser.add_argument('--dummy', action='store_true',help='An example dummy TRUE/FALSE variable (default: FALSE)')

args = parser.parse_args()

# Create (possibly redundant) pointers to files from command line arguments
SAM = args.samfile[0]
OUT = args.output[0]
maxReadLength = args.maxReadLength

# Variable initialisations
CpG_pattern = re.compile(r"[Zz]")
methylated_CpG_pattern = re.compile(r"[Z]")
unmethylated_CpG_pattern = re.compile(r"[z]")
read_counter = 0 # The number of reads processed by SAM2MS.py
CpGs_per_read = [0] * (maxReadLength/2 + 2) # The number of reads containing x CpGs (x = 0, 1, ...; the minimum number of CpGs is 0 and the maximum number of CpGs is half the read length + 1)
CpG_positions = [0] * maxReadLength # Count how often each position on the read was a CpG
methylated_CpG_positions = [0] * maxReadLength # Count how many times each position on the read was a methylated CpG
unmethylated_CpG_positions = [0] * maxReadLength # Count how many times each position on the read was an unmethylated CpG

# Function to write tab-separated outputfile
tabWriter = csv.writer(OUT, delimiter='\t', quotechar=' ', quoting=csv.QUOTE_MINIMAL)

# Process each line of SAM file
for line in SAM:
    if not line.startswith('@'): # Skip SAM header lines
        read_counter += 1 # Increment the counter
        li = line.split("\t")
        chrom = li[2]
        start = int(li[3])
        XM = li[13].split(":")[2]
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
        # If there is more than one CpG in the read then we want to process that read
        if n_CpGs > 1:
            ## Choose two of those CpGs at random
            #CpG_pair = random.sample(CpG_index, 2)
            #CpG_pair.sort() # Important to sort the result to ensure CpG_1 < CpG_2 by position
            # Choose the outermost CpG pair, i.e. the CpG with the greatest intra-pair distance
            CpG_pair = [CpG_index[1], CpG_index[-1]]
            index = array(CpG_pair)
            positions = start + index # positions are 1-based leftmost mapping positions (just like in the SAM spec)
            #print >> OUT, chrom, positions[0], positions[1], XM[index[0]], XM[index[1]]
            output=[chrom, positions[0], positions[1], XM[index[0]], XM[index[1]]]
            tabWriter.writerow(output)

# Close the file connections     
SAM.close()
OUT.close()

# Print summary statistics to standard output
print '%d reads processed' % read_counter
sys.stdout.write('x')
sys.stdout.write('\t')
sys.stdout.write('Count\n')
print '------------------'
for index in range(len(CpGs_per_read)):
    print index, '\t', CpGs_per_read[index]
print '%d (%d%%) reads contained > 1 CpG' % (sum(CpGs_per_read) - CpGs_per_read[0], 100 * (sum(CpGs_per_read) - CpGs_per_read[0])/read_counter)
print 'Position in read of methylated CpGs'
print methylated_CpG_positions
print 'Position in read of unmethylated CpGs'
print unmethylated_CpG_positions
print 'Position in read of CpGs'
print CpG_positions

