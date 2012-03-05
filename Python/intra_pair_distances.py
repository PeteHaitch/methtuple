#!/usr/bin/python
import re # Python's regular expression module
import random
from numpy import array
import argparse

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
## TODO: Add number of reads overlapping > x CpGs (%) for x = 1, 2, 3, 4,...
## TODO: Write Version 2 that looks at read content, rather than XM tag, to determine methylation status
## TODO: Write a basic version of this program implemented in C. Compare run times of python vs. C implementation
## TODO: Discuss "random choice of CpGs in a read" with Terry - naive approach results in pairs with intra-pair distance > 40 (< 30%). Prove this analytically.

# Command line passer
parser = argparse.ArgumentParser(description='Extract the methylation call at two CpGs for reads that overlap multiple CpG from a Bismark SAM file. If a read overlaps more than two CpGs, select two CpGs at random. The output of this file can be used for analysing comethylation along a read.')
parser.add_argument('samfile', type=argparse.FileType('r'), nargs=1,
                   help='The Bismark SAM file to be processed')
parser.add_argument('output', type=argparse.FileType('w'), nargs=1,
                   help='The output filename')
parser.add_argument('--ignore5', metavar = '<int>',
                  default=0,
                  help='NOT YET IMPLEMENTED ignore <int> bases from 5\' (left) end of reads (default: 0)')
parser.add_argument('--ignore3', metavar = '<int>',
                  default=0,
                  help='NOT YET IMPLEMENTED ignore <int> bases from 3\' (right) end of reads (default: 0)')
parser.add_argument('--dummy', action='store_true',help='An example dummy TRUE/FALSE variable (default: False)')

args = parser.parse_args()

# Create (possibly redundant) pointers to files
SAM = args.samfile[0]
OUT = args.output[0]

CpG_pattern = re.compile(r"[Zz]")
read_counter = 0 # The number of reads processed by SAM2MS.py
multiple_CpG_reads_counter = 0 # The number of reads containing > 1 CpG

# Process each line of SAM file
for line in SAM:
    if not line.startswith('@'): # Skip SAM header lines
        # Increment the counter
        read_counter += 1
        li = line.split("\t")
        chrom = li[2]
        start = int(li[3])
        XM = li[13].split(":")[2] ### TODO: CHECK ORDER OF FIELDS IN BISMARK 0.6.3 OUTPUT
        CpG_index = [m.start() for m in re.finditer(CpG_pattern, XM)] # Stores index of CpGs in read
        n_CpGs = len(CpG_index) # The number of CpGs the read overlaps
        # If there is more than one CpG in the read then we want to process that read
        if n_CpGs > 1:
            # Increment the counter
            multiple_CpG_reads_counter += 1
            positions = CpG_index #### TODO: Check positions are correct (i.e. 0-offset vs. 1-offset) - need to think about stranded-ness as well!
            print >> OUT, ','.join(map(str,positions))

SAM.close()
OUT.close()

print '%d reads processed' % read_counter
print '%d (%d%%) reads contained > 1 CpG' % (multiple_CpG_reads_counter, 100 * multiple_CpG_reads_counter/read_counter)

