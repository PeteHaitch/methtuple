#!/usr/bin/python
import argparse
import sys
import csv

## Summarise output of SAM2MS.py
## For each unique CpG pair count the number of reads with MM (ZZ), MU (Zz), UM (zZ) and MM (zz) methylation strings and write to file. Also tabulate strand-specific values for these counts
## Input file is the output of SAM2MS.py. SAM2MS.py can be sorted by using sortBed in BEDtools
## Output file is a tab separated file of the form: (chr, pos1, pos2, MM, MU, UM, UU, MM+, MU+, UM+, UU+, MM-, MU, UM, UU-)
## R can't handle "+" or "-" characters in header row, thus I use "_OT" = orignal top for "+" strand and "_OB" = original bottom for "-" strand sufixes.

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

## TODO: Output of this script is sortable by head -n1 summarised.out > sorted_summarised.bed; tail -n+2 summarised.out | sort -k1,1 -k2,2g -k3,3g >> sorted_summarised.bed (adapted from http://cassjohnston.wordpress.com/2011/05/10/unix-sort-bed-file/). This is a lexicographic sort.
## TODO: See if lor-distance relationship is a function of coverage

# Command line passer
parser = argparse.ArgumentParser(description='Summarise the output of SAM2MS.py. For each unique CpG pair count the number of reads with MM (ZZ), MU (Zz), UM (zZ) and MM (zz) methylation strings and write to file. Also tabulate strand-specific values for these counts. The output is unsorted.')
parser.add_argument('MS', type=argparse.FileType('r'), nargs=1,
                   help='The path to the output of SAM2MS.py ')
parser.add_argument('summarised', type=argparse.FileType('w'), nargs=1,
                   help='The filename of the summarised output')
args = parser.parse_args()

# Create (possibly redundant) pointers to files from command line arguments
IN = args.MS[0]
OUT = args.summarised[0]

# Function definitions
def makeCount():
    return  {'MM': 0, 'MU': 0, 'UM': 0, 'UU': 0, 'MM_OT': 0, 'MU_OT': 0, 'UM_OT': 0, 'UU_OT': 0, 'MM_OB': 0, 'MU_OB': 0, 'UM_OB': 0, 'UU_OB': 0} # OT = orignal top, OB = original bottom

def incrementCount(count, ms, strand, line_cnt):
    if ms == 'ZZ':
        count['MM'] += 1
        if strand == '+':
            count['MM_OT'] += 1
        elif strand == '-':
            count['MM_OB'] +=1
        else:
            exit_msg = ''.join(['Error: Invalid strand at line ', str(line)])
            sys.exit(exit_msg)
        return count
    elif ms == 'Zz':
        count['MU'] += 1
        if strand == '+':
            count['MU_OT'] += 1
        elif strand == '-':
            count['MU_OB'] +=1
        else:
            exit_msg = ''.join(['Error: Invalid strand at line ', str(line)])
            sys.exit(exit_msg)
        return count
    elif ms == 'zZ':
        count['UM'] += 1
        if strand == '+':
            count['UM_OT'] += 1
        elif strand == '-':
            count['UM_OB'] +=1
        else:
            exit_msg = ''.join(['Error: Invalid strand at line ', str(line)])
            sys.exit(exit_msg)
        return count
    elif ms == 'zz':
        count['UU'] += 1
        if strand == '+':
            count['UU_OT'] += 1
        elif strand == '-':
            count['UU_OB'] +=1
        else:
            exit_msg = ''.join(['Error: Invalid strand at line ', str(line)])
            sys.exit(exit_msg)
        return count
    else:
        exit_msg = ''.join(['Error: Invalid methylation string at line ', str(line)])
        sys.exit(exit_msg)

tabWriter = csv.writer(OUT, delimiter='\t', quotechar=' ', quoting=csv.QUOTE_MINIMAL)

# Variable initialisations
PAIRS = {} # Dictionary of CpG pair IDs with keys of form chr_pos1_pos2 and values corresponding to a makeCount object
line_cnt = 0

for line in IN:
    line_cnt += 1
    line = line.rstrip()
    line = line.split('\t')
    pair_ID = '::'.join(line[:3])
    MS = ''.join(line[3:5])
    strand = line[5]
    if not pair_ID in PAIRS:
        PAIRS[pair_ID] = makeCount() # Create count object
        PAIRS[pair_ID] = incrementCount(PAIRS[pair_ID], MS, strand, line_cnt) # Increment the count
    else:
        PAIRS[pair_ID] = incrementCount(PAIRS[pair_ID], MS, strand, line_cnt) # Increment the count
        
# Print PAIRS to a tab separated file
# Be careful of the order of items when printing dictionaries

# Need "#" character at start of header line to comply with BED-specs
header = ['#chr', 'pos1', 'pos2', 'dist' ] + sorted(PAIRS.values()[0].keys())
tabWriter.writerow(header)
for pos, count in PAIRS.iteritems():
    pos = pos.split('::')
    chrom = pos[0]
    pos1 = pos[1]
    pos2 = pos[2]
    dist = int(float(pos2) - float(pos1))
    row = [chrom, pos1, pos2, dist]
    for state, value in sorted(count.iteritems()):
        row.append(value)
    tabWriter.writerow(row)

IN.close()
OUT.close()
