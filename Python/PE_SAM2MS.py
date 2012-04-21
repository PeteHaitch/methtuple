#!/usr/bin/python
import re
import random
from numpy import array
import argparse
import sys
import csv
import pysam

# Firstly need to correct the FLAG (to properly encode the strand information) and fix the read-names using correct_Bismark_PE_SAM.py
# Input must have been corrected with correct_Bismark_PE_SAM.py!!! Only works for directional libraries.
# Only allows "outermost" CpG-pairs
# If a read1 and read2 overlap (i.e. the DNA fragment was shorter then the sum of the read lenghts) I exclude the overlapping positions from the 3' end of read2 since read2 is generally of lower quality. A better approach would be to trim bases from either read1 or read2 depending on which read had the lower quality 'overlap-bases'. The Epigenomic Roadmap recommend that "Paired-reads that have partial overlap in genome coverage should be trimmed from the 3' so as to avoid treating sequence derived from multiple passes of the same genomic DNA fragment as independent data points." (http://www.roadmapepigenomics.org/files/protocols/data/dna-methylation/MethylC-SeqStandards_FINAL.pdf). A further improvement would be a check that the overlap reports a consistent methylation call. A further complication is addressed by bad_overlap_pairs.
# CpGL and positionL always refer to the leftmost CpG in the pair while CpGR and positionR always refer to the rightmost CpG in the pair. Thus, for OT reads CpGL is from read1 while for OB reads CpGL is from read2 in readpairs where each read contains a CpG.
# TLEN is defined as the rightmost position of read2 (read1) minus the leftmost position of read1 (read2) plus 1 for OT (OB) aligned reads. Thus if TLEN < read1.alen + read2.alen the two reads overlap.

## TODO: Add a count of the number of CpGs between the CpG-pair to the output for each CpG pair
## TODO: Add flags to skip 'non-primary alignment' reads, 'failed reads', etc.
## TODO: Bismark's TLEN computation is incorrect for overlapping reads from a read-pair. Thus it should not be used for inferring that reads overlap. Rewrite any features that rely on inferring read overlaps so as not to rely on the TLEN value. A bug report has been submitted to the developer's of Bismark.
## TODO: Fix bug in calculation of overlapping reads - reads overlap if abs(TLEN) < read1.alen + read2.alen
## TODO: Fix bug when ignoring reads that have "bad" overlap.
## TODO: Add parameter to specify minimum and maximum insert sizes for a read-pair to be considered
## TODO: Write bad_overlap_reads to a separate file.
## TODO: Make sure that Case1C and Case2B point to the C in the CpG and not the G
## TODO: Test-case with reads that contain an overlap.
## TODO: Add counters for each of the cases, CpGs_per_read, CpG_positions, methylated_CpG_positions, unmethylated_CpG_positions, bad_overlap_pairs
## TODO: Add --ignore5 and --ignore3 options
## TODO: Add check if --ignoreDuplicates is set and ignore duplicates if it is set

# Command line passer
parser = argparse.ArgumentParser(description='Extract the methylation calls for a CpG pair from reads that overlap multiple CpGs from a Bismark SAM/BAM file. The output file contains the positions of the CpG pair using 1-based co-ordinates on the forward strand of the cytosine in each CpG. If a read overlaps more than two CpGs there are several ways to construct the pairs (see the --pairChoice argument). The output of this file can be used for analysing comethylation along a read.')
parser.add_argument('SAM', metavar = 'SAM',
                  help='The path to the Bismark SAM/BAM file')
parser.add_argument('output', type=argparse.FileType('w'), nargs=1,
                   help='The output filename')
args = parser.parse_args()

# Create (possibly redundant) pointers to files from command line arguments
OUT = args.output[0]

# Opens SAM/BAM file
f = pysam.Samfile(args.SAM) 

# Variable initialisations
maxReadLength = 150
CpG_pattern = re.compile(r"[Zz]")
methylated_CpG_pattern = re.compile(r"[Z]")
unmethylated_CpG_pattern = re.compile(r"[z]")
pair_counter = 0 # The number of read pairs processed by PE_SAM2MS.py
duplicate_read_counter = 0 # The number of duplicate reads ignore by PE_SAM2MS.py
CpGs_per_read = [0] * (maxReadLength/2 + 2) # The number of reads containing x CpGs (x = 0, 1, ...; the minimum number of CpGs is 0 and the maximum number of CpGs is half the read length + 1)
CpG_positions = [0] * maxReadLength # Count how often each position on the read was a CpG
methylated_CpG_positions = [0] * maxReadLength # Count how many times each position on the read was a methylated CpG
unmethylated_CpG_positions = [0] * maxReadLength # Count how many times each position on the read was an unmethylated CpG
ignore2 = 0 # The number of bases of read2 to ignore due to overlap between read1 and read2
bad_overlap_pairs = 0 # The number of read-pairs where the leftmost position of read2 (read1) is to the left of the leftmost position of read1 (read2) for reads informative for the OT (OB) strand. Such reads are ignored when computing single-strand comethylation as they can violate the assumption that posL < posR.

# Open file handles
#f = pysam.Samfile('DM_ADS-adipose_Bismark_Lister_parameters.bam')
#OUT = open('DM_ADS-adipose_Bismark_Lister_parameters.ms', 'w')
# Function to write tab-separated outputfile
tabWriter = csv.writer(OUT, delimiter='\t', quotechar=' ', quoting=csv.QUOTE_MINIMAL)

## This looping scheme properly looks at each read-pair once and only once since it only look at read1 in pairs to ensure each pair is only parsed once
## "read" refers to read1, "mate" refers to read2 in pair
f.reset()
for read in f:
    pair_counter += 1
    if read.is_proper_pair and read.is_read1 and read.tid == read.rnext and not read.is_duplicate:
        chrom = f.getrname(read.tid)
        ignore2 = 0
        pointer = f.tell() # pointer to the current position in the SAM file
        # Try to find the mate (read2) of the read-pair, if it can't be found skip this read. Return pointer to position of read1 in the SAM file
        try: 
            mate = f.mate(read) 
        except ValueError: 
            # Invalid mate (usually post-filtered) 
            continue 
        finally: 
            f.seek(pointer)
        # Rename read and mate as the more descriptive read1 and read2
        read1 = read
        read2 = mate
        read1XM = [tag[1] for tag in read1.tags if tag[0]=='XM'][0] # The XM tag for read1. Tags are stored as a Python list of 2-tuples [("TAG_NAME", "TAG_VALUE"), ...]. This assumes there is one, and one only, XM tag per read
        read2XM = [tag[1] for tag in read2.tags if tag[0]=='XM'][0] # The XM tag for read2. Tags are stored as a Python list of 2-tuples [("TAG_NAME", "TAG_VALUE"), ...]. This assumes there is one, and one only, XM tag per read
        # Case1: read-pair informative for OT, read orientation is +/-
        if (not read1.is_reverse and read2.is_reverse):
            # If read1.tlen < 0 it means the leftmost (3') position of read2 is to the left of the leftmost (5') position of read1. These reads complicate the assignment of "leftmost" CpGs in CpG-pairs, so we exclude these read-pairs.
            if read1.tlen < 0:
                bad_overlap_pairs += 1
                print 'Ignoring read-pair', read.qname, 'due to bad overlap (=', read1.tlen, ')'
                continue
            # start1 <= start2 by definition for a properly paired OT read
            start1 = read1.pos + 1 # read1.pos is 0-based in accordance with BAM file specificiations and Python standards, regardless of whether the file is in BAM or SAM format. I convert this to 1-based positions
            start2 = read2.pos + 1 # read2.pos is 0-based in accordance with BAM file specificiations and Python standards, regardless of whether the file is in BAM or SAM format. I convert this to 1-based positions
            strand = "+"
            # If the two reads overlap, i.e the rightmost position of read1 is greater than or equal to the leftmost position of read2, drop the offending positions from the leftmost (3') position of read2
            if read1.positions[-1] >= read2.positions[0]: 
                ignore2 =  int(read1.positions[-1] - read2.positions[0] + 1) # In this case, ignore2 is the number of bases to ignore from the 3' end of read2.
                read2XM = read2XM[ignore2:]
            # Identify CpGs in read1 and read2
            CpG_index1 = [m.start() for m in re.finditer(CpG_pattern, read1XM)]
            CpG_index2 = [m.start() for m in re.finditer(CpG_pattern, read2XM)]
            # CaseA: > 0 CpGs in both reads of pair    
            if (len(CpG_index1) > 0 and len(CpG_index2) > 0):
                CpGL = CpG_index1[0] # Leftmost CpG in read1
                CpGR = CpG_index2[-1] # Rightmost CpG in read2
                positionL = start1 + CpGL
                positionR = start2 + CpGR + ignore2 # +ignore2 to deal with ignored bases
                output=[chrom, positionL, positionR, read1XM[CpGL], read2XM[CpGR], strand] # Output is left-to-right, thus read1XM appears before read2XM for OT read-pairs.
                if positionL > positionR:
                    print "Error: Case1A posL > posR"
                    print read1
                    print read2
                    print output
                    break
                tabWriter.writerow(output)
            # CaseB: > 1 CpG in read1 but 0 CpGs in read2
            elif (len(CpG_index1) > 1 and len(CpG_index2) == 0):
                CpGL = CpG_index1[0]
                CpGR = CpG_index1[-1]
                positionL = start1 + CpGL
                positionR = start1 + CpGR
                output=[chrom, positionL, positionR, read1XM[CpGL], read1XM[CpGR], strand]
                if positionL > positionR:
                    print "Error: Case1B posL > posR"
                    print read1
                    print read2
                    print output
                    break
                tabWriter.writerow(output)
            # CaseC: 0 CpGs in read1 but > 1 CpGs in read2
            elif (len(CpG_index1) == 0 and len(CpG_index2) > 1):
                CpGL = CpG_index2[0]
                CpGR = CpG_index2[-1]
                positionL = start2 + CpGL
                positionR = start2 + CpGR
                output=[chrom, positionL, positionR, read2XM[CpGL], read2XM[CpGR], strand]
                if positionL > positionR:
                    print "Error: Case1C posL > posR"
                    print read1
                    print read2
                    print output
                    break
                tabWriter.writerow(output)
            # CaseD: < 2 CpGs in read-pair
            elif (len(CpG_index1) + len(CpG_index2)) < 2:
                continue
            else:
                print 'Unexpected CpG count at read', read.qname
                break
            # Case2: read-pair informative for OB
        elif (read1.is_reverse and not read2.is_reverse):
            # If read2.tlen < 0 it means the leftmost (3') position of read1 is to the left of the leftmost (5') position of read2. These reads complicate the assignment of "leftmost" CpGs in CpG-pairs, so we exclude these read-pairs
            if read2.tlen < 0:
                bad_overlap_pairs += 1
                print 'Ignoring read-pair', read.qname, 'due to bad overlap (=', read2.tlen, ')'
                continue
            # start1 >= start2 by definition for a properly paired OB read
            start1 = read1.pos + 1 # read1.pos is 0-based in accordance with BAM file specificiations and Python standards, regardless of whether the file is in BAM or SAM format. I convert this to 1-based positions
            start2 = read2.pos + 1 # read2.pos is 0-based in accordance with BAM file specificiations and Python standards, regardless of whether the file is in BAM or SAM format. I convert this to 1-based positions
            strand = "-"
            # If the two reads overlap, i.e the rightmost position of read2 is greater than or equal to the leftmost position of read1, drop the offending positions from the rightmost (3') position of read2
            if read2.positions[-1] >= read1.positions[0]:
                ignore2 = len(read2XM) - int(read2.positions[-1] - read1.positions[0] + 1) # In this case, ignore2 is the number of bases to retain from the 5' end of read2.
                read2XM = read2XM[:ignore2]
            # Identify CpGs in read1 and read2
            CpG_index1 = [m.start() for m in re.finditer(CpG_pattern, read1XM)]
            CpG_index2 = [m.start() for m in re.finditer(CpG_pattern, read2XM)]
            # CaseA: A CpG in both reads of pair
            if (len(CpG_index1) > 0 and len(CpG_index2) > 0):
                CpGL = CpG_index2[0] # Leftmost CpG in read2
                CpGR = CpG_index1[-1] # Rightmost CpG in read1
                positionL = start2 + CpGL - 1 # No need to ignore2 since 3' bases of read2 are rightmost; -1 so as to point to C on the forward strand in the CpG 
                positionR = start1 + CpGR - 1 # -1 so as to point to C on the forward strand in the CpG
                output=[chrom, positionL, positionR, read2XM[CpGL], read1XM[CpGR], strand] # Output is left-to-right, thus read2XM appears before read1XM for OB read-pairs.
                if positionL > positionR:
                    print "Error: Case2A posL > posR"
                    print read1
                    print read2
                    print output
                    break
                tabWriter.writerow(output)
                # CaseB: > 1 CpG in read1 but 0 CpGs in read2
            elif (len(CpG_index1) > 1 and len(CpG_index2) == 0):
                CpGL = CpG_index1[0] # Leftmost CpG in read1
                CpGR = CpG_index1[-1] # Rightmost CpG in read1
                positionL = start1 + CpGL - 1 # -1 so as to point to C on the forward strand in the CpG 
                positionR = start1 + CpGR - 1 # -1 so as to point to C on the forward strand in the CpG
                output=[chrom, positionL, positionR, read1XM[CpGL], read1XM[CpGR], strand]
                if positionL > positionR:
                    print "Error: Case2B posL > posR"
                    print read1
                    print read2
                    print output
                    break
                tabWriter.writerow(output)
            # CaseC: 0 CpGs in read1 but > 1 CpGs in read2
            elif (len(CpG_index1) == 0 and len(CpG_index2) > 1):
                CpGL = CpG_index2[0] # Leftmost CpG in read2
                CpGR = CpG_index2[-1] # Rightmost CpG in read2
                positionL = start2 + CpGL - 1 # -1 so as to point to C on the forward strand in the CpG 
                positionR = start2 + CpGR - 1 # -1 so as to point to C on the forward strand in the CpG
                output=[chrom, positionL, positionR, read2XM[CpGL], read2XM[CpGR], strand]
                if positionL > positionR:
                    print "Error: Case2C posL > posR"
                    print read1
                    print read2
                    print output
                    break
                tabWriter.writerow(output)
            # CaseD: < 2 CpGs in read-pair
            elif (len(CpG_index1) + len(CpG_index2)) < 2:
                continue
            else:
                print 'Unexpected CpG count at read', read.qname
                break
            # Case3: read-pair in unexpected orientation
        else:
            print 'Unexpected orientation at read', read.qname
            break

f.close()
OUT.close()
print 'Total number of read-pairs with bad overlaps = ', bad_overlap_pairs
