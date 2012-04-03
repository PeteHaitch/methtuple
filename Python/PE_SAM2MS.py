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
# CpG1 and position1 always refer to the leftmost CpG in the pair while CpG2 and position2 always refer to the rightmost CpG in the pair. Thus, for OT reads CpG1 is from read1 while for OB reads CpG1 is from read2 in readpairs where each read contains a CpG.
# CpGL and positionL always refer to the leftmost CpG in the pair while CpGR and positionR always refer to the rightmost CpG in the pair. Thus, for OT reads CpGL is from read1 while for OB reads CpGL is from read2 in readpairs where each read contains a CpG.

## TODO: Write bad_overlap_reads to a separate file.
## TODO: Test-case with reads that contain an overlap.
## TODO: Add counters for each of the cases, CpGs_per_read, CpG_positions, methylated_CpG_positions, unmethylated_CpG_positions, bad_overlap_pairs
## TODO: Add --ignore5 and --ignore3 options

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
f = pysam.Samfile('sorted_test_cases.bam')
OUT = open('sorted_test_cases.MS', 'w')
# Function to write tab-separated outputfile
tabWriter = csv.writer(OUT, delimiter='\t', quotechar=' ', quoting=csv.QUOTE_MINIMAL)

## This looping scheme properly looks at each read-pair once and only once since it only look at read1 in pairs to ensure each pair is only parsed once
## "read" refers to read1, "mate" refers to read2 in pair
f.reset()
for read in f:
    pair_counter += 1
    ## TODO: Add check if --ignoreDuplicates is set and ignore duplicates if it is set
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
        readXM = [tag[1] for tag in read.tags if tag[0]=='XM'][0] # The XM tag for read1. Tags are stored as a Python list of 2-tuples [("TAG_NAME", "TAG_VALUE"), ...]. This assumes there is one, and one only, XM tag per read
        mateXM = [tag[1] for tag in mate.tags if tag[0]=='XM'][0] # The XM tag for read2. Tags are stored as a Python list of 2-tuples [("TAG_NAME", "TAG_VALUE"), ...]. This assumes there is one, and one only, XM tag per read
        # Case1: read-pair informative for OT, read orientation is +/-
        if (not read.is_reverse and mate.is_reverse):
            # If TLEN < 0 it means the leftmost (3') position of read2 is left of the leftmost (5') position of read1. These reads complicate
            if read.tlen < 0:
                bad_overlap_pairs += 1              
            # start1 <= start2 by definition for a properly paired OT read
            start1 = read.pos + 1 # read.pos is 0-based in accordance with BAM file specificiations and Python standards, regardless of whether the file is in BAM or SAM format. I convert this to 1-based positions
            start2 = mate.pos + 1 # mate.pos is 0-based in accordance with BAM file specificiations and Python standards, regardless of whether the file is in BAM or SAM format. I convert this to 1-based positions
            strand = "+"
            # If the two reads overlap, drop the offending positions from the leftmost (3') position of read2
            if read.positions[-1] >= mate.positions[0]: 
                ignore2 =  int(read.positions[-1] - mate.positions[0] + 1)
                mateXM = mateXM[ignore2:]
            # Identify CpGs in read1 and read2
            CpG_index1 = [m.start() for m in re.finditer(CpG_pattern, readXM)]
            CpG_index2 = [m.start() for m in re.finditer(CpG_pattern, mateXM)]
            # CaseA: > 0 CpGs in both reads of pair    
            if (len(CpG_index1) > 0 and len(CpG_index2) > 0):
                CpG1 = CpG_index1[0] # Leftmost CpG in read1
                CpG2 = CpG_index2[-1] # Rightmost CpG in read2
                position1 = start1 + CpG1
                position2 = start2 + CpG2 + ignore2 # +ignore2 to deal with ignored bases
                output=[chrom, position1, position2, readXM[CpG1], mateXM[CpG2], strand] # Output is left-to-right thus read1 appears first for OT reads.
                if position1 > position2:
                    print "Error: Case1A pos1 > pos2"
                    print read
                    print mate
                    print output
                    break
                tabWriter.writerow(output)
            # CaseB: > 1 CpG in read1 but 0 CpGs in read2
            elif (len(CpG_index1) > 1 and len(CpG_index2) == 0):
                CpG1 = CpG_index1[0]
                CpG2 = CpG_index1[-1]
                position1 = start1 + CpG1
                position2 = start1 + CpG2
                output=[chrom, position1, position2, readXM[CpG1], readXM[CpG2], strand]
                if position1 > position2:
                    print "Error: Case1B pos1 > pos2"
                    print read
                    print mate
                    print output
                    break
                tabWriter.writerow(output)
            # CaseC: 0 CpGs in read1 but > 1 CpGs in read2
            elif (len(CpG_index1) == 0 and len(CpG_index2) > 1):
                CpG1 = CpG_index2[0]
                CpG2 = CpG_index2[-1]
                position1 = start2 + CpG1
                position2 = start2 + CpG2
                output=[chrom, position1, position2, mateXM[CpG1], mateXM[CpG2], strand] # Need to add ignore2 to pos2 because otherwise the position is incorrectly translated to the left
                if position1 > position2:
                    print "Error: Case1C pos1 > pos2"
                    print read
                    print mate
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
        elif (read.is_reverse and not mate.is_reverse):
            # start1 >= start2 by defintion for a properly paired OB read
            start1 = read.pos + 1 # read.pos is 0-based in accordance with BAM file specificiations and Python standards, regardless of whether the file is in BAM or SAM format. I convert this to 1-based positions
            start2 = mate.pos + 1 # mate.pos is 0-based in accordance with BAM file specificiations and Python standards, regardless of whether the file is in BAM or SAM format. I convert this to 1-based positions
            strand = "-"
            # If the two reads overlap, drop the offending positions from the leftmost (3') position of read2
            if mate.positions[-1] >= read.positions[0]:
                ignore2 = len(mateXM) - int(mate.positions[-1] - read.positions[0] + 1)
                mateXM = mateXM[:ignore2]
            # Identify CpGs in read1 and read2
            CpG_index1 = [m.start() for m in re.finditer(CpG_pattern, mateXM)]
            CpG_index2 = [m.start() for m in re.finditer(CpG_pattern, readXM)]
            # CaseA: A CpG in both reads of pair
            if (len(CpG_index1) > 0 and len(CpG_index2) > 0):
                CpG1 = CpG_index1[0] # Leftmost CpG in read2
                CpG2 = CpG_index2[-1] # Rightmost CpG in read1
                position1 = start2 + CpG1 - 1 # No need to ignore2 since 3' bases of read2 are rightmost; -1 so as to point to C on the forward strand in the CpG 
                position2 = start1 + CpG2 - 1 # -1 so as to point to C on the forward strand in the CpG
                output=[chrom, position1, position2, mateXM[CpG1], readXM[CpG2], strand] # Output is left-to-right thus read2 appears first for OB reads.
                if position1 > position2:
                    print "Error: Case2A pos1 > pos2"
                    print read
                    print mate
                    print output
                    break
                tabWriter.writerow(output)
                # CaseB: > 1 CpG in read1 but 0 CpGs in read2
            elif (len(CpG_index1) > 1 and len(CpG_index2) == 0):
                CpG1 = CpG_index1[0] # Leftmost CpG in read1
                CpG2 = CpG_index1[-1] # Rightmost CpG in read1
                position1 = start1 + CpG1 - 1 # -1 so as to point to C on the forward strand in the CpG 
                position2 = start1 + CpG2 - 1 # -1 so as to point to C on the forward strand in the CpG
                output=[chrom, position1, position2, readXM[CpG1], readXM[CpG2], strand]
                if position1 > position2:
                    print "Error: Case2B pos1 > pos2"
                    print read
                    print mate
                    print output
                    break
                tabWriter.writerow(output)
            # CaseC: 0 CpGs in read1 but > 1 CpGs in read2
            elif (len(CpG_index1) == 0 and len(CpG_index2) > 1):
                CpG1 = CpG_index2[0] # Leftmost CpG in read2
                CpG2 = CpG_index2[-1] # Rightmost CpG in read2
                position1 = start2 + CpG1 - 1 # -1 so as to point to C on the forward strand in the CpG 
                position2 = start2 + CpG2 - 1 # -1 so as to point to C on the forward strand in the CpG
                output=[chrom, position1, position2, mateXM[CpG1], mateXM[CpG2], strand] # Need to add ignore2 to pos2 because otherwise the position is incorrectly translated to the left
                if position1 > position2:
                    print "Error: Case2C pos1 > pos2"
                    print read
                    print mate
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
