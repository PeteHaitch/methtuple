#!/usr/bin/python
import re
import random
from numpy import array
import argparse
import sys
import csv
import pysam

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

### Program description ###
############################################################################################################################################################################################
# Extract co-methylation measurements at CpGs for the aligned reads from a BS-Seq experiment.
# SAM/BAM must be generated by Bismark. Other aligners will be supported in future releases. The SAM/BAM file must comply with the SAM spec (v1.4). The Bismark SAM file does not comply with these standards and must be pre-processed with correct_Bismark.py. The reason for this is that the Bismark SAM file does not correctly encode the orientation of read-pairs in the SAM flag field.
# The SAM/BAM can contain single-end reads, paired-end reads or a combination of the two types. The type of sequencing performed is determined by checking whether the 0x01 flag bit is set ("template being sequenced has multiple fragments"). SAM/BAM files containing paired-end reads must be sorted by query-name using Picard's SortSam function.
# CpGs are identified by examining the XM-tag of each read. The XM-tag is reference-based, so only reference Cpgs are examined. This means that sample-specific CpGs are ignored and reference-specific CpGs, which may not be CpGs in the sample of interest, are also examined.
# This program brings together the functionality of SAM2MS.py, QS_PE_SAM2MS.py and summarise_aggregate_methylation.py.
# There are two types of comethylation: (1) Within-fragment co-methylation (WF), which is the dependence of methylation states along a DNA read/fragment/chromosome that contains multiple Cpgs; and (2) correlation of aggregate methylation values (AM), which is the correlation of aggregate methylation values at nearby CpGs. Aggregate methylation values are based on the pileup of reads at each CpG and include beta = M/(U+M) and gamma = log((M + eps)/(U + eps)), where M = the number of methylated Cs, U = the number of unmethylated Cs and esp is a small number to ensure numerical stability of gamma, typically esp = 0.5.
# comethlyation.py does not perform the statistical analysis of comethylation, rather it outputs two files for subsequent analysis: sample.wf (within-framgent co-methylation) and sample.am (aggregate methylation), where "sample" is the name of the sample being processed.
############################################################################################################################################################################################

### Description of sample.wf ###
############################################################################################################################################################################################
# Each line corresponds to a CpG-pair and includes a count of the number of methylated Cs (M), unmethylated Cs (U) and other bases (other) for that CpG.
# CpG-pairs can be created in several ways, see comethylation.py --help for details. Each pair is defined as (CpG1, CpG2), where CpG1 is the leftmost CpG in the pair and CpG2 is the rightmost CpG in the pair (with reference to Watson-strand coordinates).
# The data for each CpG-pair are stratified by which strand the reads aligned to: Original Top (OT) aka Watson-strand aka forward-strand aka "+"-strand; Original Bottom (OB) aka Crick-strand aka reverse-strand aka "-"-strand and combined across strands.
# The strand information is encoded as a subscript, e.g. MM_OT is the count of reads informative for the OT-strand that are methylated at both CpGs in the CpG-pair, MU_OB is the count of reads informative for the OB-strand that are methylated at CpG1 and unmethylated at CpG2, UU (no subscript) is the count of reads, combined across strands, that are unmethylated at both CpGs in the CpG-pair.

############################################################################################################################################################################################

### Description of sample.am ###
############################################################################################################################################################################################
# Each line corresponds to a single (reference) CpG.
# The data for each CpG are counts of the number of methylated Cs (M), unmethylated Cs (U) and other bases (other).
# The data for each CpG are stratified by which strand the read is informative for: Original Top (OT) aka Watson-strand aka forward-strand aka "+"-strand; Original Bottom (OB) aka Crick-strand aka reverse-strand aka "-"-strand and combined across strands.
############################################################################################################################################################################################

### Explanation of XM-tag (methylation string) ###
############################################################################################################################################################################################
# . for bases not involving cytosines                       
# X for methylated C in CHG context (was protected)         
# x for not methylated C in CHG context (was converted)     
# H for methylated C in CHH context (was protected)         
# h for not methylated C in CHH context (was converted)     
# Z for methylated C in CpG context (was protected)         
# z for not methylated C in CpG context (was converted)     
############################################################################################################################################################################################

### TODOs ###
############################################################################################################################################################################################
# TODO: Define trimOverlappingPair() t
# TODO: Insert program description in arg.parse
# TODO: Define AM function
############################################################################################################################################################################################


### Command line parser ###
############################################################################################################################################################################################
parser = argparse.ArgumentParser(description='INSERT PROGRAM DESCRIPTION')
parser.add_argument('BAM', metavar = 'BAM',
                    help='The path to the SAM/BAM file')
parser.add_argument('sampleName',
                    metavar = <string>,
                    help='The name of the sample. All output files will have this prefix.')
parser.add_argument('--ignoreDuplicates',
                    action='store_true',
                    help='Ignore reads that have been flagged as PCR duplicates by Picard\'s MarkDuplicates function')
parser.add_argument('-5', '--ignore5', metavar = '<int>',
                    type = int,
                    default=0,
                    help='NOT YET IMPLEMENTED. Ignore <int> bases from 5\' (left) end of reads (default: 0). WARNING: Parameter value not sanity checked by program.')
parser.add_argument('-3', '--ignore3', metavar = '<int>',
                    type = int,
                    default=0,
                    help='NOT YET IMPLEMENTED. Ignore <int> bases from 3\' (right) end of reads (default: 0). WARNING: Parameter value not sanity checked by program.')
parser.add_argument('--pairChoice',
                    metavar = '<string>',
                    default="outermost",
                    help='Method for constructing CpG-pairs: outermost, random or all (default: outermost)')

args = parser.parse_args()
############################################################################################################################################################################################

### Open SAM/BAM file and output files ###
############################################################################################################################################################################################
BAM = pysam.Samfile(args.BAM)
WF = open(".".join([sampleName, ".wf"]), "w")
am = open(".".join([sampleName, ".am"]), "w")

############################################################################################################################################################################################

### Variable initialisation ###
#############################################################################################################################################################################################
CpG_pattern = re.compile(r"[Zz]")
methylated_CpG_pattern = re.compile(r"[Z]")
unmethylated_CpG_pattern = re.compile(r"[z]")
bad_overlap_pairs = 0 # The number of read-pairs with a bad overlap. A read-pair with a bad overlap is one where the leftmost position of read2 (read1) is to the left of the leftmost position of read1 (read2) for reads informative for the OT (OB) strand. Such reads are ignored when computing within-fragment comethylation as they can violate the assumption that posL < posR.
fragment_counter = 0 # The number of DNA fragments processed by SAM2MS.py. Both one single-end read and one paired-end read contribute a single DNA fragment.
duplicate_fragment_counter = 0 # The number of duplicate fragments ignore by SAM2MS.py.  Both one single-end read and one paired-end read contribute a single DNA fragment.
CpG_pairs = {} # Dictionary of CpG pair IDs with keys of form chr_pos1_pos2 and values corresponding to a makeCount object 
#############################################################################################################################################################################################

### Function definitions ###
#############################################################################################################################################################################################
## Create a dictionary of methylation-states for a CpG-pair
def makeWFCount():
    return  {'MM': 0, 'MU': 0, 'UM': 0, 'UU': 0, 'MM_OT': 0, 'MU_OT': 0, 'UM_OT': 0, 'UU_OT': 0, 'MM_OB': 0, 'MU_OB': 0, 'UM_OB': 0, 'UU_OB': 0}

## Increment the counts in an makeWFCount() object based on the new information from MS_read, the output of SAM2MS_SE() or SAM2MS_PE().
def incrementCount(CpG_pair, fragment_MS):
    ms = ''.join(fragment_MS[3:5]) # The CpG methylation state - ZZ, Zz, zZ or zz - for that CpG-pair from that read
    strand = fragment_MS[5]
    if ms == 'ZZ':
        CpG_pair['MM'] += 1
        if strand == 'OT':
            CpG_pair['MM_OT'] += 1
        elif strand == 'OB':
            CpG_pair['MM_OB'] +=1
        else:
            exit_msg = ''.join(['Error: Invalid strand at line ', str(line)])
            sys.exit(exit_msg)
        return CpG_pair
    elif ms == 'Zz':
        CpG_pair['MU'] += 1
        if strand == 'OT':
            CpG_pair['MU_OT'] += 1
        elif strand == 'OB':
            CpG_pair['MU_OB'] +=1
        else:
            exit_msg = ''.join(['Error: Invalid strand at line ', str(line)])
            sys.exit(exit_msg)
        return CpG_pair
    elif ms == 'zZ':
        CpG_pair['UM'] += 1
        if strand == 'OT':
            CpG_pair['UM_OT'] += 1
        elif strand == 'OB':
            CpG_pair['UM_OB'] +=1
        else:
            exit_msg = ''.join(['Error: Invalid strand at line ', str(line)])
            sys.exit(exit_msg)
        return CpG_pair
    elif ms == 'zz':
        CpG_pair['UU'] += 1
        if strand == 'OT':
            CpG_pair['UU_OT'] += 1
        elif strand == 'OB':
            CpG_pair['UU_OB'] +=1
        else:
            exit_msg = ''.join(['Error: Invalid strand at line ', str(line)])
            sys.exit(exit_msg)
        return CpG_pair
    else:
        exit_msg = ''.join(['Error: Invalid methylation string at line ', str(line)])
        sys.exit(exit_msg)

## trimOverlappingPair() trims XM-tags of reads from an overlapping read-pair. Bases in the overlapping region constitute one data point that are measured twice (once by each read in the read-pair). We only want to count this data point once for the purposes of studying comethylation. We trim the read with the lower aggregate base qualities for the overlap region to remove the overlap. We also check that each read reports the same sequence in the overlap region; if the two reads give conflicting information we exclude the read-pair from the analysis
def trimOverlappingPair(read1XM, read2XM, absTLEN): # read1XM and read2XM are the XM-tags for read1 and read2, resp. absTLEN is the absolute value of the template length (TLEN).
    n_overlap = absTLEN - (len(read1XM) + len(read2XM))
    

## SAM2MS_SE extracts, summarises and returns the methylation string (MS) information for a read mapped as single-end data
def SAM2MS_SE(read):
    chrom = BAM.getrname(read.tid)
    # Store the XM tag for the read.
    readXM = [tag[1] for tag in read.tags if tag[0] == 'XM'][0]
    start = read.pos + 1 # read.pos is 0-based in accordance with BAM file specificiations, regardless of whether the input file is SAM or BAM. I convert this to a 1-based position.
    strand = getStrand(read, "NA")
    # Identify CpGs in read
    CpG_index = [m.start() for m in re.finditer(CpG_pattern, readXM)]
    # Case A: > 1 CpG in read
    if len(CpG_index) > 1:
        if pairChoice == 'outermost':
            CpGL = CpG_index[0] # Leftmost CpG in read
            CpGR = CpG_index[-1] # Rightmost CpG in read
            positionL = start + CpGL
            positionR = start + CpGR
            if read.is_reverse: # If a read maps in the reverse orientation the Z/z characters in the XM string point to the G in the CpG - I want to point to the C in the CpG on the OT-strand so I move the position coordinates 1bp to the left
                positionL -= 1
                positionR -= 1
            output=[chrom, positionL, positionR, readXM[CpGL], readXM[CpGR], strand] # Output is left-to-right, thus read1XM appears before read2XM for OT read-pairs.
            if positionL > positionR:
                print "ERROR: Case1A posL > posR for single-end read ", read1.qname," with output ", output
                return None
            return output
        else:
            sys.exit('ERROR: Only \'--pairChoice outermost\' is implemented.')
    else: # Less than 2 CpGs in read, therefore no within-fragment comethylation measurement for this read.
        return None

## SAM2MS_PE extracts, summarises and returns the methylation string (MS) information for a read-pair mapped as paired-end data.
def SAM2MS_PE(read1, read2):
    chrom = BAM.getrname(read1.tid)
    ignore2 = 0 # Initialise ignore2 to zero. ignore2 = the number of bases of read2 to ignore due to overlap between read1 and read2 and is computed later in this function.
    # Store the XM tag for each read. Tags are stored as a Python list of 2-tuples [("TAG_NAME", "TAG_VALUE"), ...]. This assumes there is one, and one only, XM tag per read
    read1XM = [tag[1] for tag in read1.tags if tag[0] == 'XM'][0] 
    read2XM = [tag[1] for tag in read2.tags if tag[0] == 'XM'][0]
    # Case1: Read-pair orientation is +/- (read1/read2)
    if (not read1.is_reverse and read2.is_reverse):
        # If read1.tlen < 0, then the read-pair has a bad overlap
        if read1.tlen < 0:
            print 'WARNING: Ignoring read-pair ', read1.qname, ' due to bad overlap (=', read1.tlen, ')'
            return None
        # start1 <= start2 by definition for a properly paired OT read
        start1 = read1.pos + 1 # read1.pos is 0-based in accordance with BAM file specificiations, regardless of whether the input file is SAM or BAM. I convert this to a 1-based position.
        start2 = read2.pos + 1 # read2.pos is 0-based in accordance with BAM file specificiations, regardless of whether the input file is SAM or BAM. I convert this to a 1-based position.
        strand = getStrand(read1, read2)
        # If the two reads overlap, i.e the rightmost position of read1 is greater than or equal to the leftmost position of read2, drop the offending positions from the leftmost (3') position of read2
        if read1.positions[-1] >= read2.positions[0]: 
            ignore2 = int(read1.positions[-1] - read2.positions[0] + 1) # In this case, ignore2 is the number of bases to ignore from the 3' end of read2.
            read2XM = read2XM[ignore2:] # Drop the ignored positions from read2's XM tag
        # Identify CpGs in read1 and read2
        CpG_index1 = [m.start() for m in re.finditer(CpG_pattern, read1XM)]
        CpG_index2 = [m.start() for m in re.finditer(CpG_pattern, read2XM)]
        # CaseA: > 0 CpGs in both reads of pair    
        if (len(CpG_index1) > 0 and len(CpG_index2) > 0):
            if pairChoice == "outermost":
                CpGL = CpG_index1[0] # Leftmost CpG in read1
                CpGR = CpG_index2[-1] # Rightmost CpG in read2
                positionL = start1 + CpGL
                positionR = start2 + CpGR + ignore2 # +ignore2 to deal with ignored bases
                output=[chrom, positionL, positionR, read1XM[CpGL], read2XM[CpGR], strand] # Output is left-to-right, thus read1XM appears before read2XM for OT read-pairs.
                if positionL > positionR:
                    print "ERROR: Case1A posL > posR for read-pair ", read1.qname," with output ", output
                    return None
                else:
                    return output
        # CaseB: > 1 CpG in read1 but 0 CpGs in read2
        elif (len(CpG_index1) > 1 and len(CpG_index2) == 0):
            if pairChoice == "outermost":
                CpGL = CpG_index1[0]
                CpGR = CpG_index1[-1]
                positionL = start1 + CpGL
                positionR = start1 + CpGR
                output=[chrom, positionL, positionR, read1XM[CpGL], read1XM[CpGR], strand]
                if positionL > positionR:
                    print "ERROR: Case1B posL > posR for read ", read1.qname," with output ", output
                    return None
                else:
                    return output
        # CaseC: 0 CpGs in read1 but > 1 CpG in read2
        elif (len(CpG_index1) == 0 and len(CpG_index2) > 1):
            if pairChoice == "outermost":
                CpGL = CpG_index2[0]
                CpGR = CpG_index2[-1]
                positionL = start2 + CpGL
                positionR = start2 + CpGR
                output=[chrom, positionL, positionR, read2XM[CpGL], read2XM[CpGR], strand]
                if positionL > positionR:
                    print "ERROR: Case1C posL > posR for read ", read2.qname," with output ", output
                    return None
                else:
                    return output
        # CaseD: < 2 CpGs in read-pair
        elif (len(CpG_index1) + len(CpG_index2)) < 2:
            return None
        # CaseE: Unexpected CpG count in read-pair
        else:
            print 'ERROR: Unexpected CpG count for read-pair ', read1.qname
            return None
        # Case2: Read-pair orientation is -/+ (read1/read2)
    elif read1.is_reverse and not read2.is_reverse:
        # If read2.tlen < 0, then the read-pair has a bad overlap
        if read2.tlen < 0:
            print 'WARNING: Ignoring read-pair ', read2.qname, ' due to bad overlap (=', read2.tlen, ')'
            return None
        # start1 >= start2 by definition for a properly paired OB read
        start1 = read1.pos + 1 # read1.pos is 0-based in accordance with BAM file specificiations, regardless of whether the input file is SAM or BAM. I convert this to a 1-based position.
        start2 = read2.pos + 1 # read2.pos is 0-based in accordance with BAM file specificiations, regardless of whether the input file is SAM or BAM. I convert this to a 1-based position.
        strand = getStrand(read1, read2)
        # If the two reads overlap, i.e the rightmost position of read2 is greater than or equal to the leftmost position of read1, drop the offending positions from the rightmost (3') position of read2
        if read2.positions[-1] >= read1.positions[0]:
            ignore2 = len(read2XM) - int(read2.positions[-1] - read1.positions[0] + 1) # In this case, ignore2 is the number of bases to retain from the 5' end of read2.
            read2XM = read2XM[:ignore2]
        # Identify CpGs in read1 and read2
        CpG_index1 = [m.start() for m in re.finditer(CpG_pattern, read1XM)]
        CpG_index2 = [m.start() for m in re.finditer(CpG_pattern, read2XM)]
        # CaseA: A CpG in both reads of pair
        if (len(CpG_index1) > 0 and len(CpG_index2) > 0):
            if pairChoice == "outermost":
                CpGL = CpG_index2[0] # Leftmost CpG in read2
                CpGR = CpG_index1[-1] # Rightmost CpG in read1
                positionL = start2 + CpGL - 1 # No need to ignore2 since 3' bases of read2 are rightmost; -1 so as to point to C on the forward strand in the CpG 
                positionR = start1 + CpGR - 1 # -1 so as to point to C on the forward strand in the CpG
                output=[chrom, positionL, positionR, read2XM[CpGL], read1XM[CpGR], strand] # Output is left-to-right, thus read2XM appears before read1XM for OB read-pairs.
                if positionL > positionR:
                    print "ERROR: Case2A posL > posR for read-pair ", read1.qname," with output ", output
                    return None
                else:
                    return output
        # CaseB: > 1 CpG in read1 but 0 CpGs in read2
        elif (len(CpG_index1) > 1 and len(CpG_index2) == 0):
            if pairChoice == "outermost":
                CpGL = CpG_index1[0] # Leftmost CpG in read1
                CpGR = CpG_index1[-1] # Rightmost CpG in read1
                positionL = start1 + CpGL - 1 # -1 so as to point to C on the forward strand in the CpG 
                positionR = start1 + CpGR - 1 # -1 so as to point to C on the forward strand in the CpG
                output=[chrom, positionL, positionR, read1XM[CpGL], read1XM[CpGR], strand]
                if positionL > positionR:
                    print "ERROR: Case2B posL > posR for read ", read1.qname," with output ", output
                    return None
                else:
                    return output
        # CaseC: 0 CpGs in read1 but > 1 CpG in read2
        elif (len(CpG_index1) == 0 and len(CpG_index2) > 1):
            if pairChoice == "outermost":
                CpGL = CpG_index2[0] # Leftmost CpG in read2
                CpGR = CpG_index2[-1] # Rightmost CpG in read2
                positionL = start2 + CpGL - 1 # -1 so as to point to C on the forward strand in the CpG 
                positionR = start2 + CpGR - 1 # -1 so as to point to C on the forward strand in the CpG
                output=[chrom, positionL, positionR, read2XM[CpGL], read2XM[CpGR], strand]
                if positionL > positionR:
                    print "ERROR: Case2C posL > posR for read ", read2.qname," with output ", output
                    return None
                else:
                    return output
        # CaseD: < 2 CpGs in read-pair
        elif (len(CpG_index1) + len(CpG_index2)) < 2:
            return None
        else:
            print 'ERROR: Unexpected CpG count for read-pair', read1.qname
            return None
    # Case3: read-pair in unexpected orientation
    else:
        print 'ERROR: Unexpected orientation of read-pair', read1.qname
        return None

## getStrand() returns whether the read or read-pair is informative for the OT (original top, Watson) or OB (original bottom, Crick) strand.
def getStrand(read1, read2):
    if read2 == "NA": # Single-end data
        if read1.opt('XG') == 'CT': # Read aligned to CT-converted reference genome and therefore informative for the OT-strand
            strand = 'OT'
        elif read1.opt('XG') == 'GA': # Read aligned to GA-converted reference genome and therefore informative for the OB-strand
            strand = 'OB'
    else: # Paired-end data
        if read1.opt('XG') == 'CT' and read2.opt('XG') == 'CT': # Read-pair aligned to CT-converted reference genome and therefore informative for the OT-strand
            strand = 'OT'
        elif read1.opt('XG') == 'GA' and read2.opt('XG') == 'GA' : # Read-pair aligned to GA-converted reference genome and therefore informative for the OB-strand
            strand = 'OB'
        else:
            error_message = ''.join(['Read-pair is aligned to both the CT- and GA-converted reference genomes: ', read1.qname]) 
            sys.exit(error_message)
    return strand
## WFWriter writes a tab-separated output file to the filehandle WF
WFWriter = csv.writer(WF, delimiter='\t', quotechar=' ', quoting=csv.QUOTE_MINIMAL)
## AMWriter writes a tab-separated output file to the filehandle WF
AMWriter = csv.writer(AM, delimiter='\t', quotechar=' ', quoting=csv.QUOTE_MINIMAL)
## writeWF() writes the CpG_pairs object to disk as a tab-separated file.
def writeWF(CpG_pairs):
    # Create the header row
    header = ['chr', 'pos1', 'pos2', 'distance']
    example_row = CpG_pairs[CpG_pairs.keys()[1]]
    for key in sorted(example_row.iterkeys()):
        header.append(key)
    # Write the header to file
    WFWriter.writerow(header)
    # Write each CpG-pair to file
    for pair in CpG_pairs.iterkeys():
        pair_counts = []
	for count in sorted(CpG_pairs[pair].iterkeys()):
		pair_counts.append(CpG_pairs[pair][count])
        chrom = pair.rsplit(':')[0]
        pos1 = pair.rsplit(':')[1].rsplit('-')[0]
        pos2 = pair.rsplit(':')[1].rsplit('-')[1]
        dist = int(pos2) - int(pos1)
        row = [chrom, pos1, pos2, dist] + pair_counts
        WFWriter.writerow(row)
    
#############################################################################################################################################################################################
### The main program. Loops over the BAM file line-by-line (i.e. alignedRead-by-alignedRead) and extracts the XM information for each read or read-pair. ###
##############################################################################################################################################################################################
for read in BAM:
    fragment_MS = None # Reset the fragment_MS to None to erase values from previous reads/read-pairs
    if args.ignoreDuplicates and read.is_duplicate:
        continue
    else:
        if read.is_read1 and read.is_proper_pair:
            read1 = read
            continue
        elif read.is_read2 and read.is_proper_pair: 
            read2 = read
            # Check that read1 and read2 are aligned to the same chromosome and have identical read-names. If not, skip the read-pair.
            if read1.tid == read2.tid and read1.qname == read2.qname:
                fragment_MS = SAM2MS_PE(read1, read2)
            elif read1.tid != read2.tid:
                print "ERROR: Reads in pair aligned to different chromosomes: ", read1.qname, BAM.getrname(read1.tid), read2.qname, BAM.getrname(read2.tid)
                continue
            elif read1.qname != read2.qname:
                print "ERROR: The name of read1 is not identical to to that of read2 for read-pair ", read1.qname, read2.qname, "\nHas your BAM file been sorted in query-name-order with Picard's SortSam function?"
                break
            elif read.is_paired and not read.is_proper_pair:
                print "WARNING: Skipping read ", read.qname, " as it is part of an improper pair."
                continue
        elif not read.is_paired:  
            fragment_MS = SAM2MS_SE(read)
        else:
            print "Read is neither a single-end read nor part of a paired-end read. Check the SAM flag values are correctly set for read:", read.qname
            continue
    if fragment_MS is not None: 
        pair_ID = ''.join([fragment_MS[0], ':', str(fragment_MS[1]), '-', str(fragment_MS[2])])
        if not pair_ID in CpG_pairs: # CpG-pair not yet seen - create a dictionary key for it and increment its count (value)
            CpG_pairs[pair_ID] = makeWFCount()
            CpG_pairs[pair_ID] = incrementCount(CpG_pairs[pair_ID], fragment_MS) 
        else: # CpG-pair already seen - increment its count (value)
            CpG_pairs[pair_ID] = incrementCount(CpG_pairs[pair_ID], fragment_MS)

writeWF(CpG_pairs)
BAM.close()
WF.close()
AM.close()
##############################################################################################################################################################################################

# Print WF output (incomplete)
header = ['chrID']
for key in sorted(CpG_pairs.iterkeys()):
	for key2 in sorted(CpG_pairs[key].iterkeys()):
		header.append(key2)
		
vals = []
for key in sorted(CpG_pairs.iterkeys()):
	vals.append(key)
	for key2 in sorted(CpG_pairs[key].iterkeys()):
		vals.append(CpG_pairs[key][key2])

print header, '\n', vals
