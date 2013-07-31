#!/usr/bin/python

#### TODOs ####
# TODO: Make maxReadLength a (hidden) parameter option

#### Program information ####
# Peter Hickey
# 22/06/2012
# Program to check all reads of a BS-seq experiment for biases in methylation calls by read position. Such biases could be indicative of sequencing errors or problems with the bisulfite conversion. Input BAM is required to have the XM-tag (i.e. the methylation string tag produced by Bismark or my script XM_tag.py

#### Import required libraries ####
import argparse
import pysam
import re
import numpy

#### Explanation of XM-tag (methylation string) ####
############################################################################################################################################################################################
# . for bases not involving cytosines                       
# X for methylated C in CHG context (was protected)         
# x for not methylated C in CHG context (was converted)     
# H for methylated C in CHH context (was protected)         
# h for not methylated C in CHH context (was converted)     
# Z for methylated C in CpG context (was protected)         
# z for not methylated C in CpG context (was converted)
# U for methylated C in "unknown" context, e.g. CNN-context, (was protected). Non-standard Bismark XM-tag values; unique to output of XM_tag.py.
# u for not methylated C in "unknown"-context, e.g. CNN-context, (was converted). Non-standard Bismark XM-tag values; unique to output of XM_tag.py.
############################################################################################################################################################################################

#### Command line parser ####
parser = argparse.ArgumentParser(description='Check for read-position biases in methylation calls.\nWARNING: Assume maximum read-length is less than 200bp.\nWARNING: This tool is not "duplicate-read" aware.\nWARNING: Assumes there are no reads with length equal to zero.')
parser.add_argument('infile', metavar = 'in.bam',
                  help='The path to the SAM/BAM file')
parser.add_argument('sampleName', metavar = 'sampleName',
                  help='A sample ID. This is used as a prefix for all output files')
parser.add_argument('--pairedEnd',
                    action='store_true',
                    help='The data are paired-end')
args = parser.parse_args()


#### Open files ####
# SAM/BAM file to be processed
IN = pysam.Samfile(args.infile)

#### Variable initialisation ####
maxReadLength = 200 # Fairly arbitrary maximum read length. All lists are initialised to this length and unused elements are removed prior to writing the lists to file.
readLengths_OT_1 = numpy.array([0] * maxReadLength) # NB: readLengths_OT_1[i] is the count of OT-strand read1s with read-length == (i+1)
readLengths_OB_1 = numpy.array([0] * maxReadLength) # NB: readLengths_OB_1[i] is the count of OB-strand read1s with read-length == (i+1)
readLengths_OT_2 = numpy.array([0] * maxReadLength) # NB: readLengths_OT_2[i] is the count of OT-strand read2s with read-length == (i+1)
readLengths_OB_2 = numpy.array([0] * maxReadLength) # NB: readLengths_OB_2[i] is the count of OB-strand read2s with read-length == (i+1)
# List of all possible methylation calls
MS_states = ['X', 'x', 'H', 'h', 'Z', 'z', 'U', 'u']
# Lists to store counts of methylation calls at each position; e.g. MS_OT_1['X'][5] is the count of "X" calls at the position X on read1s aligned to the OT-strand
MS_OT_1 = {}
MS_OB_1 = {}
MS_OT_2 = {}
MS_OB_2 = {}
for state in MS_states:
    MS_OT_1[state] = numpy.array([0] * 200)
    MS_OB_1[state] = numpy.array([0] * 200)
    MS_OT_2[state] = numpy.array([0] * 200)
    MS_OB_2[state] = numpy.array([0] * 200)

#### Main loop - two different loops depending on whether data is single-end or paired-end ####
i = 1
if not args.pairedEnd:
    for read in IN:
        if i % 1000000 == 0:
            print 'Processed', i, 'reads...'
        if read.opt('XG') == 'CT':
            readLengths_OT_1[read.qlen - 1] += 1
            for state in MS_states:
                ms_index = [m.start() for m in re.finditer(state, read.opt('XM'))] # OT-strand read1s are ordered 5'-to-3' left-to-right
                MS_OT_1[state][ms_index] += 1 # Increment the count
        elif read.opt('XG') == 'GA':
            readLengths_OB_1[read.qlen - 1] += 1
            for state in MS_states:
                ms_index = [m.start() for m in re.finditer(state, read.opt('XM'))] # OB-strand read1s are ordered 3'-to-5' left-to-right
                ms_index = list(read.qlen - numpy.array(ms_index) - 1) # Convert the index to be 5'-to-3' left-to-right (-1 is because python arrays are 0-based)
                MS_OB_1[state][ms_index] += 1 # Increment the count
        else:
            print 'Read skipped: Undefined strand (missing XG-tag) for read', read.qname
            continue
        i += 1
#### Write output to file ####
# Trim trailing zeros from readLength arrays
    readLengths_OT_1 = numpy.trim_zeros(readLengths_OT_1, 'b')
    readLengths_OB_1 = numpy.trim_zeros(readLengths_OB_1, 'b')

# Compute how many reads are at least x bp long and convert result to float
    readLengthsCumSum_OT_1 = numpy.cumsum(readLengths_OT_1[::-1])[::-1]
    readLengthsCumSum_OB_1 = numpy.cumsum(readLengths_OB_1[::-1])[::-1]

# Write all output to file (trimming trailing zeros as necessary); output to be plotted using R
    numpy.savetxt(''.join([args.sampleName, '.RLOT1']), readLengthsCumSum_OT_1, delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.RLOB1']), readLengthsCumSum_OB_1, delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.XOT1']), MS_OT_1['X'][:len(readLengthsCumSum_OT_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.xOT1']), MS_OT_1['x'][:len(readLengthsCumSum_OT_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.HOT1']), MS_OT_1['H'][:len(readLengthsCumSum_OT_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.hOT1']), MS_OT_1['h'][:len(readLengthsCumSum_OT_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.ZOT1']), MS_OT_1['Z'][:len(readLengthsCumSum_OT_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.zOT1']), MS_OT_1['z'][:len(readLengthsCumSum_OT_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.UOT1']), MS_OT_1['U'][:len(readLengthsCumSum_OT_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.uOT1']), MS_OT_1['u'][:len(readLengthsCumSum_OT_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.XOB1']), MS_OB_1['X'][:len(readLengthsCumSum_OB_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.xOB1']), MS_OB_1['x'][:len(readLengthsCumSum_OB_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.HOB1']), MS_OB_1['H'][:len(readLengthsCumSum_OB_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.hOB1']), MS_OB_1['h'][:len(readLengthsCumSum_OB_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.ZOB1']), MS_OB_1['Z'][:len(readLengthsCumSum_OB_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.zOB1']), MS_OB_1['z'][:len(readLengthsCumSum_OB_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.UOB1']), MS_OB_1['U'][:len(readLengthsCumSum_OB_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.uOB1']), MS_OB_1['u'][:len(readLengthsCumSum_OB_1)], delimiter = '\t', fmt='%-20d')
    
elif args.pairedEnd:
    for read in IN:
        if i % 1000000 == 0:
            print 'Processed', i, 'reads...'
        if read.opt('XG') == 'CT' and read.is_read1:
            readLengths_OT_1[read.qlen - 1] += 1
            for state in MS_states:
                ms_index = [m.start() for m in re.finditer(state, read.opt('XM'))] # OT-strand read1s are ordered 5'-to-3' left-to-right
                MS_OT_1[state][ms_index] += 1 # Increment the count
        elif read.opt('XG') == 'CT' and read.is_read2:
            readLengths_OT_2[read.qlen - 1] += 1
            for state in MS_states:
                ms_index = [m.start() for m in re.finditer(state, read.opt('XM'))] # OT-strand read2s are ordered 3'-to-5' left-to-right
                ms_index = list(read.qlen - numpy.array(ms_index) - 1) # Convert the index to be 5'-to-3' left-to-right (-1 is because python arrays are 0-based)
                MS_OT_2[state][ms_index] += 1 # Increment the count
        elif read.opt('XG') == 'GA' and read.is_read1:
            readLengths_OB_1[read.qlen - 1] += 1
            for state in MS_states:
                ms_index = [m.start() for m in re.finditer(state, read.opt('XM'))] # OB-strand read1s are ordered 3'-to-5' left-to-right
                ms_index = list(read.qlen - numpy.array(ms_index) - 1) # Convert the index to be 5'-to-3' left-to-right (-1 is because python arrays are 0-based)
                MS_OB_1[state][ms_index] += 1 # Increment the count
        elif read.opt('XG') == 'GA' and read.is_read2:
            readLengths_OB_2[read.qlen - 1] += 1
            for state in MS_states:
                ms_index = [m.start() for m in re.finditer(state, read.opt('XM'))] # OB-strand read2s are ordered 5'-to-3' left-to-right
                MS_OB_2[state][ms_index] += 1 # Increment the count
        else:
            print 'Read skipped: Undefined strand (missing XG-tag) or ill-defined read1/read2 FLAG for read', read.qname
            continue
        i += 1
#### Write output to file ####
# Trim trailing zeros from readLength arrays
    readLengths_OT_1 = numpy.trim_zeros(readLengths_OT_1, 'b')
    readLengths_OB_1 = numpy.trim_zeros(readLengths_OB_1, 'b')
    readLengths_OT_2 = numpy.trim_zeros(readLengths_OT_2, 'b')
    readLengths_OB_2 = numpy.trim_zeros(readLengths_OB_2, 'b')
    
# Compute how many reads are at least x bp long (0 < x <= length_longestRead) and convert result to float; e.g. readLengthsCumSum_OT_1[i] is the count of all reads at least (i+1) nt long
    readLengthsCumSum_OT_1 = numpy.cumsum(readLengths_OT_1[::-1])[::-1]
    readLengthsCumSum_OB_1 = numpy.cumsum(readLengths_OB_1[::-1])[::-1]
    readLengthsCumSum_OT_2 = numpy.cumsum(readLengths_OT_2[::-1])[::-1]
    readLengthsCumSum_OB_2 = numpy.cumsum(readLengths_OB_2[::-1])[::-1]
    
# Write all output to file (trimming trailing zeros as necessary); output to be plotted using R
    numpy.savetxt(''.join([args.sampleName, '.RLOT1']), readLengthsCumSum_OT_1, delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.RLOB1']), readLengthsCumSum_OB_1, delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.XOT1']), MS_OT_1['X'][:len(readLengthsCumSum_OT_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.xOT1']), MS_OT_1['x'][:len(readLengthsCumSum_OT_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.HOT1']), MS_OT_1['H'][:len(readLengthsCumSum_OT_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.hOT1']), MS_OT_1['h'][:len(readLengthsCumSum_OT_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.ZOT1']), MS_OT_1['Z'][:len(readLengthsCumSum_OT_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.zOT1']), MS_OT_1['z'][:len(readLengthsCumSum_OT_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.UOT1']), MS_OT_1['U'][:len(readLengthsCumSum_OT_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.uOT1']), MS_OT_1['u'][:len(readLengthsCumSum_OT_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.XOB1']), MS_OB_1['X'][:len(readLengthsCumSum_OB_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.xOB1']), MS_OB_1['x'][:len(readLengthsCumSum_OB_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.HOB1']), MS_OB_1['H'][:len(readLengthsCumSum_OB_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.hOB1']), MS_OB_1['h'][:len(readLengthsCumSum_OB_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.ZOB1']), MS_OB_1['Z'][:len(readLengthsCumSum_OB_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.zOB1']), MS_OB_1['z'][:len(readLengthsCumSum_OB_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.UOB1']), MS_OB_1['U'][:len(readLengthsCumSum_OB_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.uOB1']), MS_OB_1['u'][:len(readLengthsCumSum_OB_1)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.RLOT2']), readLengthsCumSum_OT_2, delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.RLOB2']), readLengthsCumSum_OB_2, delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.XOT2']), MS_OT_2['X'][:len(readLengthsCumSum_OT_2)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.xOT2']), MS_OT_2['x'][:len(readLengthsCumSum_OT_2)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.HOT2']), MS_OT_2['H'][:len(readLengthsCumSum_OT_2)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.hOT2']), MS_OT_2['h'][:len(readLengthsCumSum_OT_2)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.ZOT2']), MS_OT_2['Z'][:len(readLengthsCumSum_OT_2)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.zOT2']), MS_OT_2['z'][:len(readLengthsCumSum_OT_2)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.UOT2']), MS_OT_2['U'][:len(readLengthsCumSum_OT_2)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.uOT2']), MS_OT_2['u'][:len(readLengthsCumSum_OT_2)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.XOB2']), MS_OB_2['X'][:len(readLengthsCumSum_OB_2)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.xOB2']), MS_OB_2['x'][:len(readLengthsCumSum_OB_2)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.HOB2']), MS_OB_2['H'][:len(readLengthsCumSum_OB_2)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.hOB2']), MS_OB_2['h'][:len(readLengthsCumSum_OB_2)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.ZOB2']), MS_OB_2['Z'][:len(readLengthsCumSum_OB_2)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.zOB2']), MS_OB_2['z'][:len(readLengthsCumSum_OB_2)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.UOB2']), MS_OB_2['U'][:len(readLengthsCumSum_OB_2)], delimiter = '\t', fmt='%-20d')
    numpy.savetxt(''.join([args.sampleName, '.uOB2']), MS_OB_2['u'][:len(readLengthsCumSum_OB_2)], delimiter = '\t', fmt='%-20d')
    
# Close BAM file
IN.close()
