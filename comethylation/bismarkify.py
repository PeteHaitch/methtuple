#!/usr/bin/env python
import argparse
import sys
import csv
import pysam
import warnings

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

### Explanation of XR-tag and XG-tag ###
############################################################################################################################################################################################
#   XR    XG    strand
#   CT    CT    CT
#   GA    GA    CTOB
#   GA    CT    CTOT
#   CT    GA    OB
############################################################################################################################################################################################

#### Command line passer ####
parser = argparse.ArgumentParser(description='Adds a methylation call string tag (XM), read conversion state tag (XR) and genome converstion state tag (XG) to each read of a bisulfite-sequencing SAM/BAM file and writes the updated version to a BAM file. These custom tags are based on the Bismark bilsufite mapper (http://www.bioinformatics.babraham.ac.uk/projects/download.html#bismark) and are defined in the Bismark user manual. Unmapped reads in the BAM file have their XM, XR and XG tags set to *')
parser.add_argument('infile', metavar = 'in.bam',
                  help='The path to the original SAM/BAM file.')
parser.add_argument('outfile', metavar = 'out.bam',
                  help='The path to the new SAM/BAM file.')
parser.add_argument('--reference', metavar = 'reference.fa',
                  help='The path to the reference genome FASTA file. Must be indexed.')
parser.add_argument('--aligner', metavar = '<string>',
                  help='The aligner used to create the input SAM/BAM file.')
args = parser.parse_args()

#### Open files ####
# SAM/BAM file to be processed
IN = pysam.Samfile(args.infile)
# New BAM file to be created.Add a @PG tag for bismarkify.py to the new BAM file.
header = IN.header
id = 'bismarkify.py'
vn = '0.1'
cl = ' '.join(['python bismarkify.py',  args.reference, args.aligner, args.infile, args.outfile])
XM_tag_PG = {'ID': id, 'VN': vn, 'CL': cl}
header['PG'].append(XM_tag_PG) # FIXME: Assumes that a @PG tag already exists
OUT = pysam.Samfile(args.outfile, "wb", header = header)
# Reference FASTA file. Required for looking up sequence context of cytosines.
REF = pysam.Fastafile(args.reference)

#### Function definitions ####
def get_padded_ref_seq(read, ref, IN): # 'read' is an AlignedRead object, 'ref' is a FastaFile object
    """Extract reference sequence +/- 2bp of the location of a mapped read. The various cases are to ensure that all possible cases are considered, such as when the read maps to the start (end) of a contig and therefore there are no prior (subsequent) 2 bases to extract. The result is padded with 'N' characters in the cases where there is no reference sequence at the start or end positions. In its current form this function cannot handle a read that is longer than its mapped-contig at each end, e.g. a 100bp read containing a 98bp contig internally, but these shouldn't exist and would likely not be mapped even if such short contigs did exist.
    
    Args:
        read: A pysam.AlignedRead instance.
        ref: A pysam.Fastafile instance.
        IN: A pysam.Samfile instance.

    Returns:
        A Python string containing the relevant reference sequence.
         
    """   
    try:
        refseq = ref.fetch(IN.getrname(read.tid), start = (read.pos - 2), end = (read.aend + 2)).upper() # FIXME: Process CIGAR so that correct start and end are obtained
    except ValueError:
        try:
            refseq = ''.join(['N', ref.fetch(IN.getrname(read.tid), start = (read.pos - 1), end = (read.aend + 2)).upper()]) # FIXME: Process CIGAR so that correct start and end are obtained
        except ValueError:
            try: 
                refseq = ''.join(['NN', ref.fetch(IN.getrname(read.tid), start = read.pos , end = (read.aend + 2)).upper()]) # FIXME: Process CIGAR so that correct start and end are obtained
            except ValueError:
                try:
                    refseq = ''.join([ref.fetch(IN.getrname(read.tid), start = (read.pos - 2) , end = (read.aend + 1)).upper(), 'N']) # FIXME: Process CIGAR so that correct start and end are obtained
                except ValueError:
                    try:
                        refseq = ''.join([ref.fetch(IN.getrname(read.tid), start = (read.pos - 2) , end = read.aend).upper(), 'NN']) # FIXME: Process CIGAR so that correct start and end are obtained
                    except ValueError:
                        warning_msg "Skipped read: Unable to extract reference sequence for read", read.qname
                        warnings.warn(warning_msg)
                        refseq = ''.join('N' for i in xrange(read.qlen))
    return refseq

# TODO: makeXMtag (now, make_XM_tag) needs to be completely re-written
# Get the genomic-context for every C in the read and create the XM tag for that read
# If there is a C in the read, look ahead 1 or 2 bases to determine the genomic-context, i.e. CG, CHG or CHH.
# If the reference position is a C but the mapped base is something other than a C or a T then the '.' character is recorded in the XM tag for that position
def make_XM_tag(read, refseq, strand): # read is an AlignedRead object, refseq is the output of getPaddedRefSeq, strand is a character "+" or "-" encoding which strand the read aligned to
    XM = []
    if strand == "+":
        for i in range(read.qlen): # i indexes position along the read  # TODO: Make i index more general positions in genome, e.g. i = [0, 1, 3] corresponding to a CIGAR = 2M1D1M (should work for deletions) what about insertions?
            j = i+2 # j indexes position along the padded-reference genome and is 2 ahead of i because of the padding
            if(refseq[j] == 'C'): # Cytosine on forward strand of reference sequence
                if refseq[j+1] == 'G': # CG-event
                    if read.seq[i] == 'C':
                        XM.append('Z')
                    elif read.seq[i] == 'T':
                        XM.append('z')
                    else:
                        XM.append('.')
                elif refseq[j+2] == 'G': # CHG-event
                    if read.seq[i] == 'C':
                        XM.append('X')
                    elif read.seq[i] == 'T':
                        XM.append('x')
                    else:
                        XM.append('.')
                elif refseq[j+2] != 'N': # CHH-event
                    if read.seq[i] == 'C':
                        XM.append('H')
                    elif read.seq[i] == 'T':
                        XM.append('h')
                    else:
                        XM.append('.')
                else: # Unable to determine the genomic-context of the cytosines in read - not CG-, CHG- or CHH-event, e.g. a CNN-event  
                    if read.seq[i] == 'C':
                        XM.append('U')
                    elif read.seq[i] == 'T':
                        XM.append('u')
                    else:
                        XM.append('.')
            else: # Not a cytosine
                XM.append('.')
    elif strand == '-':
        for i in range(read.qlen): # i indexes position along the read
            j = i+2 # j indexes position along the padded-reference genome and is 2 ahead of i because of the padding
            if(refseq[j] == 'G'): # Guanine on forward strand of reference sequence (therefore a cytosine on reverse strand)
                if refseq[j-1] == 'C': # CG-event
                    if read.seq[i] == 'G':
                        XM.append('Z')
                    elif read.seq[i] == 'A':
                        XM.append('z')
                    else:
                        XM.append('.')
                elif refseq[j-2] == 'C': # CHG-event
                    if read.seq[i] == 'G':
                        XM.append('X')
                    elif read.seq[i] == 'A':
                        XM.append('x')
                    else:
                        XM.append('.')
                elif refseq[j-2] != 'N': # CHH-event
                    if read.seq[i] == 'G':
                        XM.append('H')
                    elif read.seq[i] == 'A':
                        XM.append('h')
                    else:
                        XM.append('.')
                else: # Unable to determine the genomic-context of the guanine in read - not CG-, CHG- or CHH-event, e.g. a CNN-event  
                    if read.seq[i] == 'G':
                        XM.append('U')
                    elif read.seq[i] == 'A':
                        XM.append('u')
                    else:
                        XM.append('.')
            else: # Not a guanine
                XM.append('.')
                
    XM = ''.join(XM)
    return XM

# TODO: Write functions to create XR and XG tags based on the aligner used
                    
# TODO: Rewrite main loop. Keep all reads in memory (or in batches with maximum size specified by a command line argument) and then write to disk at the end instead of the current approach of writing each read to disk one-at-a-time. 
# TODO: Add counters for mCpG, CpG, mCHG, CHG, mCHH, CHH, reads that couldn't be processed, progress report. See what else Bismark reports in its mapping summary
# Main loop - create XM, XR and XG tags for each read and write the new BAM to file
for read in IN:
    if not read.is_unmapped:
        refseq = getPaddedRefSeq(read, REF, IN)
        if read.opt('XG') == 'CT':
            strand = "+"
        elif read.opt('XG') == 'GA':
            strand = "-"
        else:
            print 'Read skipped: Undefined strand (missing XG-tag) for read', read.qname
            continue
        XM = makeXMtag(read, refseq, strand)
    else:
        XM = '*'
    read.tags = read.tags + [('XM', XM)]
    OUT.write(read)

IN.close()
OUT.close()
REF.close()
