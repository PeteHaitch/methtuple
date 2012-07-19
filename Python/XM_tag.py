## TODO: Update @PG tag when making major changes to code.
## TODO: Add counters for mCpG, CpG, mCHG, CHG, mCHH, CHH, reads that couldn't be processed, progress report. See what else Bismark reports in its mapping summary

#!/usr/bin/python
import argparse
import sys
import csv
import pysam

## IMPORTANT: This script will not work properly with reads that have been soft-clipped, e.g. reads from Bowtie2. Bismark, when using Bowtie2, only allows end-to-end alignments, i.e. untrimmed and unclipped alignments, and thus avoid the difficulty caused by soft-clipped alignments.
## IMPORTANT: This script assumes there is no soft-clipping of bases since it uses the AlignedRead.seq to extract the read sequence rather than AlignedRead.query. The script may not work if AlignedRead.seq is simply changed to AlignedRead.query
## IMPORTANT: Reads aligned to the OB-strand have the XM-tag pointing to the C on the OB-strand, therefore we need to -1 (for CpGs, -2 for CHHs and CHGs) for positions to translate these to the Cs on the OT-strand.
## WARNING: Requires that the input BAM already contains a 'PG' header tag

### Explanation of XM-tag (methylation string) ###
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

# Command line passer
parser = argparse.ArgumentParser(description='Add a methylation tag XM to each read of a SAM/BAM file and write the result to BAM file. The XM tag specification is defined in the Bismark user manual.')
parser.add_argument('infile', metavar = 'in.bam',
                  help='The path to the original SAM/BAM file')
parser.add_argument('outfile', metavar = 'out.bam',
                  help='The path to the new SAM/BAM file')
parser.add_argument('reference', metavar = 'reference.fa',
                  help='The path to the reference genome FASTA file. Must be indexed.')
args = parser.parse_args()

# Open files
# SAM/BAM file to be processed
IN = pysam.Samfile(args.infile)
# New BAM file to be created. Need to make header that includes a new @PG tag for the new BAM file.
header = IN.header
id = 'XM_tag.py'
vn = '0.3'
cl = ' '.join(['python XM_tag.py', args.infile, args.outfile, args.reference])
XM_tag_PG = {'ID': id, 'VN': vn, 'CL': cl}
header['PG'].append(XM_tag_PG)
OUT = pysam.Samfile(args.outfile, "wb", header = header)
# Reference FASTA file. Required for looking up sequence context of cytosines.
REF = pysam.Fastafile(args.reference)

### Function definitions
# Extract reference sequence +/- 2bp of the location of a mapped read.
# The various cases are to ensure that all possible cases are considered, such as when the read maps to the start (end) of a contig and therefore there are no prior (subsequent) 2 bases to extract.
# The result is padded with 'N' characters in the cases where there is no reference sequence at the end positions
# In its current form this function cannot handle a read that is longer than its mapped-contig at each end, e.g. a 100bp read containing a 98bp contig internally, but these shouldn't exist and would likely not be mapped even if such short contigs did exist.
def getPaddedRefSeq(read, ref, IN): # 'read' is an AlignedRead object, 'ref' is a FastaFile object
    try:
        refseq = ref.fetch(IN.getrname(read.tid), start = (read.pos - 2), end = (read.aend + 2)).upper()
    except ValueError:
        try:
            refseq = ''.join(['N', ref.fetch(IN.getrname(read.tid), start = (read.pos - 1), end = (read.aend + 2)).upper()])
        except ValueError:
            try: 
                refseq = ''.join(['NN', ref.fetch(IN.getrname(read.tid), start = read.pos , end = (read.aend + 2)).upper()])
            except ValueError:
                try:
                    refseq = ''.join([ref.fetch(IN.getrname(read.tid), start = (read.pos - 2) , end = (read.aend + 1)).upper(), 'N'])
                except ValueError:
                    try:
                        refseq = ''.join([ref.fetch(IN.getrname(read.tid), start = (read.pos - 2) , end = read.aend).upper(), 'NN'])
                    except ValueError:
                        print "Skipped read: Unable to extract reference sequence for read", read.qname
                        refseq = ''.join('.' for i in xrange(read.qlen))
    return refseq

# Get the genomic-context for every C in the read and create the XM tag for that read
# If there is a C in the read, look ahead 1 or 2 bases to determine the genomic-context, i.e. CG, CHG or CHH.
# If the reference position is a C but the mapped base is something other than a C or a T then the '.' character is recorded in the XM tag for that position
def makeXMtag(read, refseq, strand): # read is an AlignedRead object, refseq is the output of getPaddedRefSeq, strand is a character "+" or "-" encoding which strand the read aligned to
    XM = []
    if strand == "+":
        for i in range(read.qlen): # i indexes position along the read
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
### End function definitions
                    
# Main loop - create XM tag for each read and write the new BAM to file
for read in IN:
    refseq = getPaddedRefSeq(read, REF, IN)
    if read.opt('XG') == 'CT':
        strand = "+"
    elif read.opt('XG') == 'GA':
        strand = "-"
    else:
        print 'Read skipped: Undefined strand (missing XG-tag) for read', read.qname
        continue
    XM = makeXMtag(read, refseq, strand)
    read.tags = read.tags + [('XM', XM)]
    OUT.write(read)

IN.close()
OUT.close()
REF.close()
