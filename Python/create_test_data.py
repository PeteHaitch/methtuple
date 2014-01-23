#### DESCRIPTION ####
# Extract reads to be used as test cases for comethylation.py
# Peter Hickey
# 15/01/2013

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

#### TODO ###
# TODO: Need to manually check each type of read listed on GitHub issues tracker. 

#### Import modules ####
import pysam
import warnings
import sys

#### Define BAM files that reads are to be extracted from ####
se_bam = '/home/users/lab0605/hickey/Lister_2009_BS-seq_data/H1/r1/Downloaded_mapped_reads/BAM/DM_H1_r1.bam'
pe_bam = '/home/users/lab0605/hickey/comethylation_tests/QS_ADS_n=100000.bam'

#### |-----> OT-strand with 0-4 CpGs ####
BAM = pysam.Samfile(se_bam)
test_reads = []
n_tuples = [0, 1, 2, 3, 4]
for read in BAM:
    if read.opt('XG') == 'CT' and read.opt('XR') == 'CT' and read.opt('XM').upper().count('Z') in n_tuples:
        test_reads.append(read)
        n_tuples.remove(read.opt('XM').upper().count('Z'))
    if not n_tuples:
        break

#### <-----| OB-strand with 0-4 CpGs ####
BAM = pysam.Samfile(se_bam)
n_tuples = [0, 1, 2, 3, 4]
for read in BAM:
    if read.opt('XG') == 'GA' and read.opt('XR') == 'CT' and read.opt('XM').upper().count('Z') in n_tuples:
        test_reads.append(read)
        n_tuples.remove(read.opt('XM').upper().count('Z'))
    if not n_tuples:
        break

#### |--r1--> <--r2--| OT-strand with 0-4 CpGs in r1 and 0-4 CpGs in r2####
BAM = pysam.Samfile(pe_bam)
test_reads = []
n_tuples = [(x, y) for x in [0, 1, 2, 3, 4] for y in [0, 1, 2, 3, 4]]
for read in BAM:
    if read.is_paired and read.is_read1:
        read_1 = read
        continue
    elif read.is_paired and read.is_read2:
        read_2 = read
        # Check that read_1 and read_2 are aligned to the same chromosome and have identical read-names and that readpair aligned to OT-strand
        if (read_1.tid == read_2.tid and read_1.qname == read_2.qname) and (read_1.opt('XG') == 'CT' and read_2.opt('XG') == 'CT' and read_1.opt('XR') == 'CT' and read_2.opt('XR') == 'GA'):
            # Check that readpair satisfy the conditions on the number of methylation loci per read
            if (read_1.opt('XM').upper().count('Z'), read_2.opt('XM').upper().count('Z')) in n_tuples:
                test_reads = test_reads + [read_1, read_2]
                n_tuples.remove((read_1.opt('XM').upper().count('Z'), read_2.opt('XM').upper().count('Z')))
                if not n_tuples:
                    break
                else:
                    continue
            else:
                continue
        elif read_1.tid != read_2.tid:
            warning_msg = ''.join(['Skipping readpair', read_1.qname, ' as reads aligned to different chromosomes (', BAM.getrname(read_1.tid), ' and ', BAM.getrname(read_2.tid), ')'])
            warnings.warn(warning_msg)
            continue
        elif read_1.qname != read_2.qname:
            exit_msg = "ERROR: The name of read_1 is not identical to to that of read_2 for readpair ", read_1.qname, read_2.qname, "\nPlease sort your paired-end BAM file in query-name-order with Picard's SortSam function."
            sys.exit(exit_msg)
        elif not read.is_paired:
            warning_msg = ''.join(['Skipping read', read_1.qname, ' as it is a single-end read'])
            warnings.warn(warning_msg)
            continue            
        else:
            warning_msg = ''.join(["Read is neither a single-end read nor part of a paired-end read. Check the SAM flag values are correctly set for read:", read.qname])
            warnings.warn(warning_msg)
            continue

#### |--r2--> <--r1--| OB-strand with 0-4 CpGs in r1 and 0-4 CpGs in r2####
BAM = pysam.Samfile(pe_bam)
test_reads = []
n_tuples = [(x, y) for x in [0, 1, 2, 3, 4] for y in [0, 1, 2, 3, 4]]
for read in BAM:
    if read.is_paired and read.is_read1:
        read_1 = read
        continue
    elif read.is_paired and read.is_read2:
        read_2 = read
        # Check that read_1 and read_2 are aligned to the same chromosome and have identical read-names and that readpair aligned to OB-strand
        if (read_1.tid == read_2.tid and read_1.qname == read_2.qname) and (read_1.opt('XG') == 'GA' and read_2.opt('XG') == 'GA' and read_1.opt('XR') == 'CT' and read_2.opt('XR') == 'GA'):
            # Check that readpair satisfy the conditions on the number of methylation loci per read
            if (read_1.opt('XM').upper().count('Z'), read_2.opt('XM').upper().count('Z')) in n_tuples:
                test_reads = test_reads + [read_1, read_2]
                n_tuples.remove((read_1.opt('XM').upper().count('Z'), read_2.opt('XM').upper().count('Z')))
                if not n_tuples:
                    break
                else:
                    continue
            else:
                continue
        elif read_1.tid != read_2.tid:
            warning_msg = ''.join(['Skipping readpair', read_1.qname, ' as reads aligned to different chromosomes (', BAM.getrname(read_1.tid), ' and ', BAM.getrname(read_2.tid), ')'])
            warnings.warn(warning_msg)
            continue
        elif read_1.qname != read_2.qname:
            exit_msg = "ERROR: The name of read_1 is not identical to to that of read_2 for readpair ", read_1.qname, read_2.qname, "\nPlease sort your paired-end BAM file in query-name-order with Picard's SortSam function."
            sys.exit(exit_msg)
        elif not read.is_paired:
            warning_msg = ''.join(['Skipping read', read_1.qname, ' as it is a single-end read'])
            warnings.warn(warning_msg)
            continue            
        else:
            warning_msg = ''.join(["Read is neither a single-end read nor part of a paired-end read. Check the SAM flag values are correctly set for read:", read.qname])
            warnings.warn(warning_msg)
            continue
