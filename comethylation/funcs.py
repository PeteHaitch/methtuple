from __future__ import print_function
from .mtuple import *

import re
import csv
import operator
import sys

#### Function definitions ####
def make_ignores_list(ic):
    """Make a list from a string of read positions that are to be ignored.

    Args:
        ic: A string of read positions to ignore. Multiple values should be comma-delimited and ranges can be specified by use of the hyphen, For example:

        '1-5, 80, 98-100'

        corresponds to ignoring read positions 1, 2, 3, 4, 5, 80, 98, 99, 100.

    Returns:
        A Python list of the positions to be ignored.
    """
    if ic is None:
        val = []
    else:
        val = []
        y = [x.strip() for x in ic.split(',')]
        for i in y:
            z = [x.strip() for x in i.split('-')]
            if len(z) == 2:
                val = val + list(range(int(z[0]), int(z[1]) + 1))
            elif len(z) == 1:
                val = val + [int(z[0])]
            else:
                exit_msg = ''.join(['ERROR: -ir1p and -ir2p must be comma-delimited. Ranges can be specified by use of the hyphen, e.g. \'1-5, 80, 98-100\''])
                sys.exit(exit_msg)
        if not all(isinstance(i, int) for i in val):
                exit_msg = ''.join(['ERROR: -ir1p and -ir2p must be comma-delimited. Ranges can be specified by use of the hyphen, e.g. \'1-5, 80, 98-100\''])
                sys.exit(exit_msg)
    return val

def ignore_read_pos(read, methylation_index, ignore_read_pos_list):
    """Ignore methylation loci in a read that appear in the ignore_read_pos_list. A methylation locus may be one of CpG, CHH, CHG or CNN.

    Args:
        read: A pysam.AlignedRead instance.
        methylation_index: A list of zero-based indices. Each index corresponds to the leftmost aligned position of a methylation locus in a read. For example:

        [0, 5]

        corresponds to a read with a methylation locus at the first and sixth positions of the read.
        ignore_read_pos_list: The list of read positions to be ignored.

    Returns:
        An updated version of methylation_index. Will report a warning if the FLAG does not encode whether the read is part of a paired-end or which mate of the paired-end read it is. Will report an error and call sys.exit() if the XR-tag or XG-tag is incompatible or missing.
    """
    # NOTE: Assumes that paired-end reads have FR orientation, which is always true for Bismark but might not be for other aligners
    strand = get_strand(read)
    # Single-end reads
    if not read.is_paired:
        if strand == 'OT':
            mi_updated = [mi for mi in methylation_index if mi not in ignore_read_pos_list]
        elif strand == 'OB':
            ignore_read_pos_list = [read.rlen - ic - 1 for ic in ignore_read_pos_list]
            mi_updated = [mi for mi in methylation_index if mi not in ignore_read_pos_list]

    # Paired-end reads: read_1
    elif read.is_paired and read.is_read1:
        if strand == 'OT':
            mi_updated = [mi for mi in methylation_index if mi not in ignore_read_pos_list]
        elif strand == 'OB':
            ignore_read_pos_list = [read.rlen - ic - 1 for ic in ignore_read_pos_list]
            mi_updated = [mi for mi in methylation_index if mi not in ignore_read_pos_list]
    # Paired-end reads: read_2
    elif read.is_paired and read.is_read2:
        if strand == 'OT':
            ignore_read_pos_list = [read.rlen - ic - 1 for ic in ignore_read_pos_list]
            mi_updated = [mi for mi in methylation_index if mi not in ignore_read_pos_list]
        if strand == 'OB':
            mi_updated = [mi for mi in methylation_index if mi not in ignore_read_pos_list]

    # Return updated methylation_index
    return mi_updated
        
def ignore_low_quality_bases(read, methylation_index, min_qual, phred_offset):
    """Ignore low quality bases of a read that contribute to a read's methylation_index.

    Args:
        read: A pysam.AlignedRead instance.
        methylation_index: A list of zero-based indices. Each index corresponds to the leftmost aligned position of a methylation locus in a read. For example:

        [0, 5]

        corresponds to a read with a methylation locus at the first and sixth positions of the read.
        min_qual: The minimum base quality (integer). All bases with quality < min_qual are excluded from the returned methylation_index instance.
        phred_offset: The Phred offset of the data (33 or 64).

    Returns:
        An updated version of methylation_index.
    
    """
    if (min_qual < 0) or (round(min_qual) != min_qual):
        raise ValueError("ignore_low_quality_bases: 'low_qual' must be a positive integer")
    if phred_offset != 33 and phred_offset != 64:
        raise ValueError("ignore_low_quality_bases: 'phred_offset' must be a 33 or 64")

    ignore_these_bases = []
    for i in methylation_index:
        bqual = bytearray(read.qual)
        if (bqual[i] - phred_offset) < min_qual:
            ignore_these_bases.append(i)
    return [x for x in methylation_index if x not in ignore_these_bases]

def fix_old_bismark(read):
	"""Fix the QNAME and FLAG field of a paired-end read from a SAM/BAM file generated by Bismark version < 0.8.3

	Args:
		read: A pysam.AlignedRead instance.

	Returns:
		An updated version of the read.
    
	"""
	# Strip '/1' or '/2' appended to the end of QNAMEs by Bismark version < 0.8.3. Assumes there are no forward slash characters in the QNAME field
	read.qname = read.qname.split('/')[0]
	# Fix FLAG value
	if read.flag == 67:
		read.flag = 99
	elif read.flag == 115:
		read.flag = 83
	elif read.flag == 131:
		read.flag = 147
	elif read.flag == 179:
		read.flag = 163
	else:
		exit_msg = ''.join(['ERROR: Unexpected FLAG (', str(read.flag), ') for read ', read.qname, 'Sorry, --aligner Bismark_old is unable to deal with this FLAG. Please log an issue at www.github.com/PeteHaitch/comethylation describing the error or email me at peter.hickey@gmail.com.'])
		sys.exit(exit_msg)
	return read    

def is_overlapping_sequence_identical(read_1, read_2, n_overlap, overlap_check):
    """Check whether the overlapping sequence of read_1 and read_2 passes the filter specified by overlap_check

    Args:
        read_1: A pysam.AlignedRead instance with read.is_read1 == true. Must be paired with read_2.
        read_2: A pysam.AlignedRead instance with read.is_read2 == true. Must be paired with read_1.
        n_overlap: The number of bases in the overlap of read_1 and read_2 (must be > 0)
        overlap_check: The type of check to be performed (listed by most-to-least stringent): check the entire overlapping sequence is identical (sequence), check the XM-tag is identical for the overlapping region (XM), do no check of the overlapping bases but use the read with the higher quality basecalls in the overlapping region (quality), or simply use the overlapping bases from read_1 ala bismark_methylation_extractor (bismark)

    Returns:
        True if the overlapping sequence passes the filter, False otherwise. Furthermore, if 'overlap_check = quality' or 'overlap_check = bismark' the result is always True.
    """
    if (n_overlap < 0) or (round(n_overlap) != n_overlap):
        raise ValueError("is_overlapping_sequence_identical: 'n_overlap' must be a positive integer")
    if overlap_check != 'sequence' and overlap_check != 'XM' and overlap_check != 'quality' and overlap_check != 'bismark':
        raise ValueError("is_overlapping_sequence_identical: 'overlap_check' must be one of 'sequence', 'XM', 'quality' or 'bismark'")
    strand_1 = get_strand(read_1)
    strand_2 = get_strand(read_2)
    # Readpair is informative for OT-strand
    if strand_1 == 'OT' and strand_2 == 'OT':
        if overlap_check == 'sequence':
            overlap_1 = read_1.seq[-n_overlap:]
            overlap_2 = read_2.seq[:n_overlap]
        elif overlap_check == 'XM':
            overlap_1 = read_1.opt('XM')[-n_overlap:]
            overlap_2 = read_2.opt('XM')[:n_overlap]
        elif overlap_check == 'bismark': # return True as Bismark does not actually check the overlapping sequence but rather just takes the overlap from read_1
            overlap_1 = True
            overlap_2 = True
        elif overlap_check == 'quality':
            overlap_1 = True
            overlap_2 = True
    # Readpair is informative for OB-strand
    elif strand_1 == 'OB' and strand_2 == 'OB':
        if overlap_check == 'sequence':
            overlap_1 = read_1.seq[:n_overlap]
            overlap_2 = read_2.seq[-n_overlap:]
        elif overlap_check == 'XM':
            overlap_1 = read_1.opt('XM')[:n_overlap]
            overlap_2 = read_2.opt('XM')[-n_overlap:]
        elif overlap_check == 'bismark': # return True as Bismark does not actually check the overlapping sequence but rather just takes the overlap from read_1
            overlap_1 = True
            overlap_2 = True
        elif overlap_check == 'quality':
            overlap_1 = True
            overlap_2 = True
    else:
        exit_msg = ''.join(['ERROR: The informative strands for readpair ', read_1.qname, ',  do not agree between mates. This should not happen.\nPlease log an issue at www.github.com/PeteHaitch/comethylation describing the error or email me at peter.hickey@gmail.com'])
        sys.exit(exit_msg)
    return overlap_1 == overlap_2

def does_read_contain_indel(read):
    """Check whether a read contains an insertion or deletion (InDel).

    Args:
        read: A pysam.AlignedRead instance.

    Returns:
        True if read contains an INDEL, False otherwise.
    """
    val = any([x[0] in [1, 2] for x in read.cigar]) # In pysam, the CIGAR operation for an insertion to the reference is 1 and the CIGAR operation for a deletion to the reference is 2.
    return val

def does_read_contain_complicated_cigar(read):
    """Check whether a read contains a complicated CIGAR string character, defined as anything other than a match (M; 0), insertion (I; 1) or deletion (D; 2).

    Args:
        read: A pysam.AlignedRead instance.

    Returns:
        True if read contains an complicated CIGAR string character, False otherwise.
    """
    val = any([x[0] not in [0, 1, 2] for x in read.cigar]) # In pysam, the CIGAR operation for an insertion to the reference is 1 and the CIGAR operation for a deletion to the reference is 2.
    return val

def ignore_overlapping_sequence(read_1, read_2, methylation_index_1, methylation_index_2, n_overlap, overlap_check):
    """Ignore the overlapping sequence of read_1 and read_2 from the read with the lower (sum) base qualities in the overlapping region.
       If base qualities are identical then (arbitrarily) ignore the overlapping bases from read_2.

    Args:
        read_1: A pysam.AlignedRead instance with read.is_read1 == true. Must be paired with read_2.
        read_2: A pysam.AlignedRead instance with read.is_read2 == true. Must be paired with read_1.
        methylation_index_1: A list of zero-based indices.  Each index corresponds to the leftmost aligned position of a methylation locus in read_1. For example:

        [0, 5]

        corresponds to read_1 with a methylation locus at the first and sixth positions of the read.
        methylation_index_2: As for methylation_index_1 but informative for read_2.
        n_overlap: The number of bases in the overlap (must be > 0).
        overlap_check: The type of check to be performed (listed by most-to-least stringent): check the entire overlapping sequence is identical (sequence), check the XM-tag is identical for the overlapping region (XM), do no check of the overlapping bases but use the read with the higher quality basecalls in the overlapping region (quality), or simply use the overlapping bases from read_1 ala bismark_methylation_extractor (bismark)

    Returns:
        Updated versions of methylation_index_1 and methylation_index_2.
    """
    if (n_overlap < 0) or (round(n_overlap) != n_overlap):
        raise ValueError("ignore_overlapping_sequence: 'n_overlap' must be a positive integer")
    if overlap_check != 'sequence' and overlap_check != 'XM' and overlap_check != 'quality' and overlap_check != 'bismark':
        raise ValueError("ignore_overlapping_sequence: 'overlap_check' must be one of 'sequence', 'XM', 'quality' or 'bismark'")
    strand_1 = get_strand(read_1)
    strand_2 = get_strand(read_2)
    ignore_these_bases = []
    # Readpair is informative for OT-strand
    if strand_1 == 'OT' and strand_2 == 'OT':
        bqual_1 = bytearray(read_1.qual)
        bqual_2 = bytearray(read_2.qual)
        overlap_quals_1 = sum([x for x in bqual_1[-n_overlap:]])
        overlap_quals_2 = sum([x for x in bqual_2[:n_overlap]])
        if (overlap_quals_1 >= overlap_quals_2) | (overlap_check == 'bismark'): # overlap_check == 'bismark' simply means use the overlapping sequence from read_1.
            for i in methylation_index_2:
                if i < n_overlap:
                    ignore_these_bases.append(i)
                    methylation_index_2 = [x for x in methylation_index_2 if x not in ignore_these_bases]
        else:
            for i in methylation_index_1:
                if i >= (read_1.alen - n_overlap):
                    ignore_these_bases.append(i)
                    methylation_index_1 = [x for x in methylation_index_1 if x not in ignore_these_bases]
    # Readpair is inforamtive for OB-strand
    elif strand_1 == 'OB' and strand_2 == 'OB':
        bqual_1 = bytearray(read_1.qual)
        bqual_2 = bytearray(read_2.qual)
        overlap_quals_1 = sum([x for x in bqual_1[:n_overlap]])
        overlap_quals_2 = sum([x for x in bqual_2[-n_overlap:]])
        if (overlap_quals_1 >= overlap_quals_2) | (overlap_check == 'bismark'): # overlap_check == 'bismark' simply means use the overlapping sequence from read_1.
            for i in methylation_index_2:
                if i >= (read_2.alen - n_overlap):
                    ignore_these_bases.append(i)
                    methylation_index_2 = [x for x in methylation_index_2 if x not in ignore_these_bases]
        else:
            for i in methylation_index_1:
                if i < n_overlap:
                    ignore_these_bases.append(i)
                    methylation_index_1 = [x for x in methylation_index_1 if x not in ignore_these_bases]
    else:
        exit_msg = ''.join(['ERROR: The informative strands for readpair ', read_1.qname, ',  do not agree between mates. This should not happen.\nPlease log an issue at www.github.com/PeteHaitch/comethylation describing the error or email me at peter.hickey@gmail.com'])
        sys.exit(exit_msg)
    return methylation_index_1, methylation_index_2

def extract_and_update_methylation_index_from_single_end_read(read, BAM, methylation_m_tuples, m, methylation_type, methylation_pattern, ignore_read1_pos, min_qual, phred_offset, ob_strand_offset):
    """Extracts m-tuples of methylation loci from a single-end read and adds the comethylation m-tuple to the methylation_m_tuples object.
    
    Args:
        read: An AlignedRead instance corresponding to a single-end read.
        BAM: The Samfile instance corresponding to the sample. Required in order to extract chromosome names from read.
        methylation_m_tuples: An MTuple instance.
        methylation_type: A string of the methylation type, e.g. CG for CpG methylation. Must be a valid option for the MTuple class.
        methylation_pattern: A regular expression of the methylation loci, e.g. '[Zz]' for CpG-methylation
        m: Is the "m" in "m-tuple", i.e. the size of the m-tuple. m must be an integer greater than or equal to 1. WARNING: No error or warning produced if this condition is violated.
        ignore_read1_pos: Ignore this list of read positions from each read.
        min_qual: Ignore bases with quality-score less than this value.
        phred_offset: The offset in the Phred scores. Phred33 corresponds to phred_offset = 33 and Phred64 corresponds to phred_offset 64.
        ob_strand_offset: How many bases a methylation loci on the OB-strand must be moved to the left in order to line up with the C on the OT-strand; e.g. ob_strand_offset = 1 for CpGs.
    Returns:
        methylation_m_tuples: An updated version of methylation_m_tuples.
        n_methylation_loci: The number of methylation loci extracted from the read.
    """
    # Identify methylation events in read, e.g. CpGs or CHHs. The methylation_pattern is specified by a command line argument (e.g. Z/z corresponds to CpG)
    methylation_index = [midx.start() for midx in re.finditer(methylation_pattern, read.opt('XM'))]
    # Ignore any read positions specified in ignore_read1_pos
    methylation_index = ignore_read_pos(read, methylation_index, ignore_read1_pos)
    # Ignore any positions with a base quality less than min_qual
    methylation_index = ignore_low_quality_bases(read, methylation_index, min_qual, phred_offset)
    n_methylation_loci = len(methylation_index)
    strand = get_strand(read)
    # Case A: >= m methylation loci in the read
    if n_methylation_loci >= m:
        positions = [read.pos + x + 1 for x in methylation_index] # +1 to transform from 0-based to 1-based co-ordinates.
        # If read is informative for the OB-strand then translate co-ordinate "ob_strand_offset" bases to the left so that it points to the C on the OT-strand of the methylation locus.
        if strand == 'OB':
            positions = [x - ob_strand_offset for x in positions]
        # Exit if methylation loci are incorrectly ordered
        if not positions == sorted(positions):
            exit_msg = ' '.join(["ERROR: The positions of the methylation loci are not properly ordered for single-end read", read.qname, "\n'positions' =", str(positions), '.\nPlease log an issue at www.github.com/PeteHaitch/comethylation describing the error or email me at peter.hickey@gmail.com'])
            sys.exit(exit_msg)
        # Construct each bookended methylation-loci m-tuple and add it to the methylation_m_tuple object.
        for i in range(0, len(methylation_index) - m + 1): # For a read containing k methylation loci there are (k - m + 1) m-tuples.
            # Increment count. The MTuple.increment_count() method automatically checks whether this particular m-tuple has been observed before and updated as appropriate.
            this_comethylation_pattern = ''.join([read.opt('XM')[j] for j in methylation_index[i:(i + m)]])
            this_m_tuple_positions = (BAM.getrname(read.tid),) + tuple(positions[i:(i + m)])
            methylation_m_tuples.increment_count(this_m_tuple_positions, this_comethylation_pattern, read, None)
    return methylation_m_tuples, n_methylation_loci

def extract_and_update_methylation_index_from_paired_end_reads(read_1, read_2, BAM, methylation_m_tuples, m, methylation_type, methylation_pattern, ignore_read1_pos, ignore_read2_pos, min_qual, phred_offset, ob_strand_offset, overlap_check, n_fragment_skipped_due_to_bad_overlap, FAILED_QC):
    """Extracts m-tuples of methylation loci from a readpair and adds the comethylation m-tuple to the methylation_m_tuples object.
    
    Args:
        read_1: An AlignedRead instance corresponding to read_1 of the readpair.
        read_2: An AlignedRead instance corresponding to read_2 of the readpair.
        BAM: The Samfile instance corresponding to the sample. Required in order to extract chromosome names from read.
        methylation_m_tuples: An MTuple instance.
        m: Is the "m" in "m-tuple", i.e. the size of the m-tuple. m must be an integer greater than or equal to 1. WARNING: No error or warning produced if this condition is violated.
        methylation_type: A string of the methylation type, e.g. CG for CpG methylation. Must be a valid option for the MTuple class.
        methylation_pattern: A regular expression of the methylation loci, e.g. '[Zz]' for CpG-methylation
        ignore_read1_pos: Ignore this list of positions from each read_1.
        ignore_read2_pos: Ignore this list of positions from each read_2.
        min_qual: Ignore bases with quality-score less than this value.
        phred_offset: The offset in the Phred scores. Phred33 corresponds to phred_offset = 33 and Phred64 corresponds to phred_offset 64.
        ob_strand_offset: How many bases a methylation loci on the OB-strand must be moved to the left in order to line up with the C on the OT-strand; e.g. ob_strand_offset = 1 for CpGs.
        overlap_check: The type of check to be performed (listed by most-to-least stringent): check the entire overlapping sequence is identical (sequence), check the XM-tag is identical for the overlapping region (XM), do no check of the overlapping bases but use the read with the higher quality basecalls in the overlapping region (quality), or simply use the overlapping bases from read_1 ala bismark_methylation_extractor (bismark)
        n_fragment_skipped_due_to_bad_overlap: The total number of fragments (readpairs) skipped due to the overlapping sequencing not passing the filter.
        FAILED_QC: The file object where the QNAME of readpairs that fail the overlap check are written, along with the reason the readpairs failed
    Returns:
        methylation_m_tuples: An updated version of methylation_m_tuples.
        n_methylation_loci: The number of methylation loci extracted from the read.
    """
    # Identify methylation events in read, e.g. CpGs or CHHs. The methylation_pattern is specified by a command line argument (e.g. Z/z corresponds to CpG)
    methylation_index_1 = [midx.start() for midx in re.finditer(methylation_pattern, read_1.opt('XM'))]
    methylation_index_2 = [midx.start() for midx in re.finditer(methylation_pattern, read_2.opt('XM'))]
    # Ignore any read positions specified in ignore_read1_pos or ignore_read1_pos
    methylation_index_1 = ignore_read_pos(read_1, methylation_index_1, ignore_read1_pos)
    methylation_index_2 = ignore_read_pos(read_2, methylation_index_2, ignore_read2_pos)
    # Ignore any positions with a base quality less than  min_qual
    methylation_index_1 = ignore_low_quality_bases(read_1, methylation_index_1, min_qual, phred_offset)
    methylation_index_2 = ignore_low_quality_bases(read_2, methylation_index_2, min_qual, phred_offset)
    # Check strand of each mate make sense
    strand_1 = get_strand(read_1)
    strand_2 = get_strand(read_2)

    # Check for overlapping reads from a readpair.
    # If reads overlap check whether the overlapping sequence passes the filter given by overlap_check.
    # If the overlapping sequence does not pass the filter report a warning, increment a counter and skip the readpair (by setting methylation_index_1 and methylation_index_2 to be the empty list).
    n_overlap = read_1.alen + read_2.alen - abs(read_1.tlen)
    if n_overlap > 0:
        if is_overlapping_sequence_identical(read_1, read_2, n_overlap, overlap_check):
            methylation_index_1, methylation_index_2 = ignore_overlapping_sequence(read_1, read_2, methylation_index_1, methylation_index_2, n_overlap, overlap_check)
        else:
            failed_read_msg = '\t'.join([read_1.qname, ''.join(['failed the --overlap-filter ', overlap_check, '\n'])])
            FAILED_QC.write(failed_read_msg)
            n_fragment_skipped_due_to_bad_overlap += 1
            methylation_index_1 = []
            methylation_index_2 = []
    n_methylation_loci = len(methylation_index_1) + len(methylation_index_2)
    # Only process readpair if there are at least enough CpGs to form one m-tuple.
    if n_methylation_loci >= m:
        positions_1 = [read_1.pos + x + 1 for x in methylation_index_1] # +1 to transform from 0-based to 1-based co-ordinates.
        positions_2 = [read_2.pos + x + 1 for x in methylation_index_2] # +1 to transform from 0-based to 1-based co-ordinates.
        if any(x in positions_1 for x in positions_2):
            exit_msg = ''.join(['ERROR: For readpair ', read_1.qname, ', position_1 and position_2 contain a common position. This should not happen.\nPlease log an issue at www.github.com/PeteHaitch/comethylation describing the error or email me at peter.hickey@gmail.com'])
            print(positions_1) # TODO: Remove? This is a debugging statement
            print(positions_2) # TODO: Remove? This is a debugging statement
            sys.exit(exit_msg)
        # Case 1: Readpair is informative for OT-strand
        if strand_1 == 'OT' and strand_2 == 'OT':
            # Exit if methylation loci are incorrectly ordered
            if positions_1 + positions_2 != sorted(positions_1 + positions_2):
                exit_msg = ' '.join(["ERROR: The positions of the methylation loci are not properly ordered for paired-end read", read_1.qname, ", which is informative for the OT-strand.\n'positions_1 + positions_2' =", str(positions_1 + positions_2), '\nPlease log an issue at www.github.com/PeteHaitch/comethylation describing the error or email me at peter.hickey@gmail.com'])
                sys.exit(exit_msg)
            # Firstly, create all m-tuples of methylation loci where each locus is from read_1.
            if len(methylation_index_1) >= m:
                for i in range(0, len(methylation_index_1) - m + 1): # For a read containing k methylation loci there are (k - m + 1) m-tuples.
                    # Increment count. The MTuple.increment_count() method automatically checks whether this particular m-tuple has been observed before and updated as appropriate.
                    this_comethylation_pattern = ''.join([read_1.opt('XM')[j] for j in methylation_index_1[i:(i + m)]])
                    this_m_tuple_positions = (BAM.getrname(read_1.tid),) + tuple(positions_1[i:(i + m)])
                    methylation_m_tuples.increment_count(this_m_tuple_positions, this_comethylation_pattern, read_1, read_2)
            # Secondly, create all m-tuples of methylation loci where the leftmost locus is on read_1 and the rightmost locus is on read_2
            num_shared_m_tuples = max(len(methylation_index_1) + len(methylation_index_2) - m + 1, 0) - max(len(methylation_index_1) - m + 1, 0) - max(len(methylation_index_2) - m + 1, 0) # the number of m-tuples that span read_1 and read_2
            leftmost_shared_locus_index = max(0, len(methylation_index_1) - m + 1) # The index of the leftmost locus to be part of a "shared" m-tuple. The rightmost_shared_locus_index = min(m - 2, len(methylation_index_2) - 1), however this is not required
            for i in range(0, num_shared_m_tuples):
                this_m_tuple_positions_1 = positions_1[(leftmost_shared_locus_index + i):]
                this_m_tuple_positions_2 = positions_2[:(m - len(this_m_tuple_positions_1))]
                # Exit if methylation loci are incorrectly ordered. While a similar check is performed a few lines above, this is a sanity check to make sure than nothing has gone wrong in constructing the shared m-tuples
                if this_m_tuple_positions_1 + this_m_tuple_positions_2 != sorted(this_m_tuple_positions_1 + this_m_tuple_positions_2):
                    exit_msg = ' '.join(["ERROR: The positions of the shared methylation loci are not properly ordered for paired-end read", read_1.qname, ", which is informative for the OT-strand.\n'this_m_tuple_positions_1 + this_m_tuple_positions_2' =", str(this_m_tuple_positions_1 + this_m_tuple_positions_2), '\nPlease log an issue at www.github.com/PeteHaitch/comethylation describing the error or email me at peter.hickey@gmail.com'])
                    sys.exit(exit_msg)
                # Increment count. The MTuple.increment_count() method automatically checks whether this particular m-tuple has been observed before and updated as appropriate.
                this_comethylation_pattern = ''.join([read_1.opt('XM')[j] for j in methylation_index_1[(leftmost_shared_locus_index + i):]] + [read_2.opt('XM')[j] for j in methylation_index_2[:(m - len(this_m_tuple_positions_1))]])
                this_m_tuple_positions = (BAM.getrname(read_1.tid),) + tuple(this_m_tuple_positions_1) + tuple(this_m_tuple_positions_2)
                methylation_m_tuples.increment_count(this_m_tuple_positions, this_comethylation_pattern,  read_1, read_2)
            # Finally, create all m-tuples of methylation loci where each locus is from read_2.        
            if len(methylation_index_2) >= m:
                for i in range(0, len(methylation_index_2) - m + 1): # For a read containing k methylation loci there are (k - m + 1) m-tuples.:
                    this_comethylation_pattern = ''.join([read_2.opt('XM')[j] for j in methylation_index_2[i:(i + m)]])
                    this_m_tuple_positions = (BAM.getrname(read_2.tid),) + tuple(positions_2[i:(i + m)])
                    methylation_m_tuples.increment_count(this_m_tuple_positions, this_comethylation_pattern, read_1, read_2)

        # Case 2: Readpair is informative for OB-strand
        elif strand_1 == 'OB' and strand_2 == 'OB':
            # Translate co-ordinates "ob_strand_offset" bases to the left so that it points to the C on the OT-strand of the methylation locus
            positions_1 = [x - ob_strand_offset for x in positions_1]
            positions_2 = [x - ob_strand_offset for x in positions_2]
            # Exit if methylation loci are incorrectly ordered.
            if positions_2 + positions_1 != sorted(positions_2 + positions_1): 
                exit_msg = ' '.join(["ERROR: The positions of the methylation loci are not properly ordered for paired-end read", read_1.qname, "which is informative for the OB-strand.\n'positions_2 + positions_1' =", str(positions_2 + positions_1), '\nPlease log an issue at www.github.com/PeteHaitch/comethylation describing the error or email me at peter.hickey@gmail.com'])
                sys.exit(exit_msg)
            # Firstly, create all m-tuples of methylation loci where each locus is from read_1.
            if len(methylation_index_1) >= m:
                for i in range(0, len(methylation_index_1) - m + 1): # For a read containing m methylation loci there are (m - m-tuple + 1) m-tuples.:
                    # Increment count. The MTuple.increment_count() method automatically checks whether this particular m-tuple has been observed before and updated as appropriate.
                    this_comethylation_pattern = ''.join([read_1.opt('XM')[j] for j in methylation_index_1[i:(i + m)]])
                    this_m_tuple_positions = (BAM.getrname(read_1.tid),) + tuple(positions_1[i:(i + m)])
                    methylation_m_tuples.increment_count(this_m_tuple_positions, this_comethylation_pattern,  read_1, read_2)
            # Secondly, create all m-tuples of methylation loci where the leftmost locus is on read_1 and the rightmost locus is on read_2
            num_shared_m_tuples = max(len(methylation_index_1) + len(methylation_index_2) - m + 1, 0) - max(len(methylation_index_1) - m + 1, 0) - max(len(methylation_index_2) - m + 1, 0) # the number of m-tuples that span read_1 and read_2
            leftmost_shared_locus_index = max(0, len(methylation_index_2) - m + 1) # The index of the leftmost locus to be part of a "shared" m-tuple. The rightmost_shared_locus_index = min(m - 2, len(methylation_index_1) - 1), however this is not required (m - 2 = m - 1 - 1, because Python lists are 0-indexed)
            for i in range(0, num_shared_m_tuples):
                this_m_tuple_positions_2 = positions_2[(leftmost_shared_locus_index + i):]
                this_m_tuple_positions_1 = positions_1[:(m - len(this_m_tuple_positions_2))]
                # Exit if methylation loci are incorrectly ordered. While a similar check is performed a few lines above, this is a sanity check to make sure than nothing has gone wrong in constructing the shared m-tuples
                if this_m_tuple_positions_2 + this_m_tuple_positions_1 != sorted(this_m_tuple_positions_2 + this_m_tuple_positions_1):
                    exit_msg = ' '.join(["ERROR: The positions of the shared methylation loci are not properly ordered for paired-end read", read_1.qname, "which is aligned to the OB-strand.\n'this_m_tuple_positions_2 + this_m_tuple_positions_1' =", str(this_m_tuple_positions_2 + this_m_tuple_positions_1), '\nPlease log an issue at www.github.com/PeteHaitch/comethylation describing the error or email me at peter.hickey@gmail.com'])
                    sys.exit(exit_msg)
                # Increment count. The MTuple.increment_count() method automatically checks whether this particular m-tuple has been observed before and updated as appropriate.
                this_comethylation_pattern = ''.join([read_2.opt('XM')[j] for j in methylation_index_2[(leftmost_shared_locus_index + i):]] + [read_1.opt('XM')[j] for j in methylation_index_1[:(m - len(this_m_tuple_positions_2))]])
                this_m_tuple_positions = (BAM.getrname(read_1.tid),) + tuple(this_m_tuple_positions_2) + tuple(this_m_tuple_positions_1)
                methylation_m_tuples.increment_count(this_m_tuple_positions, this_comethylation_pattern, read_1, read_2)
            # Finally, create all m-tuples of methylation loci where each locus is from read_2.        
            if len(methylation_index_2) >= m:
                for i in range(0, len(methylation_index_2) - m + 1): # For a read containing m methylation loci there are (m - m-tuple + 1) m-tuples.:
                    this_comethylation_pattern = ''.join([read_2.opt('XM')[j] for j in methylation_index_2[i:(i + m)]])
                    this_m_tuple_positions = (BAM.getrname(read_2.tid),) + tuple(positions_2[i:(i + m)])
                    methylation_m_tuples.increment_count(this_m_tuple_positions, this_comethylation_pattern, read_1, read_2)
        else:
            exit_msg = ''.join(['ERROR: The informative strands for readpair ', read_1.qname, ',  do not agree between mates. This should not happen.\nPlease log an issue at www.github.com/PeteHaitch/comethylation describing the error or email me at peter.hickey@gmail.com'])
            sys.exit(exit_msg)
    return methylation_m_tuples, n_methylation_loci, n_fragment_skipped_due_to_bad_overlap
                                  
def write_methylation_m_tuples_to_file(methylation_m_tuples, OUT):
    """Write the methylation_m_tuples instance to a tab-separated file. The m-tuples are ordered by chromosome and genomic co-ordinates.
    
    Args:
        methylation_m_tuples: An MTuple instance.
        OUT: The file to write output to.
    """
    # tab_writer writes a tab-separated output file to the filehandle OUT
    tab_writer = csv.writer(OUT, delimiter='\t', quotechar=' ', quoting=csv.QUOTE_MINIMAL)

    # Get m
    m = methylation_m_tuples.m
    # Create the header row and write to file
    header = ['chr'] + ['pos' + str(i) for i in range(1, m + 1)] + methylation_m_tuples.comethylation_patterns
    tab_writer.writerow(header)
    # Sort methylation_m_tuples.mtuples.keys() by chromosome (using methylation_m_tuples.chr_map for sort order) and then positions (pos1, pos2, ... to posm)
    for this_m_tuple in sorted(list(methylation_m_tuples.mtuples.keys()), key = lambda x: (methylation_m_tuples.chr_map[x[0]], ) + tuple(x[1:])):
        row = this_m_tuple + tuple(methylation_m_tuples.mtuples[this_m_tuple])
        tab_writer.writerow(row)

def get_strand(read):
    """
    Report whether a read is informative for the OT-strand or OB-strand. 
    Currently using a strict check that ensures the reads are in the expected orientation for the given strand. 
    See commented out lines for a less-strict version.
    Will report an error and call sys.exit() if the XR-tag or XG-tag is incompatible or missing.

    Args:
        read: A pysam.AlignedRead instance with XR-tag and XG-tag.
    Returns:
        strand: For which strand the read/readpair is informative (OT or OB)
    """
    ## Single-end
    if not read.is_paired:
        ## Check if aligned to OT- or CTOT-strand, i.e., informative for OT-strand.
        if (read.opt('XR') == 'CT' and read.opt('XG') == 'CT') or (read.opt('XR') == 'GA' and read.opt('XG') == 'CT'): 
        # if read_1.opt('XG') == 'CT'
            strand = 'OT'
        ## Else, check if aligned to OB- or CTOB-strand, i.e., informative for OB-strand.
        elif (read.opt('XR') == 'CT' and read.opt('XG') == 'GA') or (read.opt('XR') == 'GA' and read.opt('XG') == 'GA'):
        # elif read_1.opt('XG') == 'GA'
            strand = 'OB'
        ## Else, something odd about this read
        else:
            exit_msg = ''.join(['ERROR: Read ', read.qname, ' has incompatible or missing XG-tag or XR-tag. Please log an issue at www.github.com/PeteHaitch/comethylation describing the error or email me at peter.hickey@gmail.com'])
            sys.exit(exit_msg)
    ## Paired-end
    elif read.is_paired:
        if read.is_read1:
            ## Check if aligned to CT- or CTOT-strand, i.e., informative for OT-strand.
            if (read.opt('XR') == 'CT' and read.opt('XG') == 'CT') or (read.opt('XR') == 'GA' and read.opt('XG') == 'CT'):
            #if read.opt('XG') == 'CT':
                strand = 'OT'
            ## Else, check if aligned to OB- or CTOB-strand, i.e., informative for OB-strand.
            elif (read.opt('XR') == 'CT' and read.opt('XG') == 'GA') or (read.opt('XR') == 'GA' and read.opt('XG') == 'GA'):
            #elif read.opt('XG') == 'GA':
                strand = 'OB'
            ## Else, something odd about this read
            else:
                exit_msg = ''.join(['ERROR: Read ', read.qname, ' has incompatible or missing XG-tag or XR-tag. Please log an issue at www.github.com/PeteHaitch/comethylation describing the error or email me at peter.hickey@gmail.com'])
                sys.exit(exit_msg)
        elif read.is_read2:
            ## Check if aligned CT or CTOT-strand, i.e., informative for OT-strand.
            if (read.opt('XR') == 'GA' and read.opt('XG') == 'CT') or (read.opt('XR') == 'CT' and read.opt('XG') == 'CT'):
                strand = 'OT'
            ## Else, check if aligned OB- or CTOB-strand, i.e., informative for OB-strand.
            elif (read.opt('XR') == 'GA' and read.opt('XG') == 'GA') or (read.opt('XR') == 'CT' and read.opt('XG') == 'GA'):
                strand = 'OB'
            ## Else, something odd about this read
            else:
                exit_msg = ''.join(['ERROR: Read ', read.qname, ' has incompatible or missing XG-tag or XR-tag. Please log an issue at www.github.com/PeteHaitch/comethylation describing the error or email me at peter.hickey@gmail.com'])
                sys.exit(exit_msg)
    else:
        exit_msg = ''.join(['ERROR: Read ', read.qname, ' is neither a single-end read nor part of a paired-end read. Please log an issue at www.github.com/PeteHaitch/comethylation describing the error or email me at peter.hickey@gmail.com'])
    return strand

__all__ = [
    'make_ignores_list',
    'ignore_read_pos',
    'ignore_low_quality_bases',
    'fix_old_bismark',
    'is_overlapping_sequence_identical',
    'does_read_contain_indel',
    'does_read_contain_complicated_cigar',
    'ignore_overlapping_sequence',
    'extract_and_update_methylation_index_from_single_end_read',
    'extract_and_update_methylation_index_from_paired_end_reads',
    'write_methylation_m_tuples_to_file',
    'get_strand'
]
