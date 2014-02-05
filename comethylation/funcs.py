from mtuple import *

import re
import csv
import operator



#### Function definitions ####
def ignore_first_n_bases(read, methylation_index, n):
    """Ignore methylation loci occuring in the first n bases of a read. A methylation locus may be one of CpG, CHH, CHG or CNN.

    Args:
        read: A pysam.AlignedRead instance.
        methylation_index: A list of zero-based indices. Each index corresponds to the leftmost aligned position of a methylation locus in a read. For example:

        [0, 5]

        corresponds to a read with a methylation locus at the first and sixth positions of the read.
        n: The number of bases to exclude from the start of each read. The start of a read is the first *sequenced* base, not the leftmost aligned base.

    Returns:
        An updated version of methylation_index. Will report a warning if the FLAG does not encode whether the read is part of a paired-end or which mate of the paired-end read it is. Will report an error and call sys.exit() if the XR-tag or XG-tag is incompatible or missing.
    """
    ignore_these_bases = []
    # Single-end reads
    if not read.is_paired:
        # Read aligned to OT-strand |------>
        if read.opt('XG') == 'CT' and read.opt('XR') == 'CT':
            for i in methylation_index:
                if i < n:
                    ignore_these_bases.append(i)
        # Read aligned to OB-strand <------|
        elif read.opt('XG') == 'GA' and read.opt('XR') == 'CT':
            for i in methylation_index:
                if i >= (read.alen - n):
                    ignore_these_bases.append(i)
        else:
            exit_msg = ''.join(['ERROR: Read ', read.qname, ' has incompatible or missing XG-tag or XR-tag. Please log an issue at www.github.com/PeteHaitch/Comethylation describing the error or email me at peter.hickey@gmail.com'])
            sys.exit(exit_msg)
            #methylation_index = []
    # Paired-end reads: read_1
    elif read.is_paired and read.is_read1:
        # read_1 aligned to OT-strand |------>
        if read.opt('XG') == 'CT' and read.opt('XR') == 'CT':
            for i in methylation_index:
                if i < n:
                    ignore_these_bases.append(i)
        # read_1 aligned to OB-strand <------|
        elif read.opt('XG') == 'GA' and read.opt('XR') == 'CT':
            for i in methylation_index:
                if i >= (read.alen - n):
                    ignore_these_bases.append(i)
        else:
            exit_msg = ''.join(['ERROR: Read ', read.qname, ' has incompatible or missing XG-tag or XR-tag. Please log an issue at www.github.com/PeteHaitch/Comethylation describing the error or email me at peter.hickey@gmail.com'])
            sys.exit(exit_msg)
            #methylation_index = []
    # Paired-end reads: read_2
    elif read.is_paired and read.is_read2:
        # read_2 aligned to OT-strand <------|
        if read.opt('XG') == 'CT' and read.opt('XR') == 'GA':
            for i in methylation_index:
                if i >= (read.alen - n):
                    ignore_these_bases.append(i)
        # read_2 aligned to OB-strand |------>
        elif read.opt('XG') == 'GA' and read.opt('XR') == 'GA':
            for i in methylation_index:
                if i < n:
                    ignore_these_bases.append(i)
        else:
            exit_msg = ''.join(['ERROR: Read ', read.qname, ' has incompatible or missing XG-tag or XR-tag. Please log an issue at www.github.com/PeteHaitch/Comethylation describing the error or email me at peter.hickey@gmail.com'])
            sys.exit(exit_msg)
            #methylation_index = []
    # ERROR: read does not have necessary information to infer whether read is paired or whether the read is read_1 or read_2 of the readpair. This should never happen as it should have already been caught in the main loop of the program.
    else:
        exit_msg = ''.join(['ERROR: Read ', read.qname, ' is missing the 0x01, 0x40 or 0x80 FLAG bit. This should never happen. Please log an issue at www.github.com/PeteHaitch/Comethylation describing the error or email me at peter.hickey@gmail.com.'])
        sys.exit(exit_msg)
        #methylation_index = []
    return [x for x in methylation_index if x not in ignore_these_bases]

def ignore_last_n_bases(read, methylation_index, n):
    """Ignore methylation loci occuring in the last n bases of a read. A methylation locus may be one of CpG, CHH, CHG or CNN.

    Args:
        read: A pysam.AlignedRead instance.
        methylation_index: A list of zero-based indices. Each index corresponds to the leftmost aligned position of a methylation locus in a read. For example:

        [0, 5]

        corresponds to a read with a methylation locus at the first and sixth positions of the read.
        n: The number of bases to exclude from the start of each read. The start of a read is the first *sequenced* base, not the leftmost aligned base.

    Returns:
        An updated version of methylation_index. Will report a warning if the FLAG does not encode whether the read is part of a paired-end or which mate of the paired-end read it is. Will report an error and call sys.exit() if the XR-tag or XG-tag is incompatible or missing.
    """
    ignore_these_bases = []
    # Single-end reads
    if not read.is_paired:
        # Read aligned to OT-strand |------>
        if read.opt('XG') == 'CT' and read.opt('XR') == 'CT':
            for i in methylation_index:
                if i >= (read.alen - n):
                    ignore_these_bases.append(i)
        # Read aligned to OB-strand <------|
        elif read.opt('XG') == 'GA' and read.opt('XR') == 'CT':
            for i in methylation_index:
                if i < n:
                    ignore_these_bases.append(i)
        else:
            exit_msg = ''.join(['ERROR: Read ', read.qname, ' has incompatible or missing XG-tag or XR-tag. Please log an issue at www.github.com/PeteHaitch/Comethylation describing the error or email me at peter.hickey@gmail.com'])
            sys.exit(exit_msg)
            #methylation_index = []
    # Paired-end reads: read_1
    elif read.is_paired and read.is_read1:
        # read_1 aligned to OT-strand |------>
        if read.opt('XG') == 'CT' and read.opt('XR') == 'CT':
            for i in methylation_index:
                if i >= (read.alen - n):
                    ignore_these_bases.append(i)
        # read_1 aligned to OB-strand <------|
        elif read.opt('XG') == 'GA' and read.opt('XR') == 'CT':
            for i in methylation_index:
                if i < n:
                    ignore_these_bases.append(i)
        else:
            exit_msg = ''.join(['ERROR: Read ', read.qname, ' has incompatible or missing XG-tag or XR-tag. Please log an issue at www.github.com/PeteHaitch/Comethylation describing the error or email me at peter.hickey@gmail.com'])
            sys.exit(exit_msg)
            #methylation_index = []
    # Paired-end reads: read_2
    elif read.is_paired and read.is_read2:
        # read_2 aligned to OT-strand <------|
        if read.opt('XG') == 'CT' and read.opt('XR') == 'GA':
            for i in methylation_index:
                if i < n:
                    ignore_these_bases.append(i)
        # read_2 aligned to OB-strand |------>
        elif read.opt('XG') == 'GA' and read.opt('XR') == 'GA':
            for i in methylation_index:
                if i  >= (read.alen - n):
                    ignore_these_bases.append(i)
        else:
            exit_msg = ''.join(['ERROR: Read ', read.qname, ' has incompatible or missing XG-tag or XR-tag. Please log an issue at www.github.com/PeteHaitch/Comethylation describing the error or email me at peter.hickey@gmail.com'])
            sys.exit(exit_msg)
            #methylation_index = []
    # ERROR: read does not have necessary information to infer whether read is paired or whether the read is read_1 or read_2 of the readpair. This should never happen as it should have already been caught in the main loop of the program.
    else:
        exit_msg = ''.join(['ERROR: Read ', read.qname, ' is missing the 0x01, 0x40 or 0x80 FLAG bit. This should never happen. Please log an issue at www.github.com/PeteHaitch/Comethylation describing the error or email me at peter.hickey@gmail.com.'])
        sys.exit(exit_msg)
        #methylation_index = []
    return [x for x in methylation_index if x not in ignore_these_bases]
        
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
    ignore_these_bases = []
    for i in methylation_index:
        if (ord(read.qual[i]) - phred_offset) < min_qual:
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
		exit_msg = ''.join(['ERROR: Unexpected FLAG (', str(read.flag), ') for read ', read.qname, 'Sorry, --oldBismark is unable to deal with this FLAG. Please log an issue at www.github.com/PeteHaitch/Comethylation describing the error or email me at peter.hickey@gmail.com.'])
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
        True if the overlapping sequence passes the filter, False otherwise (NB: this means that readpairs that trigger the warning for having mis-specified XG- or XR-tags will also return 'False').
    """
    # Readpair aligns to OT-strand
    if read_1.opt('XG') == 'CT' and read_2.opt('XG') == 'CT' and read_1.opt('XR') == 'CT' and read_2.opt('XR') == 'GA':
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
    # Readpair aligns to OB-strand
    elif read_1.opt('XG') == 'GA' and read_2.opt('XG') == 'GA' and read_1.opt('XR') == 'CT' and read_2.opt('XR') == 'GA':
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
        warning_msg = ''.join(['XG-tags or XR-tags for readpair ', read.qname, ' are inconsistent with OT-strand or OB-strand (XG-tags = ', read_1.opt('XG'),', ', read_2.opt('XG'), '; XR-tags = ', read_1.opt('XR'), ', ', read_2.opt('XR'), ')'])
        warnings.warn(warning_msg)
        overlap_1 = True
        overlap_2 = False
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
    ignore_these_bases = []
    # Readpair aligns to OT-strand
    if read_1.opt('XG') == 'CT' and read_2.opt('XG') == 'CT' and read_1.opt('XR') == 'CT' and read_2.opt('XR') == 'GA':
        overlap_quals_1 = sum([ord(x) for x in read_1.qual[-n_overlap:]])
        overlap_quals_2 = sum([ord(x) for x in read_2.qual[-n_overlap:]])
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
    # Readpair aligns to OB-strand
    elif read_1.opt('XG') == 'GA' and read_2.opt('XG') == 'GA' and read_1.opt('XR') == 'CT' and read_2.opt('XR') == 'GA':
        overlap_quals_1 = sum([ord(x) for x in read_1.qual[:n_overlap]])
        overlap_quals_2 = sum([ord(x) for x in read_2.qual[-n_overlap:]])
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
        warning_msg = ''.join(['XG-tags or XR-tags for readpair ', read_1.qname, ' are inconsistent with OT-strand or OB-strand (XG-tags = ', read_1.opt('XG'),', ', read_2.opt('XG'), '; XR-tags = ', read_1.opt('XR'), ', ', read_2.opt('XR'), ')'])
        warnings.warn(warning_msg)
        methylation_index_1 = []
        methylation_index_2 = []
    return methylation_index_1, methylation_index_2

def extract_and_update_methylation_index_from_single_end_read(read, BAM, methylation_m_tuples, m, methylation_type, methylation_pattern, ignore_start_r1, ignore_end_r1, min_qual, phred_offset, ob_strand_offset):
    """Extracts m-tuples of methylation loci from a single-end read and adds the comethylation m-tuple to the methylation_m_tuples object.
    
    Args:
        read: An AlignedRead instance corresponding to a single-end read.
        BAM: The Samfile instance corresponding to the sample. Required in order to extract chromosome names from read.
        methylation_m_tuples: A dictionary storing all observed m-tuples of methylation events and their WithinFragmentComethylationMTuple instance.
        methylation_type: A string of the methylation type, e.g. CG for CpG methylation. Must be a valid option for the WithinFragmentComethylationMTuple class.
        methylation_pattern: A regular expression of the methylation loci, e.g. '[Zz]' for CpG-methylation
        m: Is the "m" in "m-tuple", i.e. the size of the m-tuple. m must be an integer greater than or equal to 1. WARNING: No error or warning produced if this condition is violated.
        ignore_start_r1: How many bases to ignore from start (5' end) of read.
        ignore_end_r1: How many bases to ignore from end (3' end) of read.
        min_qual: Ignore bases with quality-score less than this value.
        phred_offset: The offset in the Phred scores. Phred33 corresponds to phred_offset = 33 and Phred64 corresponds to phred_offset 64.
        ob_strand_offset: How many bases a methylation loci on the OB-strand must be moved to the left in order to line up with the C on the OT-strand; e.g. ob_strand_offset = 1 for CpGs.
    Returns:
        methylation_m_tuples: An updated version of methylation_m_tuples
        n_methylation_loci: The number of methylation loci extracted from the read.
    """
    # Identify methylation events in read, e.g. CpGs or CHHs. The methylation_pattern is specified by a command line argument (e.g. Z/z corresponds to CpG)
    methylation_index = [midx.start() for midx in re.finditer(methylation_pattern, read.opt('XM'))]
    # Ignore any start or end positions of read if required
    if ignore_start_r1 > 0:
        methylation_index = ignore_first_n_bases(read, methylation_index, ignore_start_r1)
    if ignore_end_r1 > 0:
        methylation_index = ignore_last_n_bases(read, methylation_index, ignore_end_r1)
    # Ignore any positions with a base quality less than min_qual
    methylation_index = ignore_low_quality_bases(read, methylation_index, min_qual, phred_offset)
    n_methylation_loci = len(methylation_index)
    # Case A: >= m methylation loci in the read
    if n_methylation_loci >= m:
        positions = [read.pos + x + 1 for x in methylation_index] # +1 to transform from 0-based to 1-based co-ordinates.
        # If read is aligned to OB-strand then translate co-ordinate "ob_strand_offset" bases to the left so that it points to the C on the OT-strand of the methylation locus.
        if read.opt('XG') == 'GA' and read.opt('XR') == 'CT':
            positions = [x - ob_strand_offset for x in positions]
        # Exit if methylation loci are incorrectly ordered
        if not positions == sorted(positions):
            exit_msg = ' '.join(["ERROR: The positions of the methylation loci are not properly ordered for single-end read", read.qname, "\n'positions' =", str(positions), '.\nPlease log an issue at www.github.com/PeteHaitch/Comethylation describing the error or email me at peter.hickey@gmail.com'])
            sys.exit(exit_msg)
        # Else, construct each bookended methylation-loci m-tuple and add it to the methylation_m_tuple instance.
        else:
            for i in range(0, len(methylation_index) - m + 1): # For a read containing k methylation loci there are (k - m + 1) m-tuples.
                this_m_tuple_positions = positions[i:(i + m)]
                # Create a unique ID for each m-tuple of methylation loci (of form "chromosome:position_1-position_2-...-position_M")
                m_tuple_id = ''.join([BAM.getrname(read.tid), ':', '-'.join([str(j) for j in this_m_tuple_positions])])
                # Check whether m-tuple has already been observed. If not, create a WithinFragmentMethylationMTuple instance for it and increment its count. Otherwise, just increment its count.
                if not m_tuple_id in methylation_m_tuples:
                    methylation_m_tuples[m_tuple_id] = WithinFragmentComethylationMTuple(BAM.getrname(read.tid), read.tid, m, this_m_tuple_positions, methylation_type) # read.tid acts as the chromosome_index required by WithinFragmentComethylationMTuple class
                    methylation_m_tuples[m_tuple_id].increment_count(''.join([read.opt('XM')[j] for j in methylation_index[i:(i + m)]]), read, None)
                else:
                    methylation_m_tuples[m_tuple_id].increment_count(''.join([read.opt('XM')[j] for j in methylation_index[i:(i + m)]]), read, None)
    return methylation_m_tuples, n_methylation_loci

def extract_and_update_methylation_index_from_paired_end_reads(read_1, read_2, BAM, methylation_m_tuples, m, methylation_type, methylation_pattern, ignore_start_r1, ignore_start_r2, ignore_end_r1, ignore_end_r2, min_qual, phred_offset, ob_strand_offset, overlap_check, n_fragment_skipped_due_to_bad_overlap, FAILED_QC):
    """Extracts m-tuples of methylation loci from a readpair and adds the comethylation m-tuple to the methylation_m_tuples object.
    
    Args:
        read_1: An AlignedRead instance corresponding to read_1 of the readpair.
        read_2: An AlignedRead instance corresponding to read_2 of the readpair.
        BAM: The Samfile instance corresponding to the sample. Required in order to extract chromosome names from read.
        methylation_m_tuples: A dictionary storing all observed m-tuples of methylation events and their WithinFragmentComethylationMTuple instance.
        m: Is the "m" in "m-tuple", i.e. the size of the m-tuple. m must be an integer greater than or equal to 1. WARNING: No error or warning produced if this condition is violated.
        methylation_type: A string of the methylation type, e.g. CG for CpG methylation. Must be a valid option for the WithinFragmentComethylationMTuple class.
        methylation_pattern: A regular expression of the methylation loci, e.g. '[Zz]' for CpG-methylation
        ignore_start_r1: How many bases to ignore from start (5' end) of read_1.
        ignore_start_r2: How many bases to ignore from start (5' end) of read_2.
        ignore_end_r1: How many bases to ignore from end (3' end) of read_1.
        ignore_end_r2: How many bases to ignore from end (3' end) of read_2.
        min_qual: Ignore bases with quality-score less than this value.
        phred_offset: The offset in the Phred scores. Phred33 corresponds to phred_offset = 33 and Phred64 corresponds to phred_offset 64.
        ob_strand_offset: How many bases a methylation loci on the OB-strand must be moved to the left in order to line up with the C on the OT-strand; e.g. ob_strand_offset = 1 for CpGs.
        overlap_check: The type of check to be performed (listed by most-to-least stringent): check the entire overlapping sequence is identical (sequence), check the XM-tag is identical for the overlapping region (XM), do no check of the overlapping bases but use the read with the higher quality basecalls in the overlapping region (quality), or simply use the overlapping bases from read_1 ala bismark_methylation_extractor (bismark)
        n_fragment_skipped_due_to_bad_overlap: The total number of fragments (read-pairs) skipped due to the overlapping sequencing not passing the filter.
        FAILED_QC: The file object where the QNAME of readpairs that fail the overlap check are written, along with the reason the readpairs failed
    Returns:
        methylation_m_tuples: An updated version of methylation_m_tuples
        n_methylation_loci: The number of methylation loci extracted from the read.
    """
    # Identify methylation events in read, e.g. CpGs or CHHs. The methylation_pattern is specified by a command line argument (e.g. Z/z corresponds to CpG)
    methylation_index_1 = [midx.start() for midx in re.finditer(methylation_pattern, read_1.opt('XM'))]
    methylation_index_2 = [midx.start() for midx in re.finditer(methylation_pattern, read_2.opt('XM'))]
    # Ignore any start or end positions of read_1 or read_2 if required
    if ignore_start_r1 > 0:
        methylation_index_1 = ignore_first_n_bases(read_1, methylation_index_1, ignore_start_r1)
    if ignore_end_r1 > 0:
        methylation_index_1 = ignore_last_n_bases(read_1, methylation_index_1, ignore_end_r1)
    if ignore_start_r2 > 0:
        methylation_index_2 = ignore_first_n_bases(read_2, methylation_index_2, ignore_start_r2)
    if ignore_end_r2 > 0:
        methylation_index_2 = ignore_last_n_bases(read_2, methylation_index_2, ignore_end_r2)
    # Ignore any positions with a base quality less than  min_qual
    methylation_index_1 = ignore_low_quality_bases(read_1, methylation_index_1, min_qual, phred_offset)
    methylation_index_2 = ignore_low_quality_bases(read_2, methylation_index_2, min_qual, phred_offset)

    # Check for overlapping reads from a readpair.
    # If reads overlap check whether the overlapping sequence passes the filter given by overlap_check.
    # If the overlapping sequence does not pass the filter report a warning, increment a counter and skip the readpair (by setting methylation_index_1 and methylation_index_2 to be the empty list).
    n_overlap = read_1.alen + read_2.alen - abs(read_1.tlen)
    if n_overlap > 0:
        if is_overlapping_sequence_identical(read_1, read_2, n_overlap, overlap_check):
            methylation_index_1, methylation_index_2 = ignore_overlapping_sequence(read_1, read_2, methylation_index_1, methylation_index_2, n_overlap, overlap_check)
        else:
            failed_read_msg = '\t'.join([read_1.qname, ''.join(['failed the --overlappingPairedEndFilter ', overlap_check, '\n'])])
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
            exit_msg = ''.join(['ERROR: For readpair ', read.qname, ', position_1 and position_2 contain a common position. This should not happen.\nPlease log an issue at www.github.com/PeteHaitch/Comethylation describing the error or email me at peter.hickey@gmail.com'])
            print positions_1
            print positions_2
            sys.exit(exit_msg)
        # Case 1: Readpair aligns to OT-strand
        if read_1.opt('XG') == 'CT' and read_2.opt('XG') == 'CT' and read_1.opt('XR') == 'CT' and read_2.opt('XR') == 'GA':
            # Exit if methylation loci are incorrectly ordered
            if not positions_1 + positions_2 == sorted(positions_1 + positions_2):
                exit_msg = ' '.join(["ERROR: The positions of the methylation loci are not properly ordered for paired-end read", read_1.qname, "which is aligned to the OT-strand.\n'positions_1 + positions_2' =", str(positions_1 + positions_2), '\nPlease log an issue at www.github.com/PeteHaitch/Comethylation describing the error or email me at peter.hickey@gmail.com'])
                sys.exit(exit_msg)
            else:
                # First, create all m-tuples of methylation loci where each locus is from read_1.
                if len(methylation_index_1) >= m:
                    for i in range(0, len(methylation_index_1) - m + 1): # For a read containing k methylation loci there are (k - m + 1) m-tuples.:
                        this_m_tuple_positions_1 = positions_1[i:(i + m)]
                        # Create a unique ID for each m-tuple of methylation loci (of form "chromosome:position_1-position_2-...-position_m")
                        m_tuple_id = ''.join([BAM.getrname(read_1.tid), ':', '-'.join([str(j) for j in this_m_tuple_positions_1])])
                        # Check whether m-tuple has already been observed. If not, create a WithinFragmentMethylationMTuple instance for it and increment its count. Otherwise, just increment its count.
                        if not m_tuple_id in methylation_m_tuples:
                            methylation_m_tuples[m_tuple_id] = WithinFragmentComethylationMTuple(BAM.getrname(read_1.tid), read_1.tid, m, this_m_tuple_positions_1, methylation_type) # read_1.tid acts as the chromosome_index required by WithinFragmentComethylationMTuple class
                            methylation_m_tuples[m_tuple_id].increment_count(''.join([read_1.opt('XM')[j] for j in methylation_index_1[i:(i + m)]]), read_1, read_2)
                        else:
                            methylation_m_tuples[m_tuple_id].increment_count(''.join([read_1.opt('XM')[j] for j in methylation_index_1[i:(i + m)]]), read_1, read_2)
                # Second, create all m-tuples of methylation loci where the leftmost locus is on read_1 and the rightmost locus is on read_2
                num_shared_m_tuples = max(len(methylation_index_1) + len(methylation_index_2) - m + 1, 0) - max(len(methylation_index_1) - m + 1, 0) - max(len(methylation_index_2) - m + 1, 0) # the number of m-tuples that span read_1 and read_2
                leftmost_shared_locus_index = max(0, len(methylation_index_1) - m + 1) # The index of the leftmost locus to be part of a "shared" m-tuple. The rightmost_shared_locus_index = min(m - 2, len(methylation_index_2) - 1), however this is not required
                for i in range(0, num_shared_m_tuples):
                    this_m_tuple_positions_1 = positions_1[(leftmost_shared_locus_index + i):]
                    this_m_tuple_positions_2 = positions_2[:(m - len(this_m_tuple_positions_1))]
                    # Exit if methylation loci are incorrectly ordered. While a similar check is performed a few lines above, this is a sanity check to make sure than nothing has gone wrong in constructing the shared m-tuples
                    if not this_m_tuple_positions_1 + this_m_tuple_positions_2 == sorted(this_m_tuple_positions_1 + this_m_tuple_positions_2):
                        exit_msg = ' '.join(["ERROR: The positions of the shared methylation loci are not properly ordered for paired-end read", read_1.qname, "which is aligned to the OT-strand.\n'this_m_tuple_positions_1 + this_m_tuple_positions_2' =", str(this_m_tuple_positions_1 + this_m_tuple_positions_2), '\nPlease log an issue at www.github.com/PeteHaitch/Comethylation describing the error or email me at peter.hickey@gmail.com'])
                        sys.exit(exit_msg)
                    else:
                        # Create a unique ID for each m-tuple of methylation loci (of form "chromosome:position_1-position_2-...-position_m")
                        m_tuple_id = ''.join([BAM.getrname(read_1.tid), ':', '-'.join([str(j) for j in this_m_tuple_positions_1] + [str(k) for k in this_m_tuple_positions_2])])
                    # Check whether m-tuple has already been observed. If not, create a WithinFragmentMethylationMTuple instance for it and increment its count. Otherwise, just increment its count.
                        if not m_tuple_id in methylation_m_tuples:
                            methylation_m_tuples[m_tuple_id] = WithinFragmentComethylationMTuple(BAM.getrname(read_1.tid), read_1.tid, m, this_m_tuple_positions_1 + this_m_tuple_positions_2, methylation_type) # read_1.tid acts as the chromosome_index required by WithinFragmentComethylationMTuple class
                            methylation_m_tuples[m_tuple_id].increment_count(''.join([read_1.opt('XM')[j] for j in methylation_index_1[(leftmost_shared_locus_index + i):]] + [read_2.opt('XM')[j] for j in methylation_index_2[:(m - len(this_m_tuple_positions_1))]]), read_1, read_2)
                        else:
                            methylation_m_tuples[m_tuple_id].increment_count(''.join([read_1.opt('XM')[j] for j in methylation_index_1[(leftmost_shared_locus_index + i):]] + [read_2.opt('XM')[j] for j in methylation_index_2[:(m - len(this_m_tuple_positions_1))]]), read_1, read_2)
                # Finally, create all m-tuples of methylation loci where each locus is from read_2.        
                if len(methylation_index_2) >= m:
                    for i in range(0, len(methylation_index_2) - m + 1): # For a read containing k methylation loci there are (k - m + 1) m-tuples.:
                        this_m_tuple_positions_2 = positions_2[i:(i + m)]
                        # Create a unique ID for each m-tuple of methylation loci (of form "chromosome:position_1-position_2-...-position_m")
                        m_tuple_id = ''.join([BAM.getrname(read_2.tid), ':', '-'.join([str(j) for j in this_m_tuple_positions_2])])
                        # Check whether m-tuple has already been observed. If not, create a WithinFragmentMethylationMTuple instance for it and increment its count.
                        if not m_tuple_id in methylation_m_tuples:
                            methylation_m_tuples[m_tuple_id] = WithinFragmentComethylationMTuple(BAM.getrname(read_2.tid), read_2.tid, m, this_m_tuple_positions_2, methylation_type) # read_2.tid acts as the chromosome_index required by WithinFragmentComethylationMTuple class
                            methylation_m_tuples[m_tuple_id].increment_count(''.join([read_2.opt('XM')[j] for j in methylation_index_2[i:(i + m)]]), read_1, read_2)
                        else:
                            methylation_m_tuples[m_tuple_id].increment_count(''.join([read_2.opt('XM')[j] for j in methylation_index_2[i:(i + m)]]), read_1, read_2)
        # Case 2: Readpair aligns to OB-strand
        elif read_1.opt('XG') == 'GA' and read_2.opt('XG') == 'GA' and read_1.opt('XR') == 'CT' and read_2.opt('XR') == 'GA':
            # Translate co-ordinates "ob_strand_offset" bases to the left so that it points to the C on the OT-strand of the methylation locus
            positions_1 = [x - ob_strand_offset for x in positions_1]
            positions_2 = [x - ob_strand_offset for x in positions_2]
            # Exit if methylation loci are incorrectly ordered.
            if not positions_2 + positions_1 == sorted(positions_2 + positions_1): 
                exit_msg = ' '.join(["ERROR: The positions of the methylation loci are not properly ordered for paired-end read", read_1.qname, "which is aligned to the OB-strand.\n'positions_2 + positions_1' =", str(positions_2 + positions_1), '\nPlease log an issue at www.github.com/PeteHaitch/Comethylation describing the error or email me at peter.hickey@gmail.com'])
                sys.exit(exit_msg)
            else:
                # First, create all m-tuples of methylation loci where each locus is from read_1.
                if len(methylation_index_1) >= m:
                    for i in range(0, len(methylation_index_1) - m + 1): # For a read containing m methylation loci there are (m - m-tuple + 1) m-tuples.:
                        this_m_tuple_positions_1 = positions_1[i:(i + m)]
                        # Create a unique ID for each m-tuple of methylation loci (of form "chromosome:position_1-position_2-...-position_m")
                        m_tuple_id = ''.join([BAM.getrname(read_1.tid), ':', '-'.join([str(j) for j in this_m_tuple_positions_1])])
                        # Check whether m-tuple has already been observed. If not, create a WithinFragmentMethylationMTuple instance for it and increment its count.
                        if not m_tuple_id in methylation_m_tuples:
                            methylation_m_tuples[m_tuple_id] = WithinFragmentComethylationMTuple(BAM.getrname(read_1.tid), read_1.tid, m, this_m_tuple_positions_1, methylation_type) # read_1.tid acts as the chromosome_index required by WithinFragmentComethylationMTuple class
                            methylation_m_tuples[m_tuple_id].increment_count(''.join([read_1.opt('XM')[j] for j in methylation_index_1[i:(i + m)]]), read_1, read_2)
                        else:
                            methylation_m_tuples[m_tuple_id].increment_count(''.join([read_1.opt('XM')[j] for j in methylation_index_1[i:(i + m)]]), read_1, read_2)
                # Second, create all m-tuples of methylation loci where the leftmost locus is on read_1 and the rightmost locus is on read_2
                num_shared_m_tuples = max(len(methylation_index_1) + len(methylation_index_2) - m + 1, 0) - max(len(methylation_index_1) - m + 1, 0) - max(len(methylation_index_2) - m + 1, 0) # the number of m-tuples that span read_1 and read_2
                leftmost_shared_locus_index = max(0, len(methylation_index_2) - m + 1) # The index of the leftmost locus to be part of a "shared" m-tuple. The rightmost_shared_locus_index = min(m - 2, len(methylation_index_1) - 1), however this is not required (m - 2 = m - 1 - 1, because Python lists are 0-indexed)
                for i in range(0, num_shared_m_tuples):
                    this_m_tuple_positions_2 = positions_2[(leftmost_shared_locus_index + i):]
                    this_m_tuple_positions_1 = positions_1[:(m - len(this_m_tuple_positions_2))]
                    # Exit if methylation loci are incorrectly ordered. While a similar check is performed a few lines above, this is a sanity check to make sure than nothing has gone wrong in constructing the shared m-tuples
                    if not this_m_tuple_positions_2 + this_m_tuple_positions_1 == sorted(this_m_tuple_positions_2 + this_m_tuple_positions_1):
                        exit_msg = ' '.join(["ERROR: The positions of the shared methylation loci are not properly ordered for paired-end read", read_1.qname, "which is aligned to the OB-strand.\n'this_m_tuple_positions_2 + this_m_tuple_positions_1' =", str(this_m_tuple_positions_2 + this_m_tuple_positions_1), '\nPlease log an issue at www.github.com/PeteHaitch/Comethylation describing the error or email me at peter.hickey@gmail.com'])
                        sys.exit(exit_msg)
                    else:
                        # Create a unique ID for each m-tuple of methylation loci (of form "chromosome:position_1-position_2-...-position_m")
                        m_tuple_id = ''.join([BAM.getrname(read_1.tid), ':', '-'.join([str(j) for j in this_m_tuple_positions_2] + [str(k) for k in this_m_tuple_positions_1])])
                    # Check whether m-tuple has already been observed. If not, create a WithinFragmentMethylationMTuple instance for it and increment its count.
                        if not m_tuple_id in methylation_m_tuples:
                            methylation_m_tuples[m_tuple_id] = WithinFragmentComethylationMTuple(BAM.getrname(read_1.tid), read_1.tid, m, this_m_tuple_positions_2 + this_m_tuple_positions_1, methylation_type) # read_1.tid acts as the chromosome_index required by WithinFragmentComethylationMTuple class
                            methylation_m_tuples[m_tuple_id].increment_count(''.join([read_2.opt('XM')[j] for j in methylation_index_2[(leftmost_shared_locus_index + i):]] + [read_1.opt('XM')[j] for j in methylation_index_1[:(m - len(this_m_tuple_positions_2))]]), read_1, read_2)
                        else:
                            methylation_m_tuples[m_tuple_id].increment_count(''.join([read_2.opt('XM')[j] for j in methylation_index_2[(leftmost_shared_locus_index + i):]] + [read_1.opt('XM')[j] for j in methylation_index_1[:(m - len(this_m_tuple_positions_2))]]), read_1, read_2)
                # Finally, create all m-tuples of methylation loci where each locus is from read_2.        
                if len(methylation_index_2) >= m:
                    for i in range(0, len(methylation_index_2) - m + 1): # For a read containing m methylation loci there are (m - m-tuple + 1) m-tuples.:
                        this_m_tuple_positions_2 = positions_2[i:(i + m)]
                        # Create a unique ID for each m-tuple of methylation loci (of form "chromosome:position_1-position_2-...-position_m")
                        m_tuple_id = ''.join([BAM.getrname(read_2.tid), ':', '-'.join([str(j) for j in this_m_tuple_positions_2])])
                        # Check whether m-tuple has already been observed. If not, create a WithinFragmentMethylationMTuple instance for it and increment its count.
                        if not m_tuple_id in methylation_m_tuples:
                            methylation_m_tuples[m_tuple_id] = WithinFragmentComethylationMTuple(BAM.getrname(read_2.tid), read_2.tid, m, this_m_tuple_positions_2, methylation_type) # read_2.tid acts as the chromosome_index required by WithinFragmentComethylationMTuple class
                            methylation_m_tuples[m_tuple_id].increment_count(''.join([read_2.opt('XM')[j] for j in methylation_index_2[i:(i + m)]]), read_1, read_2)
                        else:
                            methylation_m_tuples[m_tuple_id].increment_count(''.join([read_2.opt('XM')[j] for j in methylation_index_2[i:(i + m)]]), read_1, read_2)
    return methylation_m_tuples, n_methylation_loci, n_fragment_skipped_due_to_bad_overlap
                                    
def write_methylation_m_tuples_to_file(methylation_m_tuples, m, OUT):
    """Write the methylation_m_tuples instance to a tab-separated file. The m-tuples are ordered by chromosome and genomic co-ordinates.
    
    Args:
        methylation_m_tuples: A dictionary storing all observed m-tuples of methylation events and their WithinFragmentComethylationMTuple instance.
        m:
        OUT: The file to write output to.
    """
    # tab_writer writes a tab-separated output file to the filehandle OUT
    tab_writer = csv.writer(OUT, delimiter='\t', quotechar=' ', quoting=csv.QUOTE_MINIMAL)

    # Create the header row
    header = ['chr'] + ['pos' + str(i) for i in range(1, m + 1)]
    example_row = WithinFragmentComethylationMTuple('chr1', 1, m, [0] * m, 'CG').counts
    for key in sorted(example_row.iterkeys()):
        header.append(key)
    # Write the header to file
    tab_writer.writerow(header)
    # Write each m-tuple of methylation loci to file using a nested for-loop
    # (1) Loop over chromosomes by increasing chromosome_index order
    # (2) Loop over positions by increasing order
    for this_m_tuple in sorted(methylation_m_tuples.values(), key = operator.attrgetter('chromosome_index', 'positions')):
        this_m_tuple_ordered_counts = []
        for i in sorted(this_m_tuple.counts.iterkeys()):
            this_m_tuple_ordered_counts.append(this_m_tuple.counts[i])
        row = [this_m_tuple.chromosome] + this_m_tuple.positions + this_m_tuple_ordered_counts
        tab_writer.writerow(row)



__all__ = [
    'ignore_first_n_bases',
    'ignore_last_n_bases',
    'ignore_low_quality_bases',
    'fix_old_bismark',
    'is_overlapping_sequence_identical',
    'does_read_contain_indel',
    'does_read_contain_complicated_cigar',
    'ignore_overlapping_sequence',
    'extract_and_update_methylation_index_from_single_end_read',
    'extract_and_update_methylation_index_from_paired_end_reads',
    'write_methylation_m_tuples_to_file'
]
