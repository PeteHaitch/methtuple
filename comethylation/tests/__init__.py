'''unit testing code for comethylation.
'''

import unittest
import pysam
import sys
import tempfile
import os
import re
import array

from comethylation import *

from comethylation.mtuple import *
from comethylation.funcs import *

class TestIgnoreCycles(unittest.TestCase):
	'''Test the function ignore_read_pos
	'''

	def setUp(self):

		def buildOTRead1():
			'''build an example read_1 aligned to OT-strand.
			'''

			read = pysam.AlignedRead()
			read.qname = "ADS-adipose_chr1_8"
			read.seq = "AATTTTAATTTTAATTTTTGCGGTATTTTTAGTCGGTTCGTTCGTTCGGGTTTGATTTGAG"
			read.flag = 99
			read.tid = 0
			read.pos = 450
			read.mapq = 255
			read.cigar = [(0,61)]
			read.rnext = 1
			read.pnext = 512
			read.isize = 121
			read.qual = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
			read.tags = read.tags + [("XG", "CT")] + [("XM", "..hhh...hhh...hhh.z.Z....hhh.x..xZ..hxZ.hxZ.hxZ....x...hx....")] + [("XR", "CT")]
			return read

		def buildOTRead2():
			'''build an example read_2 aligned to OT-strand.
			'''

			read = pysam.AlignedRead()
			read.qname = "ADS-adipose_chr1_8"
			read.seq = "AGAATTGTGTTTCGTTTTTAGAGTATTATCGAAATTTGTGTAGAGGATAACGTAGCTTC"
			read.flag = 147
			read.tid = 0
			read.pos = 512
			read.mapq = 255
			read.cigar = [(0,59)]
			read.rnext = 1
			read.pnext = 450
			read.isize = -121
			read.qual = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
			read.tags = read.tags + [("XG", "CT")] + [("XM", "....x....h.xZ.hh..x......hh.xZ.....x....x......h..Z.x..H.xZ")] + [("XR", "GA")]
			return read

		def buildOBRead1():
			'''build an example read_1 aligned to OB-strand
			'''

			read = pysam.AlignedRead()
			read.qname = "ADS-adipose_chr1_22929891"
			read.seq = "AACGCAACTCCGCCCTCGCGATACTCTCCGAATCTATACTAAAAAAAACGCAACTCCGCCGAC"
			read.flag = 83
			read.tid = 0
			read.pos = 560
			read.mapq = 255
			read.cigar = [(0,63)]
			read.rnext = 1
			read.pnext = 492
			read.isize = -131
			read.qual = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
			read.tags = read.tags + [("XG", "GA")] + [("XM", "...Z..x....Z.....Z.Zx.h......Zxh...x.h..x.hh.h...Z.......Z..Zx.")] + [("XR", "CT")]
			return read

		def buildOBRead2():
			'''build an example read_2 aligned to OB-strand.
			'''

			read = pysam.AlignedRead()
			read.qname = "ADS-adipose_chr1_22929891"
			read.seq = "CACCCGAATCTAACCTAAAAAAAACTATACTCCGCCTTCAAAATACCACCGAAATCTATACAAAAAA"
			read.flag = 163
			read.tid = 0
			read.pos = 492
			read.mapq = 255
			read.cigar = [(0,67)]
			read.rnext = 1
			read.pnext = 560
			read.isize = 131
			read.qual = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
			read.tags = read.tags + [("XG", "GA")] + [("XM", ".z...Zxh...x....x.hh.h....x.h....Z......x.h.......Z......x.h..x.hh.")] + [("XR", "GA")]
			return read

		# Create the reads and methylation indexes (using CpGs)
		self.otr_1 = buildOTRead1()
		self.otr_2 = buildOTRead2()
		self.obr_1 = buildOBRead1()
		self.obr_2 = buildOBRead2()
		self.otm_1 = [18, 20, 33, 38, 42, 46]
		self.otm_2 = [12, 29, 50, 58]
		self.obm_1 = [3, 11, 17, 19, 29, 49, 57, 60]
		self.obm_2 = [1, 5, 33, 50]

	def test_n_0(self):
		# Shouldn't change methylation indexes
		self.assertEqual(ignore_read_pos(self.otr_1, self.otm_1, []), [18, 20, 33, 38, 42, 46])
		self.assertEqual(ignore_read_pos(self.otr_2, self.otm_2, []), [12, 29, 50, 58])
		self.assertEqual(ignore_read_pos(self.obr_1, self.obm_1, []), [3, 11, 17, 19, 29, 49, 57, 60])
		self.assertEqual(ignore_read_pos(self.obr_2, self.obm_2, []), [1, 5, 33, 50])

	def test_n_off_by_one(self):
		# Shouldn't change methylation indexes: Ignore from "start" of read
		self.assertEqual(ignore_read_pos(self.otr_1, self.otm_1, [17]), [18, 20, 33, 38, 42, 46])
		self.assertEqual(ignore_read_pos(self.otr_2, self.otm_2, [self.otr_2.qlen - 59 - 1]), [12, 29, 50, 58])
		self.assertEqual(ignore_read_pos(self.obr_1, self.obm_1, [self.obr_1.qlen - 61 - 1]), [3, 11, 17, 19, 29, 49, 57, 60])
		self.assertEqual(ignore_read_pos(self.obr_2, self.obm_2, [0]), [1, 5, 33, 50])
		# Shouldn't change methylation indexes: Ignore from "end" of read
		self.assertEqual(ignore_read_pos(self.otr_1, self.otm_1, [47]), [18, 20, 33, 38, 42, 46])
		self.assertEqual(ignore_read_pos(self.otr_2, self.otm_2, [self.otr_2.qlen - 11 - 1]), [12, 29, 50, 58])
		self.assertEqual(ignore_read_pos(self.obr_1, self.obm_1, [self.obr_1.qlen - 2 - 1]), [3, 11, 17, 19, 29, 49, 57, 60])
		self.assertEqual(ignore_read_pos(self.obr_2, self.obm_2, [51]), [1, 5, 33, 50])

	def test_ignore_one(self):
		# Should remove one element from methylation indexes in accordance with the strand/orientation of the read
		self.assertEqual(ignore_read_pos(self.otr_1, self.otm_1, [18]), [20, 33, 38, 42, 46])
		self.assertEqual(ignore_read_pos(self.otr_2, self.otm_2, [self.otr_2.qlen - 58 - 1]), [12, 29, 50])
		self.assertEqual(ignore_read_pos(self.obr_1, self.obm_1, [self.obr_1.qlen - 60 - 1]), [3, 11, 17, 19, 29, 49, 57])
		self.assertEqual(ignore_read_pos(self.obr_2, self.obm_2, [1]), [5, 33, 50])
		# Should remove one element from methylation indexes in accordance with the strand/orientation of the read
		self.assertEqual(ignore_read_pos(self.otr_1, self.otm_1, [46]), [18, 20, 33, 38, 42])
		self.assertEqual(ignore_read_pos(self.otr_2, self.otm_2, [self.otr_2.qlen - 12 - 1]), [29, 50, 58])
		self.assertEqual(ignore_read_pos(self.obr_1, self.obm_1, [self.obr_1.qlen - 3 - 1]), [11, 17, 19, 29, 49, 57, 60])
		self.assertEqual(ignore_read_pos(self.obr_2, self.obm_2, [50]), [1, 5, 33])

	def test_ignore_all(self):
		# Should remove all elements from methylation indexes (assuming read-lengths are < 10,000)
		self.assertEqual(ignore_read_pos(self.otr_1, self.otm_1, range(0, 10000)), [])
		self.assertEqual(ignore_read_pos(self.otr_2, self.otm_2, range(0, 10000)), [])
		self.assertEqual(ignore_read_pos(self.obr_1, self.obm_1, range(0, 10000)), [])
		self.assertEqual(ignore_read_pos(self.obr_2, self.obm_2, range(0, 10000)), [])

class TestIgnoreLowQualityBases(unittest.TestCase):
	'''Test the function ignore_low_quality_bases
	'''

	def setUp(self):

		def buildPhred33():
			'''build an example read_1 aligned to OT-strand.
			'''

			read = pysam.AlignedRead()
			read.qname = "ADS-adipose_chr1_8"
			read.seq = "AATTTTAATTTTAATTTTTGCGGTATTTTTAGTCGGTTCGTTCGTTCGGGTTTGATTTGAG"
			read.flag = 99
			read.tid = 0
			read.pos = 450
			read.mapq = 255
			read.cigar = [(0,61)]
			read.rnext = 1
			read.pnext = 512
			read.isize = 121
			read.qual = "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" # Phred33 = 0
			read.tags = read.tags + [("XG", "CT")] + [("XM", "..hhh...hhh...hhh.z.Z....hhh.x..xZ..hxZ.hxZ.hxZ....x...hx....")] + [("XR", "CT")]
			return read

		def buildPhred64():
			'''build an example read_1 aligned to OT-strand.
			'''

			read = pysam.AlignedRead()
			read.qname = "ADS-adipose_chr1_8"
			read.seq = "AATTTTAATTTTAATTTTTGCGGTATTTTTAGTCGGTTCGTTCGTTCGGGTTTGATTTGAG"
			read.flag = 99
			read.tid = 0
			read.pos = 450
			read.mapq = 255
			read.cigar = [(0,61)]
			read.rnext = 1
			read.pnext = 512
			read.isize = 121
			read.qual = "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" # Phred64 = 0
			read.tags = read.tags + [("XG", "CT")] + [("XM", "..hhh...hhh...hhh.z.Z....hhh.x..xZ..hxZ.hxZ.hxZ....x...hx....")] + [("XR", "CT")]
			return read


		# Create the reads and methylation indexes (using CpGs)
		self.p33 = buildPhred33()
		self.p64 = buildPhred64()
		self.p33m = [18, 20, 33, 38, 42, 46]
		self.p64m = [18, 20, 33, 38, 42, 46]

	def test_no_low_qual_filter(self):
		# Shouldn't change methylation indexes
		self.assertEqual(ignore_low_quality_bases(self.p33, self.p33m, 0, 33), [18, 20, 33, 38, 42, 46])
		self.assertEqual(ignore_low_quality_bases(self.p64, self.p64m, 0, 64), [18, 20, 33, 38, 42, 46])

	def test_all_fail_low_qual_filter(self):
		# Should remove all elements from methylation indexes 
		self.assertEqual(ignore_low_quality_bases(self.p33, self.p33m, 1, 33), [])
		self.assertEqual(ignore_low_quality_bases(self.p64, self.p64m, 1, 64), [])

	def test_some_fail_low_qual_filter(self):
		# Should remove all but the first element of the methylation index
		self.p33.qual = "!!!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" # Change one base to have Phred33 = 2
		self.p64.qual = "@@@@@@@@@@@@@@@@@@B@#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" # Change one base to have Phred64 = 2
		self.assertEqual(ignore_low_quality_bases(self.p33, self.p33m, 2, 33), [18])
		self.assertEqual(ignore_low_quality_bases(self.p64, self.p64m, 2, 64), [18])

	def test_bad_min_qual(self):
		# Should raise an exception
		self.assertRaises(ValueError, ignore_low_quality_bases, self.p33, self.p33m, -10, 33)
		self.assertRaises(ValueError, ignore_low_quality_bases, self.p33, self.p33m, -10, 64)
		self.assertRaises(ValueError, ignore_low_quality_bases, self.p33, self.p33m, 3.4, 33)
		self.assertRaises(ValueError, ignore_low_quality_bases, self.p33, self.p33m, 3.4, 64)

	def test_bad_phred_offset(self):
		self.assertRaises(ValueError, ignore_low_quality_bases, self.p33, self.p33m, 10, 34)

class TestFixOldBismark(unittest.TestCase):
	'''Test the function fix_old_bismark
	'''

	def setUp(self):

		def buildOTRead1():
			'''build an example read_1 aligned to OT-strand.
			'''

			read = pysam.AlignedRead()
			read.qname = "ADS-adipose_chr1_8"
			read.seq = "AATTTTAATTTTAATTTTTGCGGTATTTTTAGTCGGTTCGTTCGTTCGGGTTTGATTTGAG"
			read.flag = 99
			read.tid = 0
			read.pos = 450
			read.mapq = 255
			read.cigar = [(0,61)]
			read.rnext = 1
			read.pnext = 512
			read.isize = 121
			read.qual = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
			#read.tags = read.tags + [("XG", "CT")] + [("XM", "..hhh...hhh...hhh.z.Z....hhh.x..xZ..hxZ.hxZ.hxZ....x...hx....")] + [("XR", "CT")] # Not required for testing fix_old_bismark
			return read

		# Create the read
		self.read = buildOTRead1()

	def test_fix_old_bismark(self):
		self.read.flag = 67
		self.assertEqual(fix_old_bismark(self.read).flag, 99)
		self.read.flag = 115
		self.assertEqual(fix_old_bismark(self.read).flag, 83)
		self.read.flag = 131
		self.assertEqual(fix_old_bismark(self.read).flag, 147)
		self.read.flag = 179
		self.assertEqual(fix_old_bismark(self.read).flag, 163)
		self.read.flag = 2
		with self.assertRaises(SystemExit) as cm:
			fix_old_bismark(self.read)
			self.assertEqual(cm.exception.code, 1)

class TestIsOverlappingSequenceIdentical(unittest.TestCase):
	'''Test the function is_overlapping_sequence_identical
	'''

	def setUp(self):

		def buildRead1():
			'''build an example read_1 aligned to OT-strand.
			'''
			read = pysam.AlignedRead()
			read.qname = "tr"
			read.seq = "TTTTTATTATTAAAGATAGTAGTGTTTTAAGTTTAGTGTTAGAGGTATTTGTTTGTAGTCGAAGTATTTTGTTAAAGTTAGGAGGGTTTAATAAGGTTTG"
			read.flag = 99
			read.tid = 0
			read.pos = 853
			read.mapq = 255
			read.cigar = [(0,100)]
			read.rnext = 0
			read.pnext = 854
			read.isize = 100
			read.qual = "BBCFFBDEHH2AFHIGHIJFHIIIJJJJHHIIIJGIHHJJIJIJJDHIIIJIIJJHIJJJJJJJHIIJJJJJJGIGGJGGGFFHGFBACA@CCCCDCCD@"
			read.tags = read.tags + [("XG", "CT")] + [("XM", "hh..h.....x........x....hh.h....h......x.....h..x...x..x..xZ....h.h.....h.....x.......h.........h.z.")] + [("XR", "CT")]
			return read

		def buildRead2():
			'''build an example read_2 aligned to OT-strand.
			'''
			read = pysam.AlignedRead()
			read.qname = "tr"
			read.seq = "TTTTATTATTAAAGATAGTAGTGTTTTAAGTTTAGTGTTAGAGGTATTTGTTTGTAGTCGAAGTATTTTGTTAAAGTTAGGAGGGTTTAATAAGGTTTGA"
			read.flag = 147
			read.tid = 0
			read.pos = 854
			read.mapq = 255
			read.cigar = [(0,100)]
			read.rnext = 0
			read.pnext = 853
			read.isize = 100
			read.qual = "BCFFBDEHH2AFHIGHIJFHIIIJJJJHHIIIJGIHHJJIJIJJDHIIIJIIJJHIJJJJJJJHIIJJJJJJGIGGJGGGFFHGFBACA@CCCCDCCD@:"
			read.tags = read.tags + [("XG", "CT")] + [("XM", "h..h.....x........x....hh.h....h......x.....h..x...x..x..xZ....h.h.....h.....x.......h.........h.z..")] + [("XR", "GA")]
			return read

		# Create the reads
		self.read_1 = buildRead1()
		self.read_2 = buildRead2()
	
	def test_sequence(self):
		self.assertTrue(is_overlapping_sequence_identical(self.read_1, self.read_2, 99, 'sequence'))
		self.mod_read_1 = self.read_1
		self.mod_read_1.seq = ''.join([self.read_1.seq[:59], 'T', self.read_1.seq[60:]]) # Change a 'C' to a 'T'
		self.assertFalse(is_overlapping_sequence_identical(self.mod_read_1, self.read_2, 99, 'sequence'))

	def test_XM(self):
		self.assertTrue(is_overlapping_sequence_identical(self.read_1, self.read_2, 99, 'XM'))
		self.mod_read_1 = pysam.AlignedRead()
		self.mod_read_1.qname = self.read_1.qname
		self.mod_read_1.seq = self.read_1.seq
		self.mod_read_1.flag = self.read_1.flag
		self.mod_read_1.tid = self.read_1.tid
		self.mod_read_1.pos = self.read_1.pos
		self.mod_read_1.mapq = self.read_1.mapq
		self.mod_read_1.cigar = self.read_1.cigar
		self.mod_read_1.rnext = self.read_1.rnext
		self.mod_read_1.pnext = self.read_1.pnext
		self.mod_read_1.isize = self.read_1.isize
		self.mod_read_1.qual = self.read_1.qual
		self.mod_read_1.tags = self.mod_read_1.tags + [("XG", "CT")] + [("XM", "hh..h.....x........x....hh.h....h......x.....h..x...x..x..xz....h.h.....h.....x.......h.........h.z.")] + [("XR", "CT")] # Change a 'Z' to a 'z' at cycle 60
		self.assertFalse(is_overlapping_sequence_identical(self.mod_read_1, self.read_2, 99, 'XM'))

	def test_quality(self):
		self.assertTrue(is_overlapping_sequence_identical(self.read_1, self.read_2, 99, 'quality'))

	def test_bismark(self):
		self.assertTrue(is_overlapping_sequence_identical(self.read_1, self.read_2, 99, 'bismark'))

	def test_n_overlap(self):
		n_overlap = self.read_1.alen + self.read_2.alen - abs(self.read_1.tlen)
		self.assertEqual(n_overlap, self.read_1.tlen)
		self.assertEqual(n_overlap, self.read_2.tlen)

	def test_bad_n_overlap(self):
		self.assertRaises(ValueError, is_overlapping_sequence_identical, self.read_1, self.read_2, -10, 'sequence')
		self.assertRaises(ValueError, is_overlapping_sequence_identical, self.read_1, self.read_2, 3.4, 'sequence')

	def test_bad_overlap_check(self):
		self.assertRaises(ValueError, is_overlapping_sequence_identical, self.read_1, self.read_2, 10, 'apples')
		self.assertRaises(ValueError, is_overlapping_sequence_identical, self.read_1, self.read_2, 10, 'Bismark') # Should be 'bismark'

	def test_invalid_strands(self):
		self.mod_read_1 = self.read_1
		self.mod_read_1.tags = []
		self.mod_read_1.tags = 	self.mod_read_1.tags + [("XG", "GA")] + [("XM", "hh..h.....x........x....hh.h....h......x.....h..x...x..x..xZ....h.h.....h.....x.......h.........h.z.")] + [("XR", "CT")]
		with self.assertRaises(SystemExit) as cm:
			is_overlapping_sequence_identical(self.mod_read_1, self.read_2, 10, 'sequence')
			self.assertEqual(cm.exception.code, 1)
 
class TestDoesReadContainIndel(unittest.TestCase):
	'''Test the function does_read_contain_indel
	'''

	def setUp(self):

		def buildRead1():
			'''build an example read_1 aligned to OT-strand.
			'''
			read = pysam.AlignedRead()
			read.qname = "tr"
			read.seq = "TTTTTATTATTAAAGATAGTAGTGTTTTAAGTTTAGTGTTAGAGGTATTTGTTTGTAGTCGAAGTATTTTGTTAAAGTTAGGAGGGTTTAATAAGGTTTG"
			read.flag = 99
			read.tid = 0
			read.pos = 853
			read.mapq = 255
			read.cigar = [(0,100)]
			read.rnext = 0
			read.pnext = 854
			read.isize = 100
			read.qual = "BBCFFBDEHH2AFHIGHIJFHIIIJJJJHHIIIJGIHHJJIJIJJDHIIIJIIJJHIJJJJJJJHIIJJJJJJGIGGJGGGFFHGFBACA@CCCCDCCD@"
			read.tags = read.tags + [("XG", "CT")] + [("XM", "hh..h.....x........x....hh.h....h......x.....h..x...x..x..xZ....h.h.....h.....x.......h.........h.z.")] + [("XR", "CT")]
			return read

		# Create the reads
		self.read_1 = buildRead1()

	def test_no_indel(self):
		self.assertFalse(does_read_contain_indel(self.read_1))

	def test_insertion(self):
		self.read_1.cigar = [(0, 50), (1, 50)]
		self.assertTrue(does_read_contain_indel(self.read_1))

	def test_deletion(self):
		self.read_1.cigar = [(0, 50), (2, 50)]
		self.assertTrue(does_read_contain_indel(self.read_1))

class TestDoesReadContainComplicatedCigar(unittest.TestCase):
	'''Test the function does_read_contain_complicated_cigar
	'''

	def setUp(self):

		def buildRead1():
			'''build an example read_1 aligned to OT-strand.
			'''
			read = pysam.AlignedRead()
			read.qname = "tr"
			read.seq = "TTTTTATTATTAAAGATAGTAGTGTTTTAAGTTTAGTGTTAGAGGTATTTGTTTGTAGTCGAAGTATTTTGTTAAAGTTAGGAGGGTTTAATAAGGTTTG"
			read.flag = 99
			read.tid = 0
			read.pos = 853
			read.mapq = 255
			read.cigar = [(0,100)]
			read.rnext = 0
			read.pnext = 854
			read.isize = 100
			read.qual = "BBCFFBDEHH2AFHIGHIJFHIIIJJJJHHIIIJGIHHJJIJIJJDHIIIJIIJJHIJJJJJJJHIIJJJJJJGIGGJGGGFFHGFBACA@CCCCDCCD@"
			read.tags = read.tags + [("XG", "CT")] + [("XM", "hh..h.....x........x....hh.h....h......x.....h..x...x..x..xZ....h.h.....h.....x.......h.........h.z.")] + [("XR", "CT")]
			return read

		# Create the reads
		self.read_1 = buildRead1()

	def test_simple_cigar(self):
		self.assertFalse(does_read_contain_complicated_cigar(self.read_1))

	def test_insertion(self):
		self.read_1.cigar = [(0, 50), (1, 50)]
		self.assertFalse(does_read_contain_complicated_cigar(self.read_1))

	def test_deletion(self):
		self.read_1.cigar = [(0, 50), (2, 50)]
		self.assertFalse(does_read_contain_complicated_cigar(self.read_1))

	def test_ref_skip(self):
		self.read_1.cigar = [(0, 50), (3, 50)]
		self.assertTrue(does_read_contain_complicated_cigar(self.read_1))

	def test_soft_clip(self):
		self.read_1.cigar = [(0, 50), (4, 50)]
		self.assertTrue(does_read_contain_complicated_cigar(self.read_1))

	def test_hard_clip(self):
		self.read_1.cigar = [(0, 50), (5, 50)]
		self.assertTrue(does_read_contain_complicated_cigar(self.read_1))

	def test_pad(self):
		self.read_1.cigar = [(0, 50), (6, 50)]
		self.assertTrue(does_read_contain_complicated_cigar(self.read_1))

	def test_equal(self):
		self.read_1.cigar = [(0, 50), (7, 50)]
		self.assertTrue(does_read_contain_complicated_cigar(self.read_1))

	def test_diff(self):
		self.read_1.cigar = [(0, 50), (8, 50)]
		self.assertTrue(does_read_contain_complicated_cigar(self.read_1))

class TestIgnoreOverlappingSequence(unittest.TestCase):
	'''Test the function ignore_overlapping_sequence
	'''

	def setUp(self):

		def buildOTRead1():
			'''build an example read_1 aligned to OT-strand.
			'''
			read = pysam.AlignedRead()
			read.qname = "otr"
			read.seq = "TTTTTATTATTAAAGATAGTAGTGTTTTAAGTTTAGTGTTAGAGGTATTTGTTTGTAGTCGAAGTATTTTGTTAAAGTTAGGAGGGTTTAATAAGGTTTG"
			read.flag = 99
			read.tid = 0
			read.pos = 853
			read.mapq = 255
			read.cigar = [(0,100)]
			read.rnext = 0
			read.pnext = 854
			read.isize = 100
			read.qual = "BBCFFBDEHH2AFHIGHIJFHIIIJJJJHHIIIJGIHHJJIJIJJDHIIIJIIJJHIJJJJJJJHIIJJJJJJGIGGJGGGFFHGFBACA@CCCCDCCD@"
			read.tags = read.tags + [("XG", "CT")] + [("XM", "hh..h.....x........x....hh.h....h......x.....h..x...x..x..xZ....h.h.....h.....x.......h.........h.z.")] + [("XR", "CT")]
			return read

		def buildOTRead2():
			'''build an example read_2 aligned to OT-strand.
			'''
			read = pysam.AlignedRead()
			read.qname = "otr"
			read.seq = "TTTTATTATTAAAGATAGTAGTGTTTTAAGTTTAGTGTTAGAGGTATTTGTTTGTAGTCGAAGTATTTTGTTAAAGTTAGGAGGGTTTAATAAGGTTTGA"
			read.flag = 147
			read.tid = 0
			read.pos = 854
			read.mapq = 255
			read.cigar = [(0,100)]
			read.rnext = 0
			read.pnext = 853
			read.isize = 100
			read.qual = "BCFFBDEHH2AFHIGHIJFHIIIJJJJHHIIIJGIHHJJIJIJJDHIIIJIIJJHIJJJJJJJHIIJJJJJJGIGGJGGGFFHGFBACA@CCCCDCCD@:"
			read.tags = read.tags + [("XG", "CT")] + [("XM", "h..h.....x........x....hh.h....h......x.....h..x...x..x..xZ....h.h.....h.....x.......h.........h.z..")] + [("XR", "GA")]
			return read

		def buildOBRead1():
			'''build an example read_1 aligned to OB-strand
			'''
			read = pysam.AlignedRead()
			read.qname = "otb"
			read.seq = "ACGCAACTCCGCCCTCGCGATACTCTCCGAATCTATACTAAAAAAAACGCAACTCCGCCGAC"
			read.flag = 83
			read.tid = 1
			read.pos = 493
			read.mapq = 255
			read.cigar = [(0,62)]
			read.rnext = 1
			read.pnext = 492
			read.isize = 63
			read.qual = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
			read.tags = read.tags + [("XG", "GA")] + [("XM", "..Z..x....Z.....Z.Zx.h......Zxh...x.h..x.hh.h...Z.......Z..Zx.")] + [("XR", "CT")]
			return read

		def buildOBRead2():
			'''build an example read_2 aligned to OB-strand.
			'''
			read = pysam.AlignedRead()
			read.qname = "otb"
			read.seq = "AACGCAACTCCGCCCTCGCGATACTCTCCGAATCTATACTAAAAAAAACGCAACTCCGCCGA"
			read.flag = 163
			read.tid = 0
			read.pos = 492
			read.mapq = 255
			read.cigar = [(0,62)]
			read.rnext = 1
			read.pnext = 493
			read.isize = 63
			read.qual = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
			read.tags = read.tags + [("XG", "GA")] + [("XM", "...Z..x....Z.....Z.Zx.h......Zxh...x.h..x.hh.h...Z.......Z..Zx")] + [("XR", "GA")]
			return read

		# Create the reads
		self.read_1 = buildOTRead1()
		self.read_2 = buildOTRead2()

		# Create the reads and methylation indexes (using CpGs)
		self.otr_1 = buildOTRead1()
		self.otr_2 = buildOTRead2()
		self.obr_1 = buildOBRead1()
		self.obr_2 = buildOBRead2()
		self.otm_1 = [59, 98]
		self.otm_2 = [58, 97]
		self.obm_1 = [2, 10, 16, 18, 28, 48, 56, 59]
		self.obm_2 = [3, 11, 17, 19, 29, 49, 57, 60]

	def test_bismark(self):
		self.assertEqual(ignore_overlapping_sequence(self.otr_1, self.otr_2, self.otm_1, self.otm_2, 99, 'bismark'), ([59, 98], []))
		self.assertEqual(ignore_overlapping_sequence(self.obr_1, self.obr_2, self.obm_1, self.obm_2, 99, 'bismark'), ([2, 10, 16, 18, 28, 48, 56, 59], []))
		# Make read_2 have higher quality bases than read_1
		self.otr_1.qual = 'B' * len(self.otr_1.seq)
		self.otr_2.qual = 'K' * len(self.otr_2.seq)
		self.obr_1.qual = 'B' * len(self.obr_1.seq)
		self.obr_2.qual = 'K' * len(self.obr_2.seq)
		self.assertEqual(ignore_overlapping_sequence(self.otr_1, self.otr_2, self.otm_1, self.otm_2, 99, 'bismark'), ([59, 98], []))
		self.assertEqual(ignore_overlapping_sequence(self.obr_1, self.obr_2, self.obm_1, self.obm_2, 99, 'bismark'), ([2, 10, 16, 18, 28, 48, 56, 59], []))

	def test_quality(self):
		self.assertEqual(ignore_overlapping_sequence(self.otr_1, self.otr_2, self.otm_1, self.otm_2, 99, 'quality'), ([59, 98], []))
		self.assertEqual(ignore_overlapping_sequence(self.obr_1, self.obr_2, self.obm_1, self.obm_2, 99, 'quality'), ([2, 10, 16, 18, 28, 48, 56, 59], []))
		# Make read_2 have higher quality bases than read_1
		self.otr_1.qual = 'B' * len(self.otr_1.seq)
		self.otr_2.qual = 'K' * len(self.otr_2.seq)
		self.obr_1.qual = 'B' * len(self.obr_1.seq)
		self.obr_2.qual = 'K' * len(self.obr_2.seq)
		self.assertEqual(ignore_overlapping_sequence(self.otr_1, self.otr_2, self.otm_1, self.otm_2, 99, 'quality'), ([], [58, 97]))
		self.assertEqual(ignore_overlapping_sequence(self.obr_1, self.obr_2, self.obm_1, self.obm_2, 99, 'quality'), ([], [3, 11, 17, 19, 29, 49, 57, 60]))

	def test_XM(self):
		self.assertEqual(ignore_overlapping_sequence(self.otr_1, self.otr_2, self.otm_1, self.otm_2, 99, 'XM'), ([59, 98], []))
		self.assertEqual(ignore_overlapping_sequence(self.obr_1, self.obr_2, self.obm_1, self.obm_2, 99, 'XM'), ([2, 10, 16, 18, 28, 48, 56, 59], []))
		# Make read_2 have higher quality bases than read_1
		self.otr_1.qual = 'B' * len(self.otr_1.seq)
		self.otr_2.qual = 'K' * len(self.otr_2.seq)
		self.obr_1.qual = 'B' * len(self.obr_1.seq)
		self.obr_2.qual = 'K' * len(self.obr_2.seq)
		self.assertEqual(ignore_overlapping_sequence(self.otr_1, self.otr_2, self.otm_1, self.otm_2, 99, 'XM'), ([], [58, 97]))
		self.assertEqual(ignore_overlapping_sequence(self.obr_1, self.obr_2, self.obm_1, self.obm_2, 99, 'XM'), ([], [3, 11, 17, 19, 29, 49, 57, 60]))

	def test_sequence(self):
		self.assertEqual(ignore_overlapping_sequence(self.otr_1, self.otr_2, self.otm_1, self.otm_2, 99, 'sequence'), ([59, 98], []))
		# Make read_2 have higher quality bases than read_1
		self.otr_1.qual = 'B' * len(self.otr_1.seq)
		self.otr_2.qual = 'K' * len(self.otr_2.seq)
		self.obr_1.qual = 'B' * len(self.obr_1.seq)
		self.obr_2.qual = 'K' * len(self.obr_2.seq)
		self.assertEqual(ignore_overlapping_sequence(self.otr_1, self.otr_2, self.otm_1, self.otm_2, 99, 'sequence'), ([], [58, 97]))
		self.assertEqual(ignore_overlapping_sequence(self.obr_1, self.obr_2, self.obm_1, self.obm_2, 99, 'sequence'), ([], [3, 11, 17, 19, 29, 49, 57, 60]))

	def test_bad_n_overlap(self):
		# Should raise an exception
		self.assertRaises(ValueError, ignore_overlapping_sequence, self.otr_1, self.otr_2, self.otm_1, self.otm_2, -10, 'sequence')
		self.assertRaises(ValueError, ignore_overlapping_sequence, self.otr_1, self.otr_2, self.otm_1, self.otm_2, 3.4, 'sequence')

	def test_bad_overlap_check(self):
		# Should raise an exception
		self.assertRaises(ValueError, ignore_overlapping_sequence, self.otr_1, self.otr_2, self.otm_1, self.otm_2, 99, 'apples')
		self.assertRaises(ValueError, ignore_overlapping_sequence, self.otr_1, self.otr_2, self.otm_1, self.otm_2, 99, 'Bismark') # Should be 'bismark'

class TestExtractAndUpdateMethylationIndexFromSingleEndRead(unittest.TestCase):
	'''Test the function extract_and_update_methylation_index_from_single_end_read
	'''
	def setUp(self):

		def buildOTRead():
			'''build an example read aligned to OT-strand.
			'''
			read = pysam.AlignedRead()
			read.qname = "@SALK_2077_FC6295TAAXX:2:107:9396:15019#0/1"
			read.seq = "GGGGAAGGTGTTATGGAGTTTTTTACGATTTTTAGTCGTTTTCGTTTTTTTTTGTTTGTGGTTGTTGCGGTGGCGGTAGAGGAGGG"
			read.flag = 0
			read.tid = 0
			read.pos = 4536
			read.mapq = 255
			read.cigar = [(0,86)]
			read.rnext = 0
			read.pnext = 0
			read.isize = 0
			read.qual = "DBDB2;@>)@@F?EFG@GBGGGGDDBG@DGGGGEEFHHEGHHHHEFHHHHFHHHFHHHGHGBCEAA@?@?/A@>@3,.6,AA,@>="
			read.tags = read.tags + [("XG", "CT")] + [("XM", "...........h......hhhhh..Z....hhx...Z..hh.Z..hh.hh.x..hx.....x..x..Z.....Z..x.........")] + [("XR", "CT")]
			return read

		def buildOBRead():
			'''build an example read aligned to OB-strand.
			'''
			read = pysam.AlignedRead()
			read.qname = "@ECKER_1116_FC623CNAAXX:2:21:18515:1127#0/1"
			read.seq = "CTTCCTAACAAACAACTACACCACTACCTAACGCTATACCCTTCCTTTACTCTACCCACTAAAAACAATATTTATCATAAACCT"
			read.flag = 16
			read.tid = 0
			read.pos = 3334
			read.mapq = 255
			read.cigar = [(0,84)]
			read.rnext = 0
			read.pnext = 0
			read.isize = 0
			read.qual = "G7@G@BGB@GGGGGDIEEBIBA<AHEGEEEGGGDDEDFFEIIHIIGGDGGGGGGGGGGDGDBED<FAAFEGGGGGIHIFIGBDG"
			read.tags = read.tags + [("XG", "GA")] + [("XM", "......x...xh..x..x.......x...xh.Z..x.h..........h....x...z..xh.h..zx.h...h....hhh...")] + [("XR", "CT")]
			return read

		def buildBAM():
			'''build a BAM file
			'''
			header = { 'HD': {'VN': '1.0'}, 'SQ': [{'LN': 10000, 'SN': 'chr1'}, {'LN': 20000, 'SN': 'chr2'}] }
			tempfile_path = tempfile.mkstemp()[1]
			BAM = pysam.Samfile(tempfile_path, "wb", header = header)
			return BAM

		# Create the reads, BAM file and methylation m-tuples
		self.otr = buildOTRead()
		self.obr = buildOBRead()
		self.BAM = buildBAM()
		self.BAM.write(self.otr)
		self.BAM.write(self.obr)
		self.filename = self.BAM.filename
		self.BAM.close()
		self.BAM = pysam.Samfile(self.filename, 'rb')
		self.m1ot, self.nmlifot = extract_and_update_methylation_index_from_single_end_read(read = self.otr, BAM = self.BAM, methylation_m_tuples = MTuple('test', 1, 'CG', {'chr1': 0}), m = 1, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1)
		self.m2ot, self.nmlifot = extract_and_update_methylation_index_from_single_end_read(read = self.otr, BAM = self.BAM, methylation_m_tuples = MTuple('test', 2, 'CG', {'chr1': 0}), m = 2, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1)
		self.m3ot, self.nmlifot = extract_and_update_methylation_index_from_single_end_read(read = self.otr, BAM = self.BAM, methylation_m_tuples = MTuple('test', 3, 'CG', {'chr1': 0}), m = 3, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1)
		self.m4ot, self.nmlifot = extract_and_update_methylation_index_from_single_end_read(read = self.otr, BAM = self.BAM, methylation_m_tuples = MTuple('test', 4, 'CG', {'chr1': 0}), m = 4, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1)
		self.m5ot, self.nmlifot = extract_and_update_methylation_index_from_single_end_read(read = self.otr, BAM = self.BAM, methylation_m_tuples = MTuple('test', 5, 'CG', {'chr1': 0}), m = 5, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1)
		self.m1ob, self.nmlifob = extract_and_update_methylation_index_from_single_end_read(read = self.obr, BAM = self.BAM, methylation_m_tuples = MTuple('test', 1, 'CG', {'chr1': 0}), m = 1, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1)
		self.m2ob, self.nmlifob = extract_and_update_methylation_index_from_single_end_read(read = self.obr, BAM = self.BAM, methylation_m_tuples = MTuple('test', 2, 'CG', {'chr1': 0}), m = 2, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1)
		self.m3ob, self.nmlifob = extract_and_update_methylation_index_from_single_end_read(read = self.obr, BAM = self.BAM, methylation_m_tuples = MTuple('test', 3, 'CG', {'chr1': 0}), m = 3, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1)
		self.m2cgchg, self.nmlifotcgchg = extract_and_update_methylation_index_from_single_end_read(read = self.otr, BAM = self.BAM, methylation_m_tuples = MTuple('test', 2, 'CG/CHG', {'chr1': 0}), m = 2, methylation_type = 'CG/CHG', methylation_pattern = re.compile(r'[ZzXx]'), ignore_read1_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1)

	def test_correct_number_of_m_tuples(self):
		self.assertEqual(len(self.m1ot.mtuples), 5)
		self.assertEqual(len(self.m2ot.mtuples), 4)
		self.assertEqual(len(self.m3ot.mtuples), 3)
		self.assertEqual(len(self.m4ot.mtuples), 2)
		self.assertEqual(len(self.m5ot.mtuples), 1)
	 	self.assertEqual(len(self.m1ob.mtuples), 3)
	 	self.assertEqual(len(self.m2ob.mtuples), 2)
	 	self.assertEqual(len(self.m3ob.mtuples), 1)

 	def test_correct_m_tuple_ids(self):
 		self.assertItemsEqual(self.m1ot.mtuples.keys(), [('chr1', 4562), ('chr1', 4604), ('chr1', 4579), ('chr1', 4573), ('chr1', 4610)])
 		self.assertItemsEqual(self.m2ot.mtuples.keys(), [('chr1', 4579, 4604), ('chr1', 4562, 4573), ('chr1', 4604, 4610), ('chr1', 4573, 4579)])
 		self.assertItemsEqual(self.m3ot.mtuples.keys(), [('chr1', 4562, 4573, 4579), ('chr1', 4573, 4579, 4604), ('chr1', 4579, 4604, 4610)])
 		self.assertItemsEqual(self.m4ot.mtuples.keys(), [('chr1', 4562, 4573, 4579, 4604), ('chr1', 4573, 4579, 4604, 4610)])
 		self.assertItemsEqual(self.m5ot.mtuples.keys(), [('chr1', 4562, 4573, 4579, 4604, 4610)])
 		self.assertItemsEqual(self.m1ob.mtuples.keys(), [('chr1', 3400), ('chr1', 3366), ('chr1', 3391)])
 		self.assertItemsEqual(self.m2ob.mtuples.keys(), [('chr1', 3391, 3400), ('chr1', 3366, 3391)])
 		self.assertItemsEqual(self.m3ob.mtuples.keys(), [('chr1', 3366, 3391, 3400)])

	def test_correct_number_of_methylation_loci_in_fragment(self):
		self.assertEqual(self.nmlifot, 5)
		self.assertEqual(self.nmlifob, 3)

	def test_correct_methylation_type(self):
		self.assertEqual(self.m1ot.methylation_type, 'CG')
		self.assertEqual(self.m1ob.methylation_type, 'CG')
		self.assertEqual(self.m2cgchg.methylation_type, 'CG/CHG')

	def test_counts(self):
		self.assertEqual(self.m1ot.mtuples[('chr1', 4562)], array.array('i', [1, 0]))

	def tearDown(self):
		os.remove(self.BAM.filename)

class TestExtractAndUpdateMethylationIndexFromPairedEndReads(unittest.TestCase):
	'''Test the function extract_and_update_methylation_index_from_single_end_read
	'''
	def setUp(self):

		def buildOTRead1():
			'''build an example read_1 aligned to OT-strand.
			'''

			read = pysam.AlignedRead()
			read.qname = "ADS-adipose_chr1_8"
			read.seq = "AATTTTAATTTTAATTTTTGCGGTATTTTTAGTCGGTTCGTTCGTTCGGGTTTGATTTGAG"
			read.flag = 99
			read.tid = 0
			read.pos = 450
			read.mapq = 255
			read.cigar = [(0,61)]
			read.rnext = 1
			read.pnext = 512
			read.isize = 121
			read.qual = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
			read.tags = read.tags + [("XG", "CT")] + [("XM", "..hhh...hhh...hhh.z.Z....hhh.x..xZ..hxZ.hxZ.hxZ....x...hx....")] + [("XR", "CT")]
			return read

		def buildOTRead2():
			'''build an example read_2 aligned to OT-strand.
			'''

			read = pysam.AlignedRead()
			read.qname = "ADS-adipose_chr1_8"
			read.seq = "AGAATTGTGTTTCGTTTTTAGAGTATTATCGAAATTTGTGTAGAGGATAACGTAGCTTC"
			read.flag = 147
			read.tid = 0
			read.pos = 512
			read.mapq = 255
			read.cigar = [(0,59)]
			read.rnext = 1
			read.pnext = 450
			read.isize = -121
			read.qual = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
			read.tags = read.tags + [("XG", "CT")] + [("XM", "....x....h.xZ.hh..x......hh.xZ.....x....x......h..Z.x..H.xZ")] + [("XR", "GA")]
			return read

		def buildOBRead1():
			'''build an example read_1 aligned to OB-strand
			'''

			read = pysam.AlignedRead()
			read.qname = "ADS-adipose_chr1_22929891"
			read.seq = "AACGCAACTCCGCCCTCGCGATACTCTCCGAATCTATACTAAAAAAAACGCAACTCCGCCGAC"
			read.flag = 83
			read.tid = 0
			read.pos = 560
			read.mapq = 255
			read.cigar = [(0,63)]
			read.rnext = 1
			read.pnext = 492
			read.isize = -131
			read.qual = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
			read.tags = read.tags + [("XG", "GA")] + [("XM", "...Z..x....Z.....Z.Zx.h......Zxh...x.h..x.hh.h...Z.......Z..Zx.")] + [("XR", "CT")]
			return read

		def buildOBRead2():
			'''build an example read_2 aligned to OB-strand.
			'''

			read = pysam.AlignedRead()
			read.qname = "ADS-adipose_chr1_22929891"
			read.seq = "CACCCGAATCTAACCTAAAAAAAACTATACTCCGCCTTCAAAATACCACCGAAATCTATACAAAAAA"
			read.flag = 163
			read.tid = 0
			read.pos = 492
			read.mapq = 255
			read.cigar = [(0,67)]
			read.rnext = 1
			read.pnext = 560
			read.isize = 131
			read.qual = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
			read.tags = read.tags + [("XG", "GA")] + [("XM", ".z...Zxh...x....x.hh.h....x.h....Z......x.h.......Z......x.h..x.hh.")] + [("XR", "GA")]
			return read

		def buildBAM():
			'''build a BAM file
			'''

			header = { 'HD': {'VN': '1.0'}, 'SQ': [{'LN': 10000, 'SN': 'chr1'}, {'LN': 20000, 'SN': 'chr2'}] }
			tempfile_path = tempfile.mkstemp()[1]
			BAM = pysam.Samfile(tempfile_path, "wb", header = header)
			return BAM

		# Create the reads, BAM file and methylation m-tuples
		self.otr_1 = buildOTRead1()
		self.otr_2 = buildOTRead2()
		self.obr_1 = buildOBRead1()
		self.obr_2 = buildOBRead2()
		self.BAM = buildBAM()
		self.BAM.write(self.otr_1)
		self.BAM.write(self.otr_2)
		self.BAM.write(self.obr_1)
		self.BAM.write(self.obr_2)
		self.filename = self.BAM.filename
		self.BAM.close()
		self.BAM = pysam.Samfile(self.filename, 'rb')
		self.FAILED_QC = open(tempfile.mkstemp()[1], 'w')
		self.m1ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 1, 'CG', {'chr1': 0}), m = 1, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m2ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 2, 'CG', {'chr1': 0}), m = 2, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m3ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 3, 'CG', {'chr1': 0}), m = 3, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m4ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 4, 'CG', {'chr1': 0}), m = 4, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m5ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 5, 'CG', {'chr1': 0}), m = 5, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m6ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 6, 'CG', {'chr1': 0}), m = 6, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m7ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 7, 'CG', {'chr1': 0}), m = 7, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m8ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 8, 'CG', {'chr1': 0}), m = 8, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m9ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 9, 'CG', {'chr1': 0}), m = 9, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m10ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 10, 'CG', {'chr1': 0}), m = 10, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m1ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 1, 'CG', {'chr1': 0}), m = 1, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m2ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 2, 'CG', {'chr1': 0}), m = 2, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m3ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 3, 'CG', {'chr1': 0}), m = 3, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m4ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 4, 'CG', {'chr1': 0}), m = 4, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m5ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 5, 'CG', {'chr1': 0}), m = 5, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m6ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 6, 'CG', {'chr1': 0}), m = 6, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m7ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 7, 'CG', {'chr1': 0}), m = 7, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m8ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 8, 'CG', {'chr1': 0}), m = 8, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m9ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 9, 'CG', {'chr1': 0}), m = 9, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m10ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 10, 'CG', {'chr1': 0}), m = 10, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m11ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 11, 'CG', {'chr1': 0}), m = 11, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m12ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 12, 'CG', {'chr1': 0}), m = 12, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m2cgchg, self.nmlifotcgchg , self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = MTuple('test', 2, 'CG/CHG', {'chr1': 0}), m = 2, methylation_type = 'CG/CHG', methylation_pattern = re.compile(r'[Zz]'), ignore_read1_pos = [], ignore_read2_pos = [], min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)

	def test_correct_number_of_m_tuples(self):
		self.assertEqual(len(self.m1ot.mtuples), 10)
		self.assertEqual(len(self.m2ot.mtuples), 9)
		self.assertEqual(len(self.m3ot.mtuples), 8)
		self.assertEqual(len(self.m4ot.mtuples), 7)
		self.assertEqual(len(self.m5ot.mtuples), 6)
		self.assertEqual(len(self.m6ot.mtuples), 5)
		self.assertEqual(len(self.m7ot.mtuples), 4)
		self.assertEqual(len(self.m8ot.mtuples), 3)
		self.assertEqual(len(self.m9ot.mtuples), 2)
		self.assertEqual(len(self.m10ot.mtuples), 1)
		self.assertEqual(len(self.m1ob.mtuples), 12)
		self.assertEqual(len(self.m2ob.mtuples), 11)
		self.assertEqual(len(self.m3ob.mtuples), 10)
		self.assertEqual(len(self.m4ob.mtuples), 9)
		self.assertEqual(len(self.m5ob.mtuples), 8)
		self.assertEqual(len(self.m6ob.mtuples), 7)
		self.assertEqual(len(self.m7ob.mtuples), 6)
		self.assertEqual(len(self.m8ob.mtuples), 5)
		self.assertEqual(len(self.m9ob.mtuples), 4)
		self.assertEqual(len(self.m10ob.mtuples), 3)
		self.assertEqual(len(self.m11ob.mtuples), 2)
		self.assertEqual(len(self.m12ob.mtuples), 1)

	def test_correct_m_tuple_ids(self):
		self.assertItemsEqual(self.m1ot.mtuples.keys(), [('chr1', 563), ('chr1', 571), ('chr1', 525), ('chr1', 493), ('chr1', 469), ('chr1', 484), ('chr1', 489), ('chr1', 497), ('chr1', 471), ('chr1', 542)])
		self.assertItemsEqual(self.m2ot.mtuples.keys(), [('chr1', 497, 525), ('chr1', 563, 571), ('chr1', 484, 489), ('chr1', 469, 471), ('chr1', 542, 563), ('chr1', 493, 497), ('chr1', 471, 484), ('chr1', 489, 493), ('chr1', 525, 542)])
		self.assertItemsEqual(self.m3ot.mtuples.keys(), [('chr1', 493, 497, 525), ('chr1', 497, 525, 542), ('chr1', 489, 493, 497), ('chr1', 525, 542, 563), ('chr1', 471, 484, 489), ('chr1', 484, 489, 493), ('chr1', 542, 563, 571), ('chr1', 469, 471, 484)])
		self.assertItemsEqual(self.m4ot.mtuples.keys(), [('chr1', 484, 489, 493, 497), ('chr1', 525, 542, 563, 571), ('chr1', 471, 484, 489, 493), ('chr1', 493, 497, 525, 542), ('chr1', 469, 471, 484, 489), ('chr1', 497, 525, 542, 563), ('chr1', 489, 493, 497, 525)])
		self.assertItemsEqual(self.m5ot.mtuples.keys(), [('chr1', 469, 471, 484, 489, 493), ('chr1', 497, 525, 542, 563, 571), ('chr1', 493, 497, 525, 542, 563), ('chr1', 489, 493, 497, 525, 542), ('chr1', 484, 489, 493, 497, 525), ('chr1', 471, 484, 489, 493, 497)])
		self.assertItemsEqual(self.m6ot.mtuples.keys(), [('chr1', 469, 471, 484, 489, 493, 497), ('chr1', 471, 484, 489, 493, 497, 525), ('chr1', 489, 493, 497, 525, 542, 563), ('chr1', 484, 489, 493, 497, 525, 542), ('chr1', 493, 497, 525, 542, 563, 571)])
		self.assertItemsEqual(self.m7ot.mtuples.keys(), [('chr1', 469, 471, 484, 489, 493, 497, 525), ('chr1', 471, 484, 489, 493, 497, 525, 542), ('chr1', 484, 489, 493, 497, 525, 542, 563), ('chr1', 489, 493, 497, 525, 542, 563, 571)])
		self.assertItemsEqual(self.m8ot.mtuples.keys(), [('chr1', 469, 471, 484, 489, 493, 497, 525, 542), ('chr1', 484, 489, 493, 497, 525, 542, 563, 571), ('chr1', 471, 484, 489, 493, 497, 525, 542, 563)])
		self.assertItemsEqual(self.m9ot.mtuples.keys(), [('chr1', 471, 484, 489, 493, 497, 525, 542, 563, 571), ('chr1', 469, 471, 484, 489, 493, 497, 525, 542, 563)])
		self.assertItemsEqual(self.m10ot.mtuples.keys(), [('chr1', 469, 471, 484, 489, 493, 497, 525, 542, 563, 571)])
		self.assertItemsEqual(self.m1ob.mtuples.keys(), [('chr1', 563), ('chr1', 617), ('chr1', 493), ('chr1', 577), ('chr1', 525), ('chr1', 542), ('chr1', 579), ('chr1', 589), ('chr1', 571), ('chr1', 620), ('chr1', 497), ('chr1', 609)])
		self.assertItemsEqual(self.m2ob.mtuples.keys(), [('chr1', 589, 609), ('chr1', 609, 617), ('chr1', 571, 577), ('chr1', 563, 571), ('chr1', 577, 579), ('chr1', 542, 563), ('chr1', 493, 497), ('chr1', 617, 620), ('chr1', 579, 589), ('chr1', 525, 542), ('chr1', 497, 525)])
		self.assertItemsEqual(self.m3ob.mtuples.keys(), [('chr1', 571, 577, 579), ('chr1', 493, 497, 525), ('chr1', 609, 617, 620), ('chr1', 497, 525, 542), ('chr1', 589, 609, 617), ('chr1', 579, 589, 609), ('chr1', 525, 542, 563), ('chr1', 563, 571, 577), ('chr1', 542, 563, 571), ('chr1', 577, 579, 589)])
		self.assertItemsEqual(self.m4ob.mtuples.keys(), [('chr1', 577, 579, 589, 609), ('chr1', 571, 577, 579, 589), ('chr1', 563, 571, 577, 579), ('chr1', 525, 542, 563, 571), ('chr1', 542, 563, 571, 577), ('chr1', 579, 589, 609, 617), ('chr1', 493, 497, 525, 542), ('chr1', 497, 525, 542, 563), ('chr1', 589, 609, 617, 620)])
		self.assertItemsEqual(self.m5ob.mtuples.keys(), [('chr1', 563, 571, 577, 579, 589), ('chr1', 571, 577, 579, 589, 609), ('chr1', 577, 579, 589, 609, 617), ('chr1', 497, 525, 542, 563, 571), ('chr1', 493, 497, 525, 542, 563), ('chr1', 579, 589, 609, 617, 620), ('chr1', 542, 563, 571, 577, 579), ('chr1', 525, 542, 563, 571, 577)])
		self.assertItemsEqual(self.m6ob.mtuples.keys(), [('chr1', 577, 579, 589, 609, 617, 620), ('chr1', 497, 525, 542, 563, 571, 577), ('chr1', 571, 577, 579, 589, 609, 617), ('chr1', 493, 497, 525, 542, 563, 571), ('chr1', 563, 571, 577, 579, 589, 609), ('chr1', 525, 542, 563, 571, 577, 579), ('chr1', 542, 563, 571, 577, 579, 589)])
		self.assertItemsEqual(self.m7ob.mtuples.keys(), [('chr1', 525, 542, 563, 571, 577, 579, 589), ('chr1', 493, 497, 525, 542, 563, 571, 577), ('chr1', 497, 525, 542, 563, 571, 577, 579), ('chr1', 542, 563, 571, 577, 579, 589, 609), ('chr1', 571, 577, 579, 589, 609, 617, 620), ('chr1', 563, 571, 577, 579, 589, 609, 617)])
		self.assertItemsEqual(self.m8ob.mtuples.keys(), [('chr1', 493, 497, 525, 542, 563, 571, 577, 579), ('chr1', 497, 525, 542, 563, 571, 577, 579, 589), ('chr1', 563, 571, 577, 579, 589, 609, 617, 620), ('chr1', 525, 542, 563, 571, 577, 579, 589, 609), ('chr1', 542, 563, 571, 577, 579, 589, 609, 617)])
		self.assertItemsEqual(self.m9ob.mtuples.keys(), [('chr1', 497, 525, 542, 563, 571, 577, 579, 589, 609), ('chr1', 493, 497, 525, 542, 563, 571, 577, 579, 589), ('chr1', 542, 563, 571, 577, 579, 589, 609, 617, 620), ('chr1', 525, 542, 563, 571, 577, 579, 589, 609, 617)])
		self.assertItemsEqual(self.m10ob.mtuples.keys(), [('chr1', 497, 525, 542, 563, 571, 577, 579, 589, 609, 617), ('chr1', 493, 497, 525, 542, 563, 571, 577, 579, 589, 609), ('chr1', 525, 542, 563, 571, 577, 579, 589, 609, 617, 620)])
		self.assertItemsEqual(self.m11ob.mtuples.keys(), [('chr1', 497, 525, 542, 563, 571, 577, 579, 589, 609, 617, 620), ('chr1', 493, 497, 525, 542, 563, 571, 577, 579, 589, 609, 617)])
		self.assertItemsEqual(self.m12ob.mtuples.keys(), [('chr1', 493, 497, 525, 542, 563, 571, 577, 579, 589, 609, 617, 620)])

	def test_correct_number_of_methylation_loci_in_fragment(self):
		self.assertEqual(self.nmlifot, 10)
		self.assertEqual(self.nmlifob, 12)

	def test_correct_methylation_type(self):
		self.assertEqual(self.m1ot.methylation_type, 'CG')
		self.assertEqual(self.m1ob.methylation_type, 'CG')
		self.assertEqual(self.m2cgchg.methylation_type, 'CG/CHG')

	def test_counts(self):
		self.assertEqual(self.m1ot.mtuples[('chr1', 563)], array.array('i', [1, 0]))

	def tearDown(self):
		os.remove(self.BAM.filename)
		os.remove(self.FAILED_QC.name)

class TestMTuple(unittest.TestCase):
	'''Test the class MTuple and its methods
	'''

	def setUp(self):

		def buildOTRead():
			'''build an example read aligned to OT-strand.
			'''
			read = pysam.AlignedRead()
			read.qname = "@SALK_2077_FC6295TAAXX:2:107:9396:15019#0/1"
			read.seq = "GGGGAAGGTGTTATGGAGTTTTTTACGATTTTTAGTCGTTTTCGTTTTTTTTTGTTTGTGGTTGTTGCGGTGGCGGTAGAGGAGGG"
			read.flag = 0
			read.tid = 0
			read.pos = 4536
			read.mapq = 255
			read.cigar = [(0,86)]
			read.rnext = 0
			read.pnext = 0
			read.isize = 0
			read.qual = "DBDB2;@>)@@F?EFG@GBGGGGDDBG@DGGGGEEFHHEGHHHHEFHHHHFHHHFHHHGHGBCEAA@?@?/A@>@3,.6,AA,@>="
			read.tags = read.tags + [("XG", "CT")] + [("XM", "...........h......hhhhh..Z....hhx...Z..hh.Z..hh.hh.x..hx.....x..x..Z.....Z..x.........")] + [("XR", "CT")]
			return read

		def buildOBRead():
			'''build an example read aligned to OB-strand.
			'''
			read = pysam.AlignedRead()
			read.qname = "@ECKER_1116_FC623CNAAXX:2:21:18515:1127#0/1"
			read.seq = "CTTCCTAACAAACAACTACACCACTACCTAACGCTATACCCTTCCTTTACTCTACCCACTAAAAACAATATTTATCATAAACCT"
			read.flag = 16
			read.tid = 0
			read.pos = 3334
			read.mapq = 255
			read.cigar = [(0,84)]
			read.rnext = 0
			read.pnext = 0
			read.isize = 0
			read.qual = "G7@G@BGB@GGGGGDIEEBIBA<AHEGEEEGGGDDEDFFEIIHIIGGDGGGGGGGGGGDGDBED<FAAFEGGGGGIHIFIGBDG"
			read.tags = read.tags + [("XG", "GA")] + [("XM", "......x...xh..x..x.......x...xh.Z..x.h..........h....x...z..xh.h..zx.h...h....hhh...")] + [("XR", "CT")]
			return read

		def buildOTRead1():
			'''build an example read_1 aligned to OT-strand.
			'''
			read = pysam.AlignedRead()
			read.qname = "ADS-adipose_chr1_8"
			read.seq = "AATTTTAATTTTAATTTTTGCGGTATTTTTAGTCGGTTCGTTCGTTCGGGTTTGATTTGAG"
			read.flag = 99
			read.tid = 0
			read.pos = 450
			read.mapq = 255
			read.cigar = [(0,61)]
			read.rnext = 1
			read.pnext = 512
			read.isize = 121
			read.qual = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
			read.tags = read.tags + [("XG", "CT")] + [("XM", "..hhh...hhh...hhh.z.Z....hhh.x..xZ..hxZ.hxZ.hxZ....x...hx....")] + [("XR", "CT")]
			return read

		def buildOTRead2():
			'''build an example read_2 aligned to OT-strand.
			'''
			read = pysam.AlignedRead()
			read.qname = "ADS-adipose_chr1_8"
			read.seq = "AGAATTGTGTTTCGTTTTTAGAGTATTATCGAAATTTGTGTAGAGGATAACGTAGCTTC"
			read.flag = 147
			read.tid = 0
			read.pos = 512
			read.mapq = 255
			read.cigar = [(0,59)]
			read.rnext = 1
			read.pnext = 450
			read.isize = -121
			read.qual = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
			read.tags = read.tags + [("XG", "CT")] + [("XM", "....x....h.xZ.hh..x......hh.xZ.....x....x......h..Z.x..H.xZ")] + [("XR", "GA")]
			return read

		def buildOBRead1():
			'''build an example read_1 aligned to OB-strand
			'''
			read = pysam.AlignedRead()
			read.qname = "ADS-adipose_chr1_22929891"
			read.seq = "AACGCAACTCCGCCCTCGCGATACTCTCCGAATCTATACTAAAAAAAACGCAACTCCGCCGAC"
			read.flag = 83
			read.tid = 0
			read.pos = 560
			read.mapq = 255
			read.cigar = [(0,63)]
			read.rnext = 1
			read.pnext = 492
			read.isize = -131
			read.qual = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
			read.tags = read.tags + [("XG", "GA")] + [("XM", "...Z..x....Z.....Z.Zx.h......Zxh...x.h..x.hh.h...Z.......Z..Zx.")] + [("XR", "CT")]
			return read

		def buildOBRead2():
			'''build an example read_2 aligned to OB-strand.
			'''
			read = pysam.AlignedRead()
			read.qname = "ADS-adipose_chr1_22929891"
			read.seq = "CACCCGAATCTAACCTAAAAAAAACTATACTCCGCCTTCAAAATACCACCGAAATCTATACAAAAAA"
			read.flag = 163
			read.tid = 0
			read.pos = 492
			read.mapq = 255
			read.cigar = [(0,67)]
			read.rnext = 1
			read.pnext = 560
			read.isize = 131
			read.qual = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
			read.tags = read.tags + [("XG", "GA")] + [("XM", ".z...Zxh...x....x.hh.h....x.h....Z......x.h.......Z......x.h..x.hh.")] + [("XR", "GA")]
			return read

		def buildBAM():
			'''build a BAM file
			'''
			header = { 'HD': {'VN': '1.0'}, 'SQ': [{'LN': 10000, 'SN': 'chr1'}, {'LN': 20000, 'SN': 'chr2'}] }
			tempfile_path = tempfile.mkstemp()[1]
			BAM = pysam.Samfile(tempfile_path, "wb", header = header)
			return BAM

		# Create the reads, BAM file and methylation m-tuples
		# Single-end
		self.otr = buildOTRead()
		self.obr = buildOBRead()
		self.BAM = buildBAM()
		self.BAM.write(self.otr)
		self.BAM.write(self.obr)
		self.filename = self.BAM.filename
		self.BAM.close()
		self.BAM = pysam.Samfile(self.filename, 'rb')
		self.wfotr = MTuple('test', 2, 'CG', {'chr1': 0})
		self.wfotr.increment_count(('chr1', 4562, 4573), 'ZZ', self.otr, None)
		self.wfobr = MTuple('test', 2, 'CG', {'chr1': 0})
		self.wfobr.increment_count(('chr1', 3366, 3391), 'ZZ', self.otr, None)
		self.wfotrcgchg = MTuple('test', 3, 'CG/CHG', {'chr1': 0})
		self.wfotrcgchg.increment_count(('chr1', 4562, 4569, 4573), 'ZXZ', self.otr, None)

		# Paired-end
		self.otr_1 = buildOTRead1()
		self.otr_2 = buildOTRead2()
		self.obr_1 = buildOBRead1()
		self.obr_2 = buildOBRead2()
		self.BAMPE = buildBAM()
		self.BAMPE.write(self.otr_1)
		self.BAMPE.write(self.otr_2)
		self.BAMPE.write(self.obr_1)
		self.BAMPE.write(self.obr_2)
		self.filename = self.BAMPE.filename
		self.BAMPE.close()
		self.BAMPE = pysam.Samfile(self.filename, 'rb')
		self.wfotrpe = MTuple('test', 2, 'CG', {'chr1': 0})
		self.wfotrpe.increment_count(('chr1', 497, 525), 'ZZ', self.otr_1, self.otr_2)
		self.wfobrpe = MTuple('test', 2, 'CG', {'chr1': 0})
		self.wfobrpe.increment_count(('chr1', 563, 571), 'ZZ', self.obr_1, self.obr_2)
		self.wfotrpecgchg = MTuple('test', 3, 'CG/CHG', {'chr1': 0})
		self.wfotrpecgchg.increment_count(('chr1', 497, 525, 534), 'ZZX', self.otr_1, self.otr_2)

	def test_init_se(self):
		self.assertEqual(self.wfotr.chr_map, {'chr1': 0})
		self.assertEqual(self.wfotr.comethylation_patterns, ['MM', 'MU', 'UM', 'UU'])
		self.assertEqual(self.wfotr.m, 2)
		self.assertEqual(self.wfotr.mtuples.keys(), [('chr1', 4562, 4573)])
		self.assertEqual(self.wfotr.mtuples.values(), [array.array('i', [1, 0, 0, 0])])
		self.assertEqual(self.wfotr.methylation_type, 'CG')
		self.assertEqual(self.wfobr.methylation_type, 'CG')
		self.assertEqual(self.wfotrcgchg.methylation_type, 'CG/CHG')

	def test_init_pe(self):
		self.assertEqual(self.wfotrpe.methylation_type, 'CG')
		self.assertEqual(self.wfobrpe.methylation_type, 'CG')
		self.assertEqual(self.wfotrpecgchg.methylation_type, 'CG/CHG')

	def test_increment_count_ot_se(self):
		self.assertEqual(self.wfotr.mtuples[('chr1', 4562, 4573)], array.array('i', [1, 0, 0, 0]))
		self.wfotr.increment_count(('chr1', 4562, 4573), 'MM', self.otr, None)
		self.assertEqual(self.wfotr.mtuples[('chr1', 4562, 4573)], array.array('i', [2, 0, 0, 0]))
		self.wfotr.increment_count(('chr1', 4562, 4573), 'MU', self.otr, None)
		self.assertEqual(self.wfotr.mtuples[('chr1', 4562, 4573)], array.array('i', [2, 1, 0, 0]))
		self.wfotr.increment_count(('chr1', 4562, 4573), 'UM', self.otr, None)
		self.assertEqual(self.wfotr.mtuples[('chr1', 4562, 4573)], array.array('i', [2, 1, 1, 0]))
		self.wfotr.increment_count(('chr1', 4562, 4573), 'UU', self.otr, None)
		self.assertEqual(self.wfotr.mtuples[('chr1', 4562, 4573)], array.array('i', [2, 1, 1, 1]))
		self.wfotr.increment_count(('chr1', 4562, 4573), 'MM', self.otr, None)
		self.assertEqual(self.wfotr.mtuples[('chr1', 4562, 4573)], array.array('i', [3, 1, 1, 1]))

	def test_increment_count_ob_se(self):
		self.assertEqual(self.wfobr.mtuples[('chr1', 3366, 3391)], array.array('i', [1, 0, 0, 0]))
		self.wfobr.increment_count(('chr1', 3366, 3391), 'MM', self.obr, None)
		self.assertEqual(self.wfobr.mtuples[('chr1', 3366, 3391)], array.array('i', [2, 0, 0, 0]))
		self.wfobr.increment_count(('chr1', 3366, 3391), 'MU', self.obr, None)
		self.assertEqual(self.wfobr.mtuples[('chr1', 3366, 3391)], array.array('i', [2, 1, 0, 0]))
		self.wfobr.increment_count(('chr1', 3366, 3391), 'UM', self.obr, None)
		self.assertEqual(self.wfobr.mtuples[('chr1', 3366, 3391)], array.array('i', [2, 1, 1, 0]))
		self.wfobr.increment_count(('chr1', 3366, 3391), 'UU', self.obr, None)
		self.assertEqual(self.wfobr.mtuples[('chr1', 3366, 3391)], array.array('i', [2, 1, 1, 1]))
		self.wfobr.increment_count(('chr1', 3366, 3391), 'MM', self.obr, None)
		self.assertEqual(self.wfobr.mtuples[('chr1', 3366, 3391)], array.array('i', [3, 1, 1, 1]))

	def test_increment_count_ctot_se(self):
		self.otr.tags = []
		self.otr.tags = self.otr.tags + [("XG", "CT")] + [("XM", "...........h......hhhhh..Z....hhx...Z..hh.Z..hh.hh.x..hx.....x..x..Z.....Z..x.........")] + [("XR", "GA")]
		self.assertEqual(self.wfotr.mtuples[('chr1', 4562, 4573)], array.array('i', [1, 0, 0, 0]))
		self.wfotr.increment_count(('chr1', 4562, 4573), 'MM', self.otr, None)
		self.assertEqual(self.wfotr.mtuples[('chr1', 4562, 4573)], array.array('i', [2, 0, 0, 0]))
		self.wfotr.increment_count(('chr1', 4562, 4573), 'MU', self.otr, None)
		self.assertEqual(self.wfotr.mtuples[('chr1', 4562, 4573)], array.array('i', [2, 1, 0, 0]))
		self.wfotr.increment_count(('chr1', 4562, 4573), 'UM', self.otr, None)
		self.assertEqual(self.wfotr.mtuples[('chr1', 4562, 4573)], array.array('i', [2, 1, 1, 0]))
		self.wfotr.increment_count(('chr1', 4562, 4573), 'UU', self.otr, None)
		self.assertEqual(self.wfotr.mtuples[('chr1', 4562, 4573)], array.array('i', [2, 1, 1, 1]))
		self.wfotr.increment_count(('chr1', 4562, 4573), 'MM', self.otr, None)
		self.assertEqual(self.wfotr.mtuples[('chr1', 4562, 4573)], array.array('i', [3, 1, 1, 1]))

	def test_increment_count_ctob_se(self):
		self.obr.tags = []
		self.obr.tags = self.obr.tags + [("XG", "GA")] + [("XM", "......x...xh..x..x.......x...xh.Z..x.h..........h....x...z..xh.h..zx.h...h....hhh...")] + [("XR", "GA")]
		self.assertEqual(self.wfobr.mtuples[('chr1', 3366, 3391)], array.array('i', [1, 0, 0, 0]))
		self.wfobr.increment_count(('chr1', 3366, 3391), 'MM', self.obr, None)
		self.assertEqual(self.wfobr.mtuples[('chr1', 3366, 3391)], array.array('i', [2, 0, 0, 0]))
		self.wfobr.increment_count(('chr1', 3366, 3391), 'MU', self.obr, None)
		self.assertEqual(self.wfobr.mtuples[('chr1', 3366, 3391)], array.array('i', [2, 1, 0, 0]))
		self.wfobr.increment_count(('chr1', 3366, 3391), 'UM', self.obr, None)
		self.assertEqual(self.wfobr.mtuples[('chr1', 3366, 3391)], array.array('i', [2, 1, 1, 0]))
		self.wfobr.increment_count(('chr1', 3366, 3391), 'UU', self.obr, None)
		self.assertEqual(self.wfobr.mtuples[('chr1', 3366, 3391)], array.array('i', [2, 1, 1, 1]))
		self.wfobr.increment_count(('chr1', 3366, 3391), 'MM', self.obr, None)
		self.assertEqual(self.wfobr.mtuples[('chr1', 3366, 3391)], array.array('i', [3, 1, 1, 1]))

	def test_increment_count_multiple_methylation_types_se(self):
		self.assertEqual(self.wfotrcgchg.mtuples[('chr1', 4562, 4569, 4573)], array.array('i', [1, 0, 0, 0, 0, 0, 0, 0]))
		self.wfotrcgchg.increment_count(('chr1', 4562, 4569, 4573), 'MMM', self.otr, None)
		self.assertEqual(self.wfotrcgchg.mtuples[('chr1', 4562, 4569, 4573)], array.array('i', [2, 0, 0, 0, 0, 0, 0, 0]))
		self.wfotrcgchg.increment_count(('chr1', 4562, 4569, 4573), 'MMU', self.otr, None)
		self.assertEqual(self.wfotrcgchg.mtuples[('chr1', 4562, 4569, 4573)], array.array('i', [2, 1, 0, 0, 0, 0, 0, 0]))
		self.wfotrcgchg.increment_count(('chr1', 4562, 4569, 4573), 'MUM', self.otr, None)
		self.assertEqual(self.wfotrcgchg.mtuples[('chr1', 4562, 4569, 4573)], array.array('i', [2, 1, 1, 0, 0, 0, 0, 0]))
		self.wfotrcgchg.increment_count(('chr1', 4562, 4569, 4573), 'MUU', self.otr, None)
		self.assertEqual(self.wfotrcgchg.mtuples[('chr1', 4562, 4569, 4573)], array.array('i', [2, 1, 1, 1, 0, 0, 0, 0]))
		self.wfotrcgchg.increment_count(('chr1', 4562, 4569, 4573), 'UMM', self.otr, None)
		self.assertEqual(self.wfotrcgchg.mtuples[('chr1', 4562, 4569, 4573)], array.array('i', [2, 1, 1, 1, 1, 0, 0, 0]))
		self.wfotrcgchg.increment_count(('chr1', 4562, 4569, 4573), 'UMU', self.otr, None)
		self.assertEqual(self.wfotrcgchg.mtuples[('chr1', 4562, 4569, 4573)], array.array('i', [2, 1, 1, 1, 1, 1, 0, 0]))
		self.wfotrcgchg.increment_count(('chr1', 4562, 4569, 4573), 'UUM', self.otr, None)
		self.assertEqual(self.wfotrcgchg.mtuples[('chr1', 4562, 4569, 4573)], array.array('i', [2, 1, 1, 1, 1, 1, 1, 0]))
		self.wfotrcgchg.increment_count(('chr1', 4562, 4569, 4573), 'UUU', self.otr, None)
		self.assertEqual(self.wfotrcgchg.mtuples[('chr1', 4562, 4569, 4573)], array.array('i', [2, 1, 1, 1, 1, 1, 1, 1]))

	def test_increment_count_ot_pe(self):
		self.assertEqual(self.wfotrpe.mtuples[('chr1', 497, 525)], array.array('i', [1, 0, 0, 0]))
		self.wfotrpe.increment_count(('chr1', 497, 525), 'MM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.mtuples[('chr1', 497, 525)], array.array('i', [2, 0, 0, 0]))
		self.wfotrpe.increment_count(('chr1', 497, 525), 'MU', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.mtuples[('chr1', 497, 525)], array.array('i', [2, 1, 0, 0]))
		self.wfotrpe.increment_count(('chr1', 497, 525), 'UM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.mtuples[('chr1', 497, 525)], array.array('i', [2, 1, 1, 0]))
		self.wfotrpe.increment_count(('chr1', 497, 525), 'UU', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.mtuples[('chr1', 497, 525)], array.array('i', [2, 1, 1, 1]))

	def test_increment_count_ob_pe(self):
		self.assertEqual(self.wfobrpe.mtuples[('chr1', 563, 571)], array.array('i', [1, 0, 0, 0]))
		self.wfobrpe.increment_count(('chr1', 563, 571), 'MM', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.mtuples[('chr1', 563, 571)], array.array('i', [2, 0, 0, 0]))
		self.wfobrpe.increment_count(('chr1', 563, 571), 'MU', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.mtuples[('chr1', 563, 571)], array.array('i', [2, 1, 0, 0]))
		self.wfobrpe.increment_count(('chr1', 563, 571), 'UM', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.mtuples[('chr1', 563, 571)], array.array('i', [2, 1, 1, 0]))
		self.wfobrpe.increment_count(('chr1', 563, 571), 'UU', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.mtuples[('chr1', 563, 571)], array.array('i', [2, 1, 1, 1]))

	def test_increment_count_ctot_pe(self):
		self.otr_1.tags = []
		self.otr_1.tags = self.otr_1.tags + [('XG', 'CT'), ('XM', '..hhh...hhh...hhh.z.Z....hhh.x..xZ..hxZ.hxZ.hxZ....x...hx....'), ('XR', 'GA')]
		self.otr_2.tags = []
		self.otr_2.tags = self.otr_2.tags + [('XG', 'CT'), ('XM', '....x....h.xZ.hh..x......hh.xZ.....x....x......h..Z.x..H.xZ'), ('XR', 'CT')]
		self.assertEqual(self.wfotrpe.mtuples[('chr1', 497, 525)], array.array('i', [1, 0, 0, 0]))
		self.wfotrpe.increment_count(('chr1', 497, 525), 'MM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.mtuples[('chr1', 497, 525)], array.array('i', [2, 0, 0, 0]))
		self.wfotrpe.increment_count(('chr1', 497, 525), 'MU', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.mtuples[('chr1', 497, 525)], array.array('i', [2, 1, 0, 0]))
		self.wfotrpe.increment_count(('chr1', 497, 525), 'UM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.mtuples[('chr1', 497, 525)], array.array('i', [2, 1, 1, 0]))
		self.wfotrpe.increment_count(('chr1', 497, 525), 'UU', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.mtuples[('chr1', 497, 525)], array.array('i', [2, 1, 1, 1]))

	def test_increment_count_ctob_pe(self):
		self.obr_1.tags = []
		self.obr_1.tags = self.obr_1.tags + [('XG', 'GA'), ('XM', '...Z..x....Z.....Z.Zx.h......Zxh...x.h..x.hh.h...Z.......Z..Zx.'), ('XR', 'GA')]
		self.obr_2.tags = []
		self.obr_2.tags = self.obr_2.tags + [('XG', 'GA'), ('XM', '.z...Zxh...x....x.hh.h....x.h....Z......x.h.......Z......x.h..x.hh.'), ('XR', 'CT')]
		self.assertEqual(self.wfobrpe.mtuples[('chr1', 563, 571)], array.array('i', [1, 0, 0, 0]))
		self.wfobrpe.increment_count(('chr1', 563, 571), 'MM', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.mtuples[('chr1', 563, 571)], array.array('i', [2, 0, 0, 0]))
		self.wfobrpe.increment_count(('chr1', 563, 571), 'MU', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.mtuples[('chr1', 563, 571)], array.array('i', [2, 1, 0, 0]))
		self.wfobrpe.increment_count(('chr1', 563, 571), 'UM', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.mtuples[('chr1', 563, 571)], array.array('i', [2, 1, 1, 0]))
		self.wfobrpe.increment_count(('chr1', 563, 571), 'UU', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.mtuples[('chr1', 563, 571)], array.array('i', [2, 1, 1, 1]))	

	def test_increment_count_multiple_methylation_types_pe(self):
		self.assertEqual(self.wfotrpecgchg.mtuples[('chr1', 497, 525, 534)], array.array('i', [1, 0, 0, 0, 0, 0, 0, 0]))
		self.wfotrpecgchg.increment_count(('chr1', 497, 525, 534), 'MMM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpecgchg.mtuples[('chr1', 497, 525, 534)], array.array('i', [2, 0, 0, 0, 0, 0, 0, 0]))
		self.wfotrpecgchg.increment_count(('chr1', 497, 525, 534), 'MMU', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpecgchg.mtuples[('chr1', 497, 525, 534)], array.array('i', [2, 1, 0, 0, 0, 0, 0, 0]))
		self.wfotrpecgchg.increment_count(('chr1', 497, 525, 534), 'MUM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpecgchg.mtuples[('chr1', 497, 525, 534)], array.array('i', [2, 1, 1, 0, 0, 0, 0, 0]))
		self.wfotrpecgchg.increment_count(('chr1', 497, 525, 534), 'MUU', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpecgchg.mtuples[('chr1', 497, 525, 534)], array.array('i', [2, 1, 1, 1, 0, 0, 0, 0]))
		self.wfotrpecgchg.increment_count(('chr1', 497, 525, 534), 'UMM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpecgchg.mtuples[('chr1', 497, 525, 534)], array.array('i', [2, 1, 1, 1, 1, 0, 0, 0]))
		self.wfotrpecgchg.increment_count(('chr1', 497, 525, 534), 'UMU', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpecgchg.mtuples[('chr1', 497, 525, 534)], array.array('i', [2, 1, 1, 1, 1, 1, 0, 0]))
		self.wfotrpecgchg.increment_count(('chr1', 497, 525, 534), 'UUM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpecgchg.mtuples[('chr1', 497, 525, 534)], array.array('i', [2, 1, 1, 1, 1, 1, 1, 0]))
		self.wfotrpecgchg.increment_count(('chr1', 497, 525, 534), 'UUU', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpecgchg.mtuples[('chr1', 497, 525, 534)], array.array('i', [2, 1, 1, 1, 1, 1, 1, 1]))

	def test_invalid_comethylation_pattern(self):
		with self.assertRaises(SystemExit) as cm:
			self.wfotr.increment_count(('chr1', 4562, 4573), 'MMM', self.otr, None)
			self.assertEqual(cm.exception.code, 1)
		with self.assertRaises(SystemExit) as cm:
			self.wfotr.increment_count(('chr1', 4562, 4573), 'M', self.otr, None)
			self.assertEqual(cm.exception.code, 1)
		with self.assertRaises(SystemExit) as cm:
			self.wfotr.increment_count(('chr1', 4562, 4573), 'MA', self.otr, None)
			self.assertEqual(cm.exception.code, 1)
		with self.assertRaises(SystemExit) as cm:
			self.wfotr.increment_count(('chr1', 4562, 4573), 'mm', self.otr, None)
			self.assertEqual(cm.exception.code, 1)

	def test_invalid_m(self):
		self.assertRaises(ValueError, MTuple, 'test', -2, 'CG', {'chr1': 0})
		self.assertRaises(ValueError, MTuple, 'test', 2.3, 'CG', {'chr1': 0})
		self.assertRaises(ValueError, MTuple, 'test', 2.0, 'CG', {'chr1': 0})

	def test_invalid_methylation_type(self):
		self.assertRaises(ValueError, MTuple, 'test', -2, 'CT', {'chr1': 0})
		self.assertRaises(ValueError, MTuple, 'test', -2, 'CG-CHG', {'chr1': 0})
		self.assertRaises(ValueError, MTuple, 'test', -2, 'CHG/CHG', {'chr1': 0})
		self.assertRaises(ValueError, MTuple, 'test', -2, 'CHG/CG', {'chr1': 0})

	def tearDown(self):
		os.remove(self.BAM.filename)
		os.remove(self.BAMPE.filename)

class TestGetStrand(unittest.TestCase):

	def setUp(self):
		def buildOTSE():
			read = pysam.AlignedRead()
			read.qname = "SRR020138.15030048_SALK_2029:7:100:1740:1801_length=86"
			read.seq = "CGAATGTTTTTTATTATGAATGAGAGTTTGTTAAATTAGTTGGTTTTAGG"
			read.flag = 0
			read.tid = 0
			read.pos = 245746845
			read.mapq = 255
			read.cigar = [(0, 50)]
			read.rnext = 0
			read.pnext = 0
			read.isize = 0
			read.qual = "BC?BBBBBCCCCAA@BCB?AB@AB>CABB@@BB@?BB497@@B:@@B5>@"
			read.tags = read.tags + [("XM", "Z........h....h.............z.......x..x.....h....")] + [('XR', 'CT')] + [('XG', 'CT')]
			return read

		def buildOBSE():
			read = pysam.AlignedRead()
			read.qname = "SRR020138.15033460_SALK_2029:7:100:1783:2004_length=86"
			read.seq = "GTATACGCTAATTTTATAACCTAAAAATTTACTAAATTCATTAATCAAAT"
			read.flag = 16
			read.tid = 0
			read.pos = 96804459
			read.mapq = 255
			read.cigar = [(0, 50)]
			read.rnext = 0
			read.pnext = 0
			read.isize = 0
			read.qual = "CCBCACCBBCCCAACBBCAB@CCCCCCBBCCCCCCCCCBBBBCBCACBCB"
			read.tags = read.tags + [('XM', "Z.h...Z..x.....h.......hh.h......x........h....x..")] + [('XR', 'CT')] + [('XG', 'GA')]
			return read

		def buildCTOTSE():
			read = pysam.AlignedRead()
			read.qname = "SRR020138.15026483_SALK_2029:7:100:1698:1069_length=86"
			read.seq = "CCTCCATCATCATTCCTAATTTCTCCTTCCTCCCTTTCTACTTCCTCCTT"
			read.flag = 0
			read.tid = 0
			read.pos = 69393512
			read.mapq = 255
			read.cigar = [(0, 50)]
			read.rnext = 0
			read.pnext = 0
			read.isize = 0
			read.qual = "1)9<)@96'3%6@5:0=3::;;:89*;:@AA@=;A=3)2)1@>*9;-4:A"
			read.tags = read.tags + [('XM', 'H.h...h...H...HHh.................................')] + [('XR', 'GA')] + [('XG', 'CT')]
			return read

		def buildCTOBSE():
			read = pysam.AlignedRead()
			read.qname = "SRR020138.15034119_SALK_2029:7:100:1792:1889_length=86"
			read.seq = "ATGGAATGGAAAGGAAGGGAATTTAATGGAATGGAATGGAATGGAATGGA"
			read.flag = 16
			read.tid = 0
			read.pos = 10784296
			read.mapq = 255
			read.cigar = [(0, 50)]
			read.rnext = 0
			read.pnext = 0
			read.isize = 0
			read.qual = "BCBBBBCCBCAA=@?B?A?@@BCCBCCAAA=B9@B?A=9?BBB??AC@2="
			read.tags = read.tags + [('XM', '..HH...HHhh.HH...HH........HH.h.H..h.H..h.HH.h.HH.')] + [('XR', 'GA')] + [('XG', 'GA')]
			return read

		def buildOTPE():
			read_1 = pysam.AlignedRead()
			read_1.qname = "SRR400564.1684335_HAL:1133:C010EABXX:8:1108:18844:132483_length=101"
			read_1.seq = "CGAGTTCGTTTAAAGAGTAATTAGTTATTTTTGTAAGGTTTGGTTAGGGTTATAGAAGGTTTTTTTGGATGGTAATTTTGGTTGCGTTTTGTATTTGAATA"
			read_1.flag = 99
			read_1.tid = 0
			read_1.pos = 19978967
			read_1.mapq = 255
			read_1.cigar = [(0, 101)]
			read_1.rnext = 0
			read_1.pnext = 19979235
			read_1.isize = 369
			read_1.qual = "BC@FDFFFHHHHHJIHICHHHIIJJIGGJJJJJHGIF11?CGHIJIIJJHHHIIJIHHG9BFHIJFFDAC@?6;;-;AC@?DD<98@BDBDD>CDD#####"
			read_1.tags = read_1.tags + [('XM', 'Z...h.Z.h.h......h..hx..hh.hh.x..h.....x...hx....h..x......hhhhhx...........hx...x..Z...x..h.hx....h.')] + [('XR', 'CT')] + [('XG', 'CT')]
			read_2 = pysam.AlignedRead()
			read_2.qname = "SRR400564.1684335_HAL:1133:C010EABXX:8:1108:18844:132483_length=101"
			read_2.seq = "TATTGAGTTGTGTGTATGTCAGGTAAAGTTAGGTTGGTTTAAAGAGTAATTAGTTATTTTTGTAAGGTTGGGTTAGGAGAAGGCGGATTAGTTATTAATTT"
			read_2.flag = 147
			read_2.tid = 0
			read_2.pos = 19979235
			read_2.mapq = 255
			read_2.cigar = [(0, 101)]
			read_2.rnext = 0
			read_2.pnext = 19978967
			read_2.isize = -369
			read_2.qual = "EDDDDDDBBDDDDDDDCDEEEEEFFFFFHHGHIIFJJJIJJJJIJJJJJIHIJJIJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJHHHHHFFFFFCCC"
			read_2.tags = read_2.tags + [('XM', 'h.x....x......h.z.hX...h....hx...x...h.h......h..hx..hh.hh.x..h....x....hx.........Z...hx..hh.hh..hh.')] + [('XR', 'GA')] + [('XG', 'CT')]
			return read_1, read_2

		def buildOBPE():
			read_1 = pysam.AlignedRead()
			read_1.qname = "SRR400564.241291_HAL:1133:C010EABXX:8:1102:8553:52618_length=101"
			read_1.seq = "GATCACCTAAATCGAAAATTAAAAACCAACCTAACCAACACGATAAAACCCCATCTCTACTAAAATACAAAAACTAACCAAACGTAATAACAAACACCTAT"
			read_1.flag = 83
			read_1.tid = 0
			read_1.pos = 19195980
			read_1.mapq = 255
			read_1.cigar = [(0, 101)]
			read_1.rnext = 0
			read_1.pnext = 19195917
			read_1.isize = -164
			read_1.qual = "#####@:4B@;(85DDD@:;DDBBAACD?FD@HHHFIIJIIIHFJIFDJJIIFBGFGIGIGJIHFFDEJJJGIIGJIJJJJIJJHJIHHHHHHFDDDF@@@"
			read_1.tags = read_1.tags + [('XM', 'Z.......xhh..Zxh.h..h..h....x...xh.......Zx.h...............................h...xh.Z.hh.hh..x......x.')] + [('XR', 'CT')] + [('XG', 'GA')]
			read_2 = pysam.AlignedRead()
			read_2.qname = "SRR400564.241291_HAL:1133:C010EABXX:8:1102:8553:52618_length=101"
			read_2.seq = "TATAAAACAAAACATAATAACTCATACCTATAATCCCAACACTTTAAAAATCTAAAACAAACCGATCACCTAAATCGAAAATTAAAAACCAACCTAACCAA"
			read_2.flag = 163
			read_2.tid = 0
			read_2.pos = 19195917
			read_2.mapq = 255
			read_2.cigar = [(0, 101)]
			read_2.rnext = 0
			read_2.pnext = 19195980
			read_2.isize = 164
			read_2.qual = "CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJEHHIJJJJJJIJJIIJHHHFFFFFEEDEEDDDDDDDDDDDDDDDDDDDDDD"
			read_2.tags = read_2.tags + [("XM", ".x...hh..xh..z.hh.hh.....h...x........x......hhh.h...x.hh...h..Z.......xhh..Zxh.h..h..h....x...xh....")] + [('XR', 'GA')] + [('XG', 'GA')]
			return read_1, read_2

		def buildCTOTPE():
			read_1 = pysam.AlignedRead()
			read_1.qname = "SRR400564.6667900_HAL:1133:C010EABXX:8:2207:16412:102567_length=101"
			read_1.seq = "CGACCCCCCATCTATTCATCCATCCACCCCCCCCCCCACCCATCCATTAATTTATTCATCCATCCACCCACCCATCCACCATCCATTCAACCATCCATCCA"
			read_1.flag = 99
			read_1.tid = 0
			read_1.pos = 3287553
			read_1.mapq = 255
			read_1.cigar = [(0, 101)]
			read_1.rnext = 0
			read_1.pnext = 3287383
			read_1.isize = -271
			read_1.qual = "#########################################73GE@CEGECGDFAFAB@IGD?FB0HF?1JHGF@HFAJIHFJJIHGEA22222224B@:4"
			read_1.tags = read_1.tags + [('XM', 'Z...HH.HH..H....H..HH..HH..HH.HHH.HHH..HH..HH...........H..HH..HH.HHH.HHH..HH.HH..HH...H..HH..HH..HH.')] + [('XR', 'GA')] + [('XG', 'CT')]
			read_2 = pysam.AlignedRead()
			read_2.qname = "SRR400564.6667900_HAL:1133:C010EABXX:8:2207:16412:102567_length=101"
			read_2.seq = "ATCCACCATCTATTCATCCATCCGTCCACCCACCCATCCATCCATTAATTATCCATCCACCCACCCATCCACCATCCATTCATCCATCCATCCATCCATAC"
			read_2.flag = 147
			read_2.tid = 0
			read_2.pos = 3287383
			read_2.mapq = 255
			read_2.cigar = [(0, 101)]
			read_2.rnext = 0
			read_2.pnext = 3287553
			read_2.isize = 271
			read_2.qual = "CCCFFFFFHHHHHJJJJJJJJJJIIJIGHIJIJJJJJJJJJJJJJJJJJJGGCBFHCHGGCHIIJHHHFFFFECEDEDDDDE@>@CD:ACCDDD>>?::@C"
			read_2.tags = read_2.tags + [('XM', '..HH.HH..H....H..HH..XZ..HH.HHH.HHH..HH..HH.........HH..HH.HHH.HHH..HH.HH..HH...H..HH..HH..HH..HH...H')] + [('XR', 'CT')] + [('XG', 'CT')]
			return read_1, read_2

		def buildCTOBPE():
			read_1 = pysam.AlignedRead()
			read_1.qname = "SRR400564.4547217_HAL:1133:C010EABXX:8:2105:21225:192741_length=101"
			read_1.seq = "AAGGAAGGAGGGAAGGAAGGAAATAAAGAAAGGAAAAAAGGAAAGAAAGAAAAATAAAGAAATAAAGGAAGGAGGGAAGGAAGGAAAGAATGAAAGAAAGA"
			read_1.flag = 83
			read_1.tid = 0
			read_1.pos = 55291120
			read_1.mapq = 255
			read_1.cigar = [(0, 101)]
			read_1.rnext = 0
			read_1.pnext = 55291173
			read_1.isize = 154
			read_1.qual = "1:BD42222300<CGHIIIJIIJJIIJGIIIIIIIJJIIJIJIIIIIGHIIIIHHGGGFFFFFFDECCDCBBDB@BBA<?BC<0<AC?CB@:@CC@CBC:>"
			read_1.tags = read_1.tags + [('XM', 'z..Hh.HHh..Hh..Hh...h...h...h...H.h..hh..h...h..Hh...h...h.....hh......Hh...h...h..........H..hH.hhH.')] + [('XR', 'GA')] + [('XG', 'GA')]
			read_2 = pysam.AlignedRead()
			read_2.qname = "SRR400564.4547217_HAL:1133:C010EABXX:8:2105:21225:192741_length=101"
			read_2.seq = "AAAAGAAAAAGGAAAAAAGGAAAGAAAGAAAAATAAATGAAGGAGGGAAGGAAGGAAAGAAAGAAGGGAAGAAAGAAAAGTAAATGAAGGAGGGAAAGAAG"
			read_2.flag = 163
			read_2.tid = 0
			read_2.pos = 55291173
			read_2.mapq = 255
			read_2.cigar = [(0, 101)]
			read_2.rnext = 0
			read_2.pnext = 55291120
			read_2.isize = -154
			read_2.qual = "5DDA;DDEC@66DHHHHC=GJJIGJJIGJJJJIIJIHGHIHHHGFFFFBIIJHHGJJIGJJJJJIIHHJIGJJJGJJJJJJJJJIJJJHHHHHFFFDDBBB"
			read_2.tags = read_2.tags + [('XM', 'h...H.....HH......HH...H...H..........H..HH.HHH..HH..HH...H...H..HHH......H..........H..HH.HHH...H..H')] + [('XR', 'CT')] + [('XG', 'GA')]
			return read_1, read_2

		# Build reads
		self.otrse = buildOTSE()
		self.obrse = buildOBSE()
		self.ctotrse = buildCTOTSE()
		self.ctobrse = buildCTOBSE()
		self.otrpe_1, self.otrpe_2 = buildOTPE()
		self.obrpe_1, self.obrpe_2 = buildOBPE()
		self.ctotpe_1, self.ctotpe_2 = buildCTOTPE()
		self.ctobpe_1, self.ctobpe_2 = buildCTOBPE()

	def test_ot_se(self):
		self.assertEqual(get_strand(self.otrse), 'OT')

	def test_ob_se(self):
		self.assertEqual(get_strand(self.obrse), 'OB')

	def test_ctot_se(self):
		self.assertEqual(get_strand(self.ctotrse), 'OT')

	def test_ctob_se(self):
		self.assertEqual(get_strand(self.ctobrse), 'OB')

	def test_ot_pe(self):
		self.assertEqual(get_strand(self.otrpe_1), 'OT')
		self.assertEqual(get_strand(self.otrpe_2), 'OT')

	def test_ob_pe(self):
		self.assertEqual(get_strand(self.obrpe_1), 'OB')
		self.assertEqual(get_strand(self.obrpe_2), 'OB')

	def test_ctot_pe(self):
		self.assertEqual(get_strand(self.ctotpe_1), 'OT')
		self.assertEqual(get_strand(self.ctotpe_2), 'OT')

	def test_ctob_pe(self):
		self.assertEqual(get_strand(self.ctobpe_1), 'OB')
		self.assertEqual(get_strand(self.ctobpe_2), 'OB')

if __name__ == '__main__':
    unittest.main(verbosity=2)
