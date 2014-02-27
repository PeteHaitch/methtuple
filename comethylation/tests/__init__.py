'''unit testing code for comethylation.
'''

import unittest
import pysam
import sys
import tempfile
import os
import re

from comethylation import *

from comethylation.mtuple import *
from comethylation.funcs import *

class TestIgnoreFirstNBases(unittest.TestCase):
	'''Test the function ignore_first_n_bases
	'''

	def setUp(self):

		def buildOTRead1():
			'''build an example read_1 aligned to OT-strand.
			'''

			read = pysam.AlignedRead()
			read.qname = "ADS-adipose_chr1_8"
			read.seq = "AATTTTAATTTTAATTTTTGCGGTATTTTTAGTCGGTTCGTTCGTTCGGGTTTGATTTGAG"
			read.flag = 99
			read.rname = 0
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
			read.rname = 0
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
			read.rname = 0
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
			read.rname = 0
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
		self.assertEqual(ignore_first_n_bases(self.otr_1, self.otm_1, 0), [18, 20, 33, 38, 42, 46])
		self.assertEqual(ignore_first_n_bases(self.otr_2, self.otm_2, 0), [12, 29, 50, 58])
		self.assertEqual(ignore_first_n_bases(self.obr_1, self.obm_1, 0), [3, 11, 17, 19, 29, 49, 57, 60])
		self.assertEqual(ignore_first_n_bases(self.obr_2, self.obm_2, 0), [1, 5, 33, 50])

	def test_n_off_by_one(self):
		# Shouldn't change methylation indexes
		self.assertEqual(ignore_first_n_bases(self.otr_1, self.otm_1, 18), [18, 20, 33, 38, 42, 46])
		self.assertEqual(ignore_first_n_bases(self.otr_2, self.otm_2, 0), [12, 29, 50, 58])
		self.assertEqual(ignore_first_n_bases(self.obr_1, self.obm_1, 2), [3, 11, 17, 19, 29, 49, 57, 60])
		self.assertEqual(ignore_first_n_bases(self.obr_2, self.obm_2, 1), [1, 5, 33, 50])

	def test__ignore_one(self):
		# Should remove one element from methylation indexes in accordance with the strand/orientation of the read
		self.assertEqual(ignore_first_n_bases(self.otr_1, self.otm_1, 19), [20, 33, 38, 42, 46])
		self.assertEqual(ignore_first_n_bases(self.otr_2, self.otm_2, 1), [12, 29, 50])
		self.assertEqual(ignore_first_n_bases(self.obr_1, self.obm_1, 3), [3, 11, 17, 19, 29, 49, 57])
		self.assertEqual(ignore_first_n_bases(self.obr_2, self.obm_2, 2), [5, 33, 50])

	def test_bad_n(self):
		# Should raise an exception
		self.assertRaises(ValueError, ignore_first_n_bases, self.otr_1, self.otm_1, -10)
		self.assertRaises(ValueError, ignore_first_n_bases, self.otr_1, self.otm_1, 3.4)

	def test_ignore_all(self):
		# Should remove all elements from methylation indexes (assuming read-lengths are < 100,000,000)
		self.assertEqual(ignore_first_n_bases(self.otr_1, self.otm_1, 100000000), [])
		self.assertEqual(ignore_first_n_bases(self.otr_2, self.otm_2, 100000000), [])
		self.assertEqual(ignore_first_n_bases(self.obr_1, self.obm_1, 100000000), [])
		self.assertEqual(ignore_first_n_bases(self.obr_2, self.obm_2, 100000000), [])

class TestIgnoreLastNBases(unittest.TestCase):
	'''Test the function ignore_last_n_bases
	'''

	def setUp(self):

		def buildOTRead1():
			'''build an example read_1 aligned to OT-strand.
			'''

			read = pysam.AlignedRead()
			read.qname = "ADS-adipose_chr1_8"
			read.seq = "AATTTTAATTTTAATTTTTGCGGTATTTTTAGTCGGTTCGTTCGTTCGGGTTTGATTTGAG"
			read.flag = 99
			read.rname = 0
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
			read.rname = 0
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
			read.rname = 0
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
			read.rname = 0
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
		self.assertEqual(ignore_last_n_bases(self.otr_1, self.otm_1, 0), [18, 20, 33, 38, 42, 46])
		self.assertEqual(ignore_last_n_bases(self.otr_2, self.otm_2, 0), [12, 29, 50, 58])
		self.assertEqual(ignore_last_n_bases(self.obr_1, self.obm_1, 0), [3, 11, 17, 19, 29, 49, 57, 60])
		self.assertEqual(ignore_last_n_bases(self.obr_2, self.obm_2, 0), [1, 5, 33, 50])

	def test_n_off_by_one(self):
		# Shouldn't change methylation indexes
		self.assertEqual(ignore_last_n_bases(self.otr_1, self.otm_1, 14), [18, 20, 33, 38, 42, 46])
		self.assertEqual(ignore_last_n_bases(self.otr_2, self.otm_2, 12), [12, 29, 50, 58])
		self.assertEqual(ignore_last_n_bases(self.obr_1, self.obm_1, 3), [3, 11, 17, 19, 29, 49, 57, 60])
		self.assertEqual(ignore_last_n_bases(self.obr_2, self.obm_2, 16), [1, 5, 33, 50])

	def test_ignore_one(self):
		# Should remove one element from methylation indexes in accordance with the strand/orientation of the read
		self.assertEqual(ignore_last_n_bases(self.otr_1, self.otm_1, 15), [18, 20, 33, 38, 42])
		self.assertEqual(ignore_last_n_bases(self.otr_2, self.otm_2, 13), [29, 50, 58])
		self.assertEqual(ignore_last_n_bases(self.obr_1, self.obm_1, 4), [11, 17, 19, 29, 49, 57, 60])
		self.assertEqual(ignore_last_n_bases(self.obr_2, self.obm_2, 17), [1, 5, 33])

	def test_bad_n(self):
		# Should raise an exception
		self.assertRaises(ValueError, ignore_last_n_bases, self.otr_1, self.otm_1, -10)
		self.assertRaises(ValueError, ignore_last_n_bases, self.otr_1, self.otm_1, 3.4)

	def test_ignore_all(self):
		# Should remove all elements from methylation indexes (assuming read-lengths are < 100,000,000)
		self.assertEqual(ignore_last_n_bases(self.otr_1, self.otm_1, 100000000), [])
		self.assertEqual(ignore_last_n_bases(self.otr_2, self.otm_2, 100000000), [])
		self.assertEqual(ignore_last_n_bases(self.obr_1, self.obm_1, 100000000), [])
		self.assertEqual(ignore_last_n_bases(self.obr_2, self.obm_2, 100000000), [])

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
			read.rname = 0
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
			read.rname = 0
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
			read.rname = 0
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
			read.rname = 0
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
			read.rname = 0
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
		self.mod_read_1.rname = self.read_1.rname
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
		# Should raise an exception
		self.assertRaises(ValueError, is_overlapping_sequence_identical, self.read_1, self.read_2, -10, 'sequence')
		self.assertRaises(ValueError, is_overlapping_sequence_identical, self.read_1, self.read_2, 3.4, 'sequence')

	def test_bad_overlap_check(self):
		# Should raise an exception
		self.assertRaises(ValueError, is_overlapping_sequence_identical, self.read_1, self.read_2, 10, 'apples')
		self.assertRaises(ValueError, is_overlapping_sequence_identical, self.read_1, self.read_2, 10, 'Bismark') # Should be 'bismark'

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
			read.rname = 0
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
			read.rname = 0
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
			read.rname = 0
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
			read.rname = 0
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
			read.rname = 1
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
			read.rname = 0
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
			read.rname = 0
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
			read.rname = 0
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
		self.BAM.close()
		self.BAM = pysam.Samfile(self.BAM.filename, 'rb')
		self.m1ot, self.nmlifot = extract_and_update_methylation_index_from_single_end_read(read = self.otr, BAM = self.BAM, methylation_m_tuples = {}, m = 1, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1)
		self.m2ot, self.nmlifot = extract_and_update_methylation_index_from_single_end_read(read = self.otr, BAM = self.BAM, methylation_m_tuples = {}, m = 2, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1)
		self.m3ot, self.nmlifot = extract_and_update_methylation_index_from_single_end_read(read = self.otr, BAM = self.BAM, methylation_m_tuples = {}, m = 3, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1)
		self.m4ot, self.nmlifot = extract_and_update_methylation_index_from_single_end_read(read = self.otr, BAM = self.BAM, methylation_m_tuples = {}, m = 4, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1)
		self.m5ot, self.nmlifot = extract_and_update_methylation_index_from_single_end_read(read = self.otr, BAM = self.BAM, methylation_m_tuples = {}, m = 5, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1)
		self.m1ob, self.nmlifob = extract_and_update_methylation_index_from_single_end_read(read = self.obr, BAM = self.BAM, methylation_m_tuples = {}, m = 1, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1)
		self.m2ob, self.nmlifob = extract_and_update_methylation_index_from_single_end_read(read = self.obr, BAM = self.BAM, methylation_m_tuples = {}, m = 2, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1)
		self.m3ob, self.nmlifob = extract_and_update_methylation_index_from_single_end_read(read = self.obr, BAM = self.BAM, methylation_m_tuples = {}, m = 3, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1)
		self.m2cgchg, self.nmlifotcgchg = extract_and_update_methylation_index_from_single_end_read(read = self.otr, BAM = self.BAM, methylation_m_tuples = {}, m = 2, methylation_type = 'CG/CHG', methylation_pattern = re.compile(r'[ZzXx]'), ignore_start_r1 = 0, ignore_end_r1 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1)

	def test_correct_number_of_m_tuples(self):
		self.assertEqual(len(self.m1ot), 5)
		self.assertEqual(len(self.m2ot), 4)
		self.assertEqual(len(self.m3ot), 3)
		self.assertEqual(len(self.m4ot), 2)
		self.assertEqual(len(self.m5ot), 1)
		self.assertEqual(len(self.m1ob), 3)
		self.assertEqual(len(self.m2ob), 2)
		self.assertEqual(len(self.m3ob), 1)

	def test_correct_m_tuple_ids(self):
		self.assertItemsEqual(self.m1ot.keys(), ['chr1:4562', 'chr1:4604', 'chr1:4579', 'chr1:4573', 'chr1:4610'])
		self.assertItemsEqual(self.m2ot.keys(), ['chr1:4579-4604', 'chr1:4562-4573', 'chr1:4604-4610', 'chr1:4573-4579'])
		self.assertItemsEqual(self.m3ot.keys(), ['chr1:4562-4573-4579', 'chr1:4573-4579-4604', 'chr1:4579-4604-4610'])
		self.assertItemsEqual(self.m4ot.keys(), ['chr1:4562-4573-4579-4604', 'chr1:4573-4579-4604-4610'])
		self.assertItemsEqual(self.m5ot.keys(), ['chr1:4562-4573-4579-4604-4610'])
		self.assertItemsEqual(self.m1ob.keys(), ['chr1:3400', 'chr1:3366', 'chr1:3391'])
		self.assertItemsEqual(self.m2ob.keys(), ['chr1:3391-3400', 'chr1:3366-3391'])
		self.assertItemsEqual(self.m3ob.keys(), ['chr1:3366-3391-3400'])

	def test_correct_number_of_methylation_loci_in_fragment(self):
		self.assertEqual(self.nmlifot, 5)
		self.assertEqual(self.nmlifob, 3)

	def test_correct_methylation_type(self):
		self.assertEqual(self.m1ot['chr1:4562'].methylation_type, 'CG')
		self.assertEqual(self.m1ob['chr1:3366'].methylation_type, 'CG')
		self.assertEqual(self.m2cgchg['chr1:4601-4604'].methylation_type, 'CG/CHG')

	def test_counts(self):
		self.assertEqual(self.m1ot['chr1:4562'].counts, {'M': 1, 'U': 0})

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
			read.rname = 0
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
			read.rname = 0
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
			read.rname = 0
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
			read.rname = 0
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
		self.BAM.close()
		self.BAM = pysam.Samfile(self.BAM.filename, 'rb')
		self.FAILED_QC = open(tempfile.mkstemp()[1], 'w')
		self.m1ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 1, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m2ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 2, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m3ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 3, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m4ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 4, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m5ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 5, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m6ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 6, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m7ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 7, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m8ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 8, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m9ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 9, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m10ot, self.nmlifot, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 10, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m1ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 1, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m2ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 2, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m3ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 3, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m4ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 4, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m5ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 5, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m6ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 6, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m7ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 7, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m8ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 8, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m9ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 9, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m10ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 10, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m11ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 11, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m12ob, self.nmlifob, self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.obr_1, read_2 = self.obr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 12, methylation_type = 'CG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)
		self.m2cgchg, self.nmlifotcgchg , self.nfsdtbo = extract_and_update_methylation_index_from_paired_end_reads(read_1 = self.otr_1, read_2 = self.otr_2, BAM = self.BAM, methylation_m_tuples = {}, m = 2, methylation_type = 'CG/CHG', methylation_pattern = re.compile(r'[Zz]'), ignore_start_r1 = 0, ignore_end_r1 = 0, ignore_start_r2 = 0, ignore_end_r2 = 0, min_qual = 0, phred_offset = 33, ob_strand_offset = 1, overlap_check = 'bismark', n_fragment_skipped_due_to_bad_overlap = 0, FAILED_QC = self.FAILED_QC)

	def test_correct_number_of_m_tuples(self):
		self.assertEqual(len(self.m1ot), 10)
		self.assertEqual(len(self.m2ot), 9)
		self.assertEqual(len(self.m3ot), 8)
		self.assertEqual(len(self.m4ot), 7)
		self.assertEqual(len(self.m5ot), 6)
		self.assertEqual(len(self.m6ot), 5)
		self.assertEqual(len(self.m7ot), 4)
		self.assertEqual(len(self.m8ot), 3)
		self.assertEqual(len(self.m9ot), 2)
		self.assertEqual(len(self.m10ot), 1)
		self.assertEqual(len(self.m1ob), 12)
		self.assertEqual(len(self.m2ob), 11)
		self.assertEqual(len(self.m3ob), 10)
		self.assertEqual(len(self.m4ob), 9)
		self.assertEqual(len(self.m5ob), 8)
		self.assertEqual(len(self.m6ob), 7)
		self.assertEqual(len(self.m7ob), 6)
		self.assertEqual(len(self.m8ob), 5)
		self.assertEqual(len(self.m9ob), 4)
		self.assertEqual(len(self.m10ob), 3)
		self.assertEqual(len(self.m11ob), 2)
		self.assertEqual(len(self.m12ob), 1)

	def test_correct_m_tuple_ids(self):
		self.assertItemsEqual(self.m1ot.keys(), ['chr1:563', 'chr1:571', 'chr1:525', 'chr1:493', 'chr1:469', 'chr1:484', 'chr1:489', 'chr1:497', 'chr1:471', 'chr1:542'])
		self.assertItemsEqual(self.m2ot.keys(), ['chr1:497-525', 'chr1:493-497', 'chr1:489-493', 'chr1:563-571', 'chr1:542-563', 'chr1:469-471', 'chr1:484-489', 'chr1:525-542', 'chr1:471-484'])
		self.assertItemsEqual(self.m3ot.keys(), ['chr1:493-497-525', 'chr1:469-471-484', 'chr1:471-484-489', 'chr1:542-563-571', 'chr1:489-493-497', 'chr1:484-489-493', 'chr1:525-542-563', 'chr1:497-525-542'])
		self.assertItemsEqual(self.m4ot.keys(), ['chr1:497-525-542-563', 'chr1:471-484-489-493', 'chr1:469-471-484-489', 'chr1:484-489-493-497', 'chr1:525-542-563-571', 'chr1:493-497-525-542', 'chr1:489-493-497-525'])
		self.assertItemsEqual(self.m5ot.keys(), ['chr1:469-471-484-489-493', 'chr1:493-497-525-542-563', 'chr1:471-484-489-493-497', 'chr1:489-493-497-525-542', 'chr1:497-525-542-563-571', 'chr1:484-489-493-497-525'])
		self.assertItemsEqual(self.m6ot.keys(), ['chr1:471-484-489-493-497-525', 'chr1:469-471-484-489-493-497', 'chr1:489-493-497-525-542-563', 'chr1:493-497-525-542-563-571', 'chr1:484-489-493-497-525-542'])
		self.assertItemsEqual(self.m7ot.keys(), ['chr1:471-484-489-493-497-525-542', 'chr1:484-489-493-497-525-542-563', 'chr1:469-471-484-489-493-497-525', 'chr1:489-493-497-525-542-563-571'])
		self.assertItemsEqual(self.m8ot.keys(), ['chr1:484-489-493-497-525-542-563-571', 'chr1:469-471-484-489-493-497-525-542', 'chr1:471-484-489-493-497-525-542-563'])
		self.assertItemsEqual(self.m9ot.keys(), ['chr1:469-471-484-489-493-497-525-542-563', 'chr1:471-484-489-493-497-525-542-563-571'])
		self.assertItemsEqual(self.m10ot.keys(), ['chr1:469-471-484-489-493-497-525-542-563-571'])
		self.assertItemsEqual(self.m1ob.keys(), ['chr1:589', 'chr1:577', 'chr1:563', 'chr1:571', 'chr1:579', 'chr1:525', 'chr1:493', 'chr1:609', 'chr1:497', 'chr1:620', 'chr1:617', 'chr1:542'])
		self.assertItemsEqual(self.m2ob.keys(), ['chr1:497-525', 'chr1:493-497', 'chr1:617-620', 'chr1:563-571', 'chr1:542-563', 'chr1:589-609', 'chr1:579-589', 'chr1:577-579', 'chr1:571-577', 'chr1:525-542', 'chr1:609-617'])
		self.assertItemsEqual(self.m3ob.keys(), ['chr1:542-563-571', 'chr1:589-609-617', 'chr1:563-571-577', 'chr1:493-497-525', 'chr1:579-589-609', 'chr1:571-577-579', 'chr1:609-617-620', 'chr1:525-542-563', 'chr1:577-579-589', 'chr1:497-525-542'])
		self.assertItemsEqual(self.m4ob.keys(), ['chr1:497-525-542-563', 'chr1:563-571-577-579', 'chr1:493-497-525-542', 'chr1:542-563-571-577', 'chr1:525-542-563-571', 'chr1:571-577-579-589', 'chr1:577-579-589-609', 'chr1:579-589-609-617', 'chr1:589-609-617-620'])
		self.assertItemsEqual(self.m5ob.keys(), ['chr1:577-579-589-609-617', 'chr1:571-577-579-589-609', 'chr1:542-563-571-577-579', 'chr1:493-497-525-542-563', 'chr1:525-542-563-571-577', 'chr1:497-525-542-563-571', 'chr1:579-589-609-617-620', 'chr1:563-571-577-579-589'])
		self.assertItemsEqual(self.m6ob.keys(), ['chr1:577-579-589-609-617-620', 'chr1:525-542-563-571-577-579', 'chr1:493-497-525-542-563-571', 'chr1:542-563-571-577-579-589', 'chr1:497-525-542-563-571-577', 'chr1:571-577-579-589-609-617', 'chr1:563-571-577-579-589-609'])
		self.assertItemsEqual(self.m7ob.keys(), ['chr1:571-577-579-589-609-617-620', 'chr1:542-563-571-577-579-589-609', 'chr1:563-571-577-579-589-609-617', 'chr1:525-542-563-571-577-579-589', 'chr1:497-525-542-563-571-577-579', 'chr1:493-497-525-542-563-571-577'])
		self.assertItemsEqual(self.m8ob.keys(), ['chr1:563-571-577-579-589-609-617-620', 'chr1:542-563-571-577-579-589-609-617', 'chr1:525-542-563-571-577-579-589-609', 'chr1:497-525-542-563-571-577-579-589', 'chr1:493-497-525-542-563-571-577-579'])
		self.assertItemsEqual(self.m9ob.keys(), ['chr1:497-525-542-563-571-577-579-589-609', 'chr1:493-497-525-542-563-571-577-579-589', 'chr1:525-542-563-571-577-579-589-609-617', 'chr1:542-563-571-577-579-589-609-617-620'])
		self.assertItemsEqual(self.m10ob.keys(), ['chr1:525-542-563-571-577-579-589-609-617-620', 'chr1:493-497-525-542-563-571-577-579-589-609', 'chr1:497-525-542-563-571-577-579-589-609-617'])
		self.assertItemsEqual(self.m11ob.keys(), ['chr1:497-525-542-563-571-577-579-589-609-617-620', 'chr1:493-497-525-542-563-571-577-579-589-609-617'])
		self.assertItemsEqual(self.m12ob.keys(), ['chr1:493-497-525-542-563-571-577-579-589-609-617-620'])

	def test_correct_number_of_methylation_loci_in_fragment(self):
		self.assertEqual(self.nmlifot, 10)
		self.assertEqual(self.nmlifob, 12)

	def test_correct_methylation_type(self):
		self.assertEqual(self.m1ot['chr1:563'].methylation_type, 'CG')
		self.assertEqual(self.m1ob['chr1:589'].methylation_type, 'CG')
		self.assertEqual(self.m2cgchg['chr1:497-525'].methylation_type, 'CG/CHG')

	def test_counts(self):
		self.assertEqual(self.m1ot['chr1:563'].counts, {'M': 1, 'U': 0})

	def tearDown(self):
		os.remove(self.BAM.filename)
		os.remove(self.FAILED_QC.name)

class TestWithinFragmentComethylationMTuple(unittest.TestCase):
	'''Test the class WithinFragmentComethylationMTuple and its methods
	'''

	def setUp(self):

		def buildOTRead():
			'''build an example read aligned to OT-strand.
			'''
			read = pysam.AlignedRead()
			read.qname = "@SALK_2077_FC6295TAAXX:2:107:9396:15019#0/1"
			read.seq = "GGGGAAGGTGTTATGGAGTTTTTTACGATTTTTAGTCGTTTTCGTTTTTTTTTGTTTGTGGTTGTTGCGGTGGCGGTAGAGGAGGG"
			read.flag = 0
			read.rname = 0
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
			read.rname = 0
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
			read.rname = 0
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
			read.rname = 0
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
			read.rname = 0
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
			read.rname = 0
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
		self.BAM.close()
		self.BAM = pysam.Samfile(self.BAM.filename, 'rb')
		self.wfotr = WithinFragmentComethylationMTuple(self.BAM.getrname(self.otr.tid), self.otr.tid, 2, [4562, 4573], 'CG')
		self.wfobr = WithinFragmentComethylationMTuple(self.BAM.getrname(self.obr.tid), self.otr.tid, 2, [3366, 3391], 'CG')
		self.wfotrcgchg = WithinFragmentComethylationMTuple(self.BAM.getrname(self.otr.tid), self.otr.tid, 3, [4562, 4569, 4573], 'CG/CHG')
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
		self.BAMPE.close()
		self.BAMPE = pysam.Samfile(self.BAMPE.filename, 'rb')
		self.wfotrpe = WithinFragmentComethylationMTuple(self.BAMPE.getrname(self.otr_1.tid), self.otr_1.tid, 2, [497, 525], 'CG')
		self.wfobrpe = WithinFragmentComethylationMTuple(self.BAMPE.getrname(self.obr_1.tid), self.obr_1.tid, 2, [563, 571], 'CG') 
		self.wfotrpecgchg = WithinFragmentComethylationMTuple(self.BAMPE.getrname(self.otr_1.tid), self.otr_1.tid, 3, [497, 525, 534], 'CG/CHG')

	def test_init_se(self):
		self.assertEqual(self.wfotr.methylation_type, 'CG')
		self.assertEqual(self.wfotr.chromosome, self.BAM.getrname(self.otr.tid))
		self.assertEqual(self.wfotr.chromosome_index, self.otr.tid)
		self.assertEqual(self.wfotr.positions, [4562, 4573])
		self.assertEqual(self.wfotr.counts, {'UU': 0, 'UM': 0, 'MM': 0, 'MU': 0})
		self.assertEqual(self.wfobr.methylation_type, 'CG')
		self.assertEqual(self.wfobr.chromosome, self.BAM.getrname(self.obr.tid))
		self.assertEqual(self.wfobr.chromosome_index, self.obr.tid)
		self.assertEqual(self.wfobr.positions, [3366, 3391])
		self.assertEqual(self.wfobr.counts, {'UU': 0, 'UM': 0, 'MM': 0, 'MU': 0})
		self.assertEqual(self.wfotrcgchg.methylation_type, 'CG/CHG')
		self.assertEqual(self.wfotrcgchg.chromosome, self.BAM.getrname(self.otr.tid))
		self.assertEqual(self.wfotrcgchg.chromosome_index, self.otr.tid)
		self.assertEqual(self.wfotrcgchg.positions, [4562, 4569, 4573])
		self.assertEqual(self.wfotrcgchg.counts, {'MUU': 0, 'UMU': 0, 'UUU': 0, 'MMU': 0, 'UUM': 0, 'MUM': 0, 'UMM': 0, 'MMM': 0})

	def test_init_pe(self):
		self.assertEqual(self.wfotrpe.methylation_type, 'CG')
		self.assertEqual(self.wfotrpe.chromosome, self.BAM.getrname(self.otr_1.tid))
		self.assertEqual(self.wfotrpe.chromosome_index, self.otr_1.tid)
		self.assertEqual(self.wfotrpe.positions, [497, 525])
		self.assertEqual(self.wfotrpe.counts, {'UU': 0, 'UM': 0, 'MM': 0, 'MU': 0})
		self.assertEqual(self.wfobrpe.methylation_type, 'CG')
		self.assertEqual(self.wfobrpe.chromosome, self.BAM.getrname(self.obr_1.tid))
		self.assertEqual(self.wfobrpe.chromosome_index, self.obr_1.tid)
		self.assertEqual(self.wfobrpe.positions, [563, 571])
		self.assertEqual(self.wfobrpe.counts, {'UU': 0, 'UM': 0, 'MM': 0, 'MU': 0})
		self.assertEqual(self.wfotrpecgchg.methylation_type, 'CG/CHG')
		self.assertEqual(self.wfotrpecgchg.chromosome, self.BAM.getrname(self.otr_1.tid))
		self.assertEqual(self.wfotrcgchg.chromosome_index, self.otr_1.tid)
		self.assertEqual(self.wfotrpecgchg.positions, [497, 525, 534])
		self.assertEqual(self.wfotrpecgchg.counts, {'MUU': 0, 'UMU': 0, 'UUU': 0, 'MMU': 0, 'UUM': 0, 'MUM': 0, 'UMM': 0, 'MMM': 0})

	def test_increment_count_ot_se(self):
		self.wfotr.increment_count('MM', self.otr, None)
		self.assertEqual(self.wfotr.counts, {'UU': 0, 'UM': 0, 'MM': 1, 'MU': 0})
		self.wfotr.increment_count('MU', self.otr, None)
		self.assertEqual(self.wfotr.counts, {'UU': 0, 'UM': 0, 'MM': 1, 'MU': 1})
		self.wfotr.increment_count('UM', self.otr, None)
		self.assertEqual(self.wfotr.counts, {'UU': 0, 'UM': 1, 'MM': 1, 'MU': 1})
		self.wfotr.increment_count('UU', self.otr, None)
		self.assertEqual(self.wfotr.counts, {'UU': 1, 'UM': 1, 'MM': 1, 'MU': 1})
		self.wfotr.increment_count('MM', self.otr, None)
		self.assertEqual(self.wfotr.counts, {'UU': 1, 'UM': 1, 'MM': 2, 'MU': 1})
	def test_increment_count_ob_se(self):
		self.wfobr.increment_count('MM', self.obr, None)
		self.assertEqual(self.wfobr.counts, {'UU': 0, 'UM': 0, 'MM': 1, 'MU': 0})
		self.wfobr.increment_count('MU', self.obr, None)
		self.assertEqual(self.wfobr.counts, {'UU': 0, 'UM': 0, 'MM': 1, 'MU': 1})
		self.wfobr.increment_count('UM', self.obr, None)
		self.assertEqual(self.wfobr.counts, {'UU': 0, 'UM': 1, 'MM': 1, 'MU': 1})
		self.wfobr.increment_count('UU', self.obr, None)
		self.assertEqual(self.wfobr.counts, {'UU': 1, 'UM': 1, 'MM': 1, 'MU': 1})
		self.wfobr.increment_count('MM', self.obr, None)
		self.assertEqual(self.wfobr.counts, {'UU': 1, 'UM': 1, 'MM': 2, 'MU': 1})
	def test_increment_count_ctot_se(self):
		self.otr.tags = []
		self.otr.tags = self.otr.tags + [("XG", "CT")] + [("XM", "...........h......hhhhh..Z....hhx...Z..hh.Z..hh.hh.x..hx.....x..x..Z.....Z..x.........")] + [("XR", "GA")]
		self.wfotr.increment_count('MM', self.otr, None)
		self.assertEqual(self.wfotr.counts, {'UU': 0, 'UM': 0, 'MM': 1, 'MU': 0})
		self.wfotr.increment_count('MU', self.otr, None)
		self.assertEqual(self.wfotr.counts, {'UU': 0, 'UM': 0, 'MM': 1, 'MU': 1})
		self.wfotr.increment_count('UM', self.otr, None)
		self.assertEqual(self.wfotr.counts, {'UU': 0, 'UM': 1, 'MM': 1, 'MU': 1})
		self.wfotr.increment_count('UU', self.otr, None)
		self.assertEqual(self.wfotr.counts, {'UU': 1, 'UM': 1, 'MM': 1, 'MU': 1})
		self.wfotr.increment_count('MM', self.otr, None)
		self.assertEqual(self.wfotr.counts, {'UU': 1, 'UM': 1, 'MM': 2, 'MU': 1})
	def test_increment_count_ctob_se(self):
		self.obr.tags = []
		self.obr.tags = self.obr.tags + [("XG", "GA")] + [("XM", "......x...xh..x..x.......x...xh.Z..x.h..........h....x...z..xh.h..zx.h...h....hhh...")] + [("XR", "GA")]
		self.wfobr.increment_count('MM', self.obr, None)
		self.assertEqual(self.wfobr.counts, {'UU': 0, 'UM': 0, 'MM': 1, 'MU': 0})
		self.wfobr.increment_count('MU', self.obr, None)
		self.assertEqual(self.wfobr.counts, {'UU': 0, 'UM': 0, 'MM': 1, 'MU': 1})
		self.wfobr.increment_count('UM', self.obr, None)
		self.assertEqual(self.wfobr.counts, {'UU': 0, 'UM': 1, 'MM': 1, 'MU': 1})
		self.wfobr.increment_count('UU', self.obr, None)
		self.assertEqual(self.wfobr.counts, {'UU': 1, 'UM': 1, 'MM': 1, 'MU': 1})
		self.wfobr.increment_count('MM', self.obr, None)
		self.assertEqual(self.wfobr.counts, {'UU': 1, 'UM': 1, 'MM': 2, 'MU': 1})
	def test_increment_count_multiple_methylation_types_se(self):
		self.wfotrcgchg.increment_count('MMM', self.otr, None)
		self.assertEqual(self.wfotrcgchg.counts, {'MUU': 0, 'UMU': 0, 'UUU': 0, 'MMU': 0, 'UUM': 0, 'MUM': 0, 'UMM': 0, 'MMM': 1})
		self.wfotrcgchg.increment_count('MMU', self.otr, None)
		self.assertEqual(self.wfotrcgchg.counts, {'MUU': 0, 'UMU': 0, 'UUU': 0, 'MMU': 1, 'UUM': 0, 'MUM': 0, 'UMM': 0, 'MMM': 1})
		self.wfotrcgchg.increment_count('MUM', self.otr, None)
		self.assertEqual(self.wfotrcgchg.counts, {'MUU': 0, 'UMU': 0, 'UUU': 0, 'MMU': 1, 'UUM': 0, 'MUM': 1, 'UMM': 0, 'MMM': 1})
		self.wfotrcgchg.increment_count('MUU', self.otr, None)
		self.assertEqual(self.wfotrcgchg.counts, {'MUU': 1, 'UMU': 0, 'UUU': 0, 'MMU': 1, 'UUM': 0, 'MUM': 1, 'UMM': 0, 'MMM': 1})
		self.wfotrcgchg.increment_count('UMM', self.otr, None)
		self.assertEqual(self.wfotrcgchg.counts, {'MUU': 1, 'UMU': 0, 'UUU': 0, 'MMU': 1, 'UUM': 0, 'MUM': 1, 'UMM': 1, 'MMM': 1})
		self.wfotrcgchg.increment_count('UMU', self.otr, None)
		self.assertEqual(self.wfotrcgchg.counts, {'MUU': 1, 'UMU': 1, 'UUU': 0, 'MMU': 1, 'UUM': 0, 'MUM': 1, 'UMM': 1, 'MMM': 1})
		self.wfotrcgchg.increment_count('UUM', self.otr, None)
		self.assertEqual(self.wfotrcgchg.counts, {'MUU': 1, 'UMU': 1, 'UUU': 0, 'MMU': 1, 'UUM': 1, 'MUM': 1, 'UMM': 1, 'MMM': 1})
		self.wfotrcgchg.increment_count('UUU', self.otr, None)
		self.assertEqual(self.wfotrcgchg.counts, {'MUU': 1, 'UMU': 1, 'UUU': 1, 'MMU': 1, 'UUM': 1, 'MUM': 1, 'UMM': 1, 'MMM': 1})
		self.wfotrcgchg.increment_count('MMM', self.otr, None)
		self.assertEqual(self.wfotrcgchg.counts, {'MUU': 1, 'UMU': 1, 'UUU': 1, 'MMU': 1, 'UUM': 1, 'MUM': 1, 'UMM': 1, 'MMM': 2})

	def test_increment_count_ot_pe(self):
		self.wfotrpe.increment_count('MM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.counts, {'UU': 0, 'UM': 0, 'MM': 1, 'MU': 0})
		self.wfotrpe.increment_count('MU', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.counts, {'UU': 0, 'UM': 0, 'MM': 1, 'MU': 1})
		self.wfotrpe.increment_count('UM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.counts, {'UU': 0, 'UM': 1, 'MM': 1, 'MU': 1})
		self.wfotrpe.increment_count('UU', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.counts, {'UU': 1, 'UM': 1, 'MM': 1, 'MU': 1})
		self.wfotrpe.increment_count('MM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.counts, {'UU': 1, 'UM': 1, 'MM': 2, 'MU': 1})
	def test_increment_count_ob_pe(self):
		self.wfobrpe.increment_count('MM', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.counts, {'UU': 0, 'UM': 0, 'MM': 1, 'MU': 0})
		self.wfobrpe.increment_count('MU', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.counts, {'UU': 0, 'UM': 0, 'MM': 1, 'MU': 1})
		self.wfobrpe.increment_count('UM', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.counts, {'UU': 0, 'UM': 1, 'MM': 1, 'MU': 1})
		self.wfobrpe.increment_count('UU', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.counts, {'UU': 1, 'UM': 1, 'MM': 1, 'MU': 1})
		self.wfobrpe.increment_count('MM', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.counts, {'UU': 1, 'UM': 1, 'MM': 2, 'MU': 1})
	def test_increment_count_ctot_pe(self):
		self.otr_1.tags = []
		self.otr_1.tags = self.otr_1.tags + [('XG', 'CT'), ('XM', '..hhh...hhh...hhh.z.Z....hhh.x..xZ..hxZ.hxZ.hxZ....x...hx....'), ('XR', 'GA')]
		self.otr_2.tags = []
		self.otr_2.tags = self.otr_2.tags + [('XG', 'CT'), ('XM', '....x....h.xZ.hh..x......hh.xZ.....x....x......h..Z.x..H.xZ'), ('XR', 'CT')]
		self.wfotrpe.increment_count('MM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.counts, {'UU': 0, 'UM': 0, 'MM': 1, 'MU': 0})
		self.wfotrpe.increment_count('MU', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.counts, {'UU': 0, 'UM': 0, 'MM': 1, 'MU': 1})
		self.wfotrpe.increment_count('UM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.counts, {'UU': 0, 'UM': 1, 'MM': 1, 'MU': 1})
		self.wfotrpe.increment_count('UU', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.counts, {'UU': 1, 'UM': 1, 'MM': 1, 'MU': 1})
		self.wfotrpe.increment_count('MM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpe.counts, {'UU': 1, 'UM': 1, 'MM': 2, 'MU': 1})

	def test_increment_count_ctob_pe(self):
		self.obr_1.tags = []
		self.obr_1.tags = self.obr_1.tags + [('XG', 'GA'), ('XM', '...Z..x....Z.....Z.Zx.h......Zxh...x.h..x.hh.h...Z.......Z..Zx.'), ('XR', 'GA')]
		self.obr_2.tags = []
		self.obr_2.tags = self.obr_2.tags + [('XG', 'GA'), ('XM', '.z...Zxh...x....x.hh.h....x.h....Z......x.h.......Z......x.h..x.hh.'), ('XR', 'CT')]
		self.wfobrpe.increment_count('MM', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.counts, {'UU': 0, 'UM': 0, 'MM': 1, 'MU': 0})
		self.wfobrpe.increment_count('MU', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.counts, {'UU': 0, 'UM': 0, 'MM': 1, 'MU': 1})
		self.wfobrpe.increment_count('UM', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.counts, {'UU': 0, 'UM': 1, 'MM': 1, 'MU': 1})
		self.wfobrpe.increment_count('UU', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.counts, {'UU': 1, 'UM': 1, 'MM': 1, 'MU': 1})
		self.wfobrpe.increment_count('MM', self.obr_1, self.obr_2)
		self.assertEqual(self.wfobrpe.counts, {'UU': 1, 'UM': 1, 'MM': 2, 'MU': 1})		

	def test_increment_count_multiple_methylation_types_pe(self):
		self.wfotrpecgchg.increment_count('MMM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpecgchg.counts, {'MUU': 0, 'UMU': 0, 'UUU': 0, 'MMU': 0, 'UUM': 0, 'MUM': 0, 'UMM': 0, 'MMM': 1})
		self.wfotrpecgchg.increment_count('MMU', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpecgchg.counts, {'MUU': 0, 'UMU': 0, 'UUU': 0, 'MMU': 1, 'UUM': 0, 'MUM': 0, 'UMM': 0, 'MMM': 1})
		self.wfotrpecgchg.increment_count('MUM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpecgchg.counts, {'MUU': 0, 'UMU': 0, 'UUU': 0, 'MMU': 1, 'UUM': 0, 'MUM': 1, 'UMM': 0, 'MMM': 1})
		self.wfotrpecgchg.increment_count('MUU', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpecgchg.counts, {'MUU': 1, 'UMU': 0, 'UUU': 0, 'MMU': 1, 'UUM': 0, 'MUM': 1, 'UMM': 0, 'MMM': 1})
		self.wfotrpecgchg.increment_count('UMM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpecgchg.counts, {'MUU': 1, 'UMU': 0, 'UUU': 0, 'MMU': 1, 'UUM': 0, 'MUM': 1, 'UMM': 1, 'MMM': 1})
		self.wfotrpecgchg.increment_count('UMU', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpecgchg.counts, {'MUU': 1, 'UMU': 1, 'UUU': 0, 'MMU': 1, 'UUM': 0, 'MUM': 1, 'UMM': 1, 'MMM': 1})
		self.wfotrpecgchg.increment_count('UUM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpecgchg.counts, {'MUU': 1, 'UMU': 1, 'UUU': 0, 'MMU': 1, 'UUM': 1, 'MUM': 1, 'UMM': 1, 'MMM': 1})
		self.wfotrpecgchg.increment_count('UUU', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpecgchg.counts, {'MUU': 1, 'UMU': 1, 'UUU': 1, 'MMU': 1, 'UUM': 1, 'MUM': 1, 'UMM': 1, 'MMM': 1})
		self.wfotrpecgchg.increment_count('MMM', self.otr_1, self.otr_2)
		self.assertEqual(self.wfotrpecgchg.counts, {'MUU': 1, 'UMU': 1, 'UUU': 1, 'MMU': 1, 'UUM': 1, 'MUM': 1, 'UMM': 1, 'MMM': 2})

	def test_invalid_m(self):
		self.assertRaises(ValueError, WithinFragmentComethylationMTuple, self.BAM.getrname(self.otr.tid), self.otr.tid, 3, [4562, 4573], 'CG')
		self.assertRaises(ValueError, WithinFragmentComethylationMTuple, self.BAM.getrname(self.otr.tid), self.otr.tid, -3, [4562, 4573], 'CG')
		self.assertRaises(ValueError, WithinFragmentComethylationMTuple, self.BAM.getrname(self.otr.tid), self.otr.tid, 3.4, [4562, 4573], 'CG')

	def test_invalid_methylation_type(self):
		self.assertRaises(ValueError, WithinFragmentComethylationMTuple, self.BAMPE.getrname(self.otr_1.tid), self.otr_1.tid, 2, [497, 525], 'CT')
		self.assertRaises(ValueError, WithinFragmentComethylationMTuple, self.BAMPE.getrname(self.otr_1.tid), self.otr_1.tid, 2, [497, 525], 'CG-CHG')

	def test_invalid_comethylation_state(self):
		with self.assertRaises(SystemExit) as cm:
			self.wfotr.increment_count('mm', self.otr, None)
			self.assertEqual(cm.exception.code, 1)
		with self.assertRaises(SystemExit) as cm:
			self.wfotr.increment_count('MMM', self.otr, None)
			self.assertEqual(cm.exception.code, 1)

	def tearDown(self):
		os.remove(self.BAM.filename)
		os.remove(self.BAMPE.filename)

class TestGetStrand(unittest.TestCase):

	def test_ot_se(self):
		self.assertTrue(False) # TODO

	def test_ob_se(self):
		self.assertTrue(False) # TODO

	def test_ctot_se(self):
		self.assertTrue(False) # TODO

	def test_ctob_se(self):
		self.assertTrue(False) # TODO

	def test_ot_pe(self):
		self.assertTrue(False) # TODO

	def test_ob_pe(self):
		self.assertTrue(False) # TODO

	def test_ctot_pe(self):
		self.assertTrue(False) # TODO

	def test_ctob_pe(self):
		self.assertTrue(False) # TODO

	def test_remove_old_strand_check_tests_if_safe_to_do_so(self):
		self.assertTrue(False) # TODO

# FIXME: Remove?
if __name__ == '__main__':
    unittest.main()