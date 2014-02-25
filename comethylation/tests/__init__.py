'''unit testing code for comethylation.
'''

import unittest
import pysam
import sys

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
			read.rname = 1
			read.pos = 451
			read.mapq = 255
			read.cigar = [(0,61)]
			read.rnext = 1
			read.mpos = 513
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
			read.rname = 1
			read.pos = 513
			read.mapq = 255
			read.cigar = [(0,59)]
			read.rnext = 1
			read.mpos = 451
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
			read.rname = 1
			read.pos = 561
			read.mapq = 255
			read.cigar = [(0,63)]
			read.rnext = 1
			read.mpos = 493
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
			read.rname = 1
			read.pos = 493
			read.mapq = 255
			read.cigar = [(0,67)]
			read.rnext = 1
			read.mpos = 561
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
			read.rname = 1
			read.pos = 451
			read.mapq = 255
			read.cigar = [(0,61)]
			read.rnext = 1
			read.mpos = 513
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
			read.rname = 1
			read.pos = 513
			read.mapq = 255
			read.cigar = [(0,59)]
			read.rnext = 1
			read.mpos = 451
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
			read.rname = 1
			read.pos = 561
			read.mapq = 255
			read.cigar = [(0,63)]
			read.rnext = 1
			read.mpos = 493
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
			read.rname = 1
			read.pos = 493
			read.mapq = 255
			read.cigar = [(0,67)]
			read.rnext = 1
			read.mpos = 561
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
			read.rname = 1
			read.pos = 451
			read.mapq = 255
			read.cigar = [(0,61)]
			read.rnext = 1
			read.mpos = 513
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
			read.rname = 1
			read.pos = 451
			read.mapq = 255
			read.cigar = [(0,61)]
			read.rnext = 1
			read.mpos = 513
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
			read.rname = 1
			read.pos = 451
			read.mapq = 255
			read.cigar = [(0,61)]
			read.rnext = 1
			read.mpos = 513
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
			read.pos = 854
			read.mapq = 255
			read.cigar = [(0,100)]
			read.rnext = 0
			read.mpos = 855
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
			read.pos = 855
			read.mapq = 255
			read.cigar = [(0,100)]
			read.rnext = 0
			read.mpos = 854
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
		self.mod_read_1.mpos = self.read_1.mpos
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
			read.pos = 854
			read.mapq = 255
			read.cigar = [(0,100)]
			read.rnext = 0
			read.mpos = 855
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
			read.pos = 854
			read.mapq = 255
			read.cigar = [(0,100)]
			read.rnext = 0
			read.mpos = 855
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

# FIXME: Remove?
if __name__ == '__main__':
    unittest.main()