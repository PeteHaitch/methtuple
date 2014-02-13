'''unit testing code for comethylation.
'''

import unittest
import pysam

from comethylation import *

from comethylation.mtuple import *
from comethylation.funcs import *

class IgnoreFirstNBases(unittest.TestCase):
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

	def test_N0(self):
		# Shouldn't change methylation indexes
		self.assertEqual(ignore_first_n_bases(self.otr_1, self.otm_1, 0), [18, 20, 33, 38, 42, 46])
		self.assertEqual(ignore_first_n_bases(self.otr_2, self.otm_2, 0), [12, 29, 50, 58])
		self.assertEqual(ignore_first_n_bases(self.obr_1, self.obm_1, 0), [3, 11, 17, 19, 29, 49, 57, 60])
		self.assertEqual(ignore_first_n_bases(self.obr_2, self.obm_2, 0), [1, 5, 33, 50])

	def test_NOffByOne(self):
		# Shouldn't change methylation indexes
		self.assertEqual(ignore_first_n_bases(self.otr_1, self.otm_1, 18), [18, 20, 33, 38, 42, 46])
		self.assertEqual(ignore_first_n_bases(self.otr_2, self.otm_2, 0), [12, 29, 50, 58])
		self.assertEqual(ignore_first_n_bases(self.obr_1, self.obm_1, 2), [3, 11, 17, 19, 29, 49, 57, 60])
		self.assertEqual(ignore_first_n_bases(self.obr_2, self.obm_2, 1), [1, 5, 33, 50])

	def test_IgnoreOne(self):
		# Should remove one element from methylation indexes in accordance with the strand/orientation of the read
		self.assertEqual(ignore_first_n_bases(self.otr_1, self.otm_1, 19), [20, 33, 38, 42, 46])
		self.assertEqual(ignore_first_n_bases(self.otr_2, self.otm_2, 1), [12, 29, 50])
		self.assertEqual(ignore_first_n_bases(self.obr_1, self.obm_1, 3), [3, 11, 17, 19, 29, 49, 57])
		self.assertEqual(ignore_first_n_bases(self.obr_2, self.obm_2, 2), [5, 33, 50])

	def test_BadN(self):
		# Should raise an exception
		self.assertRaises(ValueError, ignore_first_n_bases, self.otr_1, self.otm_1, -10)
		self.assertRaises(ValueError, ignore_first_n_bases, self.otr_1, self.otm_1, 3.4)

	def test_IgnoreAll(self):
		# Should remove all elements from methylation indexes (assuming read-lengths are < 100,000,000)
		self.assertEqual(ignore_first_n_bases(self.otr_1, self.otm_1, 100000000), [])
		self.assertEqual(ignore_first_n_bases(self.otr_2, self.otm_2, 100000000), [])
		self.assertEqual(ignore_first_n_bases(self.obr_1, self.obm_1, 100000000), [])
		self.assertEqual(ignore_first_n_bases(self.obr_2, self.obm_2, 100000000), [])

class IgnoreLastNBases(unittest.TestCase):
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

	def test_N0(self):
		# Shouldn't change methylation indexes
		self.assertEqual(ignore_last_n_bases(self.otr_1, self.otm_1, 0), [18, 20, 33, 38, 42, 46])
		self.assertEqual(ignore_last_n_bases(self.otr_2, self.otm_2, 0), [12, 29, 50, 58])
		self.assertEqual(ignore_last_n_bases(self.obr_1, self.obm_1, 0), [3, 11, 17, 19, 29, 49, 57, 60])
		self.assertEqual(ignore_last_n_bases(self.obr_2, self.obm_2, 0), [1, 5, 33, 50])

	def test_NOffByOne(self):
		# Shouldn't change methylation indexes
		self.assertEqual(ignore_last_n_bases(self.otr_1, self.otm_1, 14), [18, 20, 33, 38, 42, 46])
		self.assertEqual(ignore_last_n_bases(self.otr_2, self.otm_2, 12), [12, 29, 50, 58])
		self.assertEqual(ignore_last_n_bases(self.obr_1, self.obm_1, 3), [3, 11, 17, 19, 29, 49, 57, 60])
		self.assertEqual(ignore_last_n_bases(self.obr_2, self.obm_2, 16), [1, 5, 33, 50])

	def test_IgnoreOne(self):
		# Should remove one element from methylation indexes in accordance with the strand/orientation of the read
		self.assertEqual(ignore_last_n_bases(self.otr_1, self.otm_1, 15), [18, 20, 33, 38, 42])
		self.assertEqual(ignore_last_n_bases(self.otr_2, self.otm_2, 13), [29, 50, 58])
		self.assertEqual(ignore_last_n_bases(self.obr_1, self.obm_1, 4), [11, 17, 19, 29, 49, 57, 60])
		self.assertEqual(ignore_last_n_bases(self.obr_2, self.obm_2, 17), [1, 5, 33])

	def test_BadN(self):
		# Should raise an exception
		self.assertRaises(ValueError, ignore_last_n_bases, self.otr_1, self.otm_1, -10)
		self.assertRaises(ValueError, ignore_last_n_bases, self.otr_1, self.otm_1, 3.4)

	def test_IgnoreAll(self):
		# Should remove all elements from methylation indexes (assuming read-lengths are < 100,000,000)
		self.assertEqual(ignore_last_n_bases(self.otr_1, self.otm_1, 100000000), [])
		self.assertEqual(ignore_last_n_bases(self.otr_2, self.otm_2, 100000000), [])
		self.assertEqual(ignore_last_n_bases(self.obr_1, self.obm_1, 100000000), [])
		self.assertEqual(ignore_last_n_bases(self.obr_2, self.obm_2, 100000000), [])


# FIXME: Remove?
if __name__ == '__main__':
    unittest.main()