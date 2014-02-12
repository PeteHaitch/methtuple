'''unit testing code for comethylation.
'''

import unittest
import pysam # Necessary?

from comethylation import *

from comethylation.mtuple import *
from comethylation.funcs import *


class SimpleTest(unittest.TestCase):
    def test(self):
        pass

class IgnoreFirstNBases(unittest.TestCase):
	'''Test the function ignore_first_n_bases
	'''

	def buildRead( self ):
		'''build an example read.
		'''

		read = pysam.AlignedRead()
		read.qname = "chr1_22929880"
		read.seq = "CCCCTAACCCTAACCCTAACCCTAACCCTCGCGATACCCTCAACCAACCCGCCCGCCCG"
		read.flag = 83
		read.rname = 1
		read.pos = 439
		read.mapq = 255
		read.cigar = [(0,59)]
		read.rnext = 1
		read.mpos=151
		read.isize=-348
		read.qual="EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
		read.tags = read.tags + [("XG", "GA")] + [("XM", "..............................Z.Zx........x..zx...Z...Z...Z")] + [("XR", "CT")]
		return read

	def test_N0(self):
		# Shouldn't change methylation_index
		read = self.buildRead()
		self.assertEqual(ignore_first_n_bases(read, [29, 31, 44, 49, 53, 57, ], 0), [29, 31, 44, 49, 53, 57, ])

	def test_N30(self):
		read = self.buildRead()
		# Should remove first element of methylation_index
		self.assertEqual(ignore_first_n_bases(read, [29, 31, 44, 49, 53, 57, ], 30), [31, 44, 49, 53, 57, ])

	def test_NNegative(self):
		read = self.buildRead()
		# Should report an error
		self.assertRaises(ignore_first_n_bases(read, [29, 31, 44, 49, 53, 57, ], -10))

# FIXME: Remove?
if __name__ == '__main__':
    unittest.main()