#!/usr/bin/python
import re
import random
from numpy import array
import argparse
import sys
import csv
import pysam

## Extract comethylation signal for reads overlapping multiple CpGs
## SAM input file is the output of Bismark (no pre-processing required)
## CpG overlap is as determined by Bismark output, i.e. only positions with Z or z in the methylation string are considered CpGs

## This program is Copyright (C) 2012, Peter Hickey (hickey@wehi.edu.au)

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.


# Order of fields in Bismark SAM file ###
 #####################################################################################
 ### (1) QNAME                                                                     ###
 ### (2) FLAG                                                                      ### 
 ### (3) RNAME                                                                     ###
 ### (4) POS                                                                       ###
 ### (5) MAPQ                                                                      ###
 ### (6) CIGAR                                                                     ###
 ### (7) RNEXT                                                                     ###
 ### (8) PNEXT                                                                     ###
 ### (9) TLEN                                                                      ###
 ### (10) SEQ                                                                      ###
 ### (11) QUAL                                                                     ###
 ### (12) NM-tag (edit distance to reference)                                      ###
 ### (13) XX-tag (base-by-base mismatches to the reference, not including indels)  ###
 ### (14) XM-tag (methylation call string)                                         ###
 ### (15) XR-tag (read conversion state for the alignment)                         ###
 ### (16) XG-tag (genome conversion state for the alignment)                       ###
######################################################################################

# Explanation of methylation string (XM) 
  #################################################################
  ### . for bases not involving cytosines                       ###
  ### X for methylated C in CHG context (was protected)         ###
  ### x for not methylated C in CHG context (was converted)     ###
  ### H for methylated C in CHH context (was protected)         ###
  ### h for not methylated C in CHH context (was converted)     ###
  ### Z for methylated C in CpG context (was protected)         ###
  ### z for not methylated C in CpG context (was converted)     ###
  #################################################################

# This program is a re-write of SAM2MS.py in order to utilise the pysam library.

