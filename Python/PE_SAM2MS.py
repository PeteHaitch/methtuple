import re
import random
from numpy import array
import argparse
import sys
import csv
import pysam

f = pysam.Samfile('sorted_bismark.bam', 'rb')

for read in f:
    if read.tid == read.rnext and read.is_proper_pair and !read.is_duplicate and !PAIR ISNT DUPLICATE  # Make sure both reads mapped to the same chromosome

for read in f: 
    if read.is_proper_pair and read.is_read1: 
      pos = f.tell() 
      try: 
        mate = f.mate(read) 
      except ValueError: 
        # Invalid mate (usually post-filtered) 
        continue 
      finally: 
        f.seek(pos)
        
      print read.qname,'read1' 
      print mate.qname,'read2' 
      #print 

