#!/usr/bin/python
import re
import random
from numpy import array
import argparse
import sys
import csv
import pysam

# Firstly need to correct the FLAG (to properly encode the strand information) and fix the read-names using correct_Bismark_PE_SAM.py
## TODO: Need to consider how to handle paired end reads that sequenced the same stretch of DNA twice. It is effectively two measurements of the same data point - perhaps throw out data point if the measurements are consistent, otherwise retain.


FLAGS = {}
for read in f:
    flag = read.flag
    if not flag in FLAGS:
        FLAGS[flag] = 1
    else:
        FLAGS[flag] +=1

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

