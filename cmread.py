#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="parses the output of cross match, cross match must be run with the flags -alignment and -tags")
parser.add_argument("input", help="cross_match output file" )
parser.add_argument("out", help="outputfile of this prgram" )
parser.add_argument('-d', action="store_true", default=False)
args = parser.parse_args()
DEBUG=args.d

import glob
import os
import sys
import re
import itertools 
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import pysam 
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO



def readCMline(line):
	line = line.split()
	# drop the tag that says ALIGNMENT
	line = line[1:]

	# check if it is a complement
	complement = False
	if("C" in line):
		complement = True

	# tells me to drop the line because it is a subset of another
	if("*" in line):
		return(None)

	# check line length
	assert len(line) == 12 + complement
	
	# move the extenshion of the sequence past the aln to match the regualr locaiton
	if(complement):
		extend = line.pop(10)
		line.append(extend)
	# a a Not complement tag to the line if NC
	if(not complement):
		line.insert(8, "NC")
	
	line = re.sub('\(|\)', '', "\t".join(line))
	return( line.split() )

#
# 
#
def readAln(lines, idx):
	header = readCMline(lines[idx])
	# ['17625', '1.81', '0.25', '0.56', 'chr1:223822113-224042658', '139281', '158695', '61851', 'C', 'AC270130.1', '117267', '97912', '81891']
	# move off of header and onto aln
	idx += 2

	# start, end, sequence alignments 
	fasta1 = ""
	fasta2 = ""
	start1 = 0
	start2 =0
	end1 =0
	end2 =0

	first = True
	while(True):
		seq1 = ("_" + lines[idx]).split()
		seq2 = ("_" + lines[idx+2]).split()
		if("Transitions" in seq1[0]):
			break
		if(first):
			first = False
			start1 = int(seq1[2])
			start2 = int(seq2[2])
		# update end
		end1 = int(seq1[4])
		end2 = int(seq2[4])
		# update seqs
		fasta1 += seq1[3]
		fasta2 += seq2[3]
		idx += 4
	
	extend = int( header[-1] ) * "-"
	fasta1= "-"*start2 + fasta1 + extend
	fasta2= "-"*start2 + fasta2 + extend

	rtn = (header,  [start1, end1, start2, end2, fasta1, fasta2] )
	return(rtn)

def read(file):
	lines = open(file).readlines()
	idx = 0
	# get past header
	while(True):
		line = lines[idx]
		idx += 1
		if( "Maximal single base matches" in line):
			break
	aln = ()
	while(idx < len(lines) ):
		line = lines[idx].split()
		# if this is true then we are in a alignment section
		#if(len(line) > 0 and line[0].isdigit() and line[1] != "matching" ):
		if(len(line)> 0 and line[0] == "ALIGNMENT"):
			aln = readAln(lines, idx)
			break
		idx += 1
	rec1 = SeqIO.SeqRecord(Seq(aln[1][4]),id=aln[0][4], description=" ".join(aln[0]))
	rec2 = SeqIO.SeqRecord(Seq(aln[1][5]),id=aln[0][9], description=" ".join(aln[0]))
	#SeqIO.write([rec1], "A."+args.out, "fasta")
	#SeqIO.write([rec2], "B."+args.out, "fasta")
	SeqIO.write([rec1,rec2], args.out, "fasta")



read(args.input)


