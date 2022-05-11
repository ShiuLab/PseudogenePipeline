
#
# This script Take the cds generated from the step 7 script and do two things:
# 1) Get rid of all in-frame stop codons and unknown codons, assuming frame 0 
#    for all sequences. The modified seq can then be used for rate calculation.
# 2) Partition sequence into before the 1st stop and after the 1st stop for
#    testing if there is selection constraint on part of the sequence.
#

import sys


ncode = {"TTT":"F","TCT":"S","TAT":"Y","TGT":"C",
		 "TTC":"F","TCC":"S","TAC":"Y","TGC":"C",
		 "TTA":"L","TCA":"S","TAA":"*","TGA":"*",
		 "TTG":"L","TCG":"S","TAG":"*","TGG":"W",			 
		 "CTT":"L","CCT":"P","CAT":"H","CGT":"R",
		 "CTC":"L","CCC":"P","CAC":"H","CGC":"R",
		 "CTA":"L","CCA":"P","CAA":"Q","CGA":"R",
		 "CTG":"L","CCG":"P","CAG":"Q","CGG":"R",
		 "ATT":"I","ACT":"T","AAT":"N","AGT":"S",
		 "ATC":"I","ACC":"T","AAC":"N","AGC":"S",
		 "ATA":"I","ACA":"T","AAA":"K","AGA":"R",
		 "ATG":"M","ACG":"T","AAG":"K","AGG":"R",
		 "GTT":"V","GCT":"A","GAT":"D","GGT":"G",
		 "GTC":"V","GCC":"A","GAC":"D","GGC":"G",
		 "GTA":"V","GCA":"A","GAA":"E","GGA":"G",
		 "GTG":"V","GCG":"A","GAG":"E","GGG":"G",			
		 "NNN":"X",""   :"/",""   :"|"}

inp = open(sys.argv[1])
inl = inp.readline()
oupA = open(sys.argv[1]+".rid_allstop","w")
oup1 = open(sys.argv[1]+".1st_half","w")
oup2 = open(sys.argv[1]+".2nd_half","w")

c = 0
while inl != "":
    # assuming seq in the next line and end with \n
	if inl[0] == ">":
		nam = inl[1:-1]
		#print nam
		if c % 1e3 == 0:
			print " %ik" % (c/1e3)
		c += 1
		seq = inp.readline()[:-1]
		firstStop = 0
		seqA = ""
		seq1 = ""
		seq2 = ""
		for j in range(0,len(seq),3):
			if seq[j:j+3] in ["TAA","TGA","TAG"]:
				#print " STOP:",j,seq[j:j+3]
				if firstStop == 0:
					firstStop = 1
			elif seq[j:j+3] not in ncode or ncode[seq[j:j+3]] == "X":
				print " UNKN:",j,seq[j:j+3]
			else:
				seqA += seq[j:j+3]
				if firstStop == 0:
					seq1 += seq[j:j+3]
				else:
					seq2 += seq[j:j+3]
		
		oupA.write(">%s\n%s\n" % (nam,seqA))
		if seq2 != "":
			oup1.write(">%s\n%s\n" % (nam,seq1))
			oup2.write(">%s\n%s\n" % (nam,seq2))
	
	inl = inp.readline()

print "Total % sequences processed" % c