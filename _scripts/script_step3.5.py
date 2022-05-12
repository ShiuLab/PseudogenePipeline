################################################################################
#
# 10/14/10 By Shinhan
# This script does two things:
# 1) Based on the pseudogene file from step3, get the coordinates of pseudogenes
#    on the subject sequences.
# 2) Generate a pair list with prot_id, genome_seq_coord_as_id (based on 1).
#
#

import sys

inp = open(sys.argv[1])
oup1= open(sys.argv[1]+".subj_coord","w")
oup2= open(sys.argv[1]+".pairs","w")

inl = inp.readline()
while inl != "":
	L = inl.split("\t")		# genome_seq_id, prot_id, genome_region, prot_region
	gr= eval(L[2])
	#print gr,
	gr= [gr[0][0],gr[-1][1]]# Get the regions out
	#print gr
	
	oup1.write("%s\t%s\t%s\n" % (L[0],gr[0],gr[1]))
	oup2.write("%s\t%s|%s-%s\n" % (L[1],L[0],gr[0],gr[1]))
	
	inl = inp.readline()


