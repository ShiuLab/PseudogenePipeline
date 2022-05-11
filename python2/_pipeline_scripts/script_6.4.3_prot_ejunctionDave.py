
#
# This is for finding the exon juctions based on exon cds cood on chr.
#

import sys

print "Process %s..." % sys.argv[1]
inp = open(sys.argv[1])
inl = inp.readlines()
cds = {} # {gene:[L1,R1,L2,R2,...]
for i in inl:
	i = i[:-1].split("\t")
	try:
		g = i[-1].split("pacid=")[1]
		if i[3] == "+":
			cL = int(i[1])
			cR = int(i[2])
		else:
			cL = int(i[2])
			cR = int(i[1])
		
		if g not in cds:
			cds[g] = [cL,cR]
		else:
			cds[g].extend([cL,cR])
			
		# same left, but this one is longer, never happen...
		#else:
		#	print " ERR:",g,cL,cR
		
	except IndexError:
		# there are features present in all.gff3 but not in all.TU.model
                print 1
		pass
	

###
# The way TIGR's gff file is set up, those with reverse ori, each feat has
# L < R so this should be changed. But for different feat of the same gene
# they are arranged from 5'(relative to TL start) to 3'. So don't need to take
# care or sorting and orientation issues.

print "Consolidate coordinates and determine junctions..."
ckeys = cds.keys()
ckeys.sort()
oup = open(sys.argv[1]+".exon_junc","w")
oup.write("Model\tcoordChr\tcoordCDS\taaJunct\n")
c = 0
for i in ckeys:
	if c % 1e4 == 0:
		print " %i x 10k" % (c/1e4)
	c += 1
	
	# gene coord (begin) at TL start
	t  = []
	c5 = cds[i][0] # 5' coordinate (relative to TL start)
	# reverse ori
	if cds[i][0] > cds[i][1]:
		for j in cds[i]:
			t.append(abs(j-c5-1))
	# forward
	else:
		for j in cds[i]:
			t.append(j-c5+1)
	
	# process junction
	# cds coordinate
	newt = []
	# junction aa coordinate
	junc = []
	for j in range(0,len(t),2):
		if j == 0:
			newL = t[j]
			newR = t[j+1]
			newt.extend([newL,newR])
			#print j,t[j],t[j+1],
		else:
			newL = newt[j-1]+1
			newR = newL + (t[j+1]-t[j])
			newt.extend([newL,newR])
			#print j,t[j],t[j+1],newL,newR,
			
		# figure out exon junction in aa coord
		if j+2 != len(t):
			# note that this will point to the aa right before the intron
			aa = newR/3
			junc.append(aa)
			#print aa
	
	t    = str(t)[1:-1].split(", ")
	newt = str(newt)[1:-1].split(", ")
	junc = str(junc)[1:-1].split(", ")
	oup.write("%s\t%s\t%s\t%s\n" % (i,
									",".join(t),
									",".join(newt),
									",".join(junc)))

print "Done!"

