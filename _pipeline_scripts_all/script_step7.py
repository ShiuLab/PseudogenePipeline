import sys,re

def finditer(p,s):
	p = re.compile(p)    # pattern with one or more "-"
	m = re.finditer(p,s) # iterative find all p
	spans = []
	for i in m:
		spans.append(i.span())
	return spans

# define univeral codes and return a dict with such
def get_nt_code():
	code = {"TTT":"F","TCT":"S","TAT":"Y","TGT":"C",
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
	return code

# return a dict with aa as key, nt code as value
def get_aa_code():
	nt_code = get_nt_code()
	code = {}
	for i in nt_code.keys():
		if not code.has_key(nt_code[i]):
			code[nt_code[i]] = [i]
		else:
			code[nt_code[i]].append(i)
	return code

# Read the pseudogene dna sequence file
inp = open(sys.argv[1])
inl = inp.readlines()
S   = {}
g   = ""
for i in inl:
	if i[0] == ">":
		g = i[1:-1]
	
	if g not in S:
		S[g] = ""
	else:
		S[g] += i[:-1]

# Read the disable_count file
inp = open(sys.argv[2])
inl = inp.readline()

g   = {}
p1  = "[-]+"  # any number of "-"
pcode = get_aa_code()
ncode = get_nt_code()
oup = open(sys.argv[2]+".ps_cds","w")
countAln = 0
countSeq = 0
while inl != "":
	if inl[0] == "#":
		# #25411562 AT3G30216 1-180:49-590 5.6e-66 0 0 1 0
		L = inl[1:-1].split(" ")
		print "\n",inl[:-1]
		ps    = L[1]
		coord = L[2].split(":")[1].split("-")
		seq1  = inp.readline()[:-1]  # prot aln seq 1: plantp
		gap1  = finditer(p1,seq1)   
		seq2  = inp.readline()[:-1]  # prot aln seq 2: intergenic
		ig    = "%s|%s-%s" % (L[1],coord[0],coord[1])
		# Account for cases when the sequences is only a subset of the ones
		# specified in the disable_count file
		if ig not in S:
			print "Seq not found:",ig
			inl = inp.readline()
			continue
		nt2   = S[ig]
		                             # intergenic DNA sequence
		countAln += 1
		
		cdna_mod = ""
		pos  = 0
		fs2  = 0   # in codon frameshift flag
		fs1  = 0   # weird case of frameshift following gap
		lastJ = "" # aa from the previous iteration
		aaN  = 0
		ignored = 0
		for j in seq2:
			print j,pos,
			if j == "-":
				print " gap"
			# frameshift between codon
			elif j == "/":
				if seq1[aaN] != "-":
					cdna_mod += "///"
				print " frameshift-/"
				fs1 = 1
			# frameshift within codon
			elif j == "|":
				if seq1[aaN] != "-":
					cdna_mod += "|||"
				print " frameshift-|"
				fs2 = 1
			elif j in pcode:
				print pcode[j]
				if fs2:
					"""	
					s = ""
					p = 0
					c = 0
					# XxXX or XxxXX or XXxX or XXxxX
					if nt2[pos]+nt2[pos+2:pos+4] in pcode[j]:
						print "  fs-|-1:",pos,nt2[pos]+nt2[pos+2:pos+4]
						s = nt2[pos]+nt2[pos+2:pos+4]
						p = 4
						c += 1
					if nt2[pos]+nt2[pos+3:pos+5] in pcode[j]:
						print "  fs-|-2:",pos,nt2[pos]+nt2[pos+3:pos+5]
						s = nt2[pos]+nt2[pos+3:pos+5]
						p = 5
						c += 1
					if nt2[pos:pos+2]+nt2[pos+3] in pcode[j]:
						print "  fs-|-3:",pos,nt2[pos:pos+2]+nt2[pos+3]
						s = nt2[pos:pos+2]+nt2[pos+3]
						p = 4
						c += 1
					if nt2[pos:pos+2]+nt2[pos+4] in pcode[j]:
						print "  fs-|-4:",pos,nt2[pos:pos+2]+nt2[pos+4]
						s = nt2[pos:pos+2]+nt2[pos+4]
						p = 5
						c += 1					
					if c != 1:
						print "WARNING: multiple alternative cases. FIX IT!"
						print L
						print pos,j,pcode[j]
						print "QUIT!"
						sys.exit(0)
					
					pos += p
					cdna_mod += s
					"""
					# simple ignore this aa
					fs2 = 0
					ignored = 1	
				# Frameshift following gap, ignore the amino acid altogether
				# for example in AT5G54206, the pep alignment:
				#  VANG--/SLVAYGRWFLTNP
				# There is no codon for S following /. Instead there are two
				# nucleotides that can't be assigned.
				#
				# Second case, AT3G30216. There is no gap before / in fact.
				# What happened here is that the codon before is "reused" to
				# generate the codon following. So fs1 flag should be applied
				# generally.
				#
 				# Case 3: 47681287 2|85230-87316 1-272:215-1974 1.7e-15 0 0 0 1
 				# There is a frameshift and somehow the immediate aa that follows
				# was ignored.
				# 
				elif fs1:
					# only search for two position for the aa right after the frame
 					# shift and also constraint that next aa should be the next 3 nt
					
					# simply ignore this aa
					fs1 = 0
					ignored = 1
				else:
					if pos >= len(nt2):
						print "WARNING: nt seq ends but not pep seq. FIX IT!"
						print L
						print seq2[aaN:]
						print "QUIT!"
						sys.exit(0)
						
					while pos < len(nt2):
						# code found
						if nt2[pos:pos+3] in pcode[j]:
							print "  y:",pos,nt2[pos:pos+3]
							if seq1[aaN] != "-":
								cdna_mod += nt2[pos:pos+3]
							pos += 3
							break
						# not found, increment position by 1
						else:
							# check if this codon has weird characters
							try:
								ncode[nt2[pos:pos+3]]
							except KeyError:
								# deal with X
								if j == "X":
									print " y:",pos,nt2[pos:pos+3],"ambiguous"
									if seq1[aaN] != "-":
										cdna_mod += nt2[pos:pos+3]
									pos += 3
									break
								else:	
									print "ERR: unknown codon,", nt2[pos:pos+3]
									print "QUIT!"
									sys.exit(0)
							print "  n:",pos,nt2[pos:pos+3]
							pos += 1
			else:
				print " unknown"
			
			lastJ = j
			aaN   += 1
		
		#if ignored:
		#	print "Frameshift problem: quit!"
		#	sys.exit(0)
		oup.write(">%s_%s|%s\n%s\n" % (L[1],coord[0],coord[1],cdna_mod))
		countSeq += 1		

	inl = inp.readline()

print " Alignment processed:",countAln
print " CDS output         :",countSeq

