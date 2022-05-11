
#
# This is for classifying if ps is retro or not
#

import sys
from bisect import bisect

ejunc = sys.argv[1] # exon junction
ginfo = sys.argv[2] # gap inforamtion
blast = sys.argv[3] # blast output, tabular format



##########################
print "Read exon junction info..."
inp = open(ejunc)
inl = inp.readlines()[1:]
J   = {} # {model:junc_list}
for i in inl:
	L = i[:-1].split("\t")
	if L[0] not in J:
		t = []
		if L[-1] != "":
			for j in L[-1].split(","):
				t.append(int(j))
		J[L[0]] = t
	#else:
	#	print "ERR:redun",L[0]

##########################
print "Read gap info..."
G = {} # {"ps pr":{prL:{prR:[psL,psR,prGap]}}}
inp = open(ginfo)
inl = inp.readline()
while inl != "":
	L = inl[:-1].split("\t")
	qs = "%s %s" % (L[0],L[1])
	qL = int(L[2])
	qR = int(L[3])
	sL = int(L[4])
	sR = int(L[5])
	g  = {} # {prGapIdx:[psL,psR]}
	if L[-1] != "":
		for j in L[-1].split(","):
			p = j.split("|")
			p[1] = p[1].split("-")
			# note that everything here is aa-based
			g[int(p[0])] = [int(p[1][0]),int(p[1][1])]
	
	if qs not in G:
		G[qs] = {sL:{sR:[qL,qR,g]}}
	elif sL not in G[qs]:
		G[qs][sL] = {sR:[qL,qR,g]}
	elif sR not in G[qs][sL]:
		G[qs][sL][sR] = [qL,qR,g]
	#else:
	#	print "Redun:",qs,sL,sR
	
	inl = inp.readline()

##########################
print "Process blast output table..."
inp = open(blast)
inl = inp.readline()
B   = {} # {"ps pr":{prL:[prR,psL,psR,ev]}}
eT  = 1e-5
lT  = 10
while inl != "":
	L = inl.split("\t")
	
	qs = "%s %s" % (L[0],L[1])
	ev = float(L[10])
	lg = int(L[3])
	qL = int(L[6])
	qR = int(L[7])
	sL = int(L[8])
	sR = int(L[9])
	
	#if qs == "Os02ps9 Os10g36984.1":
	#	print ">>>",ev,eT,lg,lT,qL,qR
		
	# note the last criteria, translation should be in the same orientation...
	if ev <= eT and lg >= lT and qL < qR:
		if qs not in B:
			B[qs] = {sL:[sR,qL,qR,ev]}
		elif sL not in B[qs]:
			B[qs][sL] = [sR,qL,qR,ev]
		# take the one with lower evalue
		else:
			#print qs, sL
			#print " new:",[sR,qL,qR,ev]
			#print " old:",B[qs][sL]
			if ev < B[qs][sL][3]:
				B[qs][sL] = [sR,qL,qR,ev]
	inl = inp.readline()

#print "--->",B["Os02ps9 Os10g36984.1"]

##########################
print "Process coordinates..."
sDOT = 20 # subjt (prot) distance/overlap threshold
qOT  = 30 # query (ps) overlap threshold (in aa)

sGT  = 10 # subj gap threshold (# of aa in query where subjt is gapped)
          # larger than this will be regarded as PRSENCE of intron
sGT2 = 2  # subj gap threshold, <= this will be regarded as ABSENCE of intron

sEGD = 5  # subj distance threshold between exon junction to gap position
endT = 10 # distance in aa of exon junction from the end of the alignment
          # there are many cases where intron sequences are aligned fairly
          # well at the end and this type of alignments tend to have several
          # fragment, so set at 10. Note that this will increase the false
          # negative rate of retro-ps identification.

R    = {} # {"ps pr":[#_align, 
		  #           [#_w_exonj, 
		  #            #_coincide_w_gap, 
		  #            # w exonj but too far from gap,
		  #            # w exonj and gap but gap too small],
		  #           [# of aln pairs w ej at aln junction
		  #            # with ps insertion >= sGT*3, 
		  #            # with ps insertion <= sGT2*3]]}

oup = open(sys.argv[3]+".ps_class","w")
oup.write("PseudoID\tProtID\tNumFrag\t" + \
		   "FragWExonJ\tExonJInQgap\tFragGapAbsOrTooFar\tFragGapTooSmall\t"+\
		   "ExonJInAlnJ\tExonJinQalngap\tAlnGapTooSmall\t[Dup,Retro]\n")

			 
for i in B:
	print "###################", i
	sLL    = B[i].keys() # subj coordL List
	sLL.sort()
	R[i]   = [0,[],[]]
	n      = i.split(" ")
	try:
                junc   = J[n[1]]
        except KeyError:
                print "key error here because of", n[1]
                continue
	
	F_DUP  = 0 # call duplicate pseudo or not
	F_RET  = 0 # call retro pseudo or not

	print "<------- Check Frag ------->"
	
	# 
	# check if each alignment prot region has exon junction, if so
	#     if it coincide with a gap > qGT, then this is a dup ps (D)
	#     elif there is no gap, this is a retro ps (R)
	#     else this is ambiguous, called fragment (F)
	countA = len(sLL)    #_align
	countE = 0		     #_w_exonj
	countG = 0           # conicide with qualified gap (>= sGT)
	countR = 0           # w exonj but too far from gap
	countS = 0           # w exonj and gap but gap too small (<= sGT2)
	# [D,R,F]
	case1 = [0,0,0]
	for j in sLL:
		sL = j
		sR,qL,qR,e = B[i][sL]
		idxL = bisect(junc,sL)
		idxR = bisect(junc,sR)
		
		print " [sL,sR]=",[sL,sR]
		print " junc:",junc
		
		# with exon_junc
		if idxL != idxR:			
			# G = {"ps pr":{prL:{prR:[psL,psR,prGap]}}}
			gap = G[i][sL][sR][2]
			gkeys = gap.keys()
			gkeys.sort()
			print " gap:",gkeys
					
			# check every exonj in this alignment
			for ej in junc[idxL:idxR]:
				print "  ej:",ej
				# exon_junc not too close to the end of the alignment
				if abs(ej-sL) >= endT and abs(ej-sR) >= endT:
					countE += 1
					# findout the insertion position for ej
					idxG = bisect(gkeys,ej)
					
					print "  gkeys,ej,idxG:",gkeys,ej,idxG
					if len(gkeys) == 0:
						print "  ->no gap but with ej"
						countR += 1
						
					else:
						# flag for those with ej_gap_dist larger than sEGD
						# needed because this would qualified in 5' and 3'
						# searches and be counted twice, which is no good.
						ej_gap_dist = 0	
						
						# need these two as well
						ejInGap  = 0
						gapSmall = 0
						
						# check the gap before iG
						if idxG-1 >= 0:
							g = gkeys[idxG-1]
							print "  ->check 5' (ej,g',gsize):",ej,g+j-1,gap[g][1]-gap[g][0]
	
							# g+j-1: adjust gap position relative to TL start
							if abs(g+j-1-ej) <= sEGD:
								if gap[g][1]-gap[g][0] >= sGT:
							   		print "  -->ej in gap"
								   	ejInGap = 1
								elif gap[g][1]-gap[g][0] <= sGT2:
									print "  -->gap too small:",gap[g][1]-gap[g][0]
									gapSmall = 1
							else:
								print "  -->EGD too large"
								ej_gap_dist = 1
						# check the gap after idxG
						if idxG != len(gkeys):
							g = gkeys[idxG]
							print "  ->check 3' (ej,g',gsize):",ej,g+j-1,gap[g][1]-gap[g][0]  	
	
							# g+j-1: adjust gap position relative to TL start
							if abs(g+j-1-ej) <= sEGD:
								if gap[g][1]-gap[g][0] >= sGT:
									print "  -->ej in gap"
							   		ejInGap = 1
							   	elif gap[g][1]-gap[g][0] <= sGT2:
									print "  -->Gap too small:",gap[g][1]-gap[g][0]
									gapSmall = 1
							else:
								print "  -->EGD too large"
								ej_gap_dist = 1
						
						if ejInGap:
							countG += 1
						if gapSmall:
							countS += 1
						if ej_gap_dist and not (ejInGap or gapSmall):
							countR += 1
				else:
					print "  ->ej too close to the end..."
					pass
		else:
			print " no enclosed ej..."
			pass
	
	print " [A:E:G:R:S]=",[countA,countE,countG,countR,countS]
	R[i][0] = countA
	R[i][1] = [countE,countG,countR,countS]
	
	# Now check pair of alignments
	print "<------- Check Aln Pair ------->"
	
	countP = 0 # aln seg pairs with ej in aln junction
	countG = 0 # the above with ps insertion >= sGT*3
	countS = 0 # the above with ps insertion <= sGT2*3
	print "<>>>> sLL:",sLL
	for j in range(len(sLL)):
		for k in range(j+1,len(sLL)):
			sL1 = sLL[j]
			sL2 = sLL[k]
			sR1,qL1,qR1,e1 = B[i][sL1]
			sR2,qL2,qR2,e2 = B[i][sL2]
			flag_ej  = 0
			flag_dup = 0
			flag_ret = 0

			print " frag1:",sL1,B[i][sL1]
			print " frag2:",sL2,B[i][sL2]
			
			# subject dist/over         query overlap
			if abs(sL2-sR1) <= sDOT and (qL2 > qR1 or qR1-qL2 < qOT):
				print " Q: sL2-sR1=",sL2-sR1
				# check between sL2 and sR1 if there is an exon junction
				print "  sL2,sR1,junc=",sL2,sR1,junc
				
				
				# check sR1
				idxJ = bisect(junc,sR1)
				if idxJ >= 1:
					print "   gapL,idxJ,5':",junc[idxJ-1]
					if abs(sR1-junc[idxJ-1]) < sEGD:
						print "   ej in, qL2-qR1=",qL2-qR1
						flag_ej = 1
						
						# note that q coords are nt-based.
						if qL2-qR1 >= sGT*3:
							print "   Dup"
							flag_dup = 1
							
						elif qL2-qR1 <= sGT2*3:
							print "   Retro"
							flag_ret = 1
							
						else:
							print "   Ambiguous"
					else:
						print "   ej out"
				
				if idxJ != len(junc):
					print "   gapL,idxJ,3':",junc[idxJ]
					if abs(sR1-junc[idxJ]) < sEGD:
						print "   ej in, qL2-qR1=",qL2-qR1
						flag_ej = 1
						
						# note that q coords are nt-based.
						if qL2-qR1 >= sGT*3:
							print "   Dup"
							flag_dup = 1
							
						elif qL2-qR1 <= sGT2*3:
							print "   Retro"
							flag_ret = 1
							
						else:
							print "   Ambiguous"
					else:
						print "   ej out"
						
				# check sL2
				idxJ = bisect(junc,sL2)
				if idxJ >= 1:
					print "   gapR,idxJ,5':",junc[idxJ-1]
					if abs(sL2-junc[idxJ-1]) < sEGD:
						print "   ej in, qL2-qR1=",qL2-qR1
						flag_ej = 1
						
						# note that q coords are nt-based.
						if qL2-qR1 >= sGT*3:
							print "   Dup"
							flag_dup = 1
							
						elif qL2-qR1 <= sGT2*3:
							print "   Retro"
							flag_ret = 1
							
						else:
							print "   Ambiguous"
					else:
						print "   ej out"
				
				if idxJ != len(junc):
					print "   gapR,idxJ,3':",junc[idxJ]
					if abs(sL2-junc[idxJ]) < sEGD:
						print "   ej in, qL2-qR1=",qL2-qR1
						flag_ej = 1
						
						# note that q coords are nt-based.
						if qL2-qR1 >= sGT*3:
							print "   Dup"
							
							flag_dup = 1
						elif qL2-qR1 <= sGT2*3:
							flag_ret = 1
							
							print "   Retro"
						else:
							print "   Ambiguous"
					else:
						print "   ej out"
				
			else:
				print " N: sL2-sR1=",sL2-sR1
				
			if flag_dup and flag_ret:
				print " AMBIGUOUS CASE"
			
			if flag_dup:
				countG += 1
			if flag_ret:
				countS += 1
			if flag_ej:
				countP += 1
			
		
	print " [P,G,S]=",[countP,countG,countS]
	R[i][2] = [countP,countG,countS]
	
	print R[i][2][0],R[i][2][1],R[i][2][2]
	fragS = R[i][1]
	pairS = R[i][2]
	
	F_DUP = F_RET  = 0
	ExonJInQgap        = fragS[1]
	ExonJunQalngap     = pairS[1]
	FragGapAbsOrTooFar = fragS[2]
	
	if ExonJInQgap > 0 or ExonJunQalngap > 0:
		F_DUP = 1
	if FragGapAbsOrTooFar > 0 and ExonJunQalngap == 0:
		F_RET = 1

	#oup.write("PseudoID\tProtID\tNumFrag\t" + \
	#	   "FragWExonJ\tExonJInQgap\tFragGapAbsOrTooFar\tFragGapTooSmall\t"+\
	#	   "ExonJInAlnJ\tExonJinQalngap\tAlnGapTooSmall\t[Dup,Retro]\n")

	oup.write("%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t[%i,%i]\n" % \
					(n[0],n[1],R[i][0],
					 fragS[0],fragS[1],fragS[2],fragS[3],
					 pairS[0],pairS[1],pairS[2],F_DUP,F_RET))

print "Done!"
	

