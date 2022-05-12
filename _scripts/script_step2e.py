import sys, pickle

if len(sys.argv) != 3:
	print("Usage: script_step3e.py blast gapL")
	print("       blast: filtered blast tabular output between prot and genome")
	print("QUIT!")
	sys.exit(0)


inp = open(sys.argv[1])
gapL= int(sys.argv[2])
inl = inp.readline()
D   = {}   # {intergSeqID:    (i)
           #	{intergL: (j)
           #       	 00      1      2     3
           #		[[intergR,protL,protR,Eval,protID],...]
           #    }
           # }
c   = 0
cS  = 0    # count # of subjects (intergenic seq)
while inl != "":
	if c % 1e4 == 0:
		print(" %i x 10k" % (c/1e4))
	c += 1
	
	L  = inl.split("\t")
	qL,qR,sL,sR = int(L[6]),int(L[7]),int(L[8]),int(L[9])
	eV = float(L[-2])
	ori = 1
	if sL > sR:    # reverse orientation
		qL,qR,sL,sR = qR,qL,sR,sL
		ori = 0
	if L[1] not in D:
		cS += 1
		D[L[1]] = {sL:[[sR,qL,qR,eV,L[0]]]}
	elif sL not in D[L[1]]:
		D[L[1]][sL] = [[sR,qL,qR,eV,L[0]]]
	else:
		D[L[1]][sL].append([sR,qL,qR,eV,L[0]])
	
	inl = inp.readline()
	
	
print("Number of intergenic sequences with hit:",cS)

def add(qC,sC,eV,oL,pr,v1,v2,v3,v4,v5,v6,v7):
	qC.append([v1,v2])
	sC.append([v3,v4])
	eV.append(v5)
	oL += v6
	pr.append(v7)

oup = open("%s_G%s.PE" % (sys.argv[1],sys.argv[2]),"w")
cP  = 0 # count # of pseudoexons
cO  = 0 # count # of pseudoexons with fragments of opposite ori
c   = 0
D2  = {}
for i in D:                   # interg seq ID
	if c % 1e3 == 0:
		print(" %i k" % (c/1e3))
	c += 1
	intL = list(D[i].keys()) # intergSeqL
	intL.sort()
	qC = []
	sC = []
	eV = []
	oL = ""
	pr = []
	for j in intL:
		#print i,k,j,D[i][j]
		for n in D[i][j]:
			if n[1] > n[2]:
				ori = "-"
			else:
				ori = "+"
			if qC == []:
				qC.append([n[1],n[2]])
				sC.append([j   ,n[0]])
				eV.append(n[3])
				oL += ori
				pr.append(n[4])
			# check if enclosed within previous entry or there is gap
			else:
				L = sC[-1] # last pair of coordinates
				#print L,[j,n[0]],
				# overlap
				if j < L[1]:
					# enclosed
					if n[0] <= L[1] or j == L[0] and n[0] >= L[1]:
						print("1",eV[-1],n[3])
						# smaller e value in new one, seach upward
						r = 0 # replace flag
						while 1:
							# either value not qualifying or not enclosed
							if eV[-1] <= n[3]:
								if r == 0:
									pass
								# enclosed and evalue is not as good
								elif n[0] <= L[1] or j == L[0] and n[0] >= L[1]:
									#print "enclosed"
									pass
								# not enclosed and e value better than the one
								# just deleted
								else:
									print("added1")
									add(qC,sC,eV,oL,pr,
										n[1],n[2],j,n[0],n[3],ori,n[4])
								break
							else:
								r = 1
								if n[0] <= L[1] or j == L[0] and n[0] >= L[1]:
									print("del",L)
									del(qC[-1])
									del(sC[-1])
									del(eV[-1])
									del(pr[-1])
									oL = oL[:-1]
									try:			
										L = sC[-1]
									except IndexError:
										add(qC,sC,eV,oL,pr,
											n[1],n[2],j,n[0],n[3],ori,n[4])
										break										
								else:
									print("added2")
									add(qC,sC,eV,oL,pr,
										n[1],n[2],j,n[0],n[3],ori,n[4])
									break									
					# not enclosed
					else:
						print("2")
						add(qC,sC,eV,oL,pr,
							n[1],n[2],j,n[0],n[3],ori,n[4])
				# does not overlap but right next to each other
				elif j-L[1] == 1:
					print("3")
					add(qC,sC,eV,oL,pr,
						n[1],n[2],j,n[0],n[3],ori,n[4])
				# does not overlap but gap <= gapL
				elif j - L[1] <= gapL:
					print("4")
					add(qC,sC,eV,oL,pr,
						n[1],n[2],j,n[0],n[3],ori,n[4])					
				# does not overlap and with gap > gapL
				else: 
					print("5")
					cP += 1
					if oL.find("+") != -1 and oL.find("-") != -1:
						cO += 1
					oup.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (i,
											str(sC),str(qC),str(eV),oL,str(pr)))
					# reset lists
					qC = [[n[1],n[2]]]
					sC = [[j   ,n[0]]]
					eV = [n[3]]
					oL = ori
					pr = [n[4]]
	if oL.find("+") != -1 and oL.find("-") != -1:
		cO += 1
	cP+=1
	# output last set									
	# intergID, protID, subjtCoord, queryCoord, evalues
	print("5")
	oup.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (i,
											str(sC),str(qC),str(eV),oL,str(pr)))


oup.close()
print("Number of pseudoexons:",cP)
print("Number of PE with opp ori:",cO)
print("Done!")
