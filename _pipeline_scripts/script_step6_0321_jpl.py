# Debug: fixed some indentation errors
# The output file name can change according to input argument of flank


import sys,re



def space(s):
	L = s.split(" ")
	tmp = []
	for i in L:
		if i != "":
			tmp.append(i)
	return tmp

def finditer(p,s):
	p = re.compile(p)	# pattern with one or more "-"
	m = re.finditer(p,s) # iterative find all p
	spans = []
	for i in m:
		spans.append(i.span())
	return spans
	
# From: http://www.d.umn.edu/~mhampton/m5233s7.html
def nw(seq1, seq2, matrix, gap = -8, match = 1, mismatch = -1):
	
	# initialize scoring and 'arrow' matrices to 0
	S = [[0 for x in range(len(seq2)+1)] for y in range(len(seq1)+1)]
	P = [[0 for x in range(len(seq2)+1)] for y in range(len(seq1)+1)]
	
	# initialize borders
	# for P (arrows), use 2 for diagonal, -1 for horizontal, and 1 for 
	# vertical moves.
	# I have tried to consistently use i for rows (vertical positions) in the 
	# score and pointer tables, and j for columns (horizontal positions).
	for i in range(len(seq1)+1):
		S[i][0] = gap*i
		P[i][0] = 1
	for j in range(len(seq2)+1):
		S[0][j] = gap*j
		P[0][j] = -1
	# fill with scores
	for i in range(1,len(seq1)+1):
		for j in range(1,len(seq2)+1):
			L1 = seq1[i-1]
			L2 = seq2[j-1]
			# treat stop or frameshift as gap
			# | is the substitute for \ because of it is also use as escape
			if L1 in ["*","/","\\","|"]:
				L1 = "-"
			if L2 in ["*","/","\\","|"]:
				L2 = "-"
			Sdg = S[i-1][j-1] + matrix[L1][L2]  # diagnoal
			Shz = S[i][j-1] + gap			   # horizontal
			Sup = S[i-1][j] + gap			   # up
			# TempS is list of the three S and their P
			TempS = [[Sdg,2],[Shz,-1],[Sup,1]]
			# Now we keep the highest score, and the associated direction (pointer)
			S[i][j], P[i][j] = max(TempS)
	# backtrace from end
	[i,j] = [len(seq1),len(seq2)]
	align1 = ''
	align2 = ''
	# does NOT traceback on multiple optimal paths.
	while [i,j] != [0,0]:
		if P[i][j] == 2:
			align1 = align1 + seq1[i-1]
			align2 = align2 + seq2[j-1]
			i = i - 1
			j = j - 1
		elif P[i][j] == -1:
			align1 = align1 + '-'
			align2 = align2 + seq2[j-1]
			j = j - 1
		else:
			align1 = align1 + seq1[i-1]
			align2 = align2 + '-'
			i = i - 1
	# the alignments have been created backwards, so we need to reverse them:
	align1 = align1[::-1]
	align2 = align2[::-1]
	# print out alignment
	return align1,align2

def read_blosum50(blosum):
	print "Read BLOSUM50 matrix..."
	# read blosum50 matrix
	inp = open(blosum)
	inl = inp.readlines()
	aa  = space(inl[0][:-1])
	b50 = {}
	inl = inl[1:]
	for i in range(len(inl)):
		#print [inl[i]]
		s = space(inl[i][:-1])
		for j in range(len(aa)):
			#print aa[j],
			if aa[j] not in b50:
				b50[aa[j]] = {aa[i]:int(s[j])}
				#print "new:",b50[aa[j]]
			else:
				b50[aa[j]][aa[i]] = int(s[j])
				#print "app:",b50[aa[j]]
	return b50

def read_swout(swout):
	print "Read the sw.out file..."
	# As of 2/16,2012, change the following so it will only allow one alignment
	# to be reported
	inp = open(swout)
	inl = inp.readline()
	pair = []
	stat = []
	seq  = []
	sw   = {} # {pair:[[stat1,seq1,stop1,stop2,frame1,frame2],...]}
				# stop1 and frame1: no alignment necessary
				# stop2 and frame2: aligned because they are too close to gaps

	while inl != "":
		if inl[0] == "#":
			pair = inl[1:-1]
			stat = []
		elif inl[0] == ">":
			stat = inl[1:-1].split(" ")
			seq  = []
		else:
			# rid of the escape character and replace it with "2"
			L = inl[:-1].split("\\")
			L = "|".join(L)
			seq.append(L)
		inl = inp.readline()

		# Modified to allow only one entry, 2/16/2012
		if stat != [] and len(seq) == 2:
			if pair not in sw:
				sw[pair] = [[stat,seq,0,0,0,0]]
			# to prevent redundant entries, 8/12, 07'
			#elif [stat,seq,0,0,0,0] not in sw[pair]:
			#	sw[pair].append([stat,seq,0,0,0,0])
	
	return sw

def calculate(segScore):
    flank = int(round(((segScore+6.32)/6.63)*(flankpro/100.0)))
    return flank

def compare_seq(swout,sw,b50):
	# Now compare the sequences
	print "Compare sequences:"
	p1   = "[-]+"  # any number of "-"
	p2   = "[*]"   # a single stop
	p3   = "[/|]"  # single "/" or "\"
	oup = open(swout+".disable_count_2"+".flank"+".%i" %flankpro,"w")
	c = 0
	swkl = len(sw.keys())
	
	print " total: %i alignments" % swkl
	for i in sw:
		if c % 1000 == 0:
			print " %i K" % (c/1000)
		c += 1
		#print "",i,sw[i]
		#if c == 100:
			#sys.exit(0)
		countSeg = 0

		# sw = {pair:[[stat1,seq1,stop1,stop2,frame1,frame2],...]}
		for j in sw[i]: # There should be only one element here... So kind of stupid
									  # to iterate... but I don't want to change it.
			# get gap list in seq1
                        p11  = "-" + p1 # any number of "-" more than 2
			gap1 = finditer(p11,j[1][0])
			for k in range(len(gap1)):
				gap1[k] = [gap1[k][0],gap1[k][1],"-"]
                                #print gap1, j[1][0]
			gapCoord = [0] # [begin,gap1L,gap1R,gap2L,gap2R....,end]
			for k in gap1:
				gapCoord.append(k[0]-1)
				gapCoord.append(k[1])
			gapCoord.append(len(j[1][0])-1)   
			
			segScore = {} # {"gap_end-gap_begin":score}
			for k in range(0,len(gapCoord),2):
				sL = gapCoord[k]      # segment L
				sR = gapCoord[k+1]    # segment R
				sS1= j[1][0][sL:sR+1] # segment sequence from query
				sS2= j[1][1][sL:sR+1] # segment sequence from subject
                                
				segScore["%i-%i" % (sL,sR+1)] = 0
				for m in range(len(sS1)):
					s1 = sS1[m]
					s2 = sS2[m]
					if "*" not in [s1,s2] and "/" not in [s1,s2] and "|" not in [s1,s2]\
                                            and "\\" not in [s1,s2]:
						#print s1,s2,b50[s1][s2]
						segScore["%i-%i" % (sL,sR+1)] += b50[s1][s2]
					#else:
					#	print s1,s2,"bad"

			#print segScore
                        gap2 = [[-2, 0,"-"]] + gap1 + [[len(j[1][0]),len(j[1][0])+2,"-"]]
                        
                        segScoreM = [] #[ [ begin, end, "score"], ...]]
                        for x in segScore:
                            contex = x.split("-")
                            segScoreM.append([int(contex[0]),int(contex[1]), str(segScore[x])])
                        segScoreM.sort()
                        
                        #print j[1][0]
                        #print j[1][1]
                        #print segScoreM
			#print gap2
			#
			# JINPENG WILL TAKE OVER FROM HERE !!!!!!!!!!!
			#
			
			# get stop list in seq2
			stop2 = finditer(p2,j[1][1])
			for k in range(len(stop2)):
				stop2[k] = [stop2[k][0],stop2[k][1],"*"]
			#print  stop2,j[1][1]
			
			# get framshift list in seq2
			fram2 = finditer(p3,j[1][1])
			for k in range(len(fram2)):
				fram2[k] = [fram2[k][0],fram2[k][1],"f"]		
			#print  fram2
			
			c1 = stop2 + fram2
			c1.sort()
			#print c1
			# combine check stop and frameshift
			for k in c1:
				#print k
				stop_ok1 = 1
				stop_ok2 = 0
				fram_ok1 = 1
				fram_ok2 = 0
                                flank_1  = 0
                                add1 = add2 = add3 = add4 = 0
				for m in gap2: 
					# gap has to be larger than the threshold (in this case it ignore the 1aa gap)
					if m[1]-m[0] > gaplen:
						# x-----xxxx
						# xxxx*xxxxx
						#print k,m                        
                                                                                   #### 012345        012345 
						if k[0] >= m[0] and k[1] <= m[1]:  #### xxx--x (3,5)  xxx--x (3,5)
							if k[2] == "*":            #### xxx*xx (3,4)  xxxx*x (4,5) should be k[1] <= m[1] or k[0] < m[1]?
								stop_ok1 = 0
							else:
								fram_ok1 = 0
							break
                                            ######################################
                                                else:
                                                                               # To difine the flank_1 as a function of flank and segment score
                                                        #flank_1  = 0       
                                                        #for y in segScoreM:
                                                            #print exon
                                                            #exL  = int(y[0])
                                                            #exR  = int(y[1])
                                                            #score = int(y[2])
                                                            #print exL, exR, segScore
                                                            
                                                            #if k[0] >= exL and k[1] <= exR: 
                                                                
                                                                #flank_1 = calculate(score)
                                                      
                                                                
                                                                
                                                         ###########################################
                                                        #      12345
                                                        # x----xxxxxx
                                                        # xxxxxxx*xxx
                                                        # minus 1 because gap end is exclusive
                                                                
                                                        if k[0] >= m[1] and k[0] < m[1]+flank:   ### the original code is k[0] >= m[0] and k[0] < m[1]-1+flank
                                                            if m[0] == -2 and m[1] ==0:             ### to check if the disable is close to the very begining
                                                                if k[2] == "*":                 ### might calculate alignment score to classify the disable as type2 or not
                                                                    stop_ok2 = 1
                                                                else:
                                                                    fram_ok2 = 1 
                                                                      
                                                            else:
                                                                seq1 = j[1][0][m[1]:k[1]+2]          # do ungapped global alignment
                                                                seq2 = j[1][1][m[0]:k[1]+2]
                                                                #print "1:",seq1,seq2
                                                                aln1,aln2 = nw(seq1,seq2,b50)
                                                                #print aln1
                                                                #print aln2
                                                                # evaluate if the stop/frameshift is sitting in a gap
                                                                # stop/frameshift position in the alignment: k[0]-m[1]
                                                                sf = k[0]-m[0]
                                                                if aln1[sf] == "-" and \
                                                                         (aln1[sf-1] == "-" or aln1[sf+1] == "-"):
                                                                        if k[2] == "*":
                                                                                stop_ok1 = 0
                                                                        else:
                                                                                fram_ok1 = 0
                                                                        break
                                                                else:
                                                                        if aln2[sf] == "*":
                                                                                stop_ok2 = 1
                                                                        else:
                                                                                fram_ok2 = 1
                                                                
                                                                #print aln2[sf],stop_ok2,fram_ok2
                                                                
                                                        #  12345
                                                        # xxxxxx---xx
                                                        # xxx*xxxxxxx
                                                        elif k[0] < m[0] and k[0]+flank >= m[0]:
                                                            if m[0] == len(j[1][0]) and m[1] == len(j[1][0]) + 2: ### to check if the disable is close to the very end
                                                                if k[2] == "*":                 ### might calculate alignment score to classify the disable as type2 or not
                                                                    stop_ok2 = 1
                                                                else:
                                                                    fram_ok2 = 1
                                                                # do ungapped alignment
                                                                # make sure left coord >= 0
                                                            else:
                                                                sL = max(k[0]-2,0)
                                                                seq1 = j[1][0][sL:m[0]]
                                                                seq2 = j[1][1][sL:m[1]]
                                                                #print "2:",seq1,seq2
                                                                aln1,aln2 = nw(seq1,seq2,b50)
                                                                #print aln1
                                                                #print aln2
                                                                sf = k[0]
                                                                if k[0]-2 >= 0:
                                                                        sf = 2
                                                                        
                                                                if aln1[sf] == "-" and \
                                                                         (aln1[sf-1] == "-" or aln1[sf+1] == "-"):
                                                                        if k[2] == "*":
                                                                                stop_ok1 = 0
                                                                        else:
                                                                                fram_ok1 = 0
                                                                        break
                                                                else:
                                                                        if aln2[sf] == "*":
                                                                                stop_ok2 = 1
                                                                        else:
                                                                                fram_ok2 = 1
				
				if k[2] == "*" and stop_ok1:
					if stop_ok2 == 0:
						add1 = 1
						sw[i][countSeg][2] += 1 # stop that are relatively far from gap
					else:
						add2 = 1
						sw[i][countSeg][3] += 1 # stop qualified after realignment
				elif k[2] != "*" and fram_ok1:
					if fram_ok2 == 0:
						add3 = 1
						sw[i][countSeg][4] += 1 # frame shift that are far from gap
					else:
						add4 = 1
						sw[i][countSeg][5] += 1 # frame shift qualified after realignment
                                #print k,add1, add2, add3, add4, flank_1
                                
                                
			# generate output
			# #plantp interg range STOP1=X STOP2=X FS1=X FS2=X
			# seq1
			# seq2
			oup.write("#%s %s %s %i %i %i %i\n" % \
							(i,sw[i][countSeg][0][-1],sw[i][countSeg][0][3],sw[i][countSeg][2],
							 sw[i][countSeg][3],sw[i][countSeg][4],sw[i][countSeg][5]))
			oup.write("%s\n%s\n\n" % (sw[i][countSeg][1][0],sw[i][countSeg][1][1]))
		
			#print sw[i][c]
			countSeg += 1

	oup.close()

#-------------------------------------------------------------------------------
if len(sys.argv) < 4:
	print "Usage: script_step6.py    BLOSUM     sw.out     flankthreshold"
	sys.exit(0)

# Set flank threshold, need stop and frameshift to be be this number of amino 
# acids away from gap with length . Default is 10% of the exon.
if len(sys.argv) == 4:
	try:
		flankpro = int(sys.argv[3])
	except ValueError:
		print "Value for flanking threshold is not an interger! Quit!"
		sys.exit(0)
        #if flankpro < 1 or flankpro > 100:
            #print "Value for flank proportion threshold should between 1 and 100"
            #sys.exit(0)
else:
	flankpro = 5

# need gap to be larger than this to be considered
gaplen = 1	 
flank = flankpro
blosum50 = sys.argv[1]
swout    = sys.argv[2]
	
b50 = read_blosum50(blosum50)
sw  = read_swout(swout)
compare_seq(swout,sw,b50)
	
	
print "Done!"
