
import os, sys
from bisect import bisect

class gff_util:

	#
	# List the contents of the first three columns
	#
	def checkfields(self,gff,tokens):
		
		# convert column indices to ints
		tokens = tokens.split(",")
		t      = []
		for i in tokens:
			t.append(int(i))
		
		# read gff
		inp = open(gff)
		inl = inp.readline()
		C = {} # content: C = {col_index:{unique_value:count}}
		c = 0
		print "Read",gff
		while inl != "":
			# ignore header
			if inl[0] == "#":
				inl = inp.readline()
				continue

			if c % 1e4 == 0:
				print " %i x 10k" % (c/1e4)
			c += 1
			L = inl.split("\t")
			try:
				for j in t:
					if j not in C:
						C[j] = {}
					elif L[j] not in C[j]:
						C[j][L[j]] = 1
					else:
						C[j][L[j]] +=1
			except IndexError:
				print "IndexError:",L
			inl = inp.readline()
		
		ckeys = C.keys()
		ckeys.sort()
		for i in ckeys:
			print "Col",i
			ukeys = C[i].keys()
			ukeys.sort()
			for j in ukeys:
				print "",j,C[i][j]
		
		print "Done!"
	
	#
	# Get the features specified in the columns
	#
	def getfeat(self,gff,col,feat):
		
		col = col.split(",")		
		feat = feat.split(",")
		
		inp = open(gff)
		oup = open(gff+"_col%s_%s" % (".".join(col),".".join(feat)),"w")

		# convert col index from str to int
		for i in range(len(col)):
			col[i] = int(col[i])

		inl = inp.readline()
		c = 0
		cQ= 0
		while inl != "":
			if inl[0] == "#":
				inl = inp.readline()
				continue
			if c % 1e4 == 0:
				print " %i x 10k" % (c/1e4)
			c += 1
			L = inl.strip().split("\t")
			for j in col:
				try:
					if L[j] in feat:
						oup.write(inl)
						cQ += 1
						break
				except IndexError:
					print "ERR:",j,[inl]
			inl = inp.readline()
							
		
		print "Total %i, qualified %i" % (c,cQ)
		print "Done!"
	
	#
	# Get the coordinate for each "name". Output: Chr, L, R, featname.
	#
	# @param target: cds, exon, or intron
	# @param brief : 0, all coords; 1, beginning and end. In the same direction
	#                as it is in the GFF.
	# @param ID    : specific kind of identifier for getting sequence name.
	#
	def getfeat2(self,gff,target="",brief=0,ID=""):
		print "MAY NOT WORK FOR ALL GFFs. CHECK CAREFULLY!!"
		
		inp = open(gff)
		inl = inp.readline()
		D   = {} # {name:[chr,ori,[L:R]]}
		while inl != "":
			L = inl.strip().split("\t")
			if L != [""]:
				c = L[0]          # chr
				f = L[2].lower()  # feature
				if f == target:
					l = int(L[3])   # left coord
					r = int(L[4])   # right coord
					o = L[6]        # orientation
					
					# Set sequence name
					n = ""
					if ID != "":
						L[-1] = L[-1].split(";")
						for j in L[-1]:
							if j.find(ID) != -1:
								n = j.split("=")[-1]
								break 
					if n == "":
						n = L[-1].split("\"")[1]# name
					if n not in D:	
						D[n] = [c,o,[l,r]]
					else:
						D[n][2].extend([l,r])
			inl = inp.readline()
		
		oup = open("%s.%s_B%i" % (gff,target,brief),"w")
		Dk = D.keys()
		Dk.sort()
		for i in Dk:
			c = D[i][0] # chr
			o = D[i][1] # ori
			d = D[i][2] # coordinate
			d.sort()
			if o == "-":
				d.reverse()
			
			oup.write("%s\t" % c)
			if brief:
				oup.write("%i\t%i\t" % (d[0],d[-1]))
			else:
				oup.write("%s\t" % str(d))
			oup.write("%s\n" % i)
			
		print "%i names" % len(Dk)
		print "Done!"

	#
	# Output: a four column file with genename,chr,L,R
	# 
	def starttostop(self,gff):
		
		# Read GFF into a dictionary
		inp = open(gff)
		inl = inp.readline()
		GFF = {}			# {featurename:{featuretype:L or R}}
		CHR = {}			# {featurename:chr}
		FEAT= ["start_codon","stop_codon"]
		while inl != "":
			L = inl.strip().split("\t")
			if L[2] in FEAT:
				C  = L[0]	# chr
				FT = L[2]	# feature type
				cL = L[3]
				cR = L[4]				
				if L[6] == "-":	# orientation
					[cL,cR] = [cR,cL]
				FN = L[8]	# feature name
				CHR[FN] = C
				if FN not in GFF:
					if FT == FEAT[0]:
						GFF[FN] = {FT:cL}
					else:
						GFF[FN] = {FT:cR}
				elif FT not in GFF[FN]:
					if FT == FEAT[0]:
						GFF[FN][FT] = cL
					else:
						GFF[FN][FT] = cR
				else:
					print "Redun:",FN,FT,cL,cR
			inl = inp.readline()
		
		oup = open(gff+".starttostop","w")
		for i in GFF:
			try:
				oup.write("%s\t%s\t%s\t%s\n" % \
								(i,CHR[i],GFF[i][FEAT[0]],GFF[i][FEAT[1]]))
			except KeyError:
				print "KeyERR:",i,GFF[i]
		print "Done!"

	#
	# Assume that the feature name field contain "protein_id" In the format:
	#
	#   gene_id "fgenesh1_pg.C_scaffold_1000003"; transcript_id "72561"; protein_id "72561"; exon_id "72561.2";
	#
	def joincds(self,gff):
		inp = open(gff)
		inl = inp.readline()
		CDS = {} # {protein_id:[L1,R1,L2,R2,....]}
		CHR = {} # {protein_id:[chr,ori]}
		while inl != "":
			L = inl.strip().split("\t")
			if L[2] == "CDS":
				cL  = int(L[3])
				cR  = int(L[4])
				PID = ""
				for j in L[8].split(";"):
					if j.find("protein_id") != -1:
						PID = j.split("\"")[1]
				if PID == "":
					print "ERR, no protein id:",L[8]
				else:
					CHR[PID] = [L[0],L[6]]
					if PID not in CDS:
						CDS[PID] = [cL,cR]
					else:
						CDS[PID].extend([cL,cR])
			
			inl = inp.readline()
		
		oup = open(gff+".joincds","w")
		for i in CDS:
			coords = CDS[i]
			CDS[i].sort()
			if CHR[i][1] == "-":
				CDS[i].reverse()
				
			oup.write("%s\t%s\t%i\t%i\n" % \
									(i,CHR[i][0],CDS[i][0],CDS[i][-1]))
		
		print "Done!"
		
	#
	# Get exon cds coordinates for each mRNA.
	#
	def get_exon_cds(self,gff,ID):
		inp = open(gff)
		oup = open(gff+".exoncds","w")
		inl = inp.readlines()
		#ECD = {} # {mRNA:{cds_L:[cds_R,phase]}}
		mid = "" # mRNA unique identifier
		countM = 0
		countC = 0
		ccds   = [] # cds coord
		phas   = [] # exon phase
		for i in inl:
			i    = i.strip().split("\t")
			ch   = i[0]
			feat = i[2]
			attr = i[8].split(";")
			fL   = i[3]
			fR   = i[4]
			orie = i[6]
			if orie == "-":
				fL,fR = fR,fL
			if feat == "mRNA":
				if mid != "":
					log = ""
					#print ccds,phas
					"""
					if phas[0] != 0:
						if orie == "+":
							ccds[0] = str(int(ccds[0])+(3-phas[0]))
						else:
							ccds[0] = str(int(ccds[0])-(3-phas[0]))
						log += "E[1]_phase:%i " % phas[0]
					if phas[-1] != 0:
						if orie == "+":
							ccds[-1] = str(int(ccds[-1]) - phas[-1])
						else:
							ccds[-1] = str(int(ccds[-1]) + phas[-1])
						log += "E[-1]_phase:%i" % phas[-1]
					"""
					oup.write("%s\t%s\t%s\n" % (",".join(ccds),mid))
					ccds = []
					phas = []
				for j in attr:
					if j.find(ID) != -1:
						mid = j.split("%s=" % ID)[-1]
					break
				oup.write("%s\t" % ch)
				#print "mRNA:",mid
				countM += 1
			
			elif feat == "CDS":
				pha  = int(i[7])
				ccds.extend([fL,fR])
				phas.append(pha)
				countC += 1
				
		oup.write("%s\t%s\n" % (",".join(ccds),mid))
	
		inp.close()
		oup.close()
		print "%i mRNA, %i exon CDS" % (countM,countC)
		print "Done!"

	#
	# Get intron based on an exon GFF file
	#
	def get_intron(self,gff):
		inp = open(gff)
		inl = inp.readlines()
		EX  = {} # Exon dict = {parent:[ch,ori,[cL1,cR1,...]]}
		print "Read gff..."
		for i in inl:
			i   = i.strip().split("\t")
			ch  = i[0]
			cL  = int(i[3])
			cR  = int(i[4])
			ori = i[6]
			ID  = i[-1].split("Parent=")[-1]
			if ID not in EX:
				EX[ID] = [ch,ori,[cL,cR]]
			else:
				EX[ID][2].extend([cL,cR])
		
		print "Generate output..."
		oup = open(gff+".intron_coord","w")
		loc = EX.keys()
		loc.sort()
		for i in loc:
			# exon coords
			[ch,ori,EC] = EX[i]
			EC.sort()
			if EC[1] == "-":
				EC.reverse()
			#print EC
			for j in range(1,len(EC)-1,2):
				iL = EC[j]
				iR = EC[j+1]
				if ori == "+":
					iL += 1
					iR -= 1
				else:
					iL -= 1
					iR += 1
				oup.write("%s\t%i\t%i\t%s\n" % (ch,iL,iR,i))
		oup.close()
	
	#
	# Take GFF with gene info only, generate an output with
	#  chr intergenicL intergenicR 5'gene|5'ori|3'gene|3'ori
	# @size    Chromosome sizes generated by FastaManager.get_sizes.
	#
	def get_intergenic(self,gff,size):
		inp = open(gff)
		inl = inp.readlines()
		GN  = {} # Gene dict = {ch:{geneL:[geneR,ori,ID]}}
		print "Read gff..."
		for i in inl:
                        if not i.startswith("#"):
                                i   = i.strip().split("\t")
                                ch  = i[0]
                                gL  = int(i[3])
                                gR  = int(i[4])
                                ori = i[6]
                                ID  = i[-1].split(";")[0].split("=")[1]
                                if ch not in GN:
                                        GN[ch] = {gL:[gR,ori,ID]}
                                elif gL not in GN[ch]:
                                        GN[ch][gL] = [gR,ori,ID]
                                # Same gL, longer gene, take longer one
                                elif GN[ch][gL][0] < gR:
                                        GN[ch][gL] = [gR,ori,ID]
				
		print "Read size file..."
		inl = open(size).readlines()
		SZ = {} # {chr:size}
		for i in inl:
                        if not i.startswith("#"):
                                [ch,sz] = i.split("\t")
                                sz = int(sz)
                                SZ[ch] = sz
		
		print "Get intergenic coord..."
		oup = open(gff+".intergenic_coord","w")
		ck = GN.keys() # chr keys
		ck.sort()
		cO = 0 # count overlap
		for i in ck:
			print i
			# Only entries in chr specified by the size file will be processed.
			if i not in SZ:
				print " chr not in size file, ignored"
				continue
				
			glk = GN[i].keys() # gene left coord keys
			glk.sort()
			#print len(glk)
			# Write the first intergenic region
			gL1 = glk[0]
			[gR1,ori1,ID1] = GN[i][gL1]
			oup.write("%s\t%i\t%i\t%s|%s|%s|%s\n" % \
																				(i,1,gL1-1,"5'end","NA",ID1,ori1))
			for j in range(0,len(glk)-1):
				#print "",glk[j]
				gL1 = glk[j]
				gL2 = glk[j+1]
				[gR1,ori1,ID1] = GN[i][gL1]
				[gR2,ori2,ID2] = GN[i][gL2]
				# check if these overlap
				if gR1 > gL2:
					print "Overlap, ignored:",[ID1,gL1,gR1],[ID2,gL2,gR2]
					cO += 1
				else:
					oup.write("%s\t%i\t%i\t%s|%s|%s|%s\n" % \
																				(i,gR1+1,gL2-1,ID1,ori1,ID2,ori2))
			# Write the last intergenic region
			gL1 = glk[-1]
			[gR1,ori1,ID1] = GN[i][gL1]
			oup.write("%s\t%i\t%i\t%s|%s|%s|%s\n" % \
																				(i,gR1+1,SZ[i],ID1,ori1,"3'end","NA"))

			
		print "%i gene pairs overlap" % cO
		oup.close()
		
	
	#
	# Get upstream regions of defined length, as long as it does not overlap with
	# flanking gene.
	# @gff     GFF containing only gene
	# @region  (up)stream, (dn)stream, or (both)
	# @rlen    Length of the regions to get
	def get_upstream(self,gff,region,rlen):
		pass
		
	#
	# Converting GFF format file to BED. Note that the description field may
	# be parsed incorrectly, please check.
	#
	def gff_to_bed(self,gff):
		print "WARNING:"
		print "  THE DESCRIPTION FIELD MAY NOT PARSE CORRECTLY, PLEASE CHECK!!"
		inp = open(gff)
		oup = open(gff+".bed","w")
		inl = inp.readlines()
		for i in inl:
			if i[0] != "#":
				L = i.strip().split("\t")
				C = L[0]	# chr
				cL= L[3]	# chr left
				cR= L[4]	# chr right
				S = L[5]	# score
				try:	
					G = L[8].split(";")[0].split("=")[1]
					oup.write("%s\t%s\t%s\t%s\t%s\n" % (C,cL,cR,G,S))
				except IndexError:
					print "Problem:",L[8]

		print "Done!"
		
	# 
	# Compare two GFF files. Will determine if any features in GFF2 overlap with
	# that in GFF1. The overlap is determine by the prox parameter.
	# 
	# @prox Overlap parameter. If it is zero (default), the features in GFF2 have
	#       to overlap with a feature in GFF1 to be included. If prox > 0, then
	#       GFF1 features that do not overlap but within that a distance <= prox
	#       will be included.
	#
	def compare_gffs(self,gffs,prox):
		
		if gffs.find(",") == -1:
			print "Need two GFF file names separated by ',', QUIT!"
			sys.exit(0)
			
		gffs = gffs.split(",")
		G1 = self.gff_to_dict(gffs[0])
		G2 = self.gff_to_dict(gffs[1])
		oup = open("%s_vs_%s" % (gffs[0],gffs[1]),"w")
		
		# Compare G2 to G1
		# for each sequence in G1
		for i in G1:
			# get a sorted list of cL from G2
			if i in G2:
				cL2 = G2[i].keys()
				cL2.sort()
				
				# for each cL in G1 in sequence i
				for j in G1[i]:
					cR1  = G1[i][j].keys()
					cR1.sort()
					
					idxL = bisect(cL2,j)
					# for each cR that have the same cL
					for k in cR1:
						idxR  = bisect(cL2,k+prox)
						ofeat = [] 		# overlapping features
						#print "G1:",j,k
						
						# just in case the G1 feature span multiple G2 features. In this
						# case, idxL < idxR.
						for m in range(idxL,idxR+1):
							#print " m:",m
							# Check overlap with the G2 feat before idxL. Only do so if there
							# can be an G2L before m. Also, only check this in the first
							# found if multiple m's are possible. Otherwise, it will be
							# redundant.
							if idxL != 0 and m == idxL:
								cL2x= cL2[m-1]
								cR2 = G2[i][cL2x].keys()
								cR2.sort()
								for n in cR2:
									#print " G2-1:",cL2x,n
									if n+prox >= k:
										# gene name, L, R, ori
										ofeat.append("%s[%i-%i,%s]" % \
																	(G2[i][cL2x][n][0],cL2x,n,G2[i][cL2x][n][2]))
						
							# Check overlap with the G2 feat right at idxL
							if idxL < len(cL2):
								cL2x= cL2[m]
								cR2 = G2[i][cL2x].keys()
								cR2.sort()
								for n in cR2:
									#print " G2:",cL2x,n
									if cL2x-prox < k:
										ofeat.append("%s[%i-%i,%s]" % \
																	(G2[i][cL2x][n][0],cL2x,n,G2[i][cL2x][n][2]))
								
							# Check overlap with the G2 feat right after idxL.
							# This is not necessary because the idxR has considered each G1
							# cR + prox. And in the for loop with m, all possible G2 L that
							# can be smaller than G1 cR + prox are considered.

						# output overlapping features for this G1 cL, cR combination
						oup.write("%s\t%s\n" % (G1[i][j][k][1].strip(), ",".join(ofeat)))
		print "Done!"

	#
	# G = {Chr:{L:{R:gff_line_list}}}
	#
	def gff_to_dict(self,gff):
		
		print "%s to dict..." % gff
		inp = open(gff)
		inl = inp.readlines()
		G   = {}
		for i in inl:
			if i[0] != "#":
				L = i.strip().split("\t")
				C = L[0]				# chr, or sequence name
				cL= int(L[3])		# L coord
				cR= int(L[4])		# R coord
				O = L[6]				# orientation
				try:						# gene name if parsable
					g = L[-1].split(";")[0].split("=")[1]
				except:
					print "ERR last field, use all:",L[-1]
					g = L[-1]
				
				if C not in G:
					G[C] = {cL:{cR:[g,i,O]}}
				elif cL not in G[C]:
					G[C][cL] = {cR:[g,i,O]}
				elif cR not in G[C][cL]:
					G[C][cL][cR] = [g,i,O]
				else:
					print "Completely redun:",C,cL,cR,L
		return G
		
		
	def help(self):
		print " -f function"
		print "    checkfields - list field contents. NEED: -gff. OPT: col"
		print "    compare_gffs - NEED: gff, OPT: prox"
		print "    getfeat - get features in columns. NEED: gff,col,feat"
		print "    getfeat2- get the coordinates for each name. NEED: gff"
		print "       target, brief, OPT: ID"
		print "    get_intron - parse an Exon GFF to get intron coord, NEED: gff"
		print "    get_intergenic - parse an gene GFF to get intergenic coord,"
		print "       NEED: gff"
		print "    get_exon_cds - Get exoncds coordinates, NEED: gff, ID"
		print "    gff_to_bed - NEED: gff"
		print "    joincds - Assume attribute format. Look at code. NEED: gff"
		print "    starttostop - get start to stop codon coord for each gene."
		print "       NEED: gff"
		print " -gff   gff file. For compare_gffs, pass two gffs separted by ','."
		print " -brief 0: all coordinates, 1: only start and stop"
		print " -col   column indices in gff to be checked. Default 0,1,2"
		print " -feat  target features, separated by ',', Default ''"
		print " -ID    Target ID in the attribute field"
		print " -prox  Proximity of features in gff2 from gff1. Default, overlap"
		print "        (0), or some positive values"
		print " -size  Chromosome sizes"
		print " -target the features, exon, cds, or intron"
		print ""
		sys.exit(0)
	
#-------------------------------------------------------------------------------
if __name__ == '__main__':

	gutil = gff_util()
	f = gff = feat = target = brief = ID = size = ""
	col = "0,1,2"
	prox = 0

	for i in range(1,len(sys.argv),2):
		if sys.argv[i] == "-f":
			f = sys.argv[i+1]
		elif sys.argv[i] == "-gff":
			gff = sys.argv[i+1]
		elif sys.argv[i] == "-col":
			col = sys.argv[i+1]
		elif sys.argv[i] == "-feat":
			feat = sys.argv[i+1]
		elif sys.argv[i] == "-size":
			size = sys.argv[i+1]
		elif sys.argv[i] == "-target":
			target = sys.argv[i+1]
		elif sys.argv[i] == "-ID":
			ID = sys.argv[i+1]
		elif sys.argv[i] == "-brief":
			brief = int(sys.argv[i+1])
		elif sys.argv[i] == "-prox":
			prox = int(sys.argv[i+1])
		else:
			print "Unknown flag:",sys.argv[i]
			gutil.help()
	
	if f == "checkfields":
		if gff == "":
			print "\nNeed GFF file\n"
			gutil.help()
		gutil.checkfields(gff,col)
			
	elif f == "compare_gffs":
		if gff == "":
			print "\nNeed two GFF file names separated by ','\n"
			gutil.help()
		gutil.compare_gffs(gff,prox)
			
	elif f == "getfeat":
		if "" in [gff,col,feat]:
			print "\nNeed GFF file, target col and features\n"
			gutil.help()
		gutil.getfeat(gff,col,feat)	
		
	elif f == "getfeat2":
		if "" in [gff,target,brief,ID]:
			print "\nNeed GFF file, target feat and brief flag\n"
			gutil.help()
		gutil.getfeat2(gff,target,brief,ID)

	elif f == "get_exon_cds":
		if "" in [gff,ID]:
			print "\nNeed GFF file and target ID in attribute\n"
			gutil.help()
		gutil.get_exon_cds(gff,ID)	

	elif f == "get_intergenic":
		if "" in [gff,size]:
			print "\nNeed GFF and chr size files.\n"
			gutil.help()
		gutil.get_intergenic(gff,size)	

	elif f == "get_intron":
		if "" in [gff]:
			print "\nNeed GFF file \n"
			gutil.help()
		gutil.get_intron(gff)	
		
	elif f == "gff_to_bed":
		if "" in [gff]:
			print "\nNeed GFF file\n"
			gutil.help()
		gutil.gff_to_bed(gff)
	
	elif f == "starttostop":
		if "" in [gff]:
			print "\nNeed GFF file\n"
			gutil.help()
		gutil.starttostop(gff)	
		
	elif f == "joincds":
		if "" in [gff]:
			print "\nNeed GFF file\n"
			gutil.help()
		gutil.joincds(gff)	
		
			
	else:
		print "\nUnknown function...\n"
		gutil.help()
	
	
