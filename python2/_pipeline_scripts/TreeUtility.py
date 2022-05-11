
##
# Program: TreeUltility
# Desc.  : Deal with some typical operation on trees
# Created: July 27, 2002
#
# 09/20,02
#  Implement rename()
# 11/15,02
#  For rename, implement codes so if distance information is missing, it will
#  still process the tree file.
##

import sys, FileUtility, types, os, string, time

class tree_utility:

	def __init__(self):
		pass


	##
	# Take a line of seq_id and break them into lines. Require some kind of
	# recognizable header or this won't work.
	# 
	# @infile  Input with one line. Normally a tree file is converted to emf by
	#	  Treeview first. The emf is displayed in Illustrator then the ids
	#	  are copied. Instead of 
	##
	def transpose(self,infile,prefix):
		pass

	##
	# This function is implemetned to simply replace the name within the tree
	# file. It is also designed to rewrite the tree in New Hampshire X format
	# as described in the "ATV" website (www.genetics.wustl.edu/eddy.atv) if
	# asked.
	#
	# @param tree  the New Hampshire style tree, such as those from  Neighbor
	# @param name  the name file in the [new][old] format, if the NHX option
	#		  is used, the new names HAVE TO have species-specific header
	#		  as the first two character. At this point, two styles are
	#		  possible:
	#		   a.XX_xxxxxxx
	#		   b.XXxxxxx
	# @param nhx   0: just replace names (default), 1: NHX tree
	# @param nflag 0: default, new name in index 0; 1: new name in index 1
	# @param call  internal method call, tree is a string, name is a dict
	##
	def rename(self,tree,name,nhx=0,nflag=0,call=0):

		# put names into ndict, old as key, new as value
		print "Read name file..."
		
		if call:
			tline = tree
			ndict = name
		else:
			inp = open(name,"r")
			inline = inp.readline()
			ndict = {}
			while inline != "":
				llist = inline.split("\t")
	
				# check for new line and carriage returns
				if llist[1][-1] == "\n":
					llist[1] = llist[1][:-1]
					if llist[1][-1] == "\r":
						llist[1] = llist[1][:-1]
			
				try:
					if nflag == 0:
						new = 0
						old = 1
					else:
						new = 1
						old = 0
						
					if ndict.has_key(llist[old]):
						print "Redundant old name:'%s'" % llist[old]
						print "QUIT WITHOUT FINISHING!"
						sys.exit(0)
					else:
						ndict[llist[old]] = llist[new]
						
				except IndexError:
					print "Name file format incorrect"
					print llist
					print "QUIT WITHOUT FINISHING!"		 
					sys.exit(0)
	
				inline = inp.readline()

			# the record for each name is bracketed in three ways:
			#  (A:dA,
			#  ,A:dA)
			#  ,A:dA,
			# the seq name is always in front of an ":"
	
			# convert tree file into a single long string
			inp = open(tree,"r")
			inline = inp.readline()
			tline = ""
			legal = ["(", ")", ",", ":", ".", "-"]
			while inline != "":
	
				tline = tline + inline[:-1]
	
				# in case that this is a DOS-generated file, with a carriage return
				if tline[-1] not in legal and type(tline[-1]) != "int":
					tline =  tline[:-1]	 
				
				inline = inp.readline()
			#print tline

		# scan tline for qualified components
		print "Scan for names and replace them in tree file..."
		c = 0
		mline = ""
		name  = ""
		full_name = 0

		name_count = 0
		name_max   = 0
		undefined = []
		while c < len(tline):
			# left bracket "("
			if tline[c] == "(" and tline[c+1] != "(":
				#print tline[c:]
				mline = mline +"("
				# find right bracket ","
				c = c +1
				while c < len(tline):
					if tline[c] != ",":
						name  = name+tline[c]
						c = c +1
					else:
						full_name = 1
						break
					
			# left bracket ","
			elif tline[c] == "," and tline[c+1] != "(":
				#print tline[c:]				
				mline = mline +","

				# find right bracket ")" or ","
				c = c +1
				while c < len(tline):

					if tline[c] != "," and tline[c] != ")":
						name  = name + tline[c]
						c = c+ 1
					else:
						full_name = 1
						break
			else:
				mline = mline +tline[c]

			# process name here and reset variable
			if full_name:

				found = 0
				nlist = name.split(":")
				if not ndict.has_key(nlist[0]):
					undefined.append(nlist[0])
					name = nlist[0]
				else:
					found = 1
					name = ndict[nlist[0]]
				
				dist = ""
				if len(nlist) > 1:
					dist = nlist[1]
				
					
				if not nhx:
					if dist == "":
						mline = mline + "%s" % name
					else:
						mline = mline + "%s:%s" % (name,dist)
				else:		   
					species = ""
					if name[2] != "_":
						species = name[:2]
					else:
						species = name[:name.find("_")]
						name	= name[name.find("_")+1:]
				
					# add full name and right_end to mline here
					if dist == "":
						mline = mline +"%s[&&NHX:S=%s]" % (name,species)
					else:
						mline = mline +"%s:%s[&&NHX:S=%s]" % (name,dist,species)
				
				# reset stuff
				name = ""
				full_name = 0
				name_count = name_count +1

			# only increment if full name is not obtained yet
			else:   
				c = c +1

		# generate output, At one point I didn't add ";" for whatever reason.
		# now I added it in.
		if call:
			oup = open("internal_call.tre","w")
		else:
			oup = open(tree+".mod","w")
		for i in mline:
			if i == "(" or i == ")" or i == ",":
				oup.write(i+"\n")
			else:
				oup.write(i)
		if mline[-1] != ";":
			oup.write(";")

		print "Total %i names processed, %i not defined" % \
			  (name_count,len(undefined))
		for i in undefined:
			print "",i
		print "Done!!"
		
	
	##
	# Encounter a strange situation where the bootstrap values for Clustal
	# generated tree are embeded in "[]". This function is to fix this so tree
	# explorer can read this.
	##	
	def move_boot(self,tree):

		inp = open(tree,"r")
		inl = inp.readline()
		oup = open(tree+".mod","w")
		found = 0
		while inl != "":
			if inl.find("[") != -1:
				found = 1
				boot = int(inl[inl.find("[")+1:inl.find("]")])
				oup.write("%i%s%s" % (boot,
									  inl[:inl.find("[")],
									  inl[inl.find("]")+1:]))
			else:
				oup.write(inl)
			inl = inp.readline()
		
		if found:
			print "File converted"
		else:
			print "No bootstrap value delimiter found"
	
	
	def treelist(self,tree):
		
		inp   = open(tree,"r")
		tree  = self.rmlb(inp.readline())
		
		if tree[-1] == ";":
			tree = tree[:-1]
		
		tlist = []
		tlist = self.parse_subtree(tree,tlist)
		print tlist
	
	def parse_subtree(self,tree,tlist):
		print "Tree:",tree 
			
		subIdx= 0
		countB = 0
		for i in tree:
			subIdx += 1
			if i == "(":
				countB += 1
			elif i == ")":
				countB -= 1
				if countB == 0:
					break		

		if subIdx == len(tree) and len(tree) > 1 and tree[0] == "(":			
			tree = tree[1:-1]
			if tree.find("(") == -1:
				tlist.extend(tree.split(","))
			else:
				self.parse_subtree(tree,tlist)
		elif subIdx == 1:
			tlist.extend(tree)
		elif tree[0] != "(":
			list1 = self.parse_subtree(tree[tree.find(",")+1:],[])
			tlist.extend([tree[:tree.find(",")],list1])
		else:
			list1 = self.parse_subtree(tree[:subIdx],[])
			list2 = self.parse_subtree(tree[subIdx+1:],[])
			tlist.extend([list1,list2])
			
		print "List:",tlist
		return tlist

	
	##
	# Convert tree to one line of string
	##
	def one_line(self,tree,dos=0,call=0):
		inp = open(tree,"r")	
		inl = inp.readline()
		tline = ""
		while inl != "":
			tline += self.rmlb(inl)
			inl = inp.readline()
			
		if call:
			return tline
		else:
			oup = open(tree+".mod","w")
			if dos:
				oup.write(tline+"\r\n")		
			else:
				oup.write(tline+"\n")
		
	
	##
	# Divide the bootstrap value by the passed factor in tree
	##
	def divide_boot(self,tree,factor):
		
		# first simplify the tree
		self.one_line(tree)
		
		# bootstrap value is always between ")" and ":"
		bvalue = ""
		inp = open(tree+".mod","r")
		inl = inp.readline()[:-1]
		inp.close()
		
		idx_rb = 0
		idx_co = 0
		outstr = ""
		got_b  = 0
		special= 0
		#print inl
		for i in range(len(inl)):
			outstr += inl[i]
			
			if inl[i] == ")":
				idx_rb = i
				#print "idx_rb:",i
			elif inl[i] == ":":
				idx_co = i
				#print "idx_co:",i
				#print inl[:idx_co]
			
			if inl[i] == ",":
				#print inl[:i+1]
				idx_rb = 0
				idx_co = 0
				got_b  = 0
				#print "reset"
			elif idx_rb > idx_co:
				special = 1
				#print "special"
			elif got_b and not special:
				pass
				#print "pass"
			elif (idx_rb > 0 and idx_rb < idx_co) or special:
				#print idx_rb,idx_co
				bvalue = inl[idx_rb+1:idx_co]
				
				try:
					bvalue = float(bvalue)
				except TypeError:
					print "bvalue format problem:",bvalue
				except ValueError:
					print "Value Error!"
					bvalue = 0
				# bvalue is modified
				bvalue = bvalue/float(factor)
				outstr = "%s%i%s" % (outstr[:outstr.rfind(")")+1],
									 bvalue,outstr[-1])
				got_b   = 1
				special = 0
				
		oup = open(tree+".mod","w")
		oup.write(outstr+"\n")
	
	#
	# Need to simplify tree first.
	#
	def get_subtree(self,tree,taxa,outname=""):
		
		inp = open(tree)
		lns = inp.readlines()
		
		if len(lns) > 1:
			print "Convert tree to oneline format..."
			lns = self.one_line(tree,0,1)
		else:
			lns = lns[0]
		
		if lns.find(":") != -1:
			print "Simplify tree..."
			lns = self.simplify(lns,"(",1,1)
		
		
		tree = lns
		#print "Full:",tree
		if tree.find(taxa) != -1:
			# The taxa has to be either "(taxa," or ",taxa)". This searching
			# scheme prevents multiple identical string to be misidentified.
			if tree.find("(%s," % taxa) != -1:
				L = tree.find("(%s," % taxa) + 1
			elif tree.find(",%s)" % taxa) != -1:
				L = tree.find(",%s)" % taxa) + 1
			else:
				print "Taxa string not found:",taxa
				sys.exit(0)
				
			R = L + len(taxa)
						
			# left bracket, search right
			if tree[L-1] == "(":
				sub  = tree[R+1:] # notice this, "," is not in sub
				# a simple situation with one taxa in sister group
				if sub[0] != "(":
					sub = sub[:sub.find(")")]
					sub = "(%s,%s);" % (taxa,sub)
				else:
					countI = 0 # count iterations
					countB = 0 # count brackets
					while countI < len(sub):
						if sub[countI] == "(":
							countB += 1
						elif sub[countI] == ")":
							countB -= 1
						if countB == 0 and countI != 0:
							break
						countI += 1
					sub = "(%s,%s);" % (taxa,sub[:countI+1])				
			# right bracket, search left
			else:
				sub  = tree[:L-1] # notice this, "," is not in sub
				if sub[-1] != ")":
					sub = sub[sub.rfind("(")+1:]
					sub = "(%s,%s);" % (sub,taxa)
				else:
					countI = -1 # count iterations
					countB = 0 # count brackets
					while countI*(-1)-1 < len(sub):
						if sub[countI] == "(":
							countB -= 1
						elif sub[countI] == ")":
							countB += 1
						if countB == 0 and countI != -1:
							break
						countI -= 1
					sub = "(%s,%s);" % (sub[countI:],taxa)				
			
			#print sub
		else:
			print "Specified taxa not found, QUIT!"
			sys.exit(0)
		
		if outname != "":
			oup = open(outname,"w")
			oup.write("%s\n" % sub)
		return sub

	#
	# Get taxa names as defined by the beginning and ending taxas
	#
	def get_taxa(self,tree,taxa1,taxa2):
		
		inp = open(tree)
		lns = inp.readlines()
		if len(lns) > 1:
			print "Convert tree to oneline format..."
			lns = self.one_line(tree,0,1)
		else:
			lns = lns[0]
			
		if lns.find(":") != -1:
			print "Simplify tree..."
			lns = self.simplify(lns,"(",1,1)
		
		if taxa1 == "":
			L = 0	
		elif lns.find(taxa1) != -1:
			L = lns.find(taxa1)
		else:
			print "%s not found! Quit!" % taxa1
			sys.exit(0)
		
		if taxa2 == "":
			R = len(lns)
		elif lns.find(taxa2) != -1:
			R = lns.find(taxa2) + len(taxa2)
		else:
			print "%s not found! Quit!" % taxa2
			sys.exit(0)
		
		# Get taxa names
		T = lns[L:R].split(",")
		#print T
		T2= []
		for i in T:
			
			if i.find("(") != -1:
				i = i.split("(")
				for j in i:
					if j != "":
						i = j
						break
			elif i.find(")") != -1:
				i = i.split(")")
				for j in i:
					if j != "":
						i = j
						break
			#print i
			T2.append(i)
	
		if taxa1 != "":
			oup = open("%s.%s_%s" % (tree,taxa1,taxa2),"w")
		else:
			oup = open("%s_taxa" % tree,"w")
		oup.write("\n".join(T2))
		
		inp.close()
		oup.close()
		print "Done!"

	#
	# Get subtrees with the organizatin (xxxx,each_name);
	# Also flag the tree if they overlap with others.
	#
	def batch_sub(self,tree,name):
	
		oup   = open(name+".subtrees","w")
		name  = futil.file_to_list(name)
		
		# name as key, a list as value with [subtree][redun_flag][redun_to]
		tdict = {}
		for i in name:
			sub = self.get_subtree(tree,i)
			if tdict.has_key(i):
				print "Redundant name:",i
			else:
				tdict[i] = [sub[1:-2],0,""] # won't have ; and outermost brac
	
		# do an all against all search, flag the entry if it is subtree to other
		# flag: 0, no overlap; 1, equal_to; 2, subtree_of
		K = tdict.keys()
		K.sort()
		for i in range(len(K)):
			if tdict[K[i]][1] == 2:
				continue
				
			for j in range(i,len(K)):
				if tdict[K[j]][1] == 2:
					continue
				
				if K[i] != K[j]:
					if tdict[K[i]][0] == tdict[K[j]][0]:
						tdict[K[j]][1] = 1
						tdict[K[j]][2] = K[i]
					elif tdict[K[i]][0].find(tdict[K[j]][0]) != -1:
						tdict[K[j]][1] = 2
						tdict[K[j]][2] = K[i]
					elif tdict[K[j]][0].find(tdict[K[i]][0]) != -1:
						tdict[K[i]][1] = 2
						tdict[K[i]][2] = K[j]
		
		for i in tdict:
			oup.write("%s\t%i\t%s\t%s\n"%(i,tdict[i][1],tdict[i][2],tdict[i][0]))
		
		print "Done!"


	#
	# Replace subtrees (such as those in the batch_sub output) with names
	# 
	# @param tree should be in oneline format.
	# @param name [new_name][subtree]. The subtree are without the outermost
	#								  brackets.
	#
	def replace_sub(self,tree,name):

		inp = open(tree,"r")
		tre = inp.readlines()
		if len(tre) > 1:
			print "Convert tree to oneline format first. QUIT!"
			sys.exit(0)
		
		tre   = tre[0]	
		ndict = futil.file_to_dict(name)
		for i in ndict:
			idx = tre.find(ndict[i])
			if  idx != -1:
				tre = "%s%s%s" % (tre[:idx-1],i,tre[idx+len(ndict[i])+1:])
			
			if tre.find(ndict[i]) != -1:
				print "Non-unique subtree?",ndict[i]
		
		oup = open(tree+".mod","w")
		oup.write(tre)
		
		print "Done!"
	
	def mega_to_nexus(self,tree):
		
		inp = open(tree,"r")
		inl = inp.readline()
		ndict = {}
		tre   = ""
		while inl != "":
			inl = self.rmlb(inl)
			if inl.find("[") != -1:
				ndict[inl[inl.find("[")+1:inl.find("]")]] = inl[inl.find("#")+1:]
			elif inl.find("#MEGATree=") != -1:
				tre = inl[inl.find("=")+1:]
			inl = inp.readline()
		
		self.rename(tre,ndict,call=1)
		print "Done!"
	
	#
	# @param tree      a file (when call=0) or a string (when call=1)
	# @param delimiter "(" (default) or "["
	# @param remove    remove ending semicolon or not, default no (0)
	#
	def simplify(self,tree,delimiter="(",remove=1,call=1):
		
		if call:
			inl = tree
		else:
			inp = open(tree,"r")	
			inl = inp.readlines()
			inl = string.joinfields(inl.split("\n"),"")
			inl = string.joinfields(inl.split("\r"),"")
		
		ostr = ""
		write = 1
		for i in inl:
			if i != ":":
				if write:
					if delimiter != "(":
						if i == "(":
							i = "["
						elif i == ")":
							i = "]"
					ostr += i
				else:
					if i == "," or i == ")":
						if i == ")" and delimiter != "(":
							i = "]"
						ostr += i
						write = 1
			else:
				write = 0
		
		if remove and ostr[-1] == ";":
			ostr = ostr[:-1]
		
		if call:
			return ostr
		else:
			oup = open(tree+".simplify","w")
			oup.write("%s\n" % ostr)

	#
	# Add '' around taxa names in tree (t)
	# 
	def addquote(self,t):
		mtree = ""
		
		LD = ["(","["]
		RD = [")","]"]
		
		for i in range(len(t)):
			if t[i] in LD and t[i+1] not in LD:   # (xxx...
				mtree += "%s'" % t[i]
			elif t[i] in RD and t[i-1] not in RD: # ...xxx)
				mtree += "'%s" % t[i]
			elif t[i] == ",":      
				if t[i-1] not in RD:  # xxx,...
					mtree += "'"
				if t[i+1] not in LD:  # ...,xxx
					mtree += "%s'" % t[i]
				else:
					mtree += t[i]
			else:
				mtree += t[i]

		# if a tree is like: (xxx,...),xxx -> an end quote is needed as well as
		# another outer brackets.
		if mtree[-1] not in RD:
			mtree += "'"
		return mtree

	#
	# Mid-point root, get clades and split them
	#
	def splitclades(self,retree,tdir):
		
		oup = open("splitclades.out","w")
		oup.write("Family\tClade_lef\tClade_right\n")
		log = open("splitclades.out.errors","w")
		for i in os.listdir(tdir):
			n = "%s/%s" % (tdir,i)
			if i[-3:] not in [".ph","tre"]:
				continue
			print "",i
			# tree in one line, rid of distance, change delimiter, add quote
			tree = self.one_line(n,call=1)
			tree = self.simplify(tree,delimiter="[",remove=1,call=1)
			tree = self.addquote(tree)
			# get number of taxa		
			sp = []
			try:
				self.print_species(tree,sp,enter=1)
			except MemoryError:
				log.write("%s\tmemory_error\n" % i)
				continue
			
			# call retree
			os.system("cp %s intree" % i)
			f = os.popen("%sretree" % retree,"w")
			f.write("Y\n") # setting correct
			f.write("M\n") # mid point 
			f.write("Q\n") # quit
			f.write("Y\n") # save yes
			f.write("R\n") # save rooted
			time.sleep(5)	
			f.close()
			#os.system("rm intree")
			os.system("mv outtree %s.mroot" % i)
			
			# read mid-point rooted tree
			tree = self.one_line(i+".mroot",call=1)	
			tree = self.simplify(tree,delimiter="[",remove=1,call=1)
			tree = self.addquote(tree)
			
			# split clades
			alist = []
			dummy = []
			dummy,alist = self.splits(tree,alist,enter=1)
			for j in alist:
				oup.write("%s\t%s\t%s\n" % (i[:-3],string.joinfields(j[0],","),
												   string.joinfields(j[1],",")))
		
		print "Done!"
	
	
	
	# Modified from:
	#  http://www.pasteur.fr/formation/infobio/python/ch13s03.html#ftn.d0e5618
	# @param tree a nested list string
	#
	def print_species(self,tree,alist,enter=1):
		if enter:
			tree = eval(tree)
		if type(tree) is types.ListType:
			for child in tree:
				self.print_species(child,alist,0)
		else:
			#print tree
			alist.append(tree)

	# Modified from:
	#  http://www.pasteur.fr/formation/infobio/python/ch13s03.html#ftn.d0e5618
	#
	# @param alist  for appending results in each iteration
	# @param enter  set to 1 in outside method call so the string can be 
	#               converted into a nested list.
	#
	def splits(self,tree,alist,enter=0):
		if enter:
			tree = eval(tree)
		if type(tree) is types.ListType:
			all_leave = []
			all_leave2= []
			for child in tree:
				child_leave,alist = self.splits(child,alist,0)
				#print "C:",child_leave,
				all_leave += child_leave
				all_leave2.append(child_leave)
			#print
			#print "A:",all_leave
			alist.append(all_leave2)
			#print "ALIST:",alist
			return all_leave,alist
		else:
			return [tree],alist
		
	# From:
	#  http://www.pasteur.fr/formation/infobio/python/ch13s03.html#ftn.d0e5618	
	def binary_splits(self,tree):
		if type(tree) is types.StringType:
			tree = eval(tree)
		if type(tree) is types.ListType:
			left = self.binary_splits(tree[0])
			right = self.binary_splits(tree[1])
			print left, right
			return left + right
		else:
			return [ tree ]
	
	#
	# Notung is assume to be in the directory where this functon is called.
	# Didn't bother to set classpath. See:
	# http://goby.compbio.cs.cmu.edu/notung/htmldocumentation18.crossref.html
	#
	# gdir  gene tree folder
	# sp    species tree
	# px    tree file postfix
	# root  run notung rooting option (1) or not (0, default)
	#
	def batch_notung(self,gdir,sp,px,root):
		
		gt = os.listdir(gdir)
		os.system("mkdir _notung")
		oup = open(sp+".notung_DL","w")
		oup.write("Family\tTaxa\tDup\tLoss\t#OG/Gene\n")
		oup2= open(sp+".notung_OG","w")
		oup2.write("Fam\tNodeTaxa\tDup(1)/Spe(0)\tLeft\tRight\n")
		for i in gt:
			if px == "" or (px != "" and px == i.split(".")[-1]):
				print "%s" % i
				# run notung
				if os.path.isfile("_notung/%s.rooting.0" % i):
					#print "here"
					pass
				else:
					self.run_notung(gdir,i,sp,root)
					
				# parse output
				if root:
					D,og = self.parse_notung("_notung/%s.rooting.0" % i,1)
				else:
					D,og = self.parse_notung("_notung/%s.reconciled" % i,1)					
				
				# problem with parsing or when I want to just parse the notung
				# results, the tree is there but there is no notung output.
				# this happens when the tree contain taxa that are not specifiy
				# in sp.tre.
				if D != "" and og != "":
					# generate D/L count output
					for j in D:
						#print j,D[j]
						#oup.write("%s\t%s\t%i\t%i\t%s\t%s\n" % \
						#				(i,j,D[j][0],D[j][1],D[j][2],D[j][3]))
						oup.write("%s\t%s\t%s\t%s\t%s\n" % \
										(i,j,str(D[j][0]),str(D[j][1]),D[j][2]))
					
					# generate orthologous group output
					for j in og:
						#print j
						for k in og[j]:
							oup2.write("%s\t%s\t%i\t%s\t%s\n" % \
								(i,j,k[0],",".join(k[1]),",".join(k[2])))
				
		print "Done!"
	
	def run_notung(self,gdir,tre,sp,root):
		
		if root:
			os.system("java -jar Notung-2.1.jar -root " +\
					  "-g %s/%s " % (gdir,tre)          +\
					  "-s %s "    % sp                  +\
				  	  "-speciestag prefix -edgeweights length >> notung.log")
			os.system("mv %s.rooting.0 _notung" % tre)
		else:
			os.system("java -jar Notung-2.1.jar -reconcile " +\
					  "-g %s/%s " % (gdir,tre)               +\
					  "-s %s "    % sp                       +\
				  	  "-speciestag prefix -edgeweights length >> notung.log")
			os.system("mv %s.reconciled _notung" % tre)
	
	#
	# Parse the tree generated by Notung
	#
	# @tree   file generated with the follwoing command
	#         java -jar Notung-2.1.jar -g test_gene.tre -s test_sp.tre -root 
	#                                 -speciestag prefix -edgeweights length
	# @batch  batch run (1) or not (0, default)
	#
	def parse_notung(self,tree,batch=0):
		origTreeName=tree
		try:
			inp = open(tree)
		except IOError:
			print " IOError:",tree
			return "",""
			
		gt  = self.rmlb(inp.readline()) # gene tree
		st  = self.rmlb(inp.readline()) # species tree
		pa  = self.rmlb(inp.readline()) # parameters
		
		
		# parse species tree
		#[&&NOTUNG-SPECIES-TREE(((Poptr1,Arath6)n56,Orysa4)n58,(Ostta1,Chlre3)n63)n64]
		def split(tre,D):
			# get node name and tre without the node name
			n   = tre[tre.rfind(")")+1:]
			tre = tre[1:tre.rfind(")")]
			#print "Enter split:",n,tre
			#sys.exit()
			# find the bifurcating ","
			c = 0
			for j in range(len(tre)):
				if j != 0 and c == 0 and tre[j] == ",":
					#print "break",j
					break
				if tre[j] == "(":
					c += 1
				elif tre[j] == ")":
					c -= 1
			#print n,[tre[:j],tre[j+1:]]
			#sys.exit()
			# put node and branch info into dict
			D[n] = [tre[:j],tre[j+1:]]
			
			# parse the subtrees further if necessary
			if tre[:j].find(",") != -1:
				D = split(tre[:j],D)
			if tre[j+1:].find(",") != -1:
				D = split(tre[j+1:],D)
                        #print ">>>",D
                        #sys.exit()
			return D
		
		# get taxa name out of [branchL,branchR] for each node
		def get_taxa(alist,merge):
			
			T1 = []
			t  = alist[0].split(",")
			for i in t:
				i = i.split("(")
				for j in i:
					j = j.split(")")
					if len(j) != 1:
						for k in j[:1]:
							if k != "" and k.find("*LOST") == -1:
								k = k.split("[")[0]
								T1.append(k)
					elif j != [""] and j[0].find("*LOST") == -1:
						j[0] = j[0].split("[")[0]
						T1.append(j[0])
			T2 = []
			t  = alist[1].split(",")
			for i in t:
				i = i.split("(")
				for j in i:
					j = j.split(")")
					if len(j) != 1:
						for k in j[:1]:
							if k != "" and k.find("*LOST") == -1:
								k = k.split("[")[0]
								T2.append(k)
					elif j != [""] and j[0].find("*LOST") == -1:
						j[0] = j[0].split("[")[0]
						T2.append(j[0])
			
			if merge:
				T1.extend(T2)
				T1.sort()
				return T1
			else:
				return [T1,T2]			
                
		st = st[st.find("("):-1]

		# Check if this is a single species tree		
		if st.find(",") == -1:
			print " single species tree"
			S = st[1:-1]
			n = len(gt.split("%s|" % S))-1
			if batch:
				return {S:["NA","NA",str(n)]},{}
			else:
				oup = open(tree+".split_stat","w")
				oup.write("Branch\tTaxa\tDup\tLoss\t#OG/Gene\n")
				oup.write("-\t%s\tNA\tNA\t%s\n" % (S,str(n)))
				sys.exit(0)			
			
		N = {}  # N = {node_name:[left_branch,right_branch]
		N = split(st,N)
		
		# put taxa into the dict
		T = {}  # taxa
		n1 = {} # same as N except added taxa
		for i in N:
			#print i
			n1[i] = N[i]
			if N[i][0].find(",") == -1:
				n1[N[i][0]] = N[i][0]
				
				if N[i][0] not in T:
					T[N[i][0]] = 0
			if N[i][1].find(",") == -1:
				n1[N[i][1]] = N[i][1]
				if N[i][1] not in T:
					T[N[i][1]] = 0
	
		# parse the subtree further to get taxa in each clade
		n2 = {}
		for i in n1:
			#print i,n1[i]			
			clade = []
			for j in T:                                
				if j in n1[i][0]:
					clade.append(j)
				if j in n1[i][1]:
					clade.append(j)	
			if clade == []:
				clade.append(i)
			n2[i] = [clade,0,0,0] # clade taxa, duplication, loss, #nodes/#gene
			
		# get the number of genes in each species
		for i in T:
			#print "-------",i
			gts = gt.split("%s|" %i)
			#for j in gts:
			#	print "",j
			T[i] = len(gts)-1
		
		# Parse gene tree
		gt = gt[gt.find("("):-1]
		G = {}
		G = split(gt,G)
		og= {} # og = {node_type:[[hasD,og1L,og1R],[hasD,og2L,og2R],...]
		for i in G:
			#print i
			
			# parse the node tag
			# eg. n5:0.01703[&&NHX:S=n48:D=N:B=0.01703]
			S = "" # node name
			D = 0  # duplication
			ntag = i[i.find("[")+1:-1].split(":")
			hasD = 0 # has duplication flag (1) = not a speciation node
			for j in ntag:
				if j[0] == "S":
					S = j[j.find("=")+1:]
				elif j[0] == "D":
					if j[-1] == "Y":
						D += 1
						hasD = 1
			#print "--->",S,D,hasD
			n2[S][1] += D
			n2[S][3] += 1 # this number plus losses minus dup is the 
			              # orthologous group number
			              # applicalbe to clades, NOT to single taxa branch
			              ###
			              # 11/18,06, should not add losses. since OG is
			              # describing the base of the taxa in question
			              # if loss, the base should not have it right before
			              # speciation anyway
			              ###
			
			# parse the branch for losses
			if G[i][0].find(",") == -1 and G[i][0].find("*LOST") != -1:
				SL = G[i][0][:G[i][0].find("*")]
				n2[SL][2] += 1
			if G[i][1].find(",") == -1 and G[i][1].find("*LOST") != -1:
				SL = G[i][1][:G[i][1].find("*")]
				n2[SL][2] += 1
			
			####
			# Getting orthologous group
			####
			# get the species in this clade
			if S in N:
				c_sp = ",".join(get_taxa(N[S],1))
			else:
				c_sp = S
			#print c_sp
			if c_sp not in og:
				og[c_sp] = []
			
			
			# get the gene names in this clade
			c_gn = get_taxa(G[i],0)
			og[c_sp].append([hasD,c_gn[0],c_gn[1]])
			#print "",hasD
			#print "",c_gn[0]
			#print "",c_gn[1]
			#if hasD:
			#	print "",i
			#	print "",G[i]

		
		# straighten out the counts in n2		
		for i in n2:
			if len(n2[i][0]) >1:
				n2[i][3] = str(n2[i][3]-n2[i][1])
				#n2[i][4] = "NA"
			else:
				#n2[i][3] = "NA"
				#n2[i][4] = str(T[i]) # T = {taxa:gene_count}
				n2[i][3] = str(T[i])
		if batch:
			bdict = {}
			for i in n2:
				sp = ",".join(n2[i][0])
				bdict[sp] = n2[i][1:]
			return bdict,og
		else:
			oup = open(tree+".split_stat","w")
			oup.write("Branch:%s\tTaxa\tDup\tLoss\t#OG/#Gene\n"%origTreeName)
			for i in n2:
				oup.write("%s\t%s\t%i\t%i\t%s\n" % \
					(i,",".join(n2[i][0]),n2[i][1],n2[i][2],n2[i][3]))
			
			oup = open(tree+".split_og","w")
			oup.write("NodeTaxa:%s\tDup(1)/Spe(0)\tLeft\tRight\n"%origTreeName)
			for i in og:
				#print i
				for j in og[i]:
					#print "",j
					oup.write("%s\t%i\t%s\t%s\n" % \
						(i,j[0],",".join(j[1]),",".join(j[2])))
		print "Done!", origTreeName
	
	def clustal_to_mega(self,tree):
		
		t = self.one_line(tree,call=1)
		t = t[:t.find("TRICHOTOMY;")]
		
		oup = open("%s.for_mega" % tree,"w")
		oup.write(t+";\n")
		oup.close()
	
	def batch_ctom(self,tree_dir):
		
		tlist = os.listdir(tree_dir)
		for i in tlist:
			print i
			self.clustal_to_mega("%s/%s" % (tree_dir,i))
		print "Done!"
	
	def rmlb(self,astr):                
		#if astr[-2:] == "\r\n":
		#	astr = astr[:-2]			
		#elif astr[-1] == "\n":
		#	astr = astr[:-1]
		astrx=astr.strip()
		return astrx
	
	def help(self):
		print " -f  Function to run:"
		print "    transpose - transpose order of taxa from one line to"
		print "       multiple lines. REQUIRES: -i"
		print "    rename - replace the names in the New Hampshire style"
		print "       tree. REQUIRES: t,n, OPTIONAL: -nhx"
		print "    move_boot - strange format generated by clustal is"
		print "       coverted by this function so TreeExplorer can read it"
		print "       treelist"
		print "    simplify - rid of anything between ':' and ',' or '('"
		print "       REQUIRES: t"
		print "    divide_boot - divide the bootstrap values by a factor"
		print "       Need: t, by"
		print "    one_line - convert tree file to one line of string."
		print "       REQUIRES: t, OPTIONAL: dos"
		print "    subtree - get subtree, NEED: t, x, OPTIONAL: o"
		print "    batch_sub - get subtree in batch, NEED: t, n"
		print "    replace_sub - replace subtree with names, NEED: t, n"
		print "    mega_to_nexus - NEED: t"
		print "    printsp - print species, NEED: t"
		print "    splits - split left and right nodes, NEED: t"
		print "    bsplits - split left and right nodes, NEED: t"
		print "    splitclades - do mid point rooting and divide clades, NEED: D"
		print "       OPTIONAL: retree"
		print "    parse_notung - parse the notung tree. NEED: t, OPT: og"
		print "    batch_notug - run notung and parse the output in batch."
		print "       NEED: D,sp, OPT: px, r"
		print "    clustal_to_mega - convert clustal tree output so mega can"
		print "       read it. Need: t"
		print "    batch_ctom - batch operation of clustal_to_mega, NEED: D"
		print "    get_taxa - NEED: t, OPT: x1,x2"
		print " -i  Input file"
		print " -t  Nexus style tree like the one generated by NJ"
		print " -n  The name file with [new_name][old_name] or the opposite"
		print "      depend on the nflag setting. For batch_sub, just a list of"
		print "      names. For replace_sub, [name][subtree]"
		print " -nhx Generate NHX style tree [1] or not [0,default]"
		print " -by the factor use to divide bootstrap values in a tree"
		print " -d  dos format [1] or not [0, default]"
		print " -g  0 as default with [new][old] name file, or 1 for the"
		print "		  opposite."
		print " -x  taxa that is the outgroup for getting subtree"
		print " -x1 x2 taxa names delineating the range for get_taxa"
		print " -o  for subtree, if a string is specified, it will be used as"
		print "	     file name, otherwise, subtree string will be printed and"
		print "	     returned"
		print " -retree dir with retree, set to be /home/shiu/bin/phylip3.66/exe"
		print " -D  directory with trees"
		print " -sp species tree"
		print " -px tree postfix, default ''"
		print " -r  root (1) or not (0, default)"
		print ""
		sys.exit(0)
		
#-------------------------------------------------------------------------------


utility = tree_utility()

if __name__ == '__main__':

	function = infile = tree = name = taxa = outname = D = sp = px = ""
	nhx	  = dos = nflag = r = 0
	factor   = 1
	futil	= FileUtility.file_util()
	retree = "/home/shiu/bin/phylip3.66/exe/"
	x1 = x2 = ""
	
	for i in range(1,len(sys.argv),2):
		if sys.argv[i] == "-i":
			infile	 = sys.argv[i+1]
		elif sys.argv[i] == "-f":	 
			function   = sys.argv[i+1]
		elif sys.argv[i] == "-t":
			tree	   = sys.argv[i+1]
		elif sys.argv[i] == "-n":		  
			name	   = sys.argv[i+1]
		elif sys.argv[i] == "-nhx":
			nhx	= int(sys.argv[i+1])
		elif sys.argv[i] == "-by":
			factor	 = int(sys.argv[i+1])
		elif sys.argv[i] == "-d":
			dos		= int(sys.argv[i+1])
		elif sys.argv[i] == "-g":
			nflag	  = int(sys.argv[i+1])
		elif sys.argv[i] == "-x":
			taxa	   = sys.argv[i+1]
		elif sys.argv[i] == "-o":
			outname	= sys.argv[i+1]
		elif sys.argv[i] == "-retree":
			retree	   = sys.argv[i+1]
		elif sys.argv[i] == "-D":
			D       = sys.argv[i+1]
		elif sys.argv[i] == "-sp":
			sp       = sys.argv[i+1]
		elif sys.argv[i] == "-px":
			px       = sys.argv[i+1]
		elif sys.argv[i] == "-r":
			r        = int(sys.argv[i+1])
		elif sys.argv[i] == "-x1":
			x1       = sys.argv[i+1]
		elif sys.argv[i] == "-x2":
			x2       = sys.argv[i+1]
		else:
			print "Unknown flag:",sys.argv[i]
			utility.help()
			
	if function == "transpose":
		if infile == "":
			print "\nNeed input file\n"
			utility.help()
		else:
			utility.tranpose(infile)
	elif function == "rename":
		if tree == "" or name == "":
			print "\nNeed tree and name files\n"
			utility.help()
		else:
			utility.rename(tree,name,nhx,nflag)
	elif function == "move_boot":
		if tree == "":
			print "\nNeed tree file\n"
			utility.help()
		else:
			utility.move_boot(tree)
	elif function == "treelist":
		if tree == "":
			print "\nNeed tree file\n"
			utility.help()
		else:
			utility.treelist(tree)
	elif function == "simplify":
		if tree == "":
			print "\nNeed tree file\n"
			utility.help()
		else:
			utility.simplify(tree)
	elif function == "divide_boot":
		if tree == "" or factor == 1:
			print "\nNeed tree file and factor specified\n"
			utility.help()
		else:
			utility.divide_boot(tree,factor)
	elif function == "one_line":
		if tree == "":
			print "\nNeed tree file\n"
			utility.help()
		else:
			utility.one_line(tree,dos)
	
	elif function == "subtree":
		if tree == "" or taxa == "":
			print "\nNeed tree file and taxa id\n"
			utility.help()
		else:
			utility.get_subtree(tree,taxa,outname)
	elif function == "batch_sub":
		if tree == "" or name == "":
			print "\nNeed tree and name files\n"
			utility.help()
		else:
			utility.batch_sub(tree,name)

	elif function == "get_taxa":
		if "" in [tree]:
			print "\nNeed tree name\n"
			utility.help()
		utility.get_taxa(tree,x1,x2)

	elif function == "replace_sub":
		if tree == "" or name == "":
			print "\nNeed tree and name files\n"
			utility.help()
		else:
			utility.replace_sub(tree,name)
	elif function == "mega_to_nexus":
		if tree == "":
			print "\nNeed tree files\n"
			utility.help()
		else:
			utility.mega_to_nexus(tree)
	elif function == "printsp":
		if tree == "":
			print "\nNeed tree files\n"
			utility.help()
		else:
			utility.print_species(tree)
	elif function == "splits":
		if tree == "":
			print "\nNeed tree files\n"
			utility.help()
		else:
			utility.splits(tree)
	elif function == "bsplits":
		if tree == "":
			print "\nNeed tree files\n"
			utility.help()
		else:
			utility.binary_splits(tree)
	elif function == "splitclades":
		if D == "":
			print "\nNeed tree dir\n"
			utility.help()
		else:
			utility.splitclades(retree,D)
	elif function == "parse_notung":
		if tree == "":
			print "\nNeed tree file\n"
			utility.help()
		utility.parse_notung(tree,0)
	elif function == "batch_notung":
		if "" in [D,sp]:
			print "\nNeed gene tree dir and species tree\n"
			utility.help()
		utility.batch_notung(D,sp,px,r)
	elif function == "clustal_to_mega":
		if "" in [tree]:
			print "\nNeed gene tree\n"
			utility.help()
		utility.clustal_to_mega(tree)
	elif function == "batch_ctom":
		if "" in [D]:
			print "\nNeed gene tree dir\n"
			utility.help()
		utility.batch_ctom(D)

	else:
		utility.help()


