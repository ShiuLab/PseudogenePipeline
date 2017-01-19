#/usr/bin/python
import os, sys, string

##
#
#
# 11/30,02
#  The file_to_dict utility is modified so new line character is dealt with. But
#  this action may impact some other utilities. Need to find out.
##

class file_util:

	# 
	# param file1  group file with [group_id][common_id]
	# param file2  [common_id][whatever]
	#
	def get_groups(self,file1,file2):
		
		# gdict = {group_id:[common_id]}
		gdict = {}
		inp = open(file1)
		inl = inp.readline()
		while inl != "":
			L = self.rmlb(inl).split("\t")
			if L[0] not in gdict:
				gdict[L[0]] = [L[1]]
			else:
				gdict[L[0]].append(L[1]) 
				
			inl = inp.readline()

		# cdict = {common_id:[whatever...]}
		cdict = {}
		inp = open(file2)
		inl = inp.readline()
		while inl != "":
			L = self.rmlb(inl).split("\t")
			if L[0] not in cdict:
				cdict[L[0]] = [L[1]]
			else:
				cdict[L[0]].append(L[1])
			inl = inp.readline()

		# output by group
		for i in gdict:
			oup = open("group_%s" % i,"w")
			# iterate common id
			for j in gdict[i]:
				# iterate whatever
				for k in cdict[j]:
					oup.write("%s\t%s\n" % (j,k))
		
		print "Done!"
	
	
	def merge_list(self,file1,file2):
		
		print "Read:",file1
		inp = open(file1)
		oup = open("%s.with.%s" % (file1,file2),"w")
		mdict = {}
		inl = inp.readline()
		while inl != "":
			mdict[self.rmlb(inl)] = 0
			oup.write(inl)
			inl = inp.readline()
		print " %i entries" % len(mdict.keys())
		
		print "\nRead:",file2
		inp = open(file2)
		inl = inp.readline()
		add = 0
		while inl != "":
			if self.rmlb(inl) not in mdict:
				add += 1
				mdict[self.rmlb(inl)] = 0
				oup.write(inl)
			inl = inp.readline()
		print " %i added" % add
		
		print "Done!"
		
	def join(self,file1,file2):
		
		print "Generate dict:", file1
		inp = open(file1)
		dict1 = {}
		inl = inp.readline()
		count1 = 0
		while inl != "":
			L = self.rmlb(inl).split("\t")
			if not dict1.has_key(L[0]):
				dict1[L[0]] = L[1]
				count1 += 1
			inl = inp.readline()
		
		print "Join with:", file2
		inp = open(file2)
		oup = open("merged","w")
		inl = inp.readline()
		count2 = 0
		countJ = 0
		while inl != "":
			L = self.rmlb(inl).split("\t")
			count2 += 1
			if dict1.has_key(L[0]):
				try:					
					L.append(dict1[L[0]])
					oup.write("%s\n" % string.joinfields(L,"\t"))
				except TypeError:
					print L
				countJ += 1
			inl = inp.readline()
		
		print "File1 entries :",count1
		print "File2 entries :",count2
		print "Joined entries:",countJ
		print "Done!"
		
	# swap the first two columns
	def swap_col(self,target):
		inp = open(target)
		oup = open(target+".swap","w")
		inl = inp.readline()
		while inl != "":
			L = self.rmlb(inl).split("\t")
			L0 = L[0]
			L1 = L[1]
			L[1] = L0
			L[0] = L1
			oup.write("%s\n" % string.joinfields(L,"\t"))
			inl = inp.readline()
	
	#
	# Delete columns
	#
	def del_col(self,inf,col):
		
		inp = open(inf)
		oup = open(inf+".dcol","w")
		inl = inp.readline()
		col = col.split(",")
		
		# delete col with larger index first, so reverse col
		tmp = []
		for i in col:
			tmp.append(int(i))
		tmp.reverse()
		col = tmp		

		while inl != "":
			L = self.rmlb(inl).split("\t")
			#print L
			if L != [""]:
				for j in col:
					#print j,L
					del(L[j])
				oup.write("%s\n" % (string.joinfields(L,"\t")))
			inl = inp.readline()
		print "Done!"
	
	#
	# Delete all lines containing the passed names
	#
	def del_line(self,target,name):		
		# generate a name dict
		ndict = self.file_to_dict(name,0)
		print "Total %i names to delete" % len(ndict.keys())
		
		# go through lines in target file
		inp = open(target,"r")
		oup = open(target+".mod","w")
		inl = inp.readline()
		countT = 0
		countR = 0
		while inl != "":
			countT += 1
			L = self.rmlb(inl).split("\t")
			found = 0
			for j in L:
				if ndict.has_key(j):
					countR += 1
					found = 1
					break
			if not found:
				oup.write(inl)
				
			inl = inp.readline()
		
		print "Total %i lines in target, %i deleted" % (countT,countR)
		print "Done!"

	#
	# Mark all lines containing the passed names
	#
	def mark_line(self,target,name):		
		# generate a name dict
		ndict = self.file_to_dict(name,0)
		print "Total %i names to mark" % len(ndict.keys())
		
		# go through lines in target file
		inp = open(target,"r")
		oup = open(target+".mod","w")
		inl = inp.readline()
		countT = 0
		countR = 0
		while inl != "":
			countT += 1
			L = self.rmlb(inl).split("\t")
			found = 0
			for j in L:
				# the following is for matching tokens in target containing
				# xxx.x but name has only xxx, special case, disable after use
				#if j.find(".") != -1:
				#	j = j.split(".")[0]
				if ndict.has_key(j):
					countR += 1
					found = 1
					break
			if not found:
				oup.write(inl)
			else:
				oup.write("%s\t*\n" % (self.rmlb(inl)))
				
			inl = inp.readline()
		
		print "Total %i lines in target, %i marked" % (countT,countR)
		print "Done!"

	
	#
	# This function scan the first two tokens
	#
	# @param target id1\tid2\twhatever
	# @param name   id1\tid2
	# @param oname  output file name
	# @param allT   all matching target lines are taken, else take the 1st
	# @param allN   all names, even if they are redundant
	#
	def twinselect(self,target,name,outfile,allT=0,allN=0):
		
		print "Note that names not in the target will NOT have output"
		
		if outfile == "":
			outfile = "%s_%s" % (target,name)
		
		# process name file
		n = []
		tmp = {}
		inp = open(name)
		inl = inp.readline()
		while inl != "":
			inl = self.rmlb(inl)
			if allN:
				n.append(inl)
			else:
				if inl not in tmp:
					tmp[inl] = 1
					n.append(inl)
			inl = inp.readline()
		
		# go through target
		tmp = {}
		inp = open(target)
		oup = open(outfile,"w")
		inl = inp.readline()
		while inl != "":
			L = self.rmlb(inl).split("\t")
			if "%s\t%s" % (L[0],L[1]) in n:
				if allT:
					oup.write(inl)
				else:
					if "%s\t%s" % (L[0],L[1]) not in tmp:
						tmp["%s\t%s" % (L[0],L[1])] = 1
						oup.write(inl)
			inl = inp.readline()
		
		print "Done!"

	##
	# Get specified columns from a tab-delimited file and order them based
	# on the order specified.
	##
	def get_column(self,table,column,delim="\t"):
		
		inp = open(table,"r")
		if "," in column:
			column = column.split(",")
			oup = open(table+".col_%s" % "-".join(column),"w")
			for i in range(len(column)):
				column[i] = int(column[i])
		else:
			oup = open("%s.col%s" % (table,column),"w")
			column = [int(column)]
			
		print "Get column(s):",column
		
		inl = inp.readline()
		
		if len(inl.split("\t")) < max(column):
			print "Column number too large, there are only %i columns. Quit!"%\
				len(inl.split("\t"))
			sys.exit(0)
		countL = 0	
		while inl != "":
			inl = self.rmlb(inl)
			llist = inl.split(delim)
			tmp = []
			try:
				for j in column:
					tmp.append(llist[j-1])
			except:
				print "Funny line:",j-1
				print ">>",llist
				print "skip!"
							
			oup.write("%s\n" % string.joinfields(tmp,"\t"))
			countL += 1
			inl = inp.readline()
		print "Total %i lines." % countL
		print "Done!"
	
	#
	# If multiple matches, just take the first
	def dredun(self,target,col):

		tdict = {}
		inp = open(target,"r")
		inl = inp.readline()
		while inl != "":
			ilist = inl.split("\t")
			# THE TARGET COLUMN SHOULD NOT HAVE ".", this is for ensembl
			# gene name
			ilist[col] = ilist[col].split(".")[0]
			key  = ilist[col]
			if not tdict.has_key(ilist[col]):					
				tdict[key] = [inl]
			else:
				tdict[key].append(inl)
			inl = inp.readline()
		
		tkeys = tdict.keys()
		oup = open(target+".mredun","w")
		for i in tkeys:
			oup.write(tdict[i][0])

	##
	# Read a file and return a dict. The return style is:
	#  0. 1st token as key, order as value
	#  1. 1st token as key, the rest as value in a single string, omit new line
	#  2. 1st token as key, the rest as a single string omitting new line, a 
	#              list as value, mutiple match to the 1st token in the file, 
	#		       the later ones are appended to the list
	#  3. 1st token as key, the rest in a list as value, omit new line
	#  4. 1st token as key, the rest in a nested list, if multiple match to the
	#		       key, later ones are appended to the list
	#  5. 2nd token as key, the 1st token as value in a list. If the same key
	#              appear again, the 1st token will be appended.
	#  6. Order as key, a list as value with [1st token, full line]
	#  7. 1st token as key, the rest as value in a dict with 2nd token as key, 
	#              the rest of the list as value.
	#
	# @param  infile a tab-delimited file with at least two tokens
	# @param  style  choose 1 or 2 
	# @return adict
	##
	def file_to_dict(self,infile,style=1):

		if not style in [0,1,2,3,4,5,6,7]:
			print "Unknown parsing style, QUIT!"
			sys.exit(0)
			
		adict  = {}
		inp    = open(infile)
		inline = inp.readline()
		order  = 1
		while inline != "":
			if inline in ["\n","\r\n"]:
				inline = inp.readline()
				continue
			if inline[-2:] == "\r\n":
				inline = inline[:-2]
			elif inline[-1] == "\n":
				inline = inline[:-1]
			
			llist = inline.split("\t")
			if len(llist)<2 and style not in [0,6]:
				print "Only one token, can't be processed with specified style!"
				return {}

			if style == 5:
				if not adict.has_key(llist[1]):
					adict[llist[1]] = [inline[:inline.find("\t")]]
				else:
					#print "Redundant:",llist[1]
					adict[llist[1]].append(inline[:inline.find("\t")])

			elif style in [1,2,3,4]:
				if not adict.has_key(llist[0]):
					if style == 1:
						adict[llist[0]] = inline[inline.find("\t")+1:]
					elif style == 2:
						adict[llist[0]] = [inline[inline.find("\t")+1:]]
					elif style == 3:
						adict[llist[0]] = llist[1:]
					elif style == 4:
						adict[llist[0]] = [llist[1:]]
				else:
					if style == 2:
						adict[llist[0]].append(inline[inline.find("\t")+1:])
					elif style == 4:
						adict[llist[0]].append(llist[1:])
			elif style == 0:
				if not adict.has_key(llist[0]):
					adict[llist[0]] = order
					order += 1
				#else:
				#	print " Redundant:",llist[0]
			elif style == 6:
				llist = inline.split("\t")
				# assuming order will not be redundant
				adict[order] = [llist[0],inline]
				order += 1
			elif style == 7:
				L = inline.split("\t")
				if L[0] not in adict:
					adict[L[0]] = {L[1]:L[2:]}
				elif L[1] not in adict[L[0]]:
					adict[L[0]][L[1]] = L[2:]
				else:
					print "Redun:",L
			inline = inp.readline()

		return adict

	##
	# Convert file to a list or nested list, all WITHOUT "\n".
	#
	# @param infile input file
	# @param style  the output format, all without "\n":
	#	       0 - default, a list of lines
	#	       1 - a nested list with elements of a line in the nested ones
	# @param delim  delimiter for style = 1, default is tab.
	##
	def file_to_list(self,infile,style=0,delim="\t"):

		inp    = open(infile,"r")
		inline = inp.readline()
		flist  = []
		while inline != "":
			# rid of empty line and newline char
			if inline == "\n":
				inline = inp.readline()
				continue
			if inline[-2:] == "\r\n":
				inline = inline[:-2]
			elif inline[-1] == "\n":
				inline = inline[:-1]
		
			if style == 0:
				flist.append(inline)
			elif style == 1:
				if inline.find(delim) == -1:
					print "Delimiter not found in inline:",delim
				else:
					flist.append(inline.split(delim))
			else:
				print "Unknown style, file is not loaded to the list..."

			inline = inp.readline()

		return flist
		
		
	def copy(self,t_dir,d_dir,name,ext):
		
		if t_dir[:-1] != "/":
			t_dir += "/"
		if d_dir[:-1] != "/":
			d_dir += "/"
		
		nlist = self.file_to_list(name)
		for i in nlist:
			i = i + ext
			os.system("cp %s%s %s" % (t_dir,i,d_dir))
			
		print "Done!"

	##
	# This function replace 1st tokens out of a file based on a passed name
	# file. A complication is that sometimes names have appended other info
	# already. To ensure names still can be found, the following convention
	# need to be adhere to:
	#
	#  <original name><space>whatever...<\t>other tokens...
	#
	# So whatever is before space is the name, between space and the first "\t"
	# is the descirption that should be added back after.
	#
	# As of 08/23,2004, this script can only deal with non-redundant target name
	# so is not useful when mulitple redundant target names are present. Modify
	# this.
	#
	# @param target file with old names: at the first token only
	# @param names  file with new names: [new name][old name]
	# @param F      force output even name is not replaced [1], or [0,default]i
	# @param tokens the token to replace, default 0 (the first)
	##
	def replace(self,target,names,outfile="",F=0,tokens="0"):
		
		if outfile == "":
			outfile = file1+".renamed"
		oup = open(outfile,"w")
				
		inp = open(names,"r")
		inl = inp.readline()
		ndict = {}
		print "Read names..."
		while inl != "":
			L = self.rmlb(inl).split("\t")
			ndict[L[1]] = L[0]
			inl = inp.readline()

		if tokens.find(",") != -1:
			tokens = tokens.split(",")
			t = []
			for i in tokens:
				t.append(int(i))
			tokens = t
		else:
			tokens = [int(tokens)]

		inp = open(target,"r")
		inl = inp.readline()
		countT = countF = 0
		print "Replace target entries..."
		while inl != "":
			L = self.rmlb(inl).split("\t")
			countT += 1
			for j in tokens:
				if ndict.has_key(L[j]):
					L[j] = ndict[L[j]]
					countF += 1
				else:
					if F:
						L[0] = "<%s>" % L[0]
			oup.write(string.joinfields(L,"\t")+"\n")
			inl = inp.readline()
		
		print "Total %i entries in target file" % countT
		print "      %i with new names" % countF
		
	
	#
	# Replace any token in a tab-delimited file that match the name dict
	# @format  name file format, 0:[new][old] (default), 1:[old][new]
	#
	def replace_all(self,target,names,format=0):
		
		if format != 1:
			ndict = self.file_to_dict(names,5)
		else:
			ndict = self.file_to_dict(names,3)		
		
		#print ndict
		inp = open(target,"r")
		oup = open(target+".replaced","w")
		inl = inp.readline()
		c = 0
		while inl != "":
			L = self.rmlb(inl).split("\t")
			for j in range(len(L)):
				if ndict.has_key(L[j]):
					c += 1
					L[j] = ndict[L[j]][0]
			ostr = string.joinfields(L,"\t")
			oup.write("%s\n" % ostr)
			inl = inp.readline()
		
		print "Token replaced in %i instances." % c
	
	# Replace any string occurred, no requirement in format
	def replace_any(self,target,names):
		ndict = self.file_to_dict(names,1)
		
		inp = open(target,"r")
		oup = open(target+".replaced","w")
		inl = inp.readline()
		c = 0
		while inl != "":
			if c % 100 == 0:
				print " %i x 100" % (c/100)
			#print [inl]
			for j in ndict:
				d = ndict[j]
				L = inl.split(j)
				if len(L) > 1:
					inl = d.join(L)
				#print j,d,len(L)
			oup.write(inl)
			inl = inp.readline()
			c += 1
		print "Done!"
				
	#
	# Replace any string, VERY DANGEROUS, PROBABLY SHOULD NOT BE USED.
	#
	def replace_abs_all(self,target,names):

		ndict = self.file_to_dict(names,5)
		
		inp = open(target,"r")
		oup = open(target+".replaced","w")
		inl = inp.readline()
		c = 0
		while inl != "":
			for j in ndict:
				inl = inl.split(j)
				inl = string.joinfields(inl,ndict[j])
			
			oup.write(inl)
			inl = inp.readline()
		
		print "Token replaced in %i instances." % c
		
	def merge_all(self,matrix1,matrix2,t1,t2,outfile):
		
		uid = {}				
		inp = open(matrix1)
		inl = inp.readline()
		while inl != "":
			L = self.rmlb(inl).split("\t")
			if L[t1] not in uid:
				uid[L[t1]] = [self.rmlb(inl),"-"]
			else:
				print "Redun1:",L[t2]			
			inl = inp.readline()
		
		inp = open(matrix2)
		inl = inp.readline()
		while inl != "":
			L = self.rmlb(inl).split("\t")
			if L[t2] not in uid:
				uid[L[t2]] = ["-",self.rmlb(inl)]
			elif uid[L[t2]][1] == "-":
				uid[L[t2]][1] = self.rmlb(inl)
			else:
				print "Redun2:",L[t2]
			inl = inp.readline()
		
		oup = open(outfile,"w")
		for i in uid:
			oup.write("%s\t%s\n" % (uid[i][0],uid[i][1]))
			
		print "Done!"
	
	#
	# Name can be redundant, suitable for small target file
	#
	def select3(self,target,names,outfile):
		
		target = self.file_to_dict(target,1)
		tokenI = 0 # assume 1st element if tab delimited
		
		inp = open(names)
		oup = open(outfile,"w")
		inl = inp.readline()
		countT = 0
		countF = 0
		while inl != "":
			L = self.rmlb(inl).split("\t")
			countT += 1
			if L[tokenI] in target:
				oup.write("%s\t%s\n" % (L[0],target[L[tokenI]]))
				countF += 1
			else:
				oup.write("%s\t-\n"  % L[0])
			inl = inp.readline()
		
		print "Total %i, found %i" % (countT,countF)
		print "Done!"
	
	#
	# Similar to select. But much simpler and faster, only search first token.
	#
	# @param strStrip  Deal with situation where there is empty space in the 1st
	#                  col.
	# 
	def select2(self,target,names,outfile,strStrip=0):
		
		if outfile == "":
			outfile = target+".selected"
		
		print "Target  :",target
		print "NameFile:",names
		print "OutFile :",outfile
		print ""
		
		# construct name dict
		print "Construct name dict..."
		names = self.file_to_dict(names,0)
		print " %i names" % len(names.keys())

		# go through target file
		print "Read target file..."
		inp = open(target)
		oup = open(outfile,"w")
		inl = inp.readline()
		countF = 0
		#print "missing:"
		while inl != "":
			L = inl.strip().split("\t")
			
			if strStrip:
				L[0] = L[0].strip()
			if names.has_key(L[0]):
				oup.write("\t".join(L)+"\n")
				countF += 1
				if countF % 1e4 == 0:
					print " %i x 10k" % (countF/1e4)
				names[L[0]] = 1
			
			inl = inp.readline()
		oup.close()
		
		print " %i found" % countF
		#for i in names:
		#	if names[i] != 1:
		#		print i
		print "\nDone!"
		
	
	
	# get selected lines
	#
	# @param force whatever is in the name has an output line, even though it
	#              is not found in the target, default: 0, don't do this.
	# @param M     index for the column to match, default 0, if 999 is passed
	#              all column are checked.
	# @param T     target column to get, default "" (get all). Or has a string
	#              passed with "," as seperator. NOT IMPLEMENTED
	# @param u     inclusive, if this is set to [1], multiple matches to
	#              the same idx in the target file will be appended. Otherwise
	#              they will be ignored.
	# @param p     allow period in gene name or not.
	#
	def select(self,target,names,outfile,force=0,M=0,T="",u=1,p=0):

		if outfile == "":
			outfile = target+".selected"
		
		# 1st token as key, rest as a member of a list
		if M == 0:
			tdict = self.file_to_dict(target,2)  
			if tdict == {}:
				print "Deal with file with only one token..."
			 	tdict = self.file_to_dict(target,0)
			 	tmp   = {}
			 	for i in tdict:
			 		tmp[i] = [i]
			 	tdict = tmp
			print len(tdict.keys())
			#for i in tdict:
			#	print i, tdict[i]
		# or, based on the specified idx
		else:
			tdict = {}
			inp = open(target,"r")
			inl = inp.readline()
			while inl != "":
				i = self.rmlb(inl)
				ilist = i.split("\t")
				
				if p:
					pass
				else:
					# THE TARGET COLUMN SHOULD NOT HAVE ".", this is for ensembl
					# gene name
					ilist[M] = ilist[M].split(".")[0]
			
				key  = ilist[M]
				if not tdict.has_key(ilist[M]):					
					del ilist[M]
					tstr = string.joinfields(ilist,"\t")
					tdict[key] = [tstr]
				else:
					del ilist[M]
					tstr = string.joinfields(ilist,"\t")
					tdict[key].append(tstr)
				inl = inp.readline()
	
		countN = 0
		countT = 0
		print "In name file but not in target file:"
		inp = open(names,"r")
		oup = open(outfile,"w")
		inl = inp.readline()
		while inl != "":
			inl = self.rmlb(inl)
			#print inl
			if tdict.has_key(inl):
				if u:
					for j in tdict[inl]:
						oup.write("%s\t%s\n" % (inl,j))
				else:
					oup.write("%s\t%s\n" % (inl,tdict[inl][0]))
			else:
				if force:
					oup.write("%s\t-\n" % inl)
				#print "",inl
				countN += 1
			countT += 1
			inl = inp.readline()
				
		print "Name file: %i names, %i not in target" % (countT,countN)

	#
	# This is originally written for compressed chr sequences. This function
	# decompress the sequence files then concatenate them.
	#
	def anneal(self,t_dir,outfile):
		
		flist  = os.listdir(t_dir)
		fnames = []
		print "Decompress:"
		for i in flist:
			print "",i
			if i[-3:] == ".gz":
				fnames.append(i[:-3])
				os.system("gunzip %s" % i)
		
		fstr = string.joinfields(fnames," ")
		print "Concatenate files..."
		os.system("cat %s > %s" % (fstr,outfile))
		
		print "Done!"
	
	#
	# delete files. It's up to the user to specify postfix, this is not very
	# smart, so "." has to be added by user.
	#
	def delete(self,names,postfix=""):
		
		ndict = self.file_to_dict(names,0)
		for i in ndict:
			os.system("rm %s" % (i+postfix))
	
	def split(self,file,by):
		inp   = open(file)
		lines = inp.readlines()
		leng  = len(lines)/by
		C = 0
		for i in range(by):
			oup  = open(file+"_%i" % (i+1),"w")
			if C == by -1:
				R = [leng*C,leng*(C+2)]
			else:
				R = [leng*C,leng*(C+1)]			
			oup.write("%s" % string.joinfields(lines[R[0]:R[1]],""))
			C = C+1
		
		print "Done!"
		
	def batch_move(self,names,src,dest,P):
		
		print "Source dir :",src
		print "Destination:",dest
		print "File ext   :",P
		
		if src[-1] == "/":
			src = src[:-1]
		if dest[-1] == "/":
			dest = dest[:-1]
		
		nlist = futil.file_to_list(names,0,"\t")
		print "Move %i files:" % len(nlist)
		countF = 1
		for i in nlist:
			if countF % 10 == 0:
				print " %i x 10" % (countF/10)
				
			if P != "":
				i = i + "." + P
			os.system("mv %s/%s %s" % (src,i,dest))
			countF += 1
		print "Done!"
		
	#
	# Exchange col in a two col file
	#
	def exchange(self,file):
		
		inp = open(file,"r")
		inl = inp.readline()
		oup = open(file+".exchange","w")
		while inl != "":
			L = self.rmlb(inl).split("\t")
			oup.write("%s\t%s\n" % (L[1],L[0]))
			inl = inp.readline()
		
		print "Done!"
		
		
	
	def join_tables(self):
		
		flist = os.listdir("./")
		flist.sort()
		colH = []*len(flist) # col header
		rowD = {}            # row dict
		print "%i files in local dir" % len(flist)
		
		print "Go through files..."
		countF = 0
		for i in flist:
			print i
			inp = open(i)
			# first line is header			
			inl = self.rmlb(inp.readline()).split("\t")
			colH.append(inl[1])
			print "",inl[1]
			inl = self.rmlb(inp.readline())
			while inl != "":
				L = self.rmlb(inl).split("\t")
				if L[0] not in rowD:
					rowD[L[0]] = ["-"]*len(flist)
				rowD[L[0]][countF] = L[1]
				inl = inp.readline()
			countF += 1
		
		print "Generate output..."
		rKeys = rowD.keys()
		rKeys.sort()
		oup = open("table_join","w")
		oup.write("\t%s\n" % (string.joinfields(colH,"\t")))
		for i in rKeys:
			oup.write("%s\t%s\n" % (i,string.joinfields(rowD[i],"\t")))
		
		print "Done!"
	
	def survey(self,infile,col):
		
		inp = open(infile)
		col = col.split(",")
		for i in range(len(col)):
			col[i] = int(col[i])
		
		inl = inp.readline()
		D = {}
		while inl != "":
			L = inl.split("\t")
			for j in col:
				if j not in D:
					D[j] = {L[j]:1}
				elif L[j] not in D[j]:
					D[j][L[j]] = 1
				else:
					D[j][L[j]] += 1
			
			inl = inp.readline()
		
		for i in D:
			print "Col:",i
			for j in D[i]:
				print "",j
		
		print "Done!"
	
	def rmlb(self,astr):
		if astr[-2:] == "\r\n":
			astr = astr[:-2]
		elif astr[-1] == "\n":
			astr = astr[:-1]
		
		return astr
	
	#
	# Covert a two column file into a matrix with 1st col as rows, 2nd col as col.
	#
	# @param infile  A two column file
	# @param orderr  Row ordering, one row name per line.
	# @param orderc  Column ordering, one column name per line.
	# 
	def list_to_matrix(self,infile,orderr,orderc):
		inp = open(infile)
		D   = {} # {row:{col:1}}
		col = {} # column name list
		inl = inp.readlines()
		for i in inl:
			[r,c] = i.strip().split("\t")[:2] # row and col
			if c not in col:
				col[c] = 1
				
			if r not in D:
				D[r] = {c:1}
			elif c not in D[r]:
				D[r][c] = 1
			else:
				print "Redun 2nd col for a given row name:",r,c
		
		# Generate output
		oup = open(infile+".matrix","w")
		if orderc == "":
			ckeys = col.keys()
			ckeys.sort()
		else:
			inl   = "".join(open(orderc).readlines())
			ckeys = inl.split("\n")
			
		oup.write("Protein\t%s\n" % "\t".join(ckeys))	# header
		if orderr == "":
			gkeys = D.keys()
			gkeys.sort()
		else:
			inl   = "".join(open(orderr).readlines())
			gkeys = inl.split("\n")
			
		for i in gkeys:
			oup.write(i)   # gene name
			for j in ckeys:
				if j in D[i]:
					oup.write("\t1")
				else:
					oup.write("\t0")
			oup.write("\n")
		
		oup.close()
		print "Done!"
		

	
	def help(self):
		print " -f function"
		print "    copy - copy a set of files. REQUIRES: -t, OPTIONAL: -d,-n,-e"
		print "    list_to_matrix - Convert a 2 col file to a matrix, NEED: i"
		print "    replace - replace the first token in the first file by a"
		print "        name file with [new][old] arrangement. REQUIRES: i,j"
		print "        OPTIONAL: o, F, tokens"
		print "    replace_all - replace any token in a tab delimited file that"
		print "    replace_any - replace any string. NEED: i,j"
		print "        is defined in name file. REQUIRES: i,j."
		print "    del_line - delete any line containing token(s) with the"
		print "        specified names. REQUIRES: i,j"
		print "    mark_line - mark any line containing tokens(s) in the name"
		print "        file. REQUIRE: i,j"
		print "    del_col - delete columns, NEED: i,c"
		print "    swap_col - swap the first two col. NEED: i"
		print "    select - get lines of target file based on passed names."
		print "        REQUIRES: i,j,o. OPTIONAL F,M,T,u,p"
		print "    select2 - simpler select. NEED: i,j,o, OPT: ss"
		print "    select3 - allow redundant j. NEED: i,j,o"
		print "    join - join two two-col files based on the 1st token."
		print "        REQUIRES: i,j"
		print "    join_tables - join all 2 column files in the working dir"
		print "        file should have [id][whatever]"
		print "    get_column - get a particular column (c) from file (i),"
		print "        will sort if c is multiple, OPTIONAL: D"
		print "    anneal - decompress gz files in a dir (d) and concatenate"
		print "    dredun - delete redundant, REQUIRES: i,c"
		print "    delete - delete files specified in a name file, REQUIRES: i"
		print "        OPTIONAL: P"
		print "    split  - NEED: i,F"
		print "    exchange - exchange columns in a two col file, REQUIRES: i"
		print "    batch_move - move a lot of files. NEED: i, t, d, OPTIONAL: p"
		print "    merge_list - merge 2 one column files together, no redun"
		print "        NEED: i,j"
		print "    merge_all - merge 2 matrix based on a particular col. NEED"
		print "        i,j,t1.t2,o"
		print "    twinselect - select based on the 1st 2 tokens. NEED: i,j,"
		print "        OPTIONAL: o,allT, allN"
		print "    get_groups - file1 [group][commonid], file2 [commonid][whatever]"
		print "        NEED: i,j"
		print "    survey - check the content of columns, NEED: i, c"
		print " -i input file 1. The target file, or name file for batch_move"
		print " -j input file 2. The name file"
		print " -t target directory, source directory"
		print " -t1 column token index for matrix 1"
		print " -t2 column token index for matrix 2"
		print " -d destination directory, default ./"
		print " -n file names. default *"
		print " -e file name extension need to add '.', default empty string"
		print " -o output file name"
		print " -orderc ordering for columns"
		print " -orderr ordering for rows"
		print " -F in select, force output of whatever specfieid in name."
		print "    default, 0, won't do this."
		print "    in replace, force output even if not replaced, defaut 0"
		print "    in split, split into this number of files"
		print " -c the column number NOTE: NOT zero-index-based, specify multi"
		print "    columns by using ',' as separator"
		print " -D delimiter"
		print " -M the column to match"
		print " -T the columns to get, separated by ',', default get all"
		print "    NOT IMPLEMENTED"
		print " -P postfix for file names"
		print " -p allow period in gene name or not"
		print " -u inclusive, multiple hits in the target index is appended"
		print " -allT all targets, default 0"
		print " -allN all names, default 0"
		print " -tokens the indices of the columns to replace, separated by ','"
		print " -m name file format, 0:[new][old] (default), 1:[old][new]"
		print " -ss For select2, strip the empty space in strings in target 1st col"
		print ""
		sys.exit(0)



#-------------------------------------------------------------------------------


if __name__ == '__main__':

	function = t_dir = ext = file1 = file2 = outfile = c = T = P = orderc = \
		orderr = ""
	d_dir  = "./"
	name   = "*"
	delim  = "\t"
	F      = M = p = allT = allN = m = ss = 0
	u      = 1
	tokens = "0"	

	futil = file_util()

	for i in range(1,len(sys.argv),2):
		if sys.argv[i] == "-f":
			function   = sys.argv[i+1]
		elif sys.argv[i] == "-t":
			t_dir      = sys.argv[i+1]
		elif sys.argv[i] == "-d":
			d_dir      = sys.argv[i+1]
		elif sys.argv[i] == "-n":
			name       = sys.argv[i+1]
		elif sys.argv[i] == "-e":
			ext        = sys.argv[i+1]
		elif sys.argv[i] == "-i":
			file1      = sys.argv[i+1]
		elif sys.argv[i] == "-j":
			file2      = sys.argv[i+1]
		elif sys.argv[i] == "-o":
			outfile    = sys.argv[i+1]
		elif sys.argv[i] == "-orderc":
			orderc     = sys.argv[i+1]
		elif sys.argv[i] == "-orderr":
			orderr     = sys.argv[i+1]
		elif sys.argv[i] == "-F":
			F          = int(sys.argv[i+1])
		elif sys.argv[i] == "-c":
			c          = sys.argv[i+1]
		elif sys.argv[i] == "-D":
			delim      = sys.argv[i+1]
		elif sys.argv[i] == "-M":
			M          = int(sys.argv[i+1])
		elif sys.argv[i] == "-T":
			T          = sys.argv[i+1]
		elif sys.argv[i] == "-P":
			P          = sys.argv[i+1]
		elif sys.argv[i] == "-u":
			u          = int(sys.argv[i+1])
		elif sys.argv[i] == "-p":
			p          = int(sys.argv[i+1])
		elif sys.argv[i] == "-allT":
			allT       = int(sys.argv[i+1])
		elif sys.argv[i] == "-allN":
			allN       = int(sys.argv[i+1])
		elif sys.argv[i] == "-tokens":
			tokens     = sys.argv[i+1]
		elif sys.argv[i] == "-t1":
			t1       = int(sys.argv[i+1])
		elif sys.argv[i] == "-t2":
			t2       = int(sys.argv[i+1])
		elif sys.argv[i] == "-m":
			m        = int(sys.argv[i+1])
		elif sys.argv[i] == "-ss":
			ss        = int(sys.argv[i+1])
		else:
			print "Unknown flag:",sys.argv[i]
			sys.exit(0)

	if function == "copy":
		if t_dir == "":
			print "\nNeed to specify target directory\n"
			futil.help()
		futil.copy(t_dir,d_dir,name,ext)
	elif function == "list_to_matrix":
		if file1 == "":
			print "\nNeed to specify input file\n"
			futil.help()
		futil.list_to_matrix(file1,orderr,orderc)
	elif function == "replace":
		if file1 == "" or file2 == "":
			print "\nNeed to specify input files\n"
			futil.help()
		futil.replace(file1,file2,outfile,F,tokens)
	elif function == "replace_any":
		if file1 == "" or file2 == "":
			print "\nNeed to specify input files\n"
			futil.help()
		futil.replace_any(file1,file2)
	elif function == "replace_all":
		if file1 == "" or file2 == "":
			print "\nNeed to specify input files\n"
			futil.help()
		futil.replace_all(file1,file2,m)
	elif function == "select":
		if file1 == "" or file2 == "" or outfile == "":
			print "\nNeed to specify input and output files\n"
			futil.help()
		futil.select(file1,file2,outfile,F,M,T,u,p)
	elif function == "select2":
		if file1 == "" or file2 == "" or outfile == "":
			print "\nNeed to specify input and output files\n"
			futil.help()
		futil.select2(file1,file2,outfile,ss)
	elif function == "select3":
		if file1 == "" or file2 == "" or outfile == "":
			print "\nNeed to specify input and output files\n"
			futil.help()
		futil.select3(file1,file2,outfile)
	elif function == "get_column":
		if file1 == "" or c == "":
			print "\nNeed to specify input file and column number\n"
			futil.help()
		futil.get_column(file1,c,delim)
	elif function == "anneal":
		if d_dir == "" or outfile == "":
			print "\nNeed to specify directory and output name\n"
			futil.help()
		futil.anneal(d_dir,outfile)
	elif function == "dredun":
		# c used to be an int. But I changed it to string for get_column. 
		# Didn't do anything to change this funciton. NEED TO CHECK!
		if file1 == "" or c == "":
			print "\nNeed to specify directory and output name\n"
			futil.help()
		futil.dredun(file1,c)
	elif function == "delete":
		if file1 == "":
			print "\nNeed to specify name file\n"
			futil.help()
		futil.delete(file1,P)
	elif function == "split":
		if file1 == "" and F == 0:
			print "\nNeed to specify file and split factor\n"
			futil.help()
		futil.split(file1,F)
	elif function == "batch_move":
		if file1 == "" and t_dir == "" or d_dir == "":
			print "\nNeed to specify name file, target and dest dir\n"
			futil.help()
		futil.batch_move(file1,t_dir,d_dir,P)	
	elif function == "exchange":
		if file1 == "":
			print "\nNeed to specify file\n"
			futil.help()
		futil.exchange(file1)	
	elif function == "del_line":
		if file1 == "":
			print "\nNeed to specify files\n"
			futil.help()
		futil.del_line(file1,file2)	
	elif function == "mark_line":
		if file1 == "":
			print "\nNeed to specify files\n"
			futil.help()
		futil.mark_line(file1,file2)	
	elif function == "del_col":
		if file1 == "" or c == "":
			print "\nNeed to specifiy files\n"
			futil.help()
		futil.del_col(file1,c)	
	elif function == "join":
		if "" in [file1,file2]:
			print "\nNeed to specify files\n"
			futil.help()
		futil.join(file1,file2)	
	elif function == "merge_list":
		if "" in [file1,file2]:
			print "\nNeed to specify files\n"
			futil.help()
		futil.merge_list(file1,file2)	
	elif function == "merge_all":
		if "" in [file1,file2,t1,t2,outfile]:
			print "\nNeed to specify files, output, and tokens\n"
			futil.help()
		futil.merge_all(file1,file2,t1,t2,outfile)	
	elif function == "twinselect":
		if "" in [file1,file2]:
			print "\nNeed to specify files\n"
			futil.help()
		futil.twinselect(file1,file2,outfile,allT,allN)
	elif function == "swap_col":
		if "" in [file1]:
			print "\tNeed Need to specify input file\n"
			futil.help()
		futil.swap_col(file1)			
	elif function == "join_tables":
		futil.join_tables()
	elif function == "get_groups":
		if "" in [file1,file2]:
			print "Need 2 input files"
			futil.help()
		futil.get_groups(file1,file2)			
	elif function == "survey":
		if "" in [file1,c]:
			print "Need input file and col numbers"
			futil.help()
		futil.survey(file1,c)			
	else:
		print "\nUnknown function...\n"
		futil.help()
