#!/usr/bin/

##
# This class is written to manipulate alignment and take care of the calculation
# of dN, dS, and dA
#
# 04/23,03
#  This class was modified many many times since maybe 4 months ago. No documen-
#  tation. Shame on myself. Anyway, rate_pair() need to be changed somewhat to
#  allow the use of full length Ks value as denominator in sliding window ana-
#  lysis.
# 05/15,03
#  Found a serious bug in rate_pair(). The increment was not correct so every
#  step will miss one codon. Fix it.
# 12/15,03
#  Problem with paml, some sequence pairs has dS = -1 in Nei & Gojobori method.
#  So the paml dN and dS are not generated.
# 08/17,04
#  Rate pair now include capability to use Pan's KaKsTools
#
##

import os, sys, FastaManager, Translation, time, TreeUtility, DrawBlocks, \
	   SVGWriter, FileUtility, string
#import AlnFunctions #Gaurav's module
from kaksToolsHacked import calculate

class aln_util:

	def __init__(self):
		pass	 

	#
	# Getting alingments and neighbor joining trees
	#
	def batch_nj(self,phydir,alndir):
	
		#phydir = "~/bin/phylip3.66/exe/" # dir for phylip executables
		#alndir = "/home/shiu/project/stress_expr_at/1_nj/_aln"
		flist = os.listdir(alndir)
		
		c = 1
		for alnfile in flist:
			
			if alnfile[-3:] != "phy":
				continue
			
			print "#########################"
			print "#  File %i out of %i" % (c,len(flist))
			c += 1
			if os.path.isfile("%s/%s.nj.tre" % (alndir,alnfile)) and \
			   os.path.getsize("%s/%s.nj.tre" % (alndir,alnfile)) == 0:
				print "   done!"
			else:
				os.system("rm %s/%s.nj.tre %s/%s.nj.out %s/%s.pd.out" % \
							(alndir,alnfile,alndir,alnfile,alndir,alnfile))
				
				# get distance
				f = os.popen("%sprotdist" % phydir,"w")
				f.write("%s\n" % alnfile)        # aln file name
				f.write("F\n")                   # rename file, have a fake outfile ready
				f.write("%s.pd.out\n" % alnfile) # output file name
				f.write("G\n")                   # gamma correction
				f.write("Y\n")                   # setting ok
				f.write("0.3162\n")              # coefficient of var, assume alpha = 10
				f.close()
				
				# get nj tree
				f = os.popen("%sneighbor" % phydir,"w")
				f.write("%s.pd.out\n" % alnfile) # distance file name
				f.write("F\n")                   # rename file, have a fake outfile ready
				f.write("%s.nj.out\n" % alnfile) # output file name
				f.write("Y\n")                   # accept default parameters
				f.write("F\n")                   # rename file, have a fake outree ready
				f.write("%s.nj.tre\n" % alnfile) # tree file name
				f.close()
		print "Done!"


	##
	# This is for getting alignments and bootstrapped phylogeny for each seq
	# group. 
	# 
	##
	def batch_tree(self,group,fasta,clus_dir,bootstrap=0,isdir=0):
		
		# create group fasta files and get group list back
		if isdir == 0:
			print "Creat group fasta files..."
			# originally setting is (fasta,group,0,1), change the
			# last to 0 so no outgroup is passed
			glist = fasta_manager.get_group_seq(fasta,group,0,0)
		else:
			print "Run fasta files in dir:",fasta
			tmp   = os.listdir(fasta)
			glist = []
			for i in tmp:
				if i.find(".fa") != -1:
					glist.append(i[:i.find(".fa")])
					

		print "\nAlign and generate phylogeny.."
		for i in glist:
			print "\n align:",i
			# generate alignment
			#print "%s/clustalw %s.fa " % (clus_dir,i) + ">> TMP.log"
			os.system("%s/clustalw %s.fa " % (clus_dir,i) + ">> TMP.log")
			
			print " tree :",i
			if bootstrap != 0:
				os.system("%s/clustalw %s.aln /BOOTSTRAP=400 /BOOTLABELS=node" % \
					(clus_dir,i) + ">> TMP.log")
			else:
				os.system("%s/clustalw %s.aln -tree" % \
					(clus_dir,i) + ">> TMP.log")
				
			#tree_util.simplify("%s.phb" % i)
			#os.system("mv %s.phb.simplify %s.tre" % (i,i))
			os.system("rm %s.dnd" % i)
		
		print "Done!"

	def batch_tree_submit(self,groups,fasta,c_dir,w_dir="./",mode=1,spflag=1,pick="",
						 aln_profile="",bootstrap=0):
		# read group info into dict
		gdict = file_util.file_to_dict(groups,2)
		
		# filter, if less than 3 taxa or has only one sp, they will be deleted
		countT = 0
		countU = 0
		oup_log = open(groups+".log","w")
		for i in gdict.keys():
			countT += 1
			same_sp = 1
			taxa_lm = 2
			sp      = ""
			for j in gdict[i]:
				if sp == "":
					sp = j[:2]
				elif sp != j[:2]:
					same_sp = 0
					break
			if len(gdict[i]) < 3:
				oup_log.write("%s\tless than 3 seq\n" % i)
				countU += 1
				del(gdict[i])			
			elif spflag and same_sp:
				oup_log.write("%s\tsame_sp\n" % i)
				countU += 1
				del(gdict[i])
		
		for i in gdict:
			"""
			print "python ~/codes/qsub2.py -s run_%s -c \"python" % i +\
				" ~/codes/AlnUtility.py -f batch_tree2 -group %s" % groups+\
				" -fasta %s -clustal %s -mode %i -sp %i" % (fasta,c_dir,mode,spflag)+\
				" -bootstrap %i -profile %s -pick %s\"" %(bootstrap,aln_profile,i)
			"""
			os.system("python ~/codes/qsub2.py -s run_%s -c \"python" % i +\
				" ~/codes/AlnUtility.py -f batch_tree2 -group %s" % groups+\
				" -fasta %s -clustal %s -mode %i -sp %i" % (fasta,c_dir,mode,spflag)+\
				" -bootstrap %i -profile %s -pick %s\"" %(bootstrap,aln_profile,i))
		
	##
	# Generate alignments based on a passed file with [group][seq_id]
	# Then generate a tree for it.
	#
	# @param mode   run alignment only [0] or both [1, default]
	# @param spflag run even if seqs are the same sp [0], not [1, default].
	#               sp is defined as the first two characters of seq names.
	# @param pick   specify what group(s) to run. Group names seperated by
	#               ",".
	# @param aln_profile name of the profile alignment, if this is not empty
	#               sequences to profile alignment mode for clustsal will
	#               be triggered. Assume the profile alignment file is in the
	#               working dir.
	##
	def batch_tree2(self,groups,fasta,c_dir,w_dir="./",mode=1,spflag=1,pick="",
						 aln_profile="",bootstrap=0):
		
		# get working directories, this is for qsub
		w_dir = os.getcwd()
		print "working_dir:",w_dir
		
		# make directories
		if not os.path.isdir("%s/seq" % w_dir):
			os.system("mkdir %s/seq"  % w_dir)
		if not os.path.isdir("%s/aln" % w_dir):	
			os.system("mkdir %s/aln"  % w_dir)
		if not os.path.isdir("%s/tre" % w_dir):
			os.system("mkdir %s/tre"  % w_dir)
		
		# check if the profile alignment file exisit
		if aln_profile != "" and \
		   not os.path.isfile("%s/%s" % (w_dir,aln_profile)):
			print "Profile %s does not exist in %s" % (aln_profile,w_dir)
			sys.exit(0)
		
		# read group info into dict
		gdict = file_util.file_to_dict(groups,2)
		
		# filter, if less than 3 taxa or has only one sp, they will be deleted
		countT = 0
		countU = 0
		oup_log = open(groups+".log","w")
		for i in gdict.keys():
			countT += 1
			same_sp = 1
			taxa_lm = 2
			sp      = ""
			for j in gdict[i]:
				if sp == "":
					sp = j[:2]
				elif sp != j[:2]:
					same_sp = 0
					break
			if len(gdict[i]) < 3:
				oup_log.write("%s\tless than 3 seq\n" % i)
				countU += 1
				del(gdict[i])			
			elif spflag and same_sp:
				oup_log.write("%s\tsame_sp\n" % i)
				countU += 1
				del(gdict[i])
		
		if pick != "":
			print "Pick groups..."
			pick  = pick.split(",")
			pdict = {}
			for i in pick:
				if gdict.has_key(i):
					pdict[i] = gdict[i]
				#else:
				#	print "Family %s not qualified" % i
			gdict = pdict
			if pdict == {}:
				print "No qualified family, quit!"
				sys.exit(0)
			else:
				print "%i families qualified out of %i picked" % \
					(len(gdict.keys()),len(pick))
			
		oup_log.close()
		
		# this is here so the variable can be reused. Meanwhile, pdict will not
		# be read unless ncessary, this is good for large-scale operation
		# where I usu get the group sequences ready so there is no need to
		# load pdict.
		pdict = {}		
		# iterating groups	
		count = 1
		print "Iterate groups..."
		for i in gdict.keys():
			print "%i %s, %i seq" % (count,i,len(gdict[i]))
			# write sequence into a file 
			if(os.path.isfile("%s/seq/%s.fa" % (w_dir,i)) and \
				os.path.getsize("%s/seq/%s.fa" % (w_dir,i)) != 0):
				print " seq exist..."
			else:
				print " get sequences..."
				pdict = fasta_manager.fasta_to_dict(fasta)
				oup = open("%s/seq/%s.fa" % (w_dir,i),"w")
				for j in gdict[i]:
					oup.write(">%s\n%s\n" % (j,pdict[j]))
				oup.close()
			# generate alignment
			if(os.path.isfile("%s/aln/%s.aln" % (w_dir,i)) and \
				os.path.getsize("%s/aln/%s.aln" % (w_dir,i)) != 0):
				print " aln exist..."
			else:
				print " align..."
				if aln_profile == "":
					os.system("%s/clustalw %s/seq/%s.fa" % \
						   		(c_dir,w_dir,i) + ">>%s/TMP_%s.log" % (w_dir,i))
				else:
					p1 = aln_profile
					#print "%sclustalw -sequences -PROFILE1=%s/%s "%(c_dir,w_dir,p1)+\
					#          "-PROFILE2=%s/seq/%s.fa "     %(w_dir,i)       +\
					#          ">> %s/TMP_%s.log" %(w_dir,i)
					os.system("%s/clustalw -sequences -PROFILE1=%s/%s "%(c_dir,w_dir,p1)+\
					          "-PROFILE2=%s/seq/%s.fa "     %(w_dir,i)       +\
					          ">> TMP_%s.log" %i)
				os.system("mv %s/seq/%s.aln %s/aln/" % (w_dir,i,w_dir))
				os.system("rm %s/seq/%s.dnd"         % (w_dir,i))
				time.sleep(2)
				
			# build tree			
			# 05/06,07 A weird bug. during tree building phase, the aln name
			# pass sometimes clustal will complain that it end with "^Y" so the
			# file is not found. There is no where that this control-Y is found
			# in the code. So I tried rid of the w_dir part of the os.system
			# parameter. And it somehow works.
			# 05/11,07 Still some trees are just not built for unknown reason.
			# Try sleep for 2 sec after alignment but before tree building. Add
			# "./" before aln. Somehow worked...??
			if mode:
				if((os.path.isfile("%s/tre/%s.ph"  % (w_dir,i))       and \
				    os.path.getsize("%s/tre/%s.ph" % (w_dir,i)) != 0) or  \
				   (os.path.isfile("%s/tre/%s.phb" % (w_dir,i))       and \
					os.path.getsize("%s/tre/%s.phb"% (w_dir,i)) != 0)):
					print " tree exist..."
				else:
					print " build tree..."
					if bootstrap == 0:
						os.system("%s/clustalw %s/aln/%s.aln -tree >> %s/TMP_%s.log" % \
								  (c_dir,w_dir,i,w_dir,i))
					else:
						os.system("%s/clustalw aln/%s.aln"       %(c_dir,i) +\
								  " /BOOTSTRAP=%i /BOOTLABELS=node"%bootstrap +\
								  " >> TMP_%s.log"                  %i)
					os.system("mv %s/aln/%s.ph* %s/tre" % (w_dir,i,w_dir))
					time.sleep(2)
			count+=1

		print "Done!"
		

	##
	# Genearte alignments based on a passed file with [group][seq_id]
	#q
	# @param groups tab delimited file with [organisms][clade_id][subtree]
	# @param fasta  fasta sequence file. This file should always have a sequence
	#               named "OUT" for outgroup sequence
	##
	def batch_align(self,groups,fasta,clus_dir,paml_dir):

		# read group info into dict
		gdict = {}
		trees = {}
		inp = open(groups,"r")
		lines = inp.readlines()
		oup = open(groups+".log","w")
		for i in lines:
			if i[-2:] == "\r\n":
				i = i[:-2]
			elif i[-1] == "\n":
				i = i[:-1]
				
			llist = i.split("\t")
			trees["%s_%s" % (llist[0],llist[2])] = llist[3]
			
			# parse tree into taxa information
			tree_str = ""
			for j in llist[3]:
				if j not in ["(",")",":",";"]:
					tree_str += j
			
			tlist = tree_str.split(",")
			gdict["%s_%s" % (llist[0],llist[2])] = tlist
			
			for j in tlist:
				oup.write("%s_%s\t%s\n" % (llist[0],llist[2],j))
		
		# filter, if less than 3 taxa or has only one sp, they will deleted
		countT = 0
		countU = 0
		for i in gdict.keys():
			countT += 1
			same_sp = 1
			taxa_lm = 2
			sp      = ""
			for j in gdict[i]:
				if sp == "":
					sp = j[:2]
				elif sp != j[:2]:
					same_sp = 0
					break
			
			if len(gdict[i]) <=2 or same_sp:
				countU += 1
				del(gdict[i])
				
		print "Total %i groups, %i not qualified" % (countT, countU)
			
		# read seq into dict
		pdict = fasta_manager.fasta_to_dict(fasta)
		
		# iterating groups
		oup_rate = open("%s.matrix" % group, "w")
		
		for i in gdict.keys():
		
			# write sequence into 2 files. 
			oup2= open("TMP2.FA","w")  # tree for rate calculation
			for j in gdict[i]:
				oup2.write(">%s\n%s\n" % (j,pdict[j][1]))
			oup2.close()
			# run clustal for TMP2
			print "%s/clustalw /INFILE=TMP2.FA /OUTFILE=TMP2.GDE /OUTPUT=GDE" % \
					   clus_dir
			os.system("%s/clustalw /INFILE=TMP2.FA /OUTFILE=TMP2.GDE /OUTPUT=GDE" % \
					   clus_dir)
			
			# save tree into temp file
			oup2 = open("TMP2.PHB","w")
			oup2.write("%s;\n" % trees[i])
			
			# convert alignment GDE to fa and then to phylip format
			inp = open("TMP2.GDE","r")
			oup = open("TMP2.GDE.FA","w")
			tlines = inp.readlines()
			for j in tlines:
				if j[0] == "%":
					oup.write(">%s" % j[1:])
				else:
					oup.write(j)
			oup.close()
			self.to_phylip("TMP2.GDE.FA","TMP2")
			
			# run paml
			os.system("%s/codeml" % paml_dir)
			
			# get the matrix part
			inp = open("TMP2.OUT","r")
			inlines = inp.readlines()
			found = 0
			rdict = {}
			nlist = []
			for j in inlines:
				if j.find("ML distances of aa seqs.") != -1:
					found = 1				
				elif found:
					if j[-2:] == "\r\n":
						j = j[:-2]
					elif j[-1] == "\n":
						j = j[:-1]
						
					jlist = j.split(" ")
					tlist = []
					nlist.append(jlist[0])
					for k in jlist[1:]:
						if k != "":
							tlist.append(k)
					rdict[jlist[0]] = tlist
			
			ilist = i.split("_")
			for j in range(len(nlist)):
				for k in rdict.keys():
					if nlist[j] != k:
						oup_rate.write("%s\t%s\t%s\t%s\t%s\n" % \
									   (ilist[0],ilist[1],nlist[j],k,rdict[k][j]))
					else:
						del(rdict[k])

			#os.system("rm -f TMP*")
			
			
	
	##
	# Back translate alignment sequences
	##
	def bt_align(self,pep,nt):
		# read fasta into dicts
		pdict = fasta_manager.fasta_to_dict(pep)
		ndict = fasta_manager.fasta_to_dict(nt)
		
		# back translate each and write into GDE
		oup = open(nt+".aligned.fa","w")
		flag = 0
		for j in pdict.keys():                        
			seq, flag = trans.back_translate2(pdict[j],ndict[j])			
			if flag:
				print "Error:",j, flag
				#print pdict[j]
				#print ""
				#print ndict[j]
				#sys.exit()
			oup.write(">%s\n%s\n" % (j,seq))					

		oup.close()		
	
	#
	# This function take a list of pairs, conduct pairwise alignment with 
	# clustal, then parse the gap information out of the alignments.
	#	
	def parse_gap(self,pairs,fasta,clus_dir):
		
		inp = open(pairs,"r")
		inlines = inp.readlines()
		pdict = fasta_manager.fasta_to_dict(fasta)

		# find -, and return a list of coordinates (alignment-based)
		def gap(astr):	
			c = 1
			gL = gR = 0
			glist = []
			for i in astr:
				if i == "-":
					if gL == 0:
						gL = gR = c
					else:
						gR += 1			
				else:
					if gL != 0:
						glist.append("%i-%i" % (gL,gR))
					gL = gR = 0
				c += 1
			if gL != 0:
				glist.append("%i-%i" % (gL,gR))
			return glist

		oup1= open(pairs+".gap","w")
		oup2= open(pairs+".seqpair","w")
		
		count = 1
		for i in inlines:
		
			# write sequences into temp fasta file and align sequences
			i = self.rmlb(i)	
			# line list, contain ids for each pair
			L = i.split("\t")
			print " %i %s-%s" % (count,L[0],L[1])
			count += 1
			
			# align
			error = self.do_clustal(L[0],L[1],pdict,"TMP",clus_dir)
			
			if error:
				print "   sequence absent"
				oup_rate.write("%s\t%s\t" % (llist[0],llist[1]))
				oup_rate.write("-\t-\tsequence_absent\n")
			else:
				# read alignment GDE file
				inp = open("TMP.GDE","r")
				tmp = inp.readlines()
				inp.close()
				tag = 0
				aln = ["",""]
				for i in tmp:
					if i[0] == "%":
						tag += 1
					if tag == 1 and i[0] != "%":
						aln[0] += self.rmlb(i)
					if tag == 2 and i[0] != "%":
						aln[1] += self.rmlb(i)
				
				qlist = gap(aln[0])
				slist = gap(aln[1])
				
				qlen  = len(pdict[L[0]][1])
				slen  = len(pdict[L[1]][1])
				
				# should be the same as sline
				alen  = len(aln[0])
				#  [q][s][qStr][qEnd][sStr][sEnd][alnL][qGaps][sGaps]
				oup1.write("%s\t%s\t1\t%s\t1\t%s\t%i\t%s\t%s\n" % \
					(L[0],L[1],qlen,slen,alen,string.joinfields(qlist,","),
											  string.joinfields(slist,",")))
				# 2 lines, each pair seperate by an empty line
				oup2.write("%s\t%s\n%s\t%s\n\n" % (L[0],aln[0],L[1],aln[1]))
			
			os.system("rm TMP*")
		
		print "Done!"
	
	##
	# Find 4x degenerate sites. The sequences passed should be aligned already.
	# 
	#
	# @param  seq1  nt seq for gene1
	# @param  seq2  nt seq for gene2
	# @return       a list of two modified sequence
	##
	def find_x4(self,seq1,seq2):
		
		#print seq1
		#print seq2
		x4   = trans.get_x4()
		code = trans.get_nt_code()
		nseq1 = nseq2 = ""
		for i in range(0,len(seq1),3):
			#print seq1[i:i+3],code[seq1[i:i+3]],
			#print seq2[i:i+3],code[seq2[i:i+3]],
			# both codons are 4fold
			if x4.has_key(seq1[i:i+3]) and x4.has_key(seq2[i:i+3]):
				# synonymous
				if x4[seq1[i:i+3]] == x4[seq2[i:i+3]]:
					nseq1 += seq1[i:i+3]
					nseq2 += seq2[i:i+3]
					#print "4x",
			#print ""
		
		#print nseq1
		#print nseq2
			
		return [nseq1,nseq2]
	
	
	##
	# Do the following:
	#  1. get each pair of orthologs from a pep fasta file based on a table
	#  2. clustalw
	#  3. back translate NT sequences based on pep alignment. Still in GDE
	#  4. parse GDE output and convert into "simple" PHYLIP format
	#     default output is TMP.GDE.phy
	#  5. run yn00 based on the control file yn00.ctl, default output is TMP.out
	#     in PAML dir.
	#  6. parse yn00 output, generate another output file with rates with
	#     [org1_geneA][org2_geneA][dN/Ds][dN][dS]
	#
	# !!!!!! FOR THIS TO WORK, THE PAML CONTROL FILE NEED TO BE COPIED TO THE 
	#        WORKING DIR!!!!!!!!
	# 
	# @param wdir    working directory, default local
	# @param pep     polypeptide sequence file to get pairs from
	# @param nt      nucleotide sequence file. The id should be exactly the same
	#                as that in pep. All sequences should also be in the same
	#                frame.
	# @param prog    program to use paml or tpl (Tzeng,Pan,Li, 2004)
	# @param pairs   tab delimited file with ortholog relationships
	# @param windows the sliding-window size for getting Ka, Ks. Default 0
	# @param step    step size of the sliding window
	# @param gdom    graphic domain file
	# @param ks      the Ks value used to calculate Ka/Ks ratio. Default is
	#                "full" which is the Ks for the whole gene, or "window" for
	#                Ks from the sequence of each window.
	# @param justka  0 [default] or 1, just show Ka.
	# @param x4      Only calculate Ka Ks for 4x sites [1] or all [0, default],
	#                will not allow window calculation if 4x = 1
	# @param debug   turn file deletion off (1). Default 0.
	##
	def rate_pair(self,pep,nt,pairs,clus_dir,prog,paml_dir,window,step,
					gdom,ks,justka,x4=0,debug=0):
		
		print "Calculate pairwise evol rate:"
		print " pep    :",pep
		print " nt     :",nt
		print " pairs  :",pairs
		print " clustal:",clus_dir
		print " 4fold  :",x4
		print " program:",prog
		if prog == "paml":
			print " paml   :",paml_dir
		if window != 0:
			print " window :",window
			print " step   :",step
			print " ks_full:",ks
		
		# read fasta into dicts
		print "Read:",pep
		pdict = fasta_manager.fasta_to_dict(pep,verbose=0)
		print "Read:",nt
		ndict = fasta_manager.fasta_to_dict(nt,verbose=0)
		mycount=0		
		if prog == "paml":
			try:
				inp = open("yn00.ctl","r")
			except IOError:
				print "Control file missing, can't do anything... QUIT!"
				sys.exit(0)
		# read and process sequences pair by pair
		inlines = []
		if pairs != "":
			inp = open(pairs,"r")
			inlines = inp.readlines()
			oup_rate = open(pairs+".rate","w")
			if window != 0:
				oup_rate2 = open(pairs+".rate_windows","w")
				
		# no pair info passed, look at all combo
		else:
			print "No pair info passed, all combination..."
			keys = pdict.keys()
			for i in range(len(keys)):
				for j in range(i+1,len(keys)):
					inlines.append("%s\t%s" % (keys[i],keys[j]))
			oup_rate = open("all_pairs.rate","w")
			if window != 0:
				oup_rate2 = open("all_pairs.rate_windows","w")
	
		errors = []
		
		# if window size is specified
		if window != 0:
			# read gdom file into dict
			if gdom == "":
				gdict = {}
			else:
				gdict = file_util.file_to_dict(gdom,1)
			Xsize = 2000
			Ysize = 300 * (len(inlines)+2)
			svg_writer = SVGWriter.svg_writer(pairs+".svg")
			svg_writer.write_header(Xsize,Ysize)

		oup_rate.write("Seq1\tSeq2\tKa\tKs\tvKa\tvKs\n")
		if window != 0:
			oup_rate2.write("Window(nt):%i, Step(nt):%i\n" % (window,step))
			oup_rate2.write("Seq1\tSeq2\tKa\tKs\tvKa\tvKs\tfull_Ks\tStartAA\n")
		
		# start to process pairs
		print "\nProcess pairs:"
		count = 1
		for i in inlines:			
			# write sequences into temp fasta file and align sequences
			i = self.rmlb(i)
				
			# line list, contain ids for each pair
			llist = i.split("\t")
			#print llist
			print " %i %s vs %s" % (count,llist[0],llist[1])
			#mycount+=1
			#if mycount==20:
                        #        sys.exit()
			
			error = self.do_clustal(llist[0],llist[1],pdict,"TMP",clus_dir)
			if error:
				print "   sequence absent"
				oup_rate.write("%s\t%s\t" % (llist[0],llist[1]))
				oup_rate.write("-\t-\t-\t-\tsequence_absent\n")
			else:
				# convert aligned sequence to fasta, then to dict again
				inp = open("TMP.GDE","r")
				oup = open("TMP.GDE.FA","w")
				tlines = inp.readlines()
				inp.close()
				for j in tlines:                                        
					if j[0] in ["%","#"]:
						oup.write(">%s" % j[1:])
					else:
						oup.write(j)
				oup.close()
				tdict = fasta_manager.fasta_to_dict("TMP.GDE.FA",verbose=0)

				#Gaurav's additions for pseudogenes. Should be commented out otherwise
				#Remove positions where one member of the pair has a gap (-)
				#n1=tdict.keys()[0]; n2=tdict.keys()[1]; s1=tdict[n1]; s2=tdict[n2]
				#ns1,ns2=AlnFunctions.removePairwiseGaps(s1,s2)
				#tdict[n1]=ns1;tdict[n2]=ns2
				#print tdict
				#sys.exit()
                                                                                               
				
				# back translate each and write into GDE
				flag = 0
				seqs = ""
				tkeys = tdict.keys()
				tkeys.sort()
				for j in tkeys:
                                        #print ">>>>", j
                                        #print tdict[j]
                                        #print "##"
                                        #print ndict[j]
                                        #print "-------------------------"
					try:
						seq, flag = trans.back_translate2(tdict[j],ndict[j])
						#print ">>>",flag
					except KeyError:
						flag = 2
						print ">>>",j,j in tdict, j in ndict
					if flag != 0:
						errors.append(j)
						break
					seqs += ">%s\n%s\n" % (j,seq)	
				
				# problem back translate sequences
				if flag != 0:
					oup_rate.write("%s\t%s\t-\t-\t-\t-\t" % (llist[0],llist[1]))
					if flag == 1:
						print "   size discrepancy between pep and nt seq"
						oup_rate.write("pep_cds_discrepancy\n")
					elif flag == 2:
						print "   sequence not found"
						oup_rate.write("seq_absent\n")
					count += 1
					continue
				
				# no problem, generate output
				else:
					S = seqs.split("\n")
					# if asking for 4fold sites
					if x4 == 1:
						newseq = self.find_x4(S[1],S[3])	
						S[1] = newseq[0]
						S[3] = newseq[1]
						# no 4x sites
						if len(S[1]) == 0:
							oup_rate.write("%s\t%s\t-\t-\t-\t-\tno_4x_site\n"%\
										(llist[0],llist[1]))
							continue
					if len(S[1]) == 0 or len(S[2]) == 0:
						oup_rate.write("%s\t%s\t-\t-\t-\t-\tno_common_site\n"%\
										(llist[0],llist[1]))
						continue
					oup = open("TMP.NT.FA","w")
					oup.write(string.joinfields(S,"\n"))
					
				
				###############################################################
				# If sequence absent or size discrepancy, won't get past here
				###############################################################
					
				oup.close()
				
				# calculate full length rate
				if prog == "paml":
					# convert FA into "simple PHYLIP"
					self.to_phylip("TMP.NT.FA","TMP")
					rate = self.do_yn00(paml_dir,llist)
				else:
					tmp  = calculate(seqs)
					rate = []
					for i in tmp:
						rate.append(str(i))                                
				
				# deal with situations where no Ka or Ks is generated.				
				prog_error = 0
				try:		
					Ka_full = float(rate[0])
					Ks_full = float(rate[1])
				except ValueError:
					prog_error = 1
				except IndexError:
					prog_error = 2
					#print rate
                                        #sys.exit()
				
				#########################
				# deal with window size
				#########################
				
				# program error
				if prog_error != 0:
					oup_rate.write("%s\t%s\t%sprog_err%i\n" % 
							(llist[0],llist[1],"-\t"*4,prog_error))
				# output full length only
				elif window == 0:
					oup_rate.write("%s\t%s\t%s\n" % (llist[0],llist[1],
										string.joinfields(rate,"\t")))
				# go through windows
				elif window > 0 and x4 == 0:
					# generate a full length ka ks output as well
					oup_rate.write("%s\t%s\t%s\n" % (llist[0],llist[1],
										string.joinfields(rate,"\t")))

					# make sure the window and step sizes are multiple of 3
					if window%3 != 0:
						window += (3-window%3)
					if step % 3 != 0:
						step   += (3-step%3)
						
					# nt sequence dict for a pair
					sdict = fasta_manager.fasta_to_dict("TMP.NT.FA")
					keys  = sdict.keys()
					leng  = len(sdict[keys[0]][1])
					idx   = 0
					while idx+window < leng:
						#print "\n>>>>> WINDOW:",(idx/3+1)
						#print sdict[keys[0]][1][idx:idx+window]
						#print sdict[keys[1]][1][idx:idx+window]
						oup_rate2.write("%s\t%s\t" % (llist[0],llist[1]))
						
						seqs = ">%s\n%s\n>%s\n%s\n" % \
								(keys[0],sdict[keys[0]][1][idx:idx+window],
								 keys[1],sdict[keys[1]][1][idx:idx+window])
						if prog == "paml":
							oup = open("TMP.NT.FA","w")
							oup.write(seqs)
							oup.close()
							self.to_phylip("TMP.NT.FA","TMP")
							rate = self.do_yn00(paml_dir,llist)
						else:
							#print seqs
							tmp  = calculate(seqs)
							rate = []
							for i in tmp:
								rate.append(str(i))
						
						# so far, the rate is empty if two sequences are exactly
						# the same, paml throws Error: DistanceF84: input err..
						# so the Ka/Ks is set to *
						if rate == []:
							oup_rate2.write("%s%s\t%i\n" % \
													("-\t"*4,Ks_full,(idx/3+1)))
						else:									
							oup_rate2.write("%s\t%f\t%i\n" % \
								(string.joinfields(rate,"\t"),Ks_full,idx/3+1))
											 
						idx += step # because idx start from 0, no need to +1
											
					#######################
					# DISABLED!!!
					"""
					# modify gdom line
					try:
						glist = self.modify_gdom([llist[0]+"\t"+gdict[llist[0]],
												  llist[1]+"\t"+gdict[llist[1]]],
												 "TMP.GDE.FA")
					except KeyError:
						glist = [llist[0]+"\terror",
								 llist[1]+"\terror"]
					# generate graphics
					rkeys = rdict.keys()
					rkeys.sort()
					print "\nKa/Ks values:"
					for k in rkeys:
						print k/3+1,"\t",rdict[k]
					draw_block.combine(glist,rdict,svg_writer,count,step)
					"""
					########################
				
				if prog == "paml":
					# Delete tmp files
					if debug == 0:
						os.system("rm TMP* 2YN* rst* rub")
					pass
					
			count += 1
			
			if window != 0:
				svg_writer.write_footer()
			
		print "Size discrepancies between Nt and Pep:",
		if error != []:
			for i in errors:
				#print "",i
                                pass
		else:
			print " none..."
		
		print "\nDone!\n"
	
	##
	# subroutine for synchronize alignments and gdom coordinates THERE ARE BUGS
	# IN THIS THING. SOME DON'T QUITE ALIGN RIGHT.
	#
	#
	# @param glist   a list of gdom lines
	# @param aligned the fasta file of polypeptide alignments
	##
	def modify_gdom(self,glist,aligned):
		
		adict = fasta_manager.fasta_to_dict(aligned)
		#print adict
		# scan for gaps
		gap = {}
		for i in adict:
			is_gap = 0
			gidx   = 0
			gap[i] = {}
			seq = adict[i][1]
			for j in range(len(seq)):
				#print seq[j]
				if seq[j] == "-":
					if not is_gap:
						is_gap = 1
						gidx   = j-1	
						gap[i][gidx] = 1
					else:
						gap[i][gidx] += 1
				else:
					is_gap = 0
		
		# modify gdom
		glist2 = []
		for i in glist:
			print "1:",i
			dlist = i.split("\t")
			seqid = dlist[0]   
			leng  = int(dlist[-1])
			leng  = self.add_gap(seqid,gap[seqid],leng)		
			
			dlist = dlist[1:-1]
			dstr  = seqid + "\t"
			print gap[seqid]
			for j in range(0,len(dlist),2):
				Nidx = int(dlist[j])
				Nidx = self.add_gap(seqid,gap[seqid],Nidx)				
				dom  = dlist[j+1][:dlist[j+1].find("|")]
				Cidx = int(dlist[j+1][dlist[j+1].find("|")+1:])
				Cidx = self.add_gap(seqid,gap[seqid],Cidx)
				dstr += "%i\t%s|%i\t" % (Nidx,dom,Cidx)
			dstr += "%i" % leng
			print "2:",dstr
			glist2.append(dstr)
		
		# generate a gap list with [seqid,start1,end1,start2,end2...]
		plist = []
		for i in gap.keys():
			gkeys = gap[i].keys()
			gkeys.sort()
			#for i in gkeys:
			
			#
			# UNFINISHED!!!!!!
			#	
		
		return glist2
		
	#
	# subroutine for modify gdom						
	#
	def add_gap(self,seqid,gapb,idx):
		
		for i in gapb.keys():
			if i < idx:
				idx += gapb[i]
		
		return idx
	
	#
	# Before 8/23,04, the values came from 2YN.* files. Realize that the SE
	# reported are not used. So starting parsing rst instead. The sample var
	# instead of std err is returned.
	#			
	def do_yn00(self,paml_dir,llist):
		
		# run paml
		#os.system("%syn00 >> run.log" % paml_dir)
		os.system("yn00 >> run.log")
		
		# parse the rst file with Ka, Ks, SEka, SEks
		inp = open("rst","r")
		tmp = inp.readline().split(" ")
		L = []
		for i in tmp:
			if i != "":
				L.append(i)
		#   			                          -6      -4  -3      -1  
		# L: [NG][Ka][Ks][w][YN:][t][kappa][omega][Ka][+-][SE][Ks][+-][SE]
		if L == [] or "YN:" not in L:
			rate = []
		else:
			#print L[-6],L[-3],L[-4],L[-1]
			vKa = float(L[-4])*float(L[-4])
			vKs = float(L[-1])*float(L[-1])
			rate = [L[-6],L[-3],str(vKa),str(vKs)]
		inp.close()
		
		return rate
	
	
	##
	# Generate pair-wise alignment with clustal, output format is set as GDE
	# will be stored in the calling directory.
	#
	# @param id1       first seq name
	# @param id2       second seq name
	# @param fdict     the dict converted from fasta
	# @param outheader the output file name, will add ".GDE" by default
	# @param clus_dir  the directory with clustalw
	##
	def do_clustal(self,id1,id2,fdict,outheader,clus_dir):
		
		error = 0
		loup = open(outheader+".FA","w")
		try:
			loup.write(">%s\n%s\n>%s\n%s\n" % (id1,fdict[id1],id2,fdict[id2]))
		except KeyError:
			if not fdict.has_key(id1):
				print "   sequence not found:",id1
			if not fdict.has_key(id2):
				print "   sequence not found:",id2
			error = 1
		loup.close()
		
		# the infile spec doesn't work with clustalw v1.83 for some reason.
		# rid of it for now
		#os.system("%s/clustalw /INFILE=%s.FA /OUTFILE=%s.GDE /OUTPUT=GDE" % \
		#		   (clus_dir,outheader,outheader))
		os.system("%sclustalw2 %s.FA /OUTFILE=%s.GDE /OUTPUT=GDE >> run.log" % \
				   (clus_dir,outheader,outheader))
		#sys.exit(0)
		return error
	
					
	##
	# Take a GDE format file and output it in simply phylip format:
	# GDE: 
	# #seq1
	# seq...
	#
	# Simply Phylip:
	#    [num]   [len]
	# seq1
	# seq...
	#
	# @param seqfile   the file to be converted
	# @param format    gde or fa [default]
	# @param outheader the output file name, will add ".PHY" by default.
	##
	def to_phylip(self,seqfile,outheader,format="fasta"):

		inp = open(seqfile,"r")
		inlines = inp.readlines()
		
		TAG = ">"
		if format == "gde":
			TAG == "%"
		
		# first, try to get these two parameters
		seq_num = 0
		seq_len = 0
		found_1st = 0
		found_2nd = 0
		c = 0
		for i in inlines:
			if i[0] == TAG:
				if not found_1st:
					seq_num += 1
					found_1st = 1
				else:
					seq_num += 1
					found_2nd = 1
			# just get length from the first sequence
			elif not found_2nd:				
				c += 1
				seq_len += len(self.rmlb(i))
				#print seq_len, c
		
		# write phylip file
		oup = open(seqfile+".PHY","w")
		oup.write("    %i    %i\n" % (seq_num,seq_len))
		for i in inlines:
			if i[0] == TAG:
				if i[0].find(" ") != -1:
					
					ilist = i[0].split(" ")
					i = ""
					for j in ilist:
						i = i + j + "_"
					i = i[:-1]
				
				oup.write(i[1:])
			else:
				oup.write(i)
		
		oup.close()
	
	##
	# This is for parsing the output from running codeml directly.
	#
	# @param paml_out codeml output
	# @param style    not implemented, for dealing with various output format
	##
	def parse_paml(self,paml_out,style=0):
		
		inp   = open(paml_out,"r")
		oup   = open(paml_out+".parse","w")
		lines = inp.readlines()
		
		S_tag   = "pairwise comparison,"
		S_found = 0
		
		taxa1 = taxa2 = dN = dS = ""
		
		for i in lines:
			
			if i[-2:] == "\r\n":
				i = i[:-2]
			elif i[-1] == "\n":
				i = i[:-1]
			
			if i[:len(S_tag)] == S_tag:
				S_found = 1
			
			if S_found and i not in ["\r\n","\n"]:
				
				# write entry
				if i.find("(") != -1:
					if taxa1 != "":
						print taxa1,taxa2,">",dN,dS
						oup.write("%s\t%s\t%s\t%s\n" % (taxa1,taxa2,dN,dS))
					taxa1 = i[i.find("(")+1:i.find(")")]
					taxa2 = i[i.rfind("(")+1:i.rfind(")")]
				elif i.find("dN/dS") != -1:
					
					# convert all "=" to space, because sometimes the key will
					# have no space from the values
					new = ""
					for j in i:
						if j == "=":
							j = " "
						new += j

					llist = new.split(" ")
					ldict = {}
					var   = ""
					for j in llist:
						if j != "":
							if var == "":
								var = j
							else:
								ldict[var] = j
								var = ""
					dN = ldict["dN"]
					dS = ldict["dS"]
	
	
	#
	# This function takes the modified crossPairs file and generate a output
	# with pairs. Before running this, the modification shoudld be done so 
	# [ and ] and ( and ) and space are gone. And members of a pair seperated by
	# ",". Pairs seperated by ",".
	#
	def get_pairs(self,clade):
		
		inp = open(clade,"r")
		oup = open(clade+".pairs","w")
		inl = inp.readline()
		while inl != "":
			inl   = self.rmlb(inl)
			ids = inl.split("\t")[-1]
			ids = ids.split(",")
			tdict = {}
			for j in ids:
				if not tdict.has_key(j):
					tdict[j] = 0
			ids = tdict.keys()
			for j in range(len(ids)):
				for k in range(j+1,len(ids)):
					oup.write("%s\t%s\n" % (ids[j],ids[k]))
			
			inl = inp.readline()
		
		print "Done!"
	
	#
	# Take the output of TreeUtility.splitclades()
	#
	# @param clade  [fam_id][clade_L][clade_R]
	#
	def get_pairs2(self,clade):
		
		inp = open(clade)
		oup = open(clade+".pairs","w")
		inp.readline()
		inl = inp.readline()
		while inl != "":
			L = self.rmlb(inl).split("\t")
			L[1] = L[1].split(",")
			L[2] = L[2].split(",")
			for j in L[1]:
				for k in L[2]:
					# only output same speices pairs
					if j[:2] == k[:2]:
						oup.write("%s\t%s\n" % (j,k))
			
			inl = inp.readline()
		
		print "Done!"
		
	#def batch_rate_pair(self,
	
	def rmlb(self,astr):
		if astr[-2:] == "\r\n":
			astr = astr[:-2]
		elif astr[-1] == "\n":
			astr = astr[:-1]		
		return astr
	
	def help(self):
		print " -f function"
		print "     rate_pair - Calculate dN, dS for sequence pairs. REQUIRES:"
		print "        pep,nt, OPTIONAL: p,paml,clutsal,window,step,gdom,"
		print "        pairs,ks,justka,x4,debug"
		print "     get_pairs - parse modified Pan xpair file. REQUIRES -clade"
		print "     get_pairs2 - parse clade file. REQUIRES -clade"
		print "     batch_align - align a bunch of sequences one group after"
		print "        another. REQUIRES:-group,-fasta. OPTIONAL:-clustal,-paml"
		print "     batch_nj - generate neighbo-joining trees"
		print "     batch_tree - generated group-based alignments then trees"
		print "        with bootstraps. REQUIRES:group,fasta, OPTIONAL:clustal"
		print "        bootstrap"
		print "     batch_tree2- generated group-based alignments then trees"
		print "        with bootstraps. REQUIRES:group,fasta, OPTIONAL:clustal,"
		print "        mode, sp, profile, bootstrap, rrange"
		print "     bt_align - back translate polypeptide alignment. REQUIRES:"
		print "        pep, nt"
		print "     to_phylip - convert Fasta to phylip format. REQUIRES: -fasta" 
		print "     parse_paml - parse the output of paml. Only for codeml at"
		print "        this point. REQUIRES: paml"
		print "     parse_gap - align a pair of sequences then get the gap info"
		print "        out. REQUIRES: pairs, fasta"
		print " -p     program to use, paml or tpl"
		print " -pep   protein sequence file"
		print " -nt    nucleotide seq file"
		print " -fasta sequence file, if isdir is set to 1, this means the"
		print "        fasta dir"
		print " -isdir passed fasta parameter is the dir of fasta files"
		print " -pairs tab delimited file defining pairwise relationship. If"
		print "        not specified, all combinations"
		print " -clade file with [famid][cladeL][cladeR]"
		print " -paml  paml folder. For parse_paml, it is the paml output."
		print " -clustal clustal folder"
		print " -group file with [group_id][seq_name]"
		print " -window sliding-window size for rate_pair, default 0. Should be"
		print "        multiples of 3"
		print " -step  step size for sliding window, default 1, should be 3N"
		print " -gdom  graphic domain text file, only required with sliding win"
		print " -ks    full (default) if full length Ks is used for calculation"
		print "        of Ka/Ks or window if values from sequence windows used"
		print " -justka 0 [default], calculate Ka/Ks; 1, just get Ka output"
		print " -mode  0, run alignment only; 1, both alignment and tree"
		print " -sp    0, run even if all seq in a group from the same sp; 1,"
		print "        won't do those only from one sp (default)"
		print " -pick  run particular groups, names seperated by ','"
		print " -bootstrap number of replicates, 0 means no bootstrap"
		print " -profile for batch_tree2, a profile alignment can be passed"
		print " -x4    calculate only 4x sites (1) or not (0, default)"
		print " -phydir dir for phylip"
		print " -alndir dir with alignment files"
		print ""
		sys.exit(0)

#
# Function call
#


if __name__ == '__main__':

	util          = aln_util()
	tree_util     = TreeUtility.tree_utility()
	file_util     = FileUtility.file_util()
	fasta_manager = FastaManager.fasta_manager()
	trans         = Translation.translate()
	draw_block    = DrawBlocks.draw_blocks()
	function = pep = nt = pairs = group = gdom = clade = pick = aln_profile = ""
	keep = window  = justka = bootstrap = isdir = x4 = debug = 0
	step = mode = sp = 1
	clustal =""
	prog    ="paml"
	paml    ="paml"
	wdir    ="./"
	ks      ="full"
	phydir = ""
	alndir = ""
	
	for i in range(1,len(sys.argv),2):
		if sys.argv[i] == "-f":
			function   = sys.argv[i+1]
		elif sys.argv[i] == "-pep":
			pep        = sys.argv[i+1]
		elif sys.argv[i] == "-nt":
			nt         = sys.argv[i+1]
		elif sys.argv[i] == "-p":
			prog       = sys.argv[i+1]
		elif sys.argv[i] == "-pairs":
			pairs      = sys.argv[i+1]
		elif sys.argv[i] == "-clade":
			clade      = sys.argv[i+1]
		elif sys.argv[i] == "-paml":
			paml       = sys.argv[i+1]
		elif sys.argv[i] == "-clustal":
			clustal    = sys.argv[i+1]
		elif sys.argv[i] == "-keep":
			keep       = int(sys.argv[i+1])
		elif sys.argv[i] == "-clustal":
			clustal    = sys.argv[i+1]
		elif sys.argv[i] == "-group":
			group      = sys.argv[i+1]
		elif sys.argv[i] == "-window":
			window     = int(sys.argv[i+1])
		elif sys.argv[i] == "-step":
			step       = int(sys.argv[i+1])
		elif sys.argv[i] == "-gdom":
			gdom       = sys.argv[i+1]
		elif sys.argv[i] == "-ks":
			ks	       = sys.argv[i+1]
		elif sys.argv[i] == "-fasta":
			fasta      = sys.argv[i+1]
		elif sys.argv[i] == "-justka":
			justka     = sys.argv[i+1]
		elif sys.argv[i] == "-mode":
			mode       = int(sys.argv[i+1])
		elif sys.argv[i] == "-sp":
			sp         = int(sys.argv[i+1])
		elif sys.argv[i] == "-pick":
			pick       = sys.argv[i+1]
		elif sys.argv[i] == "-bootstrap":
			bootstrap  = int(sys.argv[i+1])
		elif sys.argv[i] == "-isdir":
			isdir      = int(sys.argv[i+1])
		elif sys.argv[i] == "-x4":
			x4         = int(sys.argv[i+1])
		elif sys.argv[i] == "-debug":
			debug      = int(sys.argv[i+1])
		elif sys.argv[i] == "-profile":
			aln_profile= sys.argv[i+1]
		elif sys.argv[i] == "-phydir":
			phydir     = sys.argv[i+1]
		elif sys.argv[i] == "-alndir":
			alndir     = sys.argv[i+1]
		else:
			print "UNKNOWN FLAG:",sys.argv[i]
			print "Do -h to get help."
			sys.exit(0)

	if function == "rate_pair":
		if pep == "" or nt == "":
			print "\nNeed pep and nt fasta file, dirs\n"
			util.help()
		util.rate_pair(pep,nt,pairs,clustal,prog,paml,window,step,gdom,ks,
						justka,x4,debug)
	elif function == "get_pairs":
		if clade == "":
			print "\nNeed clade file\n"
			util.help()
		util.get_pairs(clade)
	elif function == "get_pairs2":
		if clade == "":
			print "\nNeed clade file\n"
			util.help()
		util.get_pairs2(clade)
	elif function == "batch_tree":
		if group == "" or fasta == "":
			print "\nNeed fasta and group files\n"
			util.help()
		util.batch_tree(group,fasta,clustal,bootstrap,isdir)
	elif function == "batch_tree2":
		if group == "" or fasta == "":
			print "\nNeed fasta and group files\n"
			util.help()
		util.batch_tree2(group,fasta,clustal,wdir,mode,sp,pick,
						 aln_profile,bootstrap)
	elif function == "batch_tree_submit":
		if group == "" or fasta == "":
			print "\nNeed fasta and group files\n"
			util.help()
		util.batch_tree_submit(group,fasta,clustal,wdir,mode,sp,pick,
						 aln_profile,bootstrap)
	elif function == "batch_align":
		if fasta == "":
			print "\nNeed fasta file\n"
			util.help()
		util.batch_align(group,fasta,clustal,paml)
	elif function == "bt_align":
		if pep == "" or nt == "":
			print "\nNeed pep and nt fasta file\n"
			util.help()
		util.bt_align(pep,nt)
	elif function == "to_phylip":
		if fasta == "":
			print "\nNeed fasta file\n"
			util.help()
		util.to_phylip(fasta,fasta)
	elif function == "parse_paml":
		if paml == "":
			print "\nNeed paml output file\n"
			util.help()
		util.parse_paml(paml)
	elif function == "parse_gap":
		if pairs == "" or fasta == "":
			print "\nNeed pairs and fasta files\n"
			util.help()
		util.parse_gap(pairs,fasta,clustal)
	elif function == "batch_nj":
		if "" in [phydir,alndir]:
			print "\nNeed dir info for phylip and alignments\n"
			util.help()
		util.batch_nj(phydir,alndir)
			
	elif function == "test":
		gdom = ["At1g79620	1	SIGNAL|27	96	LRR|24	121	LRR|23	145	LRR|23	169	LRR|23	199	LRR|23	248	LRR|22	567	TRANS|23	647	STYKc|270	980",
				"Osi007698.1	1	SIGNAL|33	101	LRR_PS|24	126	LRR_PS|23	150	LRR_PS|23	174	LRR_PS|23	204	LRR_PS|23	231	LRR_PS|21	253	LRR_PS|23	277	LRR_PS|24	302	LRR_PS|23	548	TRANS|23	619	STYKc|270	954"]
		util.modify_gdom(gdom,"TMP.GDE.FA")		
	else:
		print "\nUnknown function...\n"
		util.help()


# DEPRECATED METHODS


"""
	
	##
	# Genearte alignments based on a passed file with [group][seq_id]
	#
	# THIS IS THE OLD METHOD. Replaced by a method that take a group file
	# with [tree_id][
	#
	# @param group tab delimited file with [group_id][seq_name]
	# @param fasta fasta sequence file. This file should always have a sequence
	#              named "OUT" for outgroup sequence
	##
	def batch_align(self,group,fasta,keep,clus_dir,paml_dir):
		# read group info into dict
		gdict = {}
		inp = open(group,"r")
		lines = inp.readlines()
		for i in lines:
			if i[-2:] == "\r\n":
				i = i[:-2]
			elif i[-1] == "\n":
				i = i[:-1]
			llist = i.split("\t")
			llist[0] = llist[0].upper()
			if gdict.has_key(llist[0]):
				if llist[1] in gdict[llist[0]]:
					print " Redundant:",llist[1]
				else:
					gdict[llist[0]].append(llist[1])
			else:
				gdict[llist[0]] = [llist[1]]
		
		# filter, if less than 3 taxa or has only one sp, they will deleted		
		countT = 0
		countU = 0
		oup = open(group+".log","w")
		oup.write("excluded\tsame_sp\n")
		for i in gdict.keys():
			countT += 1
			same_sp = 1
			taxa_lm = 2
			sp      = ""
			for j in gdict[i]:
				if sp == "":
					sp = j[:2]
				elif sp != j[:2]:
					same_sp = 0
					break
			
			if len(gdict[i]) <=2 or same_sp:
				countU += 1
				oup.write("%s\t%s\n" % (i,same_sp))
				del(gdict[i])
		print "Total %i groups, %i not qualified" % (countT, countU)
			
		# read seq into dict
		pdict = fasta_manager.fasta_to_dict(pep)
			
		# iterating groups
		oup_rate = open("%s.matrix" % group, "w")
		
		for i in gdict.keys():
			
			# write sequence into 2 files. 
			oup = open("TMP1.FA","w") # with outgroup, tree for use later
			oup2= open("TMP2.FA","w")  # tree for rate calculation
			for j in gdict[i]:
				oup.write( ">%s\n%s\n" % (j,pdict[j][1]))
				oup2.write(">%s\n%s\n" % (j,pdict[j][1]))
			oup.write(">OUT\n%s\n" % (pdict["OUT"][1]))
			oup.close()
			oup2.close()
			
			# run clustal for TMP1
			os.system("%s/clustalw /INFILE=TMP1.FA /OUTFILE=TMP1.GDE /OUTPUT=GDE" % \
					   clus_dir)
			
			os.system("%s/clustalw /INFILE=TMP1.GDE /BOOTSTRAP=100 /BOOTLABELS=node" % clus_dir)
			tree_util.simplify("TMP1.PHB")
			os.system("mv TMP1.PHB.simplify %s.TRE" % i)
			
			# run clustal for TMP2
			os.system("%s/clustalw /INFILE=TMP2.FA /OUTFILE=TMP2.GDE /OUTPUT=GDE" % \
					   clus_dir)
			
			# can't generate bootstrap tree for 3 taxa
			if len(gdict[i]) == 3:
				oup_tree = open("TMP2.PHB.simplify","w")
				oup_tree.write("(%s,%s,%s);" % \
							   (gdict[i][0],gdict[i][1],gdict[i][2]))
			else:
				os.system("%s/clustalw /INFILE=TMP2.GDE /BOOTSTRAP=100" % \
						  clus_dir)
				tree_util.simplify("TMP2.PHB")
			
			# convert GDE to fa
			inp = open("TMP2.GDE","r")
			oup = open("TMP2.GDE.FA","w")
			tlines = inp.readlines()
			for j in tlines:
				if j[0] == "%":
					oup.write(">%s" % j[1:])
				else:
					oup.write(j)
			oup.close()
			
			# run paml
			self.to_phylip("TMP2.GDE.FA","TMP2")
			os.system("%s/codeml" % paml_dir)
			
			# get the matrix part
			inp = open("TMP2.OUT","r")
			inlines = inp.readlines()
			found = 0
			rdict = {}
			nlist = []
			for j in inlines:
				if j.find("ML distances of aa seqs.") != -1:
					found = 1				
				elif found:
					if j[-2:] == "\r\n":
						j = j[:-2]
					elif j[-1] == "\n":
						j = j[:-1]
						
					jlist = j.split(" ")
					tlist = []
					nlist.append(jlist[0])
					for k in jlist[1:]:
						if k != "":
							tlist.append(k)
					rdict[jlist[0]] = tlist
			
			for j in range(len(nlist)):
				print j
				for k in rdict.keys():
					if nlist[j] != k:
						oup_rate.write("%s\t%s\t%s\n" % (nlist[j],k,rdict[k][j]))
					else:
						del(rdict[k])

			#os.system("rm -f TMP*")
"""
