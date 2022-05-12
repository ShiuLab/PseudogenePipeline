
#!/usr/local/bin/python

#-------------------------------------------------------------------------------
# Project:      pHOG
# Class:	ParseBlast.py
# Desc.:    This class is responsible for talking to the database for cluster
#	   related operations
# Modification history
# 25.05,02
#  Created
# 03,06,02
#  Problem with inconsistent accessions in sputnik and in blast database. 
#  Sputnik delete "." and "-" in accessions then store in pep_sputnik. This
#  created a hugh problem in clustering. Implement check_acc() to check the con-
#  sistency of the parsed file. This checking mechanism shold be implemented
#  in a planned full blast parser.
# 07.06,02
#  Implement get_cluster_scores()
# 15/06,02
#  Implement apply_threshold() and modify get_cluster_scores() so it will
#  calculate the log value unlike that of tribe-parse.pl. The output of
#  get_cluster_score() will be:
#	     [seq1]\t[seq2]\t[evalue]\t[-log10(evalue)]
# 10.07,02
#  Finish get_score_for_npc().
# 18.07,02
#  Memory problem for large matrices in get_score_for_npc(). With > 1M scores,
#  ca. 750MB memory is taken.
# 25,07,02
#  Realize that I got the name wrong. It should be SPC not NPC.
# 02.08,02
#  Implement get_selected_scores().
# 08/17,02
#  Change back to US-based date system. Implement extract_cds()
# 09/03.02
#  Bug in extract_cds(). For mutiple matching area between a query and a subj,
#  the different matching areas are merged and the coordinates are wrong. Fixed.
# 09/08,02
#  Bug in extract_cds(). There are cases where matches can be found in both
#  orientation. In this situation, implement a mechanism that will take only
#  the orientation of the majority
# 09/10,02
#  Modify parse_blast_file() and symmetrify() to allow the use of percentID
#  or percentSIM.
# 09/18,02
#  Some names are longer than 10 char. Need to truncate them and generate an
#  output for a list of modified genes.
# 09/23,02
#  Some self score is less than 100%, somehow I decided that the normalization
#  in symmetrify should use the smaller self score of the pair. Don't understand
#  the logic behind that... The result is that, some negative scores is present 
#  due to that after normalization. To eliminate this artifact and to be 
#  conservative in the estimate, the score used should be the higher one.
# 10/10,02
#  Implement parse_alignment
# 10/19,02
#  Modify symmetrify, for a pair m,n, the score for (m,n) is found but not (n,m)
#  in this case, (n,m) will be set to (m,n) instead of 0.
# 01/15,03
#  How do parse_table deal with situations with multiple matches? Actually, it
#  is symmetrify who make the decision. It takes the top score.
# 02/25,03
#  Bugs in parse_align. Query line disappears. End of line char problem. Fixed.
#  Has to make sure that get_qualified still works.
# 06/30,03
#  Problem with symmetrify. The paires that do not have score need to be 0.
# 08/08,03
#  Implement match_fam
# 08/27,03
#  In symmetrifym cases exist where self score is smaller than the default 
#  cutoff. The self scores were discarded by mistake. Fix this.
# 11/12,03
#  Need to deal with WashU blast output for the first time. Accommondate only
#  in extract_cds at this point. It seems the major difference is the spacing
#  between Sbjct: or Query: to the L coordinates.
# 06/10,04
#  A bug in parse_table. The evalue part is not done. And both evalue and ident
#  -based filterings doesn't take care of multiple hit problem (ie multiple
#  a-b hits). Dealt with.
# 07/09/09
#  Gaurav found couple bugs in Chain.
#-------------------------------------------------------------------------------

import sys, math, SingleLinkage, FastaManager, FileUtility, string, os

print("here")

class parser:

	def __init__(self, dbtask="", config=""):
		self.dbtask = dbtask
		self.config = config	    

	
	def for_mega(self,bdir,fasta):
		
		if bdir[-1] != "/":
			bdir += "/"
		
		print("1. blast...")
		#os.system("%sformatdb -i %s" % (bdir,fasta))
		#os.system("%sblastall -p blastp -d %s -i %s " % (bdir,fasta,fasta) + \
		#	 	  "-o %s.out -F F -v 0 -b 5000 -m 8"  % fasta)
		print("2. deal with blast score...")
		self.parse_table("%s.out" % fasta,"evalue",T=1,wself=1) 	  
		self.symmetrify("%s.out_T1_evalue.parse"    % fasta)
		self.mega_score("%s.out_T1_evalue_c0h0.sym" % fasta)
		
		
	##
	# This function checks if subjects are more closely related to members of
	# a particular family. Comparison is done based on evalue.
	#
	# @param blast  blast output, tabular format
	# @param matrix family assignment with [seq_id][fam]
	# @param target the target family
	# @param unknown the ones with unknown groups are query[1, default] or
	#               subj[0]
	##
	def match_fam(self,blast,matrix,target,unknown=0):
		
		# read matrix into dict
		matrix = futil.file_to_dict(matrix,1)

		inp = open(blast,"r")
		inl = inp.readline()
		# first get top matches, subj as key, [query,pid,evalue] as value
		# if unknown is set at 0. Other wise, query as key.
		sdict = {}
		print("Get top matches...")
		while inl != "":
			llist  = inl.split("\t")
			query  = llist[0]
			subj   = llist[1]
			pid    = float(llist[2])
			evalue = llist[-2]
			changed = 0
			
			if unknown:
				if(query not in sdict or 
				   (query in sdict and pid > sdict[query][1])):
					sdict[query] = [subj,pid,evalue]
					changed = 1			
			else:
				if(subj not in sdict or 
				   (subj in sdict and pid > sdict[subj][1])):
					sdict[subj] = [query,pid,evalue]
					changed = 1
			inl = inp.readline()
		
		# check query group and output subj, query, percentID, evalue
		print("Match targets and generate output...")
		oup = open(blast+"."+target,"w")
		for i in list(sdict.keys()):
			#try:
			if matrix[sdict[i][0]] == target:
				oup.write("%s\t%s\t%i\t%s\n" % \
							(i,sdict[i][0],sdict[i][1],sdict[i][2]))
			#except KeyError:
			#	print "Are you sure you set in unknown flag right??"
			#	print "QUIT!!"
			#	sys.exit(0)
		
		

	##
	# Format e values for Super-Paramagnetic Clustering (SPC). The input file has
	# to conform to the following format:
	#
	#  seq1<space>coord1L<->coord1R<\t>seq2<space>coord2L<->coord2R<\t>e_value
	#
	# Three output files will be generated. The first is ".spc" or ".spc.h", de-
	# pends on whether homogenization is done or notm. It looks the same as input
	# except that:
	#
	#  1. evalue is replaced by normalized, symmetrified score
	#  2. the score is stored in matrix-like format, where score for 1-2 is stor-
	#     ed but not 2-1.
	#  3. The seq name is relaced by index numbers
	#
	# The second file is: seq_name<\t>index
	#
	# The third file is log file with running parameter and a list for sequences
	# not included in the score file due to the rules applied in the codes.
	#
	# @param outbase    output base name
	# @param score_dict from the FastaOutIO.get_scores(). This will be based on
	#		   dissimilarity matrix.
	# @param per_id     number of top score VALUES to retrieve for each gene,
	#		   defailt is 0, which will retrieve all scores. It is
	#		   expected that some will have the same scores.
	# @param matrix     score output in what kind of matrix
	#		    0 - dissimilarity matrix, similar from 0 to dissimilar 1
	#		    1 - similarity matrix, similar from 1 to dissimilar 0
	# @param homogenize to make the scores more linearly distributed. At default
	#		    0 - no homogenization
	#		    1 - homogenized by applying square of the value
	##      
	def get_scores_for_spc(self,score_file,outbase,per_id=10,homogenize=0):

		print("\nGenerate score file for Non-parametric clustering:")

		##
		# Set default values
		##
		MAX_SCORE = 200
		MIN_SCORE = 0
		
		if outbase == "":
			outbase = score_file

		oup3 = open(outbase+"_p%ih%i.log" % (per_id,homogenize),"w")

		oup3.write("ParseBlast.py get_score_for_spc")
		oup3.write(" score_file= %s\n" % score_file)
		oup3.write(" outbase   = %s\n" % outbase)
		oup3.write(" per_id    = %i\n" % per_id)
		oup3.write(" homogenize= %i\n\n" % homogenize)


		print("Convert e value, apply threshold...\n")
		# convert e value into -log10(e)
		inp = open(score_file,"r")

		# idx1+" "+idx2 as key, [len1, len2, -log10(e)] as value
		score_dict = {}

		# idx as key, included in final output [1] or not [0] as value
		name_dict    = {}
		inline   = inp.readline()
		while inline != "":	     
			llist = inline.split("\t")
			idx1  = llist[0]
			len1  = int(idx1[idx1.find("-")+1:]) - \
					int(idx1[idx1.find(" ")+1:idx1.find("-")]) + 1
			idx2  = llist[1]
			len2  = int(idx2[idx2.find("-")+1:]) - \
					int(idx2[idx2.find(" ")+1:idx2.find("-")]) + 1

			if idx1 not in name_dict:
				name_dict[idx1] = 0
			if idx2 not in name_dict:
				name_dict[idx2] = 0					     

			# convert scores to -log10
			score = llist[2][:-1]
			if score[0] == "e":
				try:
					score = -math.log10(float("1"+score))
				except OverflowError:
					print("OVERFLOW:",idx1,idx2,score)
					inline = inp.readline()
					continue
			elif score  == "0.0": 
				score = MAX_SCORE
			else:
				try:
					score = -math.log10(float(score))
				except OverflowError:
					print("OVERFLOW:",idx1,idx2,score)
					inline = inp.readline()
					continue
				except ValueError:
					print("ValueError:",idx1,idx2,score)
					inline = inp.readline()
					continue

			# apply threshold
			if score < MIN_SCORE:
				inline = inp.readline()
				continue

			# homogenize score to its square root
			if homogenize:
				score_dict[idx1+"_"+idx2] = [len1,len2,math.sqrt(score)]
			else:
				score_dict[idx1+"_"+idx2] = [len1,len2,score]
					
			inline = inp.readline()

		print("Symmetrify and normalize scores...\n")
		# symmetrify and normalize score with self score and len of longer entry
		# notice that tscore has different format from that of score_dict
		tscores = {}
		keys = list(score_dict.keys())
		for i in keys:
			idx1    = i[:i.find("_")]
			idx2    = i[i.find("_")+1:]

			# this key may not exist because it is deleted after processing
			# idx2_idx1 previously, if so, skip
			if idx1+"_"+idx2 in score_dict:
				score12 = score_dict[i][2]
			else:
				continue

			# if reciprocal score is not present, this entry is ignored
			if idx2+"_"+idx1 in score_dict:
				len1    = score_dict[i][0]
				len2    = score_dict[i][1]
				score21 = score_dict[idx2+"_"+idx1][2]
			else:
				continue

			# symmetrify
			score = (score12+score21)/2

			# normalize
			if len1 >= len2:
				selfL = score_dict[idx1+"_"+idx1][2]
			else:
				selfL = score_dict[idx2+"_"+idx2][2]
				
			score = score/selfL

			# delete score_dict element, but don't delete self scores
			if idx1 != idx2:
				del score_dict[i]
				del score_dict[idx2+"_"+idx1]
			
			if score > 1:
				print("Score bigger than 1:",idx1,idx2,score)
				score = 1

			# construct transformed score dict
			if idx1 not in tscores:
				tscores[idx1] = [[idx2],[score]]
			else:
				if not idx2 in tscores[idx1][0]:
					tscores[idx1][0].append(idx2)
					tscores[idx1][1].append(score)

			if idx2 not in tscores:
				tscores[idx2] = [[idx1],[score]]
			else:
				if not idx1 in tscores[idx2][0]:
					tscores[idx2][0].append(idx1)
					tscores[idx2][1].append(score)			  

		del score_dict
	
		# reconstruct tscores based on per_id if it is bigger than 0. The self
		# scores will not be included, which is always 1.
		if per_id != 0:

			print("Get top %i scores...\n" % per_id)
			for i in list(tscores.keys()):    # match for each index
				#print i
				
				# map scores, score as key, index of score in the list as value
				map = {}
				for j in range(len(tscores[i][1])):
					if tscores[i][1][j] not in map:
						map[tscores[i][1][j]] = [j]
					else:
						map[tscores[i][1][j]].append(j)

				#print " map_dict:",map

				order = list(map.keys())
				order.sort()
				order.reverse()

				#print " sorted score:",order

				# rank -> order index -> order element -> map key -> map value:
				#  tscore[0] index
				#  tscore[1] index

				ilist = []
				slist = []		      
				for j in range(1,per_id+1):
					if j<len(order):
						
						score   = order[j]
						members = map[score]

						for k in members:
							ilist.append(tscores[i][0][k])
							slist.append(score)
					else:
						break

				tscores[i] = [ilist,slist]


			# check reciprocity, build index dict, and generate output
			print("Check reciprocity, generate output...\n")
			idx_dict   = {}
			count = 0
			
			oup1 = open(outbase+"_p%ih%i.idx" % (per_id,homogenize),"w")
			oup2 = open(outbase+"_p%ih%i.spc" % (per_id,homogenize),"w")

			pairs = {}
			for i in list(tscores.keys()):
				for j in range(len(tscores[i][0])):

					index = tscores[i][0][j]
					score = tscores[i][1][j]

					# check if tscores[index] has i as a score
					if i in tscores[index][0]:
						
						if i not in idx_dict:
							idx_dict[i]  = count
							name_dict[i] = 1
							oup1.write("%s\t%i\n" % (i,count))      
							count = count +1

						if index not in idx_dict:
							idx_dict[index]  = count
							name_dict[index] = 1
							oup1.write("%s\t%i\n" % (index,count))
							count = count+1

						# Output as full name
						#oup2.write("%s\t%s\t%f\n" % (i,
						#							index,
						#							score))

						# Output as index number and as a matrix
						if(i+"_"+index not in pairs and
						   index+"_"+i not in pairs):
							oup2.write("%i\t%i\t%f\n" % (idx_dict[i],
														 idx_dict[index],
														 score))
							pairs[i+"_"+index] = 1
					else:
						continue

			del tscores

			oup3.write("Singletons:\n")
			count = 0
			for i in list(name_dict.keys()):
				if name_dict[i] == 0:
					oup3.write(" "+i+"\n")
					count = count+1
			oup3.write("Total = %i" % count)
			
			print("Done!\n")
			sys.exit(0)
		

	##
	# Take the output of parse_blast_file() or parse_blast_db(), perform some
	# operation to get MCL input file. Notice that NO normalization is done over
	# here.
	#
	# @param outbase    output base name
	# @param score_file from the parse_db() or parse_file.
	# @param cutoff     evalue cutoff
	# @param homogenize to make the scores more linearly distributed. At default
	#		    0 - no homogenization
	#		    1 - homogenized by applying square of the value
	##
	def get_scores_for_mcl(self,score_file,outbase,cutoff=1,homogenize=0):

		if outbase == "":
			outbase = score_file

		inp = open(score_file,"r")
		oup = open(outbase+".idx","w")

		# construct index
		idx = {}
		c = 0
		o_scores = {}

		print("Construct idx and score dict...")
		inline = inp.readline()
		countS = 0
		while inline != "":
			if countS % 100000 == 0:
				print(" %i x 100k" % (countS/100000))
			countS += 1
			lnlist = inline.split("\t")
			# check number of column
			if len(lnlist) != 3:
				print("Score format problem: should be [query][subj][-log(e)]")
				print("Quit!")
				sys.exit(0)
			
			if lnlist[0] not in idx:
				idx[lnlist[0]] = c
				oup.write("%i\t%s\n" % (c,lnlist[0]))
				c = c+1
			if lnlist[1] not in idx:
				idx[lnlist[1]] = c
				oup.write("%i\t%s\n" % (c,lnlist[1]))
				c = c+1				 

			if lnlist[0] not in o_scores:
				# apply cutoff
				if float(lnlist[2]) >= cutoff:
					o_scores[lnlist[0]] = {lnlist[1]:float(lnlist[2])}
			else:
				# apply cutoff
				if float(lnlist[2]) >= cutoff:
					# more than one scores, take the largest one
					if lnlist[1] in o_scores[lnlist[0]]:
						if float(lnlist[2]) > o_scores[lnlist[0]][lnlist[1]]:
							o_scores[lnlist[0]][lnlist[1]] = float(lnlist[2])
					else:
						o_scores[lnlist[0]][lnlist[1]] = float(lnlist[2])

			inline = inp.readline()

		##
		# SET DIMENSION
		##
		D = c
		
		inp.close()
		oup.close()
		
		# symmetrify scores
		s_scores = {}
		okeys = list(o_scores.keys())
		print("Symmetrify scores...")
		countS = 0
		for i in okeys:
			if countS % 1e4 == 0:
				print(" %i x 10k" % (countS/1e4))
			countS+= 1
			oikeys = list(o_scores[i].keys())
			for j in oikeys:
				#print "",idx[j] 
				# check backward score, if absent, this match will be discarded
				score = 0

				if j in o_scores and i in o_scores[j]:
					if o_scores[i][j] == o_scores[j][i]:
						score = o_scores[i][j]
					else:
						score = (o_scores[i][j]+o_scores[j][i])/2.0

					del o_scores[j][i]
					if i != j:
						del o_scores[i][j]
				
				# apply cutoff, homogenize if necessary
				if score !=  0:				 
					if homogenize:
						score = math.sqrt(score)

					if idx[i] in s_scores:
						if idx[j] not in s_scores[idx[i]]:
							#print "  add1:",idx[j],"to",idx[i],score
							s_scores[idx[i]][idx[j]] = score
					else:
						#print "  add2:",idx[j],"to",idx[i],score
						s_scores[idx[i]] = {idx[j]:score}
						

					if idx[j] in s_scores:
						if idx[i] not in s_scores[idx[j]]:
							#print "  add3:",idx[i],"to",idx[j],score
							s_scores[idx[j]][idx[i]] = score
					else:
						#print "  add4:",idx[i],"to",idx[j],score
						s_scores[idx[j]] = {idx[i]:score}
				#else:
				#       print "  below_cutoff"
		o_scores = {}
		
		#for i in s_scores:
		#	print i
		#	print s_scores[i]
		# write mcl file
		print("Generate output file...")
		oup = open(outbase+".mci","w")
		oup.write("(mclheader\nmcltype matrix\ndimensions %ix%i\n)\n" % (D,D))
		oup.write("(mclmatrix\nbegin\n")
		
		for i in list(s_scores.keys()):			       
			oup.write("%i " % i)
			for j in s_scores[i]:
				s = str(s_scores[i][j])
				oup.write("%i:%s " % (j,s[:s.find(".")+4]))
				#scores = str(s_scores[i][1][j])
				#out_str = out_str + "%i:%s " % (s_scores[i][0][j],
				#								scores[:scores.find(".")+4])
			oup.write("$\n")

		# add singletons
		for i in range(D):
			if i not in s_scores:
				score = 200.0
				if homogenize:
					score = math.sqrt(score)
				oup.write("%i %i:%f $\n" % (i,i,score))

		oup.write(")")
		
		print("Close outputstream...")
		oup.close()
		
				  

	##
	# Common operation on scores:
	#  1. symmetrification - so Sab = Sba
	#  2. cutoff	   - apply cutoff
	#  3. normalization    - normalize by smaller self
	#  4. homogenization   - apply sqare on normalized -log(E) value
	#
	# Assume the score is in the last token of the line
	# 
	# @param  score_file   file generated by parse_blast_file or parse_blast_db
	# @output xxx_cxhx.sym file name include the cutoff value and homogenization
	#		      option. Contains [seq1][seq2][symmed score]
	##
	def symmetrify(self,score_file,outbase="",cutoff=0,homogenize=0):
		
		oup_log = open(score_file+"_sym.log","w")
		
		print("Construct score dict...")
		o_scores = {}
		inp = open(score_file,"r")
		inline = inp.readline()

		c = 0
		ndict = {}
		print("Processed lines:")
		oup_log.write("More than one scores:\n")
		while inline != "":
			if c % 100000 == 0:
				print(" %ix100k" % (c/100000))
			c += 1
				
			# [seq1][seq2][E][-log(E)] or [seq1][seq2][-][percentID or SIM]
			lnlist = inline.split("\t")

			# assuming score is in the last token
			lnlist[-1] = self.rmlb(lnlist[-1])
			
			# store id into dict
			if lnlist[0] not in ndict:
				ndict[lnlist[0]] = 0
			if lnlist[1] not in ndict:
				ndict[lnlist[1]] = 0
			
			if lnlist[0] not in o_scores:
				o_scores[lnlist[0]] = {lnlist[1]:float(lnlist[-1])}
			else:
				# more than one scores, take the largest one
				if lnlist[1] in o_scores[lnlist[0]]:
				
					oup_log.write(" %s,%s %f -> %s\n" % \
										(lnlist[0],
										 lnlist[1],
										 o_scores[lnlist[0]][lnlist[1]],
										 lnlist[-1]))
					if float(lnlist[-1]) > o_scores[lnlist[0]][lnlist[1]]:
						o_scores[lnlist[0]][lnlist[1]] = float(lnlist[-1])											   
				else:
					o_scores[lnlist[0]][lnlist[1]] = float(lnlist[-1])
			inline = inp.readline()
		
		print("Apply cutoff...")
		okeys = list(o_scores.keys())
		selfC = 0
		oup_log.write("\nSelf score below cutoff at %i:\n" % cutoff)
		for i in okeys:
			jkeys = list(o_scores[i].keys())
			for j in jkeys:
				if o_scores[i][j] < cutoff:
					# for log
					if i == j:      
						oup_log.write(" %s\t%f\n" % (i,o_scores[i][i]))
						selfC += 1
					o_scores[i][j] = cutoff

		print(" %i self score below cutoff %i" % (selfC,cutoff))
		
		print("Symmetrify and normalize scores...")
		s_scores = {}
		print("Process total %i taxa:" % len(okeys))
		c = 0
		oup_log.write("\nNot in score dict, most likely not in BLAST file:\n")
		for i in okeys:		 
			if c % 10 == 0:
				print("",c+1)
			c += 1
			oikeys = list(o_scores[i].keys())
			for j in oikeys:
				score = 0
				if j in o_scores:
					if i in o_scores[j]:
						# at the moment, score is -log(E), except if absent
						if o_scores[i][j] == o_scores[j][i]:
							score = o_scores[i][j]
						else:
							score = (o_scores[i][j]+o_scores[j][i])/2.0
					else:
						o_scores[j][i] = o_scores[i][j]
						score = o_scores[i][j]
				t_score = score
				
				# check if scores are missing, if so, add self score as cutoff
				if i in o_scores:
					if i in o_scores[i]:
						pass
					else:
						print("Self score missing:",i)
						print("Did you parse the table -wself 1? QUIT!")
						sys.exit(0)
				else:
					print(" score dict misses:",i)
					oup_log.write(" %i" % i)
					o_scores[i] = {i:cutoff}
					#print "Quit!"
					#sys.exit(0)
				
				if j in o_scores:
					if j in o_scores[i]:
						pass
					else:
						print("Self score missing:",j)
						print("Did you parse the table -wself 1? QUIT!")
						sys.exit(0)
				else:
					print(" score dict misses:",[j])
					oup_log.write(" %s" % j)
					o_scores[j] = {j:cutoff}
					#print "Quit!"
					#sys.exit(0)				
					
				# normalize, make sure not divided by zero
				if o_scores[i][i] >= o_scores[j][j]:
					if o_scores[i][i] == 0:
						score = 1.0
					else:
						score = 1.0 - score/o_scores[i][i]
				else:
					if o_scores[j][j] == 0:
						score = 1.0
					else:
						score = 1.0 - score/o_scores[j][j]
					
				# homogenize by 
				if homogenize > 0 and score != 0:
					if homogenize == 1:
						score = score*score      # multiplication
					elif homogenize == 2:
						score = math.sqrt(score) # square root

				if i in s_scores:
					if not j in s_scores[i][0]:
						s_scores[i][0].append(j)
						s_scores[i][1].append(score)
				else:
					s_scores[i] = [[j],[score]]

				if j in s_scores:
					if not i in s_scores[j][0]:
						s_scores[j][0].append(i)
						s_scores[j][1].append(score)
				else:
					s_scores[j] = [[i],[score]]

				#print i,j,t_score,score,o_scores[i][i],o_scores[j][j]

		if outbase == "":
			outbase = score_file[:score_file.rfind(".")]
		outbase = "%s_c%ih%i" % (outbase,cutoff,homogenize)
		
		oup = open(outbase+".sym","w")
		print("Generate output: %s" % (outbase+".sym"))
		for i in list(s_scores.keys()):
			for j in range(len(s_scores[i][0])):
				oup.write("%s\t%s\t%f\n" % (i,
											s_scores[i][0][j],
											s_scores[i][1][j]))							
		print("Closing output stream...")
		oup.close()
		
	
	
	def sym2(self,score):
		
		# Not done!
		
		inp = open(score)
		inl = inp.readline()
		S   = {}
		while inl != "":
			L = inl.split("\t")
			if L[0] not in S:
				if L[1] not in S:
					S[L[0]] = {L[1]:float(L[-1])}
					S[L[1]] = {L[0]:float(L[-1])}
				elif L[0] not in S[L[1]]:
					S[L[0]] = {L[1]:float(L[-1])}
					S[L[1]] = {L[0]:float(L[-1])}
				else:
					S[L[0]][L[1]] = (S[L[0]][L[1]]+	float(L[-1]))/2.0
					S[L[1]][L[0]] = (S[L[1]][L[0]]+	float(L[-1]))/2.0
									
				
			elif L[1] not in S[L[0]]:
				S[L[0]][L[1]] = float(L[-1])
			
			
			if L[1] not in S:
				S[L[1]] = {L[0]:float(L[-1])}
			elif L[0] not in S[L[1]]:
				S[L[1]][L[0]] = float(L[-1])
			
			
			
			inl = inp.readline()
							

	##
	# This will get the score for use in the neighbor program from Phylip. The 
	# format is:
	#
	#	number of taxa
	# seq1   self  S12  S13  S14...
	# seq2   S12   self S23  S24...
	# ...
	#
	# The score matrix is similarity matrix. Self is 0. So things need to be
	# normalized somehow. The scores also need to be symmetrified. The method
	# symmetrify will generate score suitable.
	#
	# Only the first 6 char of the score will be taken so the precision to: 0.0001
	#
	# @param  sym_score  symmetrified score file generated by symmetrify
	# @output xxx.neighbor  As explained above.
	##
	def get_scores_for_neighbor(self,sym_score):

		# generate a nested dict for score
		# generate a dict for index
		
		scores = {}
		index  = {}	     
		inp = open(sym_score,"r")
		inline = inp.readline()
		print("Generate score dict...")
		c = 0
		while inline != "":
			if c % 100000 == 0:
				print(" %ix100k" % (c/100000))
			c += 1

			llist = inline.split("\t")
			llist[2] = llist[2][:6]
			
			if llist[0] not in index:
				index[llist[0]] = 1
			if llist[1] not in index:
				index[llist[1]] = 1			     

			if llist[0] not in scores:
				scores[llist[0]] = {llist[1]:llist[2]}
			else:
				scores[llist[0]][llist[1]] = llist[2]
					   
			if llist[1] not in scores:
				scores[llist[1]] = {llist[0]:llist[2]}
			else:
				scores[llist[1]][llist[0]] = llist[2]			   
			inline = inp.readline()

		oup  = open(sym_score+".neigh","w")
		sline = ""	      
		ikeys = list(index.keys())
		oup.write("     %i\n" % len(ikeys))
		print("Write scores...")
		for i in list(index.keys()):
			#print i
			count = 0
			for j in ikeys:
				count = count + 1
				if j in scores[i]:
					sline = sline + scores[i][j] + "  "
				else:
					sline = sline + "1.0000  "		      
			if len(i)>10:
				print("The identifier is longer than 10 characters.")
				print("Should run index_names\nExit!")
				sys.exit(0)
			else:
				oup.write(i+(10-len(i))*" "+sline[:-2]+"\n")
			sline = ""
		

	##
	# This will generate a score matrix for MEGA. The format is like:
	#
	# #mega
	# ! Title: something;
	# ! Format DataType=distance; 
	# 
	# #taxa1
	# #taxa2
	# #...
	# 
	#  d12 d13 d14 ...
	#      d23 d24 ... 
	#	      ...
	#
	# Input file should have 3 or more columns:
	#     Seq1<\t>Seq2<\t>....<\t>Distance
	# Distance has to be between 0 and 1.
	##
	def mega_score(self,sym_score):

		# generate a nested dict for score
		# generate a dict for index
		
		scores = {}
		index  = {}	     
		inp = open(sym_score,"r")
		inline = inp.readline()
		print("Generate score dict...")
		c = 0
		while inline != "":
			if c % 100000 == 0:
				print(" %ix100k" % (c/100000))
			c += 1

			llist = self.rmlb(inline).split("\t")
			#llist[2] = llist[2][:6]
			
			if llist[0] not in index:
				index[llist[0]] = 1
			if llist[1] not in index:
				index[llist[1]] = 1			     

			if llist[0] not in scores:
				scores[llist[0]] = {llist[1]:llist[-1]}
			else:
				scores[llist[0]][llist[1]] = llist[-1]
					   
			if llist[1] not in scores:
				scores[llist[1]] = {llist[0]:llist[-1]}
			else:
				scores[llist[1]][llist[0]] = llist[-1]			   
			inline = inp.readline()

		oup  = open(sym_score+".meg","w")
		oup.write("#mega\n!Title %s;\n"%sym_score + \
			  "!format DataType=distance DataFormat=upperright;\n\n")
		
		print("Write taxa...")
		ikeys = list(index.keys())
		ikeys.sort()
		for i in ikeys:
			oup.write("#%s\n" % i)
		oup.write("\n")
		
		# don't need self, only upper half matrix
		print("Write scores...")
		sline = ""	      
		for i in range(len(ikeys)):
			#print i
			for j in range(i+1,len(ikeys)):
				if ikeys[j] in scores[ikeys[i]]:
					sline = sline + scores[ikeys[i]][ikeys[j]] + "  "
				else:
					sline = sline + "1.0000  "		      

			oup.write(" "+sline[:-2]+"\n")
			sline = ""
		
		print("Close output stream...")
		oup.close()
		   

	#
	# Take symmetrified score file and generate a score matrix.
	#
	def score_matrix(self, score):
		
		print("Read scores into a dictionary...")
		D = {} #{id1:{id2:score}}
		inp = open(score)
		inl = inp.readline()
		while inl != "":
			[id1,id2,S] = inl.strip().split("\t")
			if id1 not in D:
				D[id1] = {id2:S}
			elif id2 not in D[id1]:
				D[id1][id2] = S
			if id2 not in D:
				D[id2] = {id1:S}
			elif id1 not in D[id2]:
				D[id2][id1] = S
			inl = inp.readline()
		
		print("Write matrix...")
		ids = list(D.keys())
		ids.sort()
		oup = open(score+".matrix","w")
		oup.write("\t%s\n" % "\t".join(ids))
		for i in ids:
			oup.write("%s" % i)
			for j in ids:
				if i == j:            # Same sequence
					oup.write("\t0")    
				elif j not in D[i]:   # Score not there due to thresholding
					oup.write("\t1")
				else:
					oup.write("\t%s" % D[i][j])
			oup.write("\n")
		
		
	
	##
	# This will parse the scores out of the alignment part of blast output. The
	# output from this program is:
	#
	#   [seq1][seq2][E][-log(E)]
	#
	# @param blastout  output file from blast
	# @param T	 threshold in the form of E value or percentage
	# @param target    evalue [default], percentID, percentSIM
	##
	def parse_blast_file(self,blastout,T,target):
		
		inp = open(blastout,"r")
		oup = open(blastout+"_T%i_%s.parse" % (T,target),"w")
		inline = inp.readline()
		
		query = sbjct = ""
		evalue = 10.0
		percent = 0
		match   = "-"
		got_e = 0
		got_p = 0
		#countQ = 0
		
		print("Query processed:")
		while inline != "":
			if inline[-2:] == "\r\n":
				inline = inline[:-2]
			elif inline[-1] == "\n":
				inline = inline[:-1]
			
			if len(inline) > 7 and inline[:7] == "Query= ":
				query = inline[7:]
				print(" ",query)

			elif len(inline) > 1 and inline[0] == ">":
				sbjct = inline[1:]
			
			# parse evalue
			elif(len(inline) > 6 and target == "evalue" and \
				 inline[:6] == " Score"):
				evalue = inline[inline.find("Expect = ")+9:]				  
				if evalue[0] == "e":
					evalue = "1"+evalue
				ef = float(evalue)						      
				# will not deal with evalue >= T
				if float(evalue) >= T:
					inline = inp.readline()
					continue
				if evalue == "0.0":
					evalue = "1e-200"		       
				got_e = 1
			
			elif((target == "percentID" or target == "percentSIM") and \
				 inline.find("Identities =") != -1):		    

				lnlist = inline.split(",")

				# the following line is only good for polypeptide blast
				#if len(lnlist) < 2:
				#       print "Missing percentID or percentSIM:",lnlist
				#       inline = inp.readline()
				#       continue
					
				if target == "percentID":
					match   = lnlist[0][lnlist[0].find("=")+1:\
										lnlist[0].find("(")]
					match   = match[match.find(" ")+1:match.rfind(" ")]
					percent = lnlist[0][lnlist[0].find("(")+1:\
										lnlist[0].find("%")]
				elif target == "percentSIM":
					percent = lnlist[1][lnlist[1].find("(")+1:\
										lnlist[1].find("%")]

				if float(percent) <= T:
					inline = inp.readline()
					continue

				got_p = 1
				
			# generate output accordingly
			if target == "evalue" and got_e:
				try:
					convert = -math.log10(float(evalue))
					oup.write("%s\t%s\t%s\t%f\n" % (query,sbjct,evalue,convert))
				except ValueError:
					print("Value can't be transformed:",evalue)
				got_e = 0
	
			elif(target == "percentID" or target == "percentSIM") and got_p:	
				oup.write("%s\t%s\t%s\t%s\n" % (query,sbjct,match,percent))
				got_p = 0
				
			inline = inp.readline()

		inp.close()
		oup.close()

	##
	# This is aimed at parsing the tabular output with large pieces of DNA
	# as BLAST subjects, such as chromosome. Match to query is parsed to get
	# the coord and orientation. The coords are merged between overlapping
	# subject regions and an output is generated with the format:
	#    [chr][id][ori][L][R]
	#
	# @param blast   blast output with -m 8 or 9
	# @param feature the sequence features used as query sequences
	##
	def merge_match(self,blast,feature,mtype="blast_match"):
		
		
		# scan for subj
		inp = open(blast,"r")
		inline = inp.readline()
		slist = []
		while inline != "":
			llist = inline.split("\t")
			if llist[1] not in slist:
				slist.append(llist[1])
			inline = inp.readline()
		
		# 1. Construct a dict for each subj. With the L coord as key, 
		#    [R,ori] as values 
		# 2. Sort keys then check for overlap
		oup = open(blast+".merged.coord","w")
		print("Process subj:")
		for i in slist:
			print("",i)
			
			cdict  = {}
			inp    = open(blast,"r")
			inline = inp.readline()
			while inline != "":
				llist = inline.split("\t")
				# only process the lines of a particular subj
				if i == llist[1]:
					if int(llist[8]) < int(llist[9]):
						sL  = int(llist[8])
						sR  = int(llist[9])
						ori = "W"
					else:
						sR  = int(llist[8])
						sL  = int(llist[9])
						ori = "C"
					
					# if two entries have exactly the same L coord, the longer 
					# one is taken
					if sL in cdict:
						if sR > cdict[sL][0]:
							cdict[sL] = [sR,ori]
					else:
						cdict[sL] = [sR,ori]
				
				inline = inp.readline() 
			
			
			# now merge coords
			ckeys = list(cdict.keys())
			ckeys.sort()
			coord_list = []
			last = -1
			not_added = 0
			for j in ckeys:
				#print "coord:",i
				# the first entry, so set the last examined to it
				if last == -1:
					last = j

				# not the first, check if overlap
				else:
					# the right coord of last is further downstream of the
					# current left, expand the last, delete current
					if cdict[last][0] > j:
						#print "overlap"
						# check ori, just print a warning
						if cdict[last][1] != cdict[j][1]:
							print("  Discrepancy: %s/%s,%s <-> %s/%s,%s" % \
									 (cdict[last][1],last,cdict[last][0],
									  cdict[j][1],i,cdict[j][0]))
						cdict[last][0] = cdict[j][0]
						
						# check if i is the last one, if so, since its coord
						# is incoporated into last, output last
						if ckeys[-1] == j:
							coord_list.append([cdict[last][1],last,\
											   cdict[last][0]])
						
					# the last entry is disjunct from the current, store it
					else:
						coord_list.append([cdict[last][1],last,cdict[last][0]])
						if ckeys[-1] == j:
							coord_list.append([cdict[j][1],j,cdict[j][0]])
						else:
							last = j
			
			# write lines		
			for j in coord_list:
				# subj,feature_id,ori,L,R
				fid = "%s_%s_%i_%i" % (feature,j[0],j[1],j[2])
				oup.write("%s\t%s\t%s\t%s\t%s\t%s\n" % \
						  (i,fid,mtype,j[0],j[1],j[2]))
			
			
	##
	# Parse the blast tabular output from invoking the option -m 9
	#
	# Header format:
	#
	# # BLASTP 2.2.4 [Aug-26-2002]
	# # Query: ensmusp00000042149 S_TKc 77-343
	# # Database: test
	#	        0	      1	          2	          3		            4
	# # Fields: Query id, Subject id, % identity, alignment length, mismatches,
	#	   5	         6	       7       8	 9       10
	#	   gap openings, q. start, q. end, s. start, s. end, e-value,
	#	   11
	#	   bit score
	#
	# @param target specified the kind of stuff to get can be:
	#		    evalue [query][subj][e-value][-log(e)]
	#           convert [query][subj][-log(e)]
	#		    percentID [aln len][% id]
	#		    or specfy the tokens separated by ","
	# @param T         threshold filter for:
	#		    evalue if -log(e) is smaller, it will be discarded
	#		    percentID if ID is smaller, it will be discarded
	# @param verbose   also output qL, qR, sL, sR
	# @param wself     default is no [0], but for clustering, self score should
	#                  be included.   
	# @param lenT      % length threshold. Should match the q or s sequence for
	#                  at least certain length
	# @param QorS      apply lenT based on query or subj length. Q [0], S [1]
	## 
	def parse_table(self,blast,target,T,verbose=0,wself=0,lenT="",fasta="",QorS=0):

		print("Parse blast table:")
		print(" Table name:",blast)
		print(" Target    :",target)
		print(" Threshold :",T)
		print(" Verbose   :",verbose)
		print(" Wself     :",wself)
		print("")
		
		inp = open(blast,"r")
		oup = open(blast+"_T%i_%s.parse" % (T,target),"w")
		
		sizes  = {}
		minLen = 100
		if lenT != "":
			lenT = float(lenT)/100.0
			print("lenT      :",lenT)
			print("Fasta     :",fasta)
			print("QorS      :",QorS)
			if fasta == "":
				print("If you set lenT, you need the sequence file!")
				sys.exit(0)
			else:
				sizes = fmanager.get_sizes(fasta,1)

		inline = inp.readline()
		c = 0
		countW = 0
		sdict  = {} # store the idx of those with scores taken
		while inline != "":
			if c%10000 == 0:
				print(" %i x10k" % (c/10000))
			c += 1  
			if inline[0] != "#":
				llist = inline.split("\t")			
				out_str = ""
				qualify = 1
				# has % length threshold passed
				if lenT != "":
					qMatch = float(llist[7]) - float(llist[6])
					sMatch = float(llist[9]) - float(llist[8])
					qLen   = float(sizes[llist[0]])
					sLen   = float(sizes[llist[1]])				
					# too short
					if qLen < minLen or sLen < minLen:
						qualify = 0			
					# matching parts too short
					if ((QorS == 0 and qMatch/qLen < lenT) or
	 			        (QorS == 1 and sMatch/sLen < lenT)):
						qualify = 0
					
				if target in ["evalue","convert"] and qualify:
					evalue = llist[10]
					if evalue[0] == "e":
						evalue = "1"+evalue
					if evalue == "0.0":
						evalue = "1e-200"
					try:
						ef = -math.log10(float(evalue))
					except OverflowError:
						print("Overflow:",evalue)
						sys.exit(0)
					# will regard evalue >= T as T
					#print llist[0],llist[1],float(evalue),ef,T
					if ef > 0 and ef > T and \
						"%s%s" % (llist[0],llist[1]) not in sdict:
						#print llist[0],llist[1],ef,T
						out_str = "%s\t%f\t%s\t%s\t%s\t%s" % \
								 (evalue,T,llist[6],llist[7],llist[8],llist[9])
						sdict["%s%s" % (llist[0],llist[1])] = 1
					else:
						inline = inp.readline()
						continue
						
					if verbose:
						out_str = "%s\t%f\t%s\t%s\t%s\t%s" % \
								 (evalue,ef,llist[6],llist[7],llist[8],llist[9])
					elif target == "convert":
						out_str = "%f" % ef
					else:
						out_str = "%s\t%f" % (evalue,ef)
					
				elif target == "percentID" and qualify:
					ident = float(llist[2])
					# apply threshold
					if ident < T or "%s%s" % (llist[0],llist[1]) in sdict:
						inline = inp.readline()
						continue
					match = int(float(llist[3])*ident/100)
					sdict["%s%s" % (llist[0],llist[1])] = 1
						
					if verbose:
						out_str = "%i/%s\t%f\t%s\t%s\t%s\t%s" % (match,llist[3],ident,llist[6],llist[7],llist[8],llist[9])
					else:	
						out_str = "%i/%s\t%f" % (match,llist[3],ident)
				elif qualify:
					print("HERE")
					tlist = target.split(",")
					for j in tlist:
						out_str = out_str + llist[int(j)] + "\t"
					out_str = out_str[:-1]
				
				if ((not wself and llist[0] != llist[1]) or wself) and qualify:
					oup.write("%s\t%s\t%s\n" % (llist[0],llist[1],out_str))
					countW += 1
					
			inline = inp.readline()

		print("Total %i pairs, %i qualified and written" % (c,countW))
		
		

	##
	# Parse blast table output based on passed thresholds. Output is a summary
	# of sequences with: [query][qualified_subj_count][subj1][subj2]...
	#
	# @param blast blast tabular output
	# @param E     E value threshold
	# @param H     length threshold
	# @param I     identity threshold
	##
	def parse_table2(self,blast,E,I,H):
		
		print("Check inconsistency, E, I defaults are set to zero")
		
		inp = open(blast,"r")
		inl = inp.readline()
		qdict = {}
		print("Parse blast output:")
		while inl != "":
			if inl[0] != "#":
				L = inl.split("\t")
				try:
					evalue = float(L[-2])
				except IndexError:
					print("Evalue ERR:",L)
					sys.exit(0)
				if evalue == 0:
					evalue = 200
				else:
					evalue = -math.log(evalue,10)
				ident  = float(L[2])
				length = int(L[3])
				
				if evalue >= E and ident >= I and length >= H:
					if L[0] in qdict:
						if L[1] not in qdict[L[0]]:
							qdict[L[0]][0] += 1
							qdict[L[0]][1].append(L[1])
					else:
						#print "",L[0]
						qdict[L[0]] = [1,[L[1]]]
			inl = inp.readline()

		print("Generate output...")
		oup = open(blast+".E%iI%iL%i.parsed" % (int(E),int(I),H),"w")
		for i in qdict:
			oup.write("%s\t%i\t%s\n" % (i,qdict[i][0],
									  string.joinfields(qdict[i][1],",")))
		


		
	#
	# Simply method, no storage involved, parse blast output of style 8.
	# Output as [query][subjt][-log(e)]
	#
	def parse_table3(self,blast,E,I):
		inp = open(blast,"r")		
		oup = open("%s_E%i_I%i.out" % (blast,E,I),"w")
		inl = inp.readline()
		print("Parse blast output:")
		c = 0
		countE = 0
		while inl != "":
			if c%10000 == 0:
				print(" %i x10k" % (c/10000))
			c += 1
			L = inl.split("\t")
			evalue = L[-2]
			try:
				ident  = float(L[2])
			except ValueError:
				# problem with line format
				print("ERR:",L)
				inl = inp.readline()
				continue
			if evalue[0] == "e":
				evalue = "1"+evalue
			if evalue == "0.0":
				evalue = "1e-200"
			try:
				ef = -math.log10(float(evalue))
			except OverflowError:
				print("Overflow:",evalue)
				sys.exit(0)
			
			if ef >= E and ident > float(I):
				oup.write("%s\t%s\t%f\n" % (L[0],L[1],ef))
			else:
				#print "%s\t%s\t%f\t%i" % (L[0],L[1],ef,ident)
				countE += 1
			inl = inp.readline()
		
		print("Total %i scores, %i eliminated" % (c,countE))
			
		
	
	
	def parse_table4(self,blast,E,I,H):
		inp = open(blast)
		oup = open(blast+".E%iI%iL%i.parse4" % (int(E),int(I),H),"w")
		inl = inp.readline()
		countT = countQ = 0
		print("BLAST file :", blast)
		print("Output file:", blast+".E%iI%iL%i.parse4" % (int(E),int(I),H))
		c = 0
		while inl != "":
			if inl[0] != "#":
				L = inl.split("\t")
				if c % 1e5 == 0:
					print(" %i x 100k" % (c/1e5))
				c += 1
				try:
					evalue = float(L[-2])
				except IndexError:
					print(L)
					sys.exit(0)
				if evalue == 0:
					evalue = 400
				else:
					evalue = -math.log(evalue,10)
				ident  = float(L[2])
				length = int(L[3])
				
				if evalue >= E and ident >= I and length >= H:
					oup.write(inl)
					countQ += 1
				countT += 1
			inl = inp.readline()
		
		print("Total %i lines, %i qualified" % (countT,countQ))
		
		
		
	#
	# For entries with -log(e) > T, only one of their score is kept. The rest
	# are thrown out. This function is useful for clustering purpose. No reason
	# to cluster very similar sequences. But be careful with the T setting.
	# Tho sequence id for those which got thrown out are stored for adding
	# back later.
	#
	def merge(self,score,T):
		inp = open(score,"r")
		oup = open("%s_T%i.merged" % (score,T),"w")
		oup2= open("%s_T%i.merged.idx" % (score,T),"w")
		inl = inp.readline()
		print("Merge blast score:")
		c = 0
		
		# the seq id of the kept score as key, a dict as value with
		# the seq id thrown out as key, "" as value
		mdict = {}
		# the seq id thrown out as key, "" as value
		tdict = {} 
		while inl != "":
			if c%10000 == 0:
				print(" %i x10k" % (c/10000))
			c += 1  
			L = inl.split("\t")
			if float(L[-1]) >= T:
				pass
				# MOT FINISHED. 
			
			inl = inp.readline()
		
		
	##
	# This will extract blast score out of database, 1000 at a time, store in a
	# file for use by get_scores_for_mcl(), in the format:
	#
	# id_1 <\t> id_2 <\t> evalue <\t> -log10(evalue)
	#
	# eg. At1g01080  vs  At3g16380 with evalue 1e-08
	#     At1g01080 <\t> At3g16380 <\t> 1e-08 <\t> 8
	#
	##
	def parse_blast_db(self,outbase,T):
		
		# get number of records
		qTuple = self.dbtask.select("SELECT count(*) FROM %s" % \
									self.config["table1"])

		print("Parse blast score from %i queries:" % qTuple[0][0])

		if outbase == "":
			outbase = self.config["table1"]
		
		oup = open(outbase+".parse","w")
		
		has_more = 1
		c = 0

		INCREMENT = 100
		
		while has_more:
			print(" ",c*INCREMENT+1,"-",(c+1)*INCREMENT)

			qTuple = self.dbtask.select(\
				"SELECT data FROM %s WHERE " % self.config["table1"] +\
				"sputnik_ref >= %i AND " % (c*INCREMENT+1) +\
				"sputnik_ref <  %i" % ((c+1)*INCREMENT))

			if qTuple == []:
				has_more = 0
				print(" done...")
			else:

				for j in qTuple:
		
					llist = j[0].split("\n")
					query = sbjct = ""
					evalue = 10.0
					got_e = 0
					for k in llist:
						if len(k) > 7 and k[:7] == "Query= ":
							query = k[7:]
						if len(k) > 1 and k[0] == ">":
							sbjct = k[1:]
						if len(k) > 6 and k[:6] == " Score":
							evalue = k[k.find("Expect = ")+9:]
							
							if evalue[0] == "e":
								evalue = "1"+evalue
							ef = float(evalue)
							
							# will not deal with evalue >= T
							if float(evalue) >= T:
								continue
							if evalue == "0.0":
								evalue = "1e-200"					       
							got_e = 1

						if got_e:
							#print query,sbjct,evalue
							convert = math.fabs(math.log10(float(evalue)))
							oup.write("%s\t%s\t%f\n" % \
									  (query,sbjct,convert))
							got_e = 0
							
			# increment while loop
			c = c+1

		oup.close()

	##
	# The blast output generated using sputnik is problematic. The query
	# sequence, if with "." or "-", the name will be changed to into one with
	# out which is inconsistent with the database it blast against. The problem
	# is then exactly the same squence will not only represent twice but also
	# have asymmetrical scores.
	#
	# @param acc_file    the original fasta flatfile with correct accession
	# @param parsed_file the file generated by parse_db() or parse_file(). The
	#		    second accession should have no problem.
	##
	def check_acc(self,acc_file,parsed_file):

		# constructing acc_dict, sputnik_variant as key, real acc as value
		print("\nConstructing acc_dict...\n")
		inp = open(acc_file,"r")
		inline = inp.readline()

		acc_dict = {}
		while inline != "":
			if inline[0] == ">":

				c_acc = inline[1:-1]

				# constructing accessions formated by sputnik
				# only modified one will be in the dictionary
				i_acc = ""		      
				if "." in c_acc:
					i_acc = c_acc[:c_acc.find(".")] + c_acc[c_acc.find(".")+1:]
					
				elif "-" in c_acc:
					i_acc = c_acc[:c_acc.find("-")] + c_acc[c_acc.find("-")+1:]

				acc_dict[i_acc] = c_acc
				
			inline = inp.readline()

		# compare acc_dict against parsed_file, regenerate output
		print("Checking parsed file...\n")
		inp = open(parsed_file,"r")
		oup = open(parsed_file+".chk","w")
		inline = inp.readline()

		count = 0
		while inline != "":
			
			b_point = inline.find("\t")
			
			acc1 = inline[:b_point]
			if acc1 in acc_dict:
				oup.write(acc_dict[acc1]+inline[b_point:])
			else:
				oup.write(inline)
			
			inline = inp.readline()
			count = count+1
			if count % 100000 == 0:
				print(" %ik" % (count/1000))

	##
	# Check if query or subject is missing in the BLAST output
	## 
	def check_missing(self,fasta,blast,qors):
		
		fnames = fmanager.get_names(fasta,1)
		FD = {}
		for i in fnames:
			FD[i] = 0
		
		inp = open(blast)
		inl = inp.readline()
		while inl != "":
			L = inl.split("\t")
			if L[qors] not in FD:
				print("Should not happen:",L[qors])
			else:
				FD[L[qors]] = 1			
			inl = inp.readline()
		
		oup = open(blast+".missing.%s" % qors,"w")
		for i in FD:
			if FD[i] == 0:
				oup.write("%s\n" % i)
						

	##
	# Get scores from a parsed score file based on a gene list. This function
	# is useful when the BLAST was run on a set of putatives and a confirmed
	# set from these is consolidated. The scores for this confirmed set can
	# be extracted.
	#
	# @param seq_list   file with seq_id
	# @param score_file file with scores, output of the parse function here
	## 
	def get_selected(self,seq_list,score_file):
		
		print("Get scores based on: %s" % seq_list)
		# construct seq name dict
		ndict = {}
		inp = open(seq_list,"r")
		inline = inp.readline()
		while inline != "":
			if self.rmlb(inline) not in ndict:
				ndict[self.rmlb(inline)] = 1
			else:
				print("Redundant:",self.rmlb(inline))
			inline = inp.readline()
				
		# scan through score file
		inp = open(score_file,"r")
		oup = open(seq_list+".scores","w")
		inline = inp.readline()
		while inline != "":
			llist = inline.split("\t")
			if llist[0] in ndict or llist[1] in ndict:
				oup.write(inline)
			inline = inp.readline()

		
	
	##
	# This function gets the subject part of the alignment. But not that simple
	# though. Since most likely there will be multiple matches, the best thing
	# to do would be doing a multiple sequence alignment for each subject and
	# determine the consensus. All sequences of one subj are most likely the
	# same ori for the EST contig this module is designed for. In case it is 
	# not, some screwy situation WILL arise. This is NOT dealt with yet.
	#
	# Now, in order to make this application useful for genomic sequences. It is
	# allowed to have more than one stretch for the same subject sequence. 
	# 
	# Two outputs are generated:
	#  1. "_ext.fa" - all sequences, regardless of the length
	#  2. "_assm.fa" - the assembled one
	#
	# @param stype  sequence in protein [0] or nucleotide [1,default] coords
	# @param wu    wash u blast or not, default 0 (not).
	# @param ignorestop don't care about stop [1] or exclude it [0,defualt]
	##
	def extract_cds(self,blast,stype,T,wu,ignorestop=0):
	
		if stype:
			INC = 3
		else:
			INC = 1
		
		print("Start extract coding sequences")

		#
		# this part gets all subj sequences out of the blast output file
		#
		inp = open(blast,"r")
		oup = open(blast+"_ext.fa","w")
		inline = inp.readline()
		subj = seq = L = R = ""
		hasL = 0
		cdict = {}
		
		countSubj = 0
		while inline != "":
			if inline[:7] == "Query= ":
				print("Query:",inline[7:-1])
				pass
			# the need to find "score =" is due to the presence of multiple
			# matching areas that are disjuct.
			elif inline[0] == ">" or inline.find("Score =") != -1:
				# if sequence is not empty and ">" or "score =" is encountered
				# this signify the start of a new entry, so put the preivious
				# one into dict
				if seq != "":
					if L == '':
						print("Is this WU blast? Need to set the flag!")
						sys.exit()
					# rid of non-alpha characters
					if seq.find("-") != -1:
						segments = seq.split("-")
						seq = ""
						for i in segments:
							seq = seq + i
					if not ignorestop and seq.find("*") != -1:
						segments = seq.split("*")
						seq = ""
						for i in segments:
							seq = seq + i				   
					
					# store subj into dict and write sequence
					if subj not in cdict:
						countSubj = countSubj+1
						
						if int(L)<int(R):
							cdict[subj] = [[int(L),int(R)]]
						else:
							cdict[subj] = [[-int(L),-int(R)]]
					else:
						if int(L)<int(R):
							cdict[subj].append([int(L),int(R)])
						else:
							cdict[subj].append([-int(L),-int(R)])	   
					oup.write(">"+subj+"_"+L+"_"+R+"\n"+seq+"\n")
					
					# reset
					seq = L = R = ""					
					# only reset subj if new subj is present
					if inline[0] == ">":
						subj = ""										       
					hasL = 0
									
				if inline[0] == ">":
					subj = self.rmlb(inline)[1:]
							
			elif inline.find("Sbjct") != -1:
				llist = inline.split(" ")
				tmp = []
				for j in llist:
					if j != "":
						tmp.append(j)
				llist = tmp				
				
				if not hasL:
					# if wu blast, L coord is the second non-empty element.
					if wu:
						print("NOT DEALT WITH")
						sys.exit(0)
						for j in llist[1:]:
							if j != "":
								L = j
								break
					else:
						L = llist[1]
					R = self.rmlb(llist[-1])
					hasL = 1
				else:
					R = self.rmlb(llist[-1])   # keep getting R   
				#print [L,R]
				
				# this is for situation where the LEFT coord is sticked to seq
				# no space in between			   
				alphaIndex = 0
				try:
					if llist[-2][0].isdigit():
						for i in range(len(llist[-2])):
							if llist[-2][i].isalpha():
								alphaIndex = i
								break
				# there are some situations where the sbjct line contain only
				# coord and no sequence, these lines will be ignored.
				except IndexError:
					pass
				
				seq = seq + llist[-2][alphaIndex:]
								
			inline = inp.readline()	 
		
		# rid of non-alpha characters in last sequence
		if seq.find("-") != -1:
			segments = seq.split("-")
			seq = ""
			for i in segments:
				seq = seq + i
		if not ignorestop and seq.find("*") != -1:
			print("STOP")
			segments = seq.split("*")
			seq = ""
			for i in segments:
				seq = seq + i				   
		
		# store subj into dict and write sequence for last
		if subj not in cdict:
			countSubj = countSubj+1
			
			if int(L)<int(R):
				cdict[subj] = [[int(L),int(R)]]
			else:
				cdict[subj] = [[-int(L),-int(R)]]
		else:
			if int(L)<int(R):
				cdict[subj].append([int(L),int(R)])
			else:
				cdict[subj].append([-int(L),-int(R)])
		
		oup.write(">"+subj+"_"+L+"_"+R+"\n"+seq+"\n")
		
		# Separate + and - orientations
		tdict = {}
		for i in list(cdict.keys()):
			countW = 0
			countC = 0
			wlist = []
			clist = []
			for j in range(len(cdict[i])):
				if abs(cdict[i][j][0]) < abs(cdict[i][j][1]):
					countW += 1
					wlist.append(cdict[i][j])
				elif abs(cdict[i][j][0]) > abs(cdict[i][j][1]):
					countC += 1
					clist.append(cdict[i][j])
			if countW > 0:		
				if i+"_0" in tdict:
					tdict[i+"_0"].append(wlist)
				else:
					tdict[i+"_0"] = wlist
			if countC > 0:
				if i+"_1" in tdict:		
					tdict[i+"_1"].append(clist)
				else:
					tdict[i+"_1"] = clist
		
		cdict = tdict
		#
		# consolidate the coordinates for sequences in each subj
		#
		print("\nConsolidate coordinates:")
		for i in list(cdict.keys()):
			
			#print "Elements:",len(cdict[i])
			idxL = idxR = 0
			mdict = {} # idx of element which is eliminated		 
			
			# find entries encompassed by other sequences
			for j in range(len(cdict[i])):
				
				#print j
				if j in mdict:
					#print " skipped1, in mdict already"
					continue
					
				for k in range(j+1,len(cdict[i])):
					
					if k in mdict:
						#print " skipped2, in mdict already"
						continue		
					#print "",cdict[i][j][0],":",cdict[i][j][1],"vs",\
					#		cdict[i][k][0],":",cdict[i][k][1]
					if cdict[i][j][0] < cdict[i][k][0]:      # jL<kL
						if cdict[i][j][1] >= cdict[i][k][1]: # and jR>=kR, j(k)
							mdict[k] = 1
							#print "  j(k)"
						else:
							pass
							#print "  j/k"
					elif cdict[i][j][0] > cdict[i][k][0]:    # jL>kL
						if cdict[i][j][1] <= cdict[i][k][1]: # and jR<=kR, k(j)
							mdict[j] = 1 
							#print "  k(j)"
							break					   
						else:
							pass
							#print "  j/k"
					else:							         # jL = kL
						if cdict[i][j][1] > cdict[i][k][1]:  # and jR>kR, j(k)
							mdict[k] = 1
							#print "  j(k)"
						elif cdict[i][j][1]<cdict[i][k][1]:  # and jR<kR, k(j)
							mdict[j] = 1
							#print "  k(j)"
							break   
						else:						    # and jR=kR, k=j
							mdict[k] = 1    
							#print "  k=j"  
			
			#print cdict[i]
			#print mdict
			# clean up cdict
			tlist = []
			for j in range(len(cdict[i])):
				if j not in mdict:
					tlist.append(cdict[i][j])
			cdict[i] = tlist
			#print cdict[i]
			
			# order elements in cdict
			# L coord as key, element idx as value
			Ldict = {}
			for j in range(len(cdict[i])):
				Ldict[cdict[i][j][0]] = j
			keys = list(Ldict.keys())
			keys.sort()
			clist = []
			for j in keys:
				clist.append(cdict[i][Ldict[j]])
			
			# define contribution of each sequence: L,R,contributeL
			for j in range(len(clist)):
				if j == 0:
					clist[j].append(clist[j][0])
				else:
					clist[j].append(clist[j-1][1]+INC)		      
			cdict[i] = clist
			
			#print i,clist
		
		#
		# load and assemble sequences
		#
		print("\nAssemble sequences:")
		
		# load sequences into memory... this will be a problem...
		inp = open(blast+"_ext.fa","r")
		oup = open(blast+"_assm.fa","w")
		inline = inp.readline()
		sdict = {}

		countAssemblee = 0
		while inline != "":
			if inline[0] == ">":			    
				inline = inline[1:]
				#print "",inline[:-1]
				llist  = inline.split("_")
				#print llist
				subj   = ""
				L      = int(llist[-2])
				R      = int(llist[-1])
				for j in llist[:-2]:
					subj = subj + j + "_"
				subj = subj[:-1]
				
				# make sure the name has ori info
				if L < R:
					subj += "_0"
				else:
					subj += "_1"
						
				# check each cdict element and add sequence stretch
				for j in range(len(cdict[subj])):
					
					# make sure that this subj is an element in cdict[j]
					# also find out its index
					if(L!= abs(cdict[subj][j][0]) or \
					   R!= abs(cdict[subj][j][1])):
						continue
					
					countAssemblee = countAssemblee + 1
					seq = inp.readline()[:-1]
					#print "1>"
					#print cdict[subj][j][0],cdict[subj][j][2]
					#print seq
					if cdict[subj][j][0] < cdict[subj][j][2]:
						seq = seq[(cdict[subj][j][2] - \
								   cdict[subj][j][1] + 1)/INC - 1:]
						#print "High: %i-%i" % (cdict[subj][j][2],cdict[subj][j][1])
						#print seq
							
					# the contributed L coord is smaller due to the
					# presence of un-matched area within subj, add a
					# special character "^" 
					elif cdict[subj][j][0] > cdict[subj][j][2]:
						seq = "(%i-%i)" % (abs(cdict[subj][j][2]),
										   abs(cdict[subj][j][0])) + seq
						#print "Low"
						#print seq
					elif cdict[subj][j][0] == cdict[subj][j][2]:
						seq = seq
					# this is a very lame fix... I notice that some fused entry
					# has one extra residue in sequences not going through the
					# above two if-elif statements. This fix works if there
					# are sequences after this. But if this is the only segment,
					# then I will lose one residue.
					else:
						seq = seq[:-1]
						#print "Lower"
						#print seq
									
						
					cdict[subj][j].append(seq)
																
			inline = inp.readline()


		# output assembled sequence
		print("\nOutput assembled sequences...")
		
		countAssembled = 0
		for i in list(cdict.keys()):
			#print cdict[i]
			L = abs(cdict[i][0][0])
			R = abs(cdict[i][0][1])
			seq = ""
			slist = []
			for j in range(len(cdict[i])):
				idx = cdict[i][j][3].find(")")
				disrupt = 0
				if idx != -1:
					gap = cdict[i][j][3][1:idx].split("-")
					if abs(int(gap[0])-int(gap[1])) > T:
						disrupt = 1
						
				if disrupt:					
					slist.append([i,L,R,seq])
					countAssembled = countAssembled+ 1
					# reset
					L = abs(cdict[i][j][0])
					R = abs(cdict[i][j][1])
					# rid of the connection info in front
					seq = cdict[i][j][3][idx+1:]
				else:
					seq = seq + cdict[i][j][3]
					R = abs(cdict[i][j][1])
			
			# the last sequence
			slist.append([i,L,R,seq])
			countAssembled = countAssembled+ 1
			
			for j in slist:
				oup.write(">%s_%i_%i\n%s\n" % (j[0],j[1],j[2],j[3]))
		
		print("\nTotal %i unique subjects"		   % countSubj)
		print("      %i sequences to be assembled" % countAssemblee)
		print("      %i assembled sequences"	% countAssembled)
		print("Done!\n")

	##
	# Output a non-redundant list of matching subjects. The code for parsing
	# descriptors ARE NOT IMPLEMENTED YET. Just did it.
	#
	# @param desc    boolean variable with 1 = parse subj from descriptor lines
	#		         or 0 = parse subj from alignment [default]
	# @param style   flat[default] or tabular
	# @param ttype   threshold stype, evalue or identity
	# @param T       evalue threshold, inclusive. Only work for tabular format.
	##
	def get_subj(self,blast,desc,style,ttype,T):

		print("Get matching subjects:")
		print(" Blast output:",blast)
		print(" Desc flag   :",desc)
		print(" Style       :",style)
		
		if T != 0 and style != 9:
			print("WARNING: T can't be apply to style 1 at this point!")
		
		inp = open(blast,"r")
		oup = open(blast+".gi","w")
		inline = inp.readline()
		sdict = {}
		count = 0
		while inline != "":
			# parse alignment
			if not desc:
				# typical alignment format
				if style == 1:
					if inline.find("Query=") != -1:
						print(inline[6:-1])
					elif inline[0] == ">":
						if inline[1:-1] not in sdict:
							count = count+ 1
							sdict[inline[1:-1]] = 0
							oup.write(inline[1:])
				# tabular format
				elif style == 9:
					if inline[0] != "#":
						llist = inline.split("\t")
						try:
							if llist[1] not in sdict:
								if T != 0:
									if ttype == "evalue":
										if float(llist[-2]) <= T:
											sdict[llist[1]] = float(llist[-2])
											count = count+ 1
											oup.write(llist[1]+"\n")
									else:
										if float(llist[2])  >= T:
											sdict[llist[1]] = float(llist[-2])
											count = count+ 1
											oup.write(llist[1]+"\n")
								else:
									sdict[llist[1]] = float(llist[2])
									count = count+ 1
									oup.write(llist[1]+"\n")
						except IndexError:
							print("ERROR-line format:",[inline])
				else:
					print("Unknown style type, quit!")
					sys.exit(0)
			# parse descriptor lines
			else:
				if inline.find("Sequences producing") != -1:
					inp.readline() # rid of the empty line
					inline = inp.readline()
					while inline != "" and inline not in ["","\n","\r\n"]:
						llist  = inline.split(" ")
						firstS = 0
						for k in range(len(llist)):
							if llist[k] == "":
								if firstS == 0:
									firstS = 1
								# consecutive space, put together strings before
								else:
									key = ""
									for m in llist[:k-1]:
										key = key + m + " "
									key = key[:-1]
									if key in sdict:
										sdict[key] += 1
									else:
										count += 1
										sdict[key] = 1
									break
						inline = inp.readline()
					
			inline = inp.readline()
			
		if desc:
			for j in list(sdict.keys()):
				#print [j,sdict[j]]
				#oup.write("%s\t%i\n" % (j,sdict[j]))
				oup.write(j+"\n")
				
		print("Total %i non-redundant matching subjects" % count)
		

	##
	# Replace the names with indices
	#
	# @param  score  the score file to be converted with the following format:
	#		[name1][name2]...[score]
	#		Notice that the score will be in the last token.
	# @output .mod   the modified score file with indices
	# @output .name  the index to name relationship
	##
	def index_names(self,score):
		
		print("Start converting names to indices:")
		inp  = open(score,"r")
		oup1 = open(score+".mod","w")
		oup2 = open(score+".name","w")
		
		inline = inp.readline()
		
		# name as key, index as value
		ndict = {}
		count = 0
		c = 0
		while inline != "":
			if c%10000 == 0:
				print(" %ix10k" % (c/10000))
			c += 1
			llist = inline.split("\t")
			if llist[0] not in ndict:
				ndict[llist[0]] = count
				oup2.write(llist[0]+"\t%i" % count+"\n")
				count = count + 1
			if llist[1] not in ndict:
				ndict[llist[1]] = count
				oup2.write(llist[1]+"\t%i" % count+"\n")
				count = count+ 1

			# the last token will have "\n" so it is not added here
			oup1.write("%i\t%i\t%s" % \
				   (ndict[llist[0]],ndict[llist[1]],llist[-1]))
					
			inline = inp.readline()

		print("Total %i score pairs with %i unique indices" % (c,count))
		


	##
	# Rename id within fasta file
	#
	# @param fasta
	# @param name  name file in [new][old] format. Tab-delimited.
	##
	def rename(self,fasta,name):

		# put names to a dict, old as key, new as value
		inp    = open(name,"r")
		inline = inp.readline()
		ndict  = {}
		while inline != "":
			L = inline[:-1].split("\t")
			if L[1] in ndict:
				print("Redundant names")
			else:
				ndict[L[1]] = L[0]
			inline = inp.readline()

		# scan fasta descriptor line
		inp    = open(fasta,"r")
		inline = inp.readline()
		countF = 0  # found
		countN = 0  # not found
		while inline != "":
			if inline[0] == ">":
				if inline[1:-1] in ndict:
					oup.write(">%s\n" % ndict[inline[1:-1]])
					countF = countF+1
				else:
					oup.write(inline)
					countN = countN+1
			else:
				oup.write(inline)
			inline = inp.readline()

		print("Found %s, not found %s" % (countF,countN))
		
	
	##
	# Specifically parse the m = 0 output into m = 8 format
	#
	# 0         1           2           3                 4
	# Query id, Subject id, % identity, alignment length, mismatches, 
	# 5             6         7       8         9       10       11
	# gap openings, q. start, q. end, s. start, s. end, e-value, bit score
	#
	# @param flag  generate screen output or not
	#
	##
	def parse_align2(self,blast,flag=1,gap=0):

		QUERY_END = "  Database:"
		QUERY_END2= "Lambda"
		ID_TAG    = " Score ="
		ID_STR    = "("
		ID_END    = "%"
		
		inp = open(blast,"r")
		oup = open(blast+".mod","w")
		inl = inp.readline()
		subj = query = ""

		olist = [""]*12 
		qL = qR = sL = sR = 0
		countL = 0
		if flag:
			print("Start parsing %s" % blast)
		while inl != "":
			if flag and countL % 100000 == 0:
				print(" %ix100k" % (countL/100000))
			countL += 1
			if inl[-2:] == "\r\n":
				inl = inl[:-2]
			elif inl[-1] == "\n":
				inl = inl[:-1]
			
			# get query
			if inl.find("Query=") != -1:
				query = inl[inl.find("=")+2:]
				# for bl2seq, no query name will be given, so here it is.
				if query == "":
					query = "query"
			# get subj
			elif inl.find(">") != -1:
				subj = inl[1:]				
			# parse score and expect	
			elif inl.find("Score =") != -1:
				
				# OUTPUT point, so every stretch has output
				if olist[0] != "":
					olist[6:10] = [qL,qR,sL,sR]
					# make sure all items in olist are strings 
					olist = [str(item) for item in olist]
					#print olist     
					oup.write("%s\n" % (string.joinfields(olist,"\t")))					
				
				# reset
				olist = [""]*12
				olist[0] = query
				olist[1] = subj
				qL = qR = sL = sR = 0
								
				# rid of space
				ilist = inl.split(" ")
				tmp = []
				for j in ilist:
					if j != "":
						tmp.append(j)
				ilist = tmp
				
				# expect is the last non-empty element, make sure it doesn't
				# start with e
				if ilist[-1][0] == "e":
					ilist[-1] = "1"+ilist[-1]
				olist[10] = ilist[-1]
				
				# score is the second non-empty element
				olist[11] = ilist[2]
				
			# parse % id, alignment length, mismatch, gap				
			elif inl.find("Identities =") != -1:
				ilist = inl.split(" ")
				# rid of empty elements
				tmp = []
				for j in ilist:
					if j != "":
						tmp.append(j)
				ilist = tmp
				# alignment length and mismatch
				mismatch = ilist[2].split("/")
				olist[3] = mismatch[1]
				mismatch = int(mismatch[1]) - int(mismatch[0])
				olist[4] = "%i" % mismatch
				# id
				olist[2] = ilist[3][ilist[3].find("(")+1:ilist[3].find("%")]
				# gap length, not really gap openings! 
				if ilist[-4] != "Gaps":
					olist[5] = "0"
				else:
					olist[5] = ilist[-2][:ilist[-2].find("/")]
			
			# parse coords				
			elif inl.find("Query:") != -1 or inl.find("Sbjct:") != -1:
				ilist = inl.split(" ")
				# L is the second non-empty, R is the last non-empty
				c = 0
				L = R = ""
				tmp = []
				for j in ilist:
					if j != "":
						tmp.append(j)
				ilist = tmp
				
				if ilist[0] == "Query:":
					if qL == 0:
						qL = ilist[1]
					qR = ilist[-1]
				elif ilist[0] == "Sbjct:":
					if sL == 0:				
						sL = ilist[1]
					sR = ilist[-1]
		
			inl = inp.readline()
		
		# write the last query
		olist[6:10] = [qL,qR,sL,sR]
		#print olist
		# make sure all items in olist are strings
		olist = [str(item) for item in olist]
		oup.write("%s\n" % (string.joinfields(olist,"\t")))
		
		if flag:
			print("Total %i subject lines. Done!" % countL)
		
	
	## [NLP, 9-26-2014]
	# Same as parse_align2 but modifed for Blast+ syntax
	# (i.e. no colons in the Query/Sbjct lines)
	#
	# @param flag  generate screen output or not
	##
	def parse_align2_blastPlus(self,blast,flag=1,gap=0):

		QUERY_END = "  Database:"
		QUERY_END2= "Lambda"
		ID_TAG    = " Score ="
		ID_STR    = "("
		ID_END    = "%"
		
		inp = open(blast,"r")

		# check if the format is correct
		if not inp.readline().startswith("TBLASTN "):
			print("ERR: anticipated TBLASTN with output format 0.")
			sys.exit(0)

		oup = open(blast+".mod","w")
		inl = inp.readline()
		subj = query = ""

		olist = [""]*12 
		qL = qR = sL = sR = 0
		countL = 0
		if flag:
			print("Start parsing %s" % blast)
		while inl != "":
			if flag and countL % 100000 == 0:
				print(" %ix100k" % (countL/100000))
			countL += 1
			if inl[-2:] == "\r\n":
				inl = inl[:-2]
			elif inl[-1] == "\n":
				inl = inl[:-1]
			
			# get query
			if inl.find("Query=") != -1:
				query = inl[inl.find("=")+2:]
				# for bl2seq, no query name will be given, so here it is.
				if query == "":
					query = "query"
			# get subj
			elif inl.find(">") != -1:
				subj = inl[1:]				
			# parse score and expect	
			elif inl.find("Score =") != -1:
				
				# OUTPUT point, so every stretch has output
				if olist[0] != "":
					olist[6:10] = [qL,qR,sL,sR]
					# make sure all items in olist are strings 
					olist = [str(item) for item in olist]
					#print olist     
					oup.write("%s\n" % ("\t".join(olist)))					
				
				# reset
				olist = [""]*12
				olist[0] = query.strip(" ")
				olist[1] = subj.strip(" ")
				qL = qR = sL = sR = 0
								
				# rid of space
				ilist = inl.split(" ")
				tmp = []
				for j in ilist:
					if j != "":
						tmp.append(j)
				ilist = tmp				
	
				# In blastplus, method is the end of this line
				# Instead, index expect value from front (7th value)
				if ilist[7][0] == "e":
						ilist[7] = "1"+ilist[7]
				olist[10] = ilist[7].strip(",")
				
				# score is the second non-empty element
				olist[11] = ilist[2]
				
			# parse % id, alignment length, mismatch, gap				
			elif inl.find("Identities =") != -1:
				ilist = inl.split(" ")
				# rid of empty elements
				tmp = []
				for j in ilist:
					if j != "":
						tmp.append(j)
				ilist = tmp
				# alignment length and mismatch
				mismatch = ilist[2].split("/")
				olist[3] = mismatch[1]
				# id
				olist[2] = ilist[3][ilist[3].find("(")+1:ilist[3].find("%")]
				# gap length, not really gap openings! 
				if ilist[-4] != "Gaps":
					olist[5] = "0"
				else:
					olist[5] = ilist[-2][:ilist[-2].find("/")]
				mismatch = int(mismatch[1]) - int(mismatch[0]) - int(olist[5])
				olist[4] = "%i" % mismatch

			
			# parse coords				
			elif inl.find("Query") != -1 or inl.find("Sbjct") != -1:
				ilist = inl.split(" ")
				# L is the second non-empty, R is the last non-empty
				c = 0
				L = R = ""
				tmp = []
				for j in ilist:
					if j != "":
						tmp.append(j)
				ilist = tmp
				
				if ilist[0] == "Query":
					if qL == 0:
						qL = ilist[1]
					qR = ilist[-1]
				elif ilist[0] == "Sbjct":
					if sL == 0:				
						sL = ilist[1]
					sR = ilist[-1]
		
			inl = inp.readline()
		
		# write the last query
		olist[6:10] = [qL,qR,sL,sR]
		#print olist
		# make sure all items in olist are strings and remove trailing spaces
		olist = [str(item.strip()) for item in olist]
		oup.write("%s\n" % ("\t".join(olist)))
		
		if flag:
			print("Total %i subject lines. Done!" % countL)
		
	## [NLP, 9-26-2014]
	# Same as parse_align2 but modifed for blastall syntax
	# (i.e. )
	#
	# @param flag  generate screen output or not
	##
	def parse_align2_blastall(self,blast,flag=1,gap=0):

		QUERY_END = "  Database:"
		QUERY_END2= "Lambda"
		ID_TAG    = " Score ="
		ID_STR    = "("
		ID_END    = "%"
		
		inp = open(blast,"r")
		oup = open(blast+".mod","w")
		inl = inp.readline()
		subj = query = ""

		olist = [""]*12 
		qL = qR = sL = sR = 0
		countL = 0
		if flag:
			print("Start parsing %s" % blast)
		while inl != "":
			if flag and countL % 100000 == 0:
				print(" %ix100k" % (countL/100000))
			countL += 1
			if inl[-2:] == "\r\n":
				inl = inl[:-2]
			elif inl[-1] == "\n":
				inl = inl[:-1]
			
			# get query
			if inl.find("Query=") != -1:
				query = inl[inl.find("=")+2:]
				# for bl2seq, no query name will be given, so here it is.
				if query == "":
					query = "query"
			# get subj
			elif inl.find(">") != -1:
				subj = inl[1:]				
			# parse score and expect	
			elif inl.find("Score =") != -1:
				
				# OUTPUT point, so every stretch has output
				if olist[0] != "":
					olist[6:10] = [qL,qR,sL,sR]
					# make sure all items in olist are strings 
					olist = [str(item) for item in olist]
					#print olist     
					oup.write("%s\n" % (string.joinfields(olist,"\t")))					
				
				# reset
				olist = [""]*12
				olist[0] = query.strip(" ")
				olist[1] = subj.strip(" ")
				qL = qR = sL = sR = 0
								
				# rid of space
				ilist = inl.split(" ")
				tmp = []
				for j in ilist:
					if j != "":
						tmp.append(j)
				ilist = tmp
				
				# In blastall, method is the end of this line
				# Instead, index expect value from front (7th value)
				if ilist[7][0] == "e":
					ilist[7] = "1"+ilist[7]
				olist[10] = ilist[7].strip(",")
				
				# score is the second non-empty element
				olist[11] = ilist[2]
				
			# parse % id, alignment length, mismatch, gap				
			elif inl.find("Identities =") != -1:
				ilist = inl.split(" ")
				# rid of empty elements
				tmp = []
				for j in ilist:
					if j != "":
						tmp.append(j)
				ilist = tmp
				# alignment length and mismatch
				mismatch = ilist[2].split("/")
				olist[3] = mismatch[1]
				# id
				olist[2] = ilist[3][ilist[3].find("(")+1:ilist[3].find("%")]
				# gap length, not really gap openings! 
				if ilist[-4] != "Gaps":
					olist[5] = "0"
				else:
					olist[5] = ilist[-2][:ilist[-2].find("/")]
				mismatch = int(mismatch[1]) - int(mismatch[0]) - int(olist[5])
				olist[4] = "%i" % mismatch

			
			# parse coords				
			elif inl.find("Query:") != -1 or inl.find("Sbjct:") != -1:
				ilist = inl.split(" ")
				# L is the second non-empty, R is the last non-empty
				c = 0
				L = R = ""
				tmp = []
				for j in ilist:
					if j != "":
						tmp.append(j)
				ilist = tmp
				
				if ilist[0] == "Query:":
					if qL == 0:
						qL = ilist[1]
					qR = ilist[-1]
				elif ilist[0] == "Sbjct:":
					if sL == 0:				
						sL = ilist[1]
					sR = ilist[-1]
		
			inl = inp.readline()
		
		# write the last query
		olist[6:10] = [qL,qR,sL,sR]
		#print olist
		# make sure all items in olist are strings
		olist = [str(item).strip(" ") for item in olist]
		oup.write("%s\n" % (string.joinfields(olist,"\t")))
		
		if flag:
			print("Total %i subject lines. Done!" % countL)
		
	
	#
	# Parse the alignments of a particular ID
	#
	def parse_align3(self,blast,seq_id):
		
		print("\nGet %s from %s:" % (blast,seq_id))
		inp = open(blast,"r")
		oup = open(seq_id+".out","w")
		inl = inp.readline()
		olist = []
		findT = 0
		count = 0
		query = ""
		while inl != "":
			if inl.find("Query=") != -1:
				if olist != []:
					oup.write("%s\n" % string.joinfields(olist,""))
					olist = []
					findT = 0
				else:
					olist.append(inl+"\n")
				query = self.rmlb(inl)
			
			if inl.find(">") != -1:
				if self.rmlb(inl)[1:] == seq_id:
					#print "",query
					findT = 1
					count += 1
					olist.append(inl)
				else:
					findT = 0
			elif findT:
				olist.append(inl)
			inl = inp.readline()
		
		print(" found %i fimes." % count)

	##
	# Parse the alignment parts and output the gap information:
	#  [query][subj][qStart][qEnd][sStart][sEnd][aln_len][qGaps][sGaps]
	# 
	# Gap format: gap1L-gap1R,gap2L-gap2R,...
	# 
	# Also output the alignment
	##
	def parse_gap(self,blast):
		
		inp = open(blast,"r")
		oup1= open(blast+".gap","w")
		oup2= open(blast+".seqpair","w")

		# Query= Os_9634.m04513
		#         (917 letters)
		#
		# >Os_9630.m00627
		#           Length = 908
		# 
		#  Score = 1526 bits (3950), Expect = 0.0
		#  Identities = 770/920 (83%), Positives = 821/920 (89%), Gaps = 15/920 (1%)
		# 
		# Query: 1   MRLSSSSGSVLPAQAASPEAVEEQKCLNSELWHACAGPLVSLPAVGSRVVYFPQGHSEQV 60
		#            M+LS S+G V     + PE  EEQKCLNSELWHACAGPLVSLPAVGSRVVYFPQGHSEQV
		# Sbjct: 1   MKLSPSAGGVSDQPPSPPEVAEEQKCLNSELWHACAGPLVSLPAVGSRVVYFPQGHSEQV 60
		#
		#  Score = ....
		#
		#
		Qtag  = "Query= "
		Qend1 = "  Database:"
		Qend2 = "Lambda"
		Stag  = ">"
		Qline = "Query: "
		Sline = "Sbjct: "
		Rtag  = "Score ="
		
		# Revised gap method, for example:
		#   1  2345 678
		# Q X--XXXX-XXX
		# S XXXX-XXX--X
		#   1234 567  8
		# Should return:
		# qlist = [[1,2-3],[5,7-7]]
		# slist = [[4,2-4],[7,6-7]]
		def gap(q,s):
			# print q
			# print s
			qL = []
			sL = []
			if len(q) == len(s):
				qI        = 0 # keep track of query string non-gap count
				sI        = 0 # keep track of subjt string non-gap count
				qgap      = 0 # query gap position
				sgap      = 0 # subject gap position
				qgap_size = 0
				sgap_size = 0

				for i in range(len(q)):
					# output qgap info
					if q[i] != "-":
						if qgap_size != 0:
							#print " qgap: %i,%i-%i" % (qI,sI-qgap_size+1,sI)
							qL.append("%i|%i-%i" % (qI,sI-qgap_size+1,sI))
							qgap_size = 0

					# output sgap info
					if s[i] != "-":
						if sgap_size != 0:
							#print " sgap: %i,%i-%i" % (sI,qI-sgap_size+1,qI)
							sL.append("%i|%i-%i" % (sI,qI-sgap_size+1,qI))
							sgap_size = 0

					# This part need to be separated from the gap location
					# output part because situation where the gap is right next
					# to each other like (which I don't recall happening...):
					#  XX--XXXX
					#  XXXX--XX
					if q[i] != "-":
						qI += 1
					else:
						qgap_size +=1
				
					if s[i] != "-":	
						sI += 1
					else:
						sgap_size +=1
											
					#print "Q:",qI,qgap_size
					#print "S:",sI,sgap_size		
			else:
				print("ERR:query and subjt strings different length")
				print(qstr)
				print(sstr)
				sys.exit(0)
			
			return qL,sL
				
		# find -, and return a list of coordinates (alignment-based)
		"""
		def gap(astr):	
			c = 1
			gL = gR = 0
			glist = []
			print [astr]
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
		"""
		
		print("Find gaps...")
		inl = inp.readline()
		query = subjt = queryNew = subjtNew = qline = sline = ""
		qC = []
		sC = []
		while inl != "":
			inl = self.rmlb(inl)
			#print [inl]
			# new query
			if inl.find(Qtag) != -1:
				if query == "":
					query = inl[len(Qtag):]
				queryNew = inl[len(Qtag):]
				#print "Q:",[query,queryNew]
			# new subject								
			elif inl.find(Stag) != -1:
				if subjt == "":
					subjt = inl[1:]
				subjtNew = inl[1:]
				#print "S:",[subjt,subjtNew]
			# new score segment
			elif inl.find(Rtag) != -1:
				#print [query,subjt]
				if sline != "":
					
					# process gap info					
					#qlist = gap(qline)
					#slist = gap(sline)
					qlist,slist = gap(qline,sline)
					
					# should be the same as sline
					alen  = len(qline)
					#  [q][s][qStr][qEnd][sStr][sEnd][alnL][qGaps][sGaps]
					#print ">>>>>write1:",query,subjt,qC[0],qC[-1],sC[0],sC[-1]
					oup1.write("%s\t%s\t%s\t%s\t%s\t%s\t%i\t%s\t%s\n" % \
						(query,subjt,qC[0],qC[-1],sC[0],sC[-1],alen,
						 string.joinfields(qlist,","),
						 string.joinfields(slist,",")))
					# 2 lines, each pair seperate by an empty line
					oup2.write("%s\t%s\n%s\t%s\n\n" % (query,qline,subjt,sline))

					qline = sline = ""
					query = queryNew
					subjt = subjtNew
					qC = []
					sC = []

			elif inl.find(Qline) != -1:
				L = inl.split(" ")
				qline += L[-2]
				#print "Q: append:",L[1],L[-1]
				qC.append(L[1])
				qC.append(L[-1])
			elif inl.find(Sline) != -1:
				L = inl.split(" ")
				sline += L[-2]
				sC.append(L[1])
				sC.append(L[-1])
			inl = inp.readline()

		# write the last entry
		qlist,slist = gap(qline,sline)
		alen  = len(qline)
		#print "writeE:",query,subjt,qC[0],qC[-1],sC[0],sC[-1]
		oup1.write("%s\t%s\t%s\t%s\t%s\t%s\t%i\t%s\t%s\n" % \
			(query,subjt,qC[0],qC[-1],sC[0],sC[-1],alen,
			 string.joinfields(qlist,","),
			 string.joinfields(slist,",")))
		oup2.write("%s\t%s\n%s\t%s\n\n" % (query,qline,subjt,sline))

	
		


	##
	# Parse the alignment part of the file with a threshold filter.
	# 
	# @param blast  the output from blastall in default format (m = 0)
	# @param T      percentID cutoff
	# @param wself  include self (1) or not (0, default)
	# @param format include alignment lines [0,default] or only unmatched
	#	       regions [1]
	# @param style  involve no nt [0,default], or do [1]
	# @param query  
	##
	def parse_align(self,blast,T,wself,format,style,query=""):

		#print "Parse alignment:",blast
		inp = open(blast,"r")
		oup = open(blast+".log","w")

		oup.write("BLAST output:"+blast+"\n")
		oup.write("Threshold   :"+str(T)+"\n")
		oup.write("Include self:"+str(wself)+"\n")
		
		QUERY_STR = "Query="
		QUERY_END = "  Database:"
		QUERY_END2= "Lambda"
		SUBJ_STR  = ">"
		ID_TAG    = " Score ="
		ID_STR    = "("
		ID_END    = "%"

		inline  = inp.readline()
		outlines= []
		subj = "" 
		id   = ""
		w    = 0
		skip = 0
		#query_end = 0
		
		print("Start parsing...")
		while inline != "":
			
			if inline[-2:] == "\r\n":
				inline = inline[:-2]
			elif inline[-1] == "\n":
				inline = inline[:-1]
			
			# find query
			if inline.find(QUERY_STR) != -1:
				
				# new query, write output for the previous query
				if outlines != []:
					for i in outlines:
						oup.write(i+"\n")
					
					outlines = []
					#query_end = 0
					
				# set query
				
				# This following if else statement is inactivated because I am
				# not sure what this is trying to do.
				#if not style:
				#       query = inline[inline.find(" ")+1:]
				#       outlines.append("%s\n%s" % ("-"*40,inline))
				#else:
				#       outlines.append("%s\nQuery= %s\n" % ("-"*40,query))
				
				query = inline[inline.find(" ")+1:]
				outlines.append("\n%s\nQuery= %s" % ("-"*40,query))
				#print "",query
							
			# if it is not query, find subject start"
			elif inline.find(SUBJ_STR) != -1:
				subj = inline[1:]
				# see if include self or not
				if not wself:
					if subj == query:
						skip = 1
					else:
						skip = 0
					#print "",skip		
				w = 0
			
			if not skip:
				# find identity and verify if it is above threshold
				if inline.find(ID_TAG) != -1:
					inline = inp.readline()	
					if inline[-2:] == "\r\n":
						inline = inline[:-2]
					elif inline[-1] == "\n":
						inline = inline[:-1]
							 
					id = int(inline[inline.find(ID_STR)+1:inline.find(ID_END)])
					if id > T:
						F = " "
						if id == 100:
							F = "y"
						outlines.append("\n[%s] SUBJ:%s(%i" % (F,subj,id) + "%)\n")
						#print "subj:%s(%i" % (subj,id) + "%)"
						inline = inp.readline()
						w = 1
					else:
						w = 0
				elif(inline.find(QUERY_END) != -1 or
					 inline.find(QUERY_END2)!= -1):
					
					#query_end = 1 
					#print "---------IN QUERY END----------"
					#outlines = outlines[:-1]
					inline = inp.readline()
					continue
				
				# subj match above threshold
				if w == 1:
					# rid of query and subj alignment lines, and reformat match
					# line

					if len(inline) > 0:					
						if inline[-2:] == "\r\n":
							inline = inline[:-2]
						elif inline[-1] == "\n":
							inline = inline[:-1]
					if len(inline)>7:			     
						if inline[:6] == "Query:":
							# the second last token is the align part
							llist   = inline.split(" ")
							aln_idx = inline.find(llist[-2])
							if not format:
								outlines.append(inline[7:])

							# get match line
							inline = inp.readline()	
								
							if inline[-2:] == "\r\n":
								inline = inline[:-2]
							elif inline[-1] == "\n":
								inline = inline[:-1]
							
							inline = inline[aln_idx:]
							convert = ""
							for j in inline:
								if j not in [" ","+"]:
									convert = convert + "."
								elif j == " ":
									convert = convert + "X"
								elif j == "+":
									convert = convert + "+"
							
							if not format:
								outlines.append(" "*(aln_idx-7)+convert)						
							else:
								outlines.append(convert)
							
							# rid of the subj align line
							if not format:
								inline = inp.readline()
								outlines.append(inline[7:])							
							else:
								inp.readline()
			
			inline = inp.readline()
		
		# write the last query
		if outlines != []:
			for j in outlines:
				oup.write(j+"\n")
					
		print("\nParse_align done!\n")
	
	#
	# Generate output like:
	# #Query\tSubj\t%ident\tqL-qR|sL-sR\tEval
	# seq1
	# seq2
	# 
	def parse_align4(self,blast):
		QUERY_END = "  Database:"
		QUERY_END2= "Lambda"
		ID_TAG    = " Score ="
		ID_STR    = "("
		ID_END    = "%"
		
		inp = open(blast,"r")
		oup = open(blast+".align_seq","w")
		inl = inp.readline()
		subj = query = ""
		
		# 0 1 2  3   4        5   6  7  8  9  10 11    12   13
		# Q S ID LEN MISMATCH GAP QL QR SL SR E  SCORE SEQL SEQR
		OL = [""]*14 
		qL = qR = sL = sR = 0
		countL = 0
		while inl != "":
			
			inl = inl.strip()
			
			# get query
			if inl.find("Query=") != -1:
				query = inl[inl.find("=")+2:]
				# for bl2seq, no query name will be given, so here it is.
				if query == "":
					query = "query"
			# get subj
			elif inl.find(">") != -1:
				subj = inl[1:]				
			# parse score and expect	
			elif inl.find("Score =") != -1:
				# OUTPUT point, so every stretch has output
				if OL[0] != "":
					if countL % 1e3 == 0:
						print(" %i k" % (countL/1e3))
					countL += 1
					OL[6:10] = [qL,qR,sL,sR]
					oup.write("#%s %s %s-%s|%s-%s %s\n" % \
								(OL[0],OL[1],OL[6],OL[7],OL[8],OL[9],OL[10]))
					oup.write("%s\n%s\n" % (OL[12],OL[13]))
				
				# reset
				OL = [""]*14
				OL[0] = query
				OL[1] = subj
				qL = qR = sL = sR = 0
								
				# rid of space
				ilist = inl.split(" ")
				tmp = []
				for j in ilist:
					if j != "":
						tmp.append(j)
				ilist = tmp
				
				# expect is the last non-empty element, make sure it doesn't
				# start with e
				if ilist[-1][0] == "e":
					ilist[-1] = "1"+ilist[-1]
				OL[10] = ilist[-1]
				
				# score is the second non-empty element
				OL[11] = ilist[2]
				
			# parse % id, alignment length, mismatch, gap				
			elif inl.find("Identities =") != -1:
				ilist = inl.split(" ")
				# rid of empty elements
				tmp = []
				for j in ilist:
					if j != "":
						tmp.append(j)
				ilist = tmp
				# alignment length and mismatch
				mismatch = ilist[2].split("/")
				OL[3] = mismatch[1]
				mismatch = int(mismatch[1]) - int(mismatch[0])
				OL[4] = "%i" % mismatch
				# id
				OL[2] = ilist[3][ilist[3].find("(")+1:ilist[3].find("%")]
				# gap length, not really gap openings! 
				if ilist[-4] != "Gaps":
					OL[5] = "0"
				else:
					OL[5] = ilist[-2][:ilist[-2].find("/")]
			
			# parse coords				
			elif inl.find("Query:") != -1 or inl.find("Sbjct:") != -1:
				#Query:1   xxxxxxx 57
				ilist = inl.split(" ")
				# L is the second non-empty, R is the last non-empty
				c = 0
				L = R = ""
				tmp = []
				for j in ilist:
					if j != "":
						tmp.append(j)
				ilist = tmp
				if ilist[0] == "Query:":
					if qL == 0:
						qL = ilist[1]
					qR = ilist[-1]
					OL[12] += ilist[2]
				elif ilist[0] == "Sbjct:":
					if sL == 0:				
						sL = ilist[1]
					sR = ilist[-1]
					OL[13] += ilist[2]
		
			inl = inp.readline()
		
		# write the last query
		OL[6:10] = [qL,qR,sL,sR]
		#print OL
		oup.write("#%s %s %s-%s|%s-%s %s\n" % \
							(OL[0],OL[1],OL[6],OL[7],OL[8],OL[9],OL[10]))
		oup.write("%s\n%s\n" % (OL[12],OL[13]))
		countL += 1
		
		print("Total %i pairs. Done!" % countL)		

	def refine(self,log):
	
		lenT = 50 # length threshold, set at 50 nt
		inp = open(log,"r")
		oup = open(log+".refined","w")
		inl = inp.readline()
		# read into a dict: subjt as key, queries as value in a list
		sdict = {}
		print("Start refinement...")
		while inl != "":
			inl = self.rmlb(inl)
		
			if inl.find("Query=") != -1:
				query = inl[len("Query=")+1:]
			elif inl.find("[") != -1:
				subjt = inl[inl.find(":")+1:inl.find("(")]
				pid   = int(inl[inl.find("(")+1:inl.find(")")-1])
				# found subj, keep reading till empty line come up
				inp.readline()
				leng = 0
				inl = inp.readline()
				inl = self.rmlb(inl)
				while inl != "":
					leng += len(inl)
					inl = inp.readline()
					try:
						inl = self.rmlb(inl)
					except IndexError:
						break
				
				if subjt in sdict:
					sdict[subjt].append([pid,leng,query])
				else:
					sdict[subjt] = [[pid,leng,query]]
			inl = inp.readline()

		# modify sdict so it is ranked
		tdict = {}
		print("Rank...")
		for i in sdict:
			qlist = sdict[i]
			lens = {}        # length as key, pid as value
			qs   = {}        # pid and length as key, query as value
			for j in qlist:
				if j[1] in lens:
					lens[j[1]].append(j[0])
				else:
					lens[j[1]] = [j[0]]
				if "%i %i" % (j[0],j[1]) in qs:
					qs["%i %i" % (j[0],j[1])].append(j[2])
				else:
					qs["%i %i" % (j[0],j[1])] = [j[2]]
			lkeys = list(lens.keys())
			lkeys.sort()
			lkeys.reverse()
			qlist = []
			for j in lkeys:
				if len(lens[j]) > 1:
					lens[j].sort()
					lens[j].reverse()
					for k in lens[j]:
						if len(qs["%i %i" % (k,j)]) > 1:
							for l in qs["%i %i" % (k,j)]:
								qlist.append([k,j,l])
						else:
							qlist.append([k,j,qs["%i %i" % (k,j)][0]])
				else:
					if len(qs["%i %i" % (lens[j][0],j)]) > 1:
						for l in qs["%i %i" % (lens[j][0],j)]:
							qlist.append([lens[j][0],j,l])
					else:
						qlist.append([lens[j][0],j,qs["%i %i" % (lens[j][0],j)][0]])
			tdict[i] = qlist
		sdict = tdict			

		# now put the ranking info back to the log file.
		print("Generate output...")
		inp = open(log,"r")
		oup = open(log+".refined","w")
		inl = inp.readline()
		while inl != "":
			inl = self.rmlb(inl)
			tag = " "
			if inl.find("Query=") != -1:
				query = inl[len("Query=")+1:]
				ostr = inl
				print("",query)
			elif inl.find("[") != -1:
				subjt = inl[inl.find(":")+1:inl.find("(")]
				pid   = int(inl[inl.find("(")+1:inl.find(")")-1])
				
				if len(sdict[subjt]) == 1:  # unique hit
					if sdict[subjt][0][1] > lenT:
						tag = "u"
					else:
						tag = "-"
				else:
					c = 0
					for j in sdict[subjt]:
						if pid == j[0] and query == j[2]:
							if j[1] > lenT:
								tag = str(c)
							else:
								tag = "-"
							break
						c += 1
						
				if tag == "-": # get rid of these below threshold stuff
					inp.readline()
					inl = inp.readline()
					inl = self.rmlb(inl)
					while inl != "":
						inl = inp.readline()
						try:
							inl = self.rmlb(inl)
						except IndexError:
							break
				else:
					ostr = "[%s%s" % (tag,inl[2:])
			else:
				ostr = inl
			
			if tag != "-":
				oup.write(ostr+"\n")
			inl = inp.readline()
			

	##
	# The output generated by parse_align is manually inspected. The "qualified"
	# sequences, in this case things that are alternative spliced, will have a
	# "y" in the bracket within the text output.A dict with query as key, a list
	# of qualified subjects as values will be used for single linkage
	#
	# @param log    output generated by parse_align
	# @param priority treat the query set with higher priority, disregard
	#               length considerations.
	# @param qtag   qualified tags
	##
	def get_qualified(self,log,fasta,qtag="y",priority=0,):

		# read file into a dict
		print("Start get_qualified:")
		print("     log:",log)
		print("  format:",format)
		print("   fasta:",fasta)
		print("priority:",priority)
		print("    qtag:",qtag,"\n")

		print("Read log file, parse qualified entries...")
		inp = open(log,"r")
		qdict = {}
		
		inline = inp.readline()
		query = subj = ""
		while inline != "":
			# find query, space as delimiter
			if inline.find("Query=") != -1:
				query = inline[inline.find(" ")+1:-1]
				if query in qdict:
					print("This shouldn't happen")
				else:
					qdict[query] = []
			elif inline[0] == "[":
				if inline[1] in qtag.split(","):
					subj = inline[inline.find(":")+1:inline.find("(")]
					if subj not in qdict[query]:
						# check if sp header is the same...
						if subj.find("_") != -1 and query.find("_") != -1:
							if subj[:subj.find("_")] != query[:query.find("_")]:
								pass
							else:
								qdict[query].append(subj)
						else:		
							qdict[query].append(subj)
			inline = inp.readline()
		# single linkage
		print("Single linkage...")
		clusters = link.get_relations2(qdict)
			
		# get sequence sizes
		print("Get sequence size...")
		fmanager = FastaManager.fasta_manager()
		sdict     = fmanager.get_sizes(fasta,1)
		
		# generate output
		print("Generate output...")
		oup = open(log+".cluster","w")
		for i in clusters:
			#print "",i
			outline = ""
			if len(i) > 1:
				# add entries one to one to a line, and check their sizes.
				max_id = i[0]
				
				for j in i:
					# priority and not a query seq
					if priority and j not in qdict:
						outline += "%s(%i)\t" % (j,sdict[j])
						continue
					# no priority set or set and is a query seq
					else:
						if max_id != j and (sdict[j] > sdict[max_id]):
							max_id = j
							
					outline += "%s(%i)\t" % (j,sdict[j])
								
				# format 1
				#oup.write("%s\t%s\n" % (max_id,outline[:-1]))	
				# format 2
				# rid of the last \t when splitted
				#print outline
				for j in outline.split("\t")[:-1]:
					if j[:j.find("(")] != max_id:
						#print " >>%s\t%s\n" % (j[:j.find("(")],max_id)
						oup.write("%s\t%s\n" % (j[:j.find("(")],max_id))
				oup.write("%s\t-\n" % max_id)
			else:
				oup.write("%s\t-\n" % i[0])

	#
	# This is similar for get_qualified3. But just get qualified entries. No
	# clustering is done
	#
	# param qOnly   [0]: look into both query and subject sequences, 
	#               [1]: only query
	#               [2]: only subject 
	# 
	def get_qualified4(self,blast,fasta,eT,idenT=95,matchL=150,lengT=0.9,
						qOnly=0):
		
		# get size dict
		print("Get sequence sizes...")
		sizes = fmanager.get_sizes(fasta,1)
		inp = open(blast,"r")
		inl = inp.readline()
		print("Parse blast output:")
		print(" eT    :",eT)
		print(" idenT :",idenT)
		print(" matchL:",matchL)
		print(" lengT :",lengT)
		countQ = 0
		#oupL = open(blast+".err_log","w")
		c = 0
		oupL = open(blast+"_E%iI%iL%iP%iQ%i.qlines" % \
									(eT,idenT,matchL,int(lengT*100),qOnly),"w")
		while inl != "":
			if c > 0 and c%10000 == 0:
				print(" %i x 10k" % (c/10000))
			c += 1
			
			if inl[0] == "#":
				inl = inp.readline()
				continue
			
			L = inl.split("\t")
			#tmp = open("temp", "w")
			#tmp.write(L)
			#sys.exit()
			I = float(L[2])  # identity
			M = int(L[3])    # match length
			E = L[-2]        # evalue
			if E[0] == "e":
				E = "1"+E
			if E == "0.0":
				E = "1e-200"
			
			try:
				E = math.fabs(math.log10(float(E)))
			except ValueError:
				print("Malformed Evalue:",E)
				inl = inp.readline()
				continue
			
			# whether evaluate the alignment lengths of both query and subject
			if qOnly == 1:
				try:
					N = sizes[L[0]]
				except KeyError:
					print("%s not in the size library" % L[0])
					print(len(L))
					print(L[:1])
					sys.exit(0)
			elif qOnly == 2:
				N = sizes[L[1]]
			else:
				# smaller sequence length
				try:             
					N = min(sizes[L[0]],sizes[L[1]]) 
				except:
					if L[0] in sizes:
						N = sizes[L[0]]
					else:
						N = sizes[L[1]]
			
			#print "I,M,N,P,E"
			#print I,M,N,float(M)/float(N),E
			if L[0] != L[1] and I >= float(idenT) and M >= matchL and \
				float(M)/float(N) >= lengT and E >= eT:
				oupL.write(inl)
				countQ += 1
			#else:
			#	print L[0],L[1],I,M,L[-2],E
			inl = inp.readline()
		
		print("%i total, %i qualified" % (c,countQ))
		
		
		
	#
	# This is similar to get_qualified. But the blast input is tabular. And the
	# thresholds can be passed.
	#
	# @param blast   blast output in tabular format
	# @param fasta   sequence file for sequence size
	# @param idenT   identity threshold
	# @param matchL  match length threshold
	# @param lengT   threshold for the proportion between match length and the
	#                shorter sequence
	#
	def get_qualified3(self,blast,fasta,idenT=95,matchL=150,lengT=0.9):
		
		# get size dict
		print("Get sequence sizes...")
		sizes = fmanager.get_sizes(fasta,1)
		
		# construct a dict with: query as key, qualified entries in a list as
		# value. Qualified: >= idenT, >= lentT for the shorter seq.
		qdict = {}
		inp = open(blast,"r")
		inl = inp.readline()
		print("Parse blast output:")
		print(" idenT :",idenT)
		print(" matchL:",matchL)
		print(" lengT :",lengT)
		countQ = 0
		#oupL = open(blast+".err_log","w")
		c = 0
		oupL = open(blast+".qlines","w")
		while inl != "":
			L = inl.split("\t")
			I = float(L[2])                  # identity
			M = int(L[3])                    # match length
			try:
				N = min(sizes[L[0]],sizes[L[1]]) # smaller sequence length
			except:
				#oupL.write("%s,%s\tsize_not_found\n" % (L[0],L[1]))
				if L[0] in sizes:
					N = sizes[L[0]]
				else:
					N = sizes[L[1]]
			
			if c > 0 and c%10000 == 0:
				print(" %i x 10k" % (c/10000))
			c += 1
			if L[0] != L[1] and I >= float(idenT) and M >= matchL and \
				float(M)/float(N) >= lengT:
				#print " %i -> %s %s(%i)" % \
				#		(M,L[0],L[1],sizes[L[1]])
				oupL.write(inl)
				if L[0] in qdict:
					if L[1] not in qdict[L[0]]:
						countQ += 1
						qdict[L[0]].append(L[1])
				else:
					countQ += 1
					qdict[L[0]] = [L[1]]
			inl = inp.readline()
		print("",countQ,"qualified pairs.")
		
		# single linkage
		print("Single linkage...")
		clusters = link.get_relations(qdict,isdict=1)
			
		# generate output
		print("Generate output...")
		oup = open(blast+".cluster","w")
		oup.write("Thresholds: idenT=%i, matchL=%i, lengT=%f\n" % \
					(idenT,matchL,lengT))
		for i in clusters:
			#print "",i
			outline = ""
			if len(i) > 1:
				#print i
				# add entries one to one to a line, and check their sizes.
				for j in i:
					if j in sizes:
						max_id = j
						break
						
				for j in i:
					if j in sizes and max_id != j and \
						(sizes[j] > sizes[max_id]):
						max_id = j
					#outline += "%s(%i)\t" % (j,sizes[j])
								
				# format 1
				#oup.write("%s\t%s\n" % (max_id,outline[:-1]))	
				# format 2
				# rid of the last \t when splitted
				for j in outline.split("\t")[:-1]:
					if j[:j.find("(")] != max_id:
						#print " >>%s\t%s\n" % (j[:j.find("(")],max_id)
						oup.write("%s\t%s\n" % (j[:j.find("(")],max_id))
				oup.write("%s\t-\n" % max_id)
			else:
				oup.write("%s\t-\n" % i[0])
		
		

	
	##
	# This is only for getting qualified subjects for query sequences.
	##
	def get_qualified2(self,log,qtag):
		
		inp = open(log,"r")
		inline = inp.readline()
		
		query = ""
		subjt = ""
		oup   = open(log+".qlist","w")
		written = {}
		while inline != "":
			if inline[-2:] == "\r\n":
				inline = inline[:-2]
			elif inline[-1] == "\n":
				inline = inline[:-1]
				
			if inline.find("Query=") != -1:
				query = inline[inline.find(" ")+1:]
			elif inline.find("SUBJ:") != -1:
				if inline[1] in qtag.split(","):
					subjt = inline[inline.find(":")+1:inline.find("(")]
					if "%s %s" % (query,subjt) in written:
						pass
					else:
						oup.write("%s\t%s\t%s\n" % (query,subjt,inline[1]))
						written["%s %s" % (query,subjt)] = 1
			
			inline = inp.readline()
			
		
		
	

	##
	# This only works for blast table format. For two species, notice that
	# the BLAST search should be conducted for using A+B as both query and
	# subjects.
	# 
	# Before 10/14,2003, the size threshold is based on a comparison of
	# full length sequences, instead of the matching area of the gene. But come
	# to think about it, what is better?? The second one makes more sense.
	# 
	# Before 12/1,03, a bug in wself set to 0, suppose to find only matches
	# between organisms. When fasta is not set, there will be problem. fix it.
	# 
	# Before 05/10,05', the sizeP setting was not useful because of a setting
	# problem. 
	#
	# Before 08/05,08', the shorterM/shorterL give an int which is not good.
	#
	# sizeP  proportion of the alignment length divided by the smaller sequence
	#        of a reciprocal best hit pair. Default 0.8
	##
	def get_reciprocal(self,blast,target,fasta,sizeP,wself,oflag=0):
		
		print("Blast  :",blast)
		print("Target :",target)
		print("Fasta  :",fasta)
		print("sizeP  :",sizeP)
		print("wself  :",wself)
		
		inp = open(blast,"r")
		inl = inp.readline()
		
		if target == "percentSIM":
			print("Can't do percentSIM here, do percentID")
			tidx = 2
		elif target == "percentID":
			tidx = 2
		else:
			tidx = 10
		
		bdict = {}
		while inl != "":
			if inl[-2:] == "\r\n":
				inl = inl[:-2]
			elif inl[-1] == "\n":
				inl = inl[:-1]
				
			ilist = inl.split("\t")
			if ilist[0] != ilist[1]:
				score = float(ilist[tidx])
				
				# notice the match length is defined as:
				# alignment leng - gap openings
				match_leng = float(ilist[3])-float(ilist[5])
												
				# Both query and subject will be verified and put into dict
				# accordingly.
				# for query
				if ilist[0] in bdict:
					if bdict[ilist[0]][1] > score:
						bdict[ilist[0]] = [ilist[1],score,match_leng,inl]
				else:
					bdict[ilist[0]] = [ilist[1],score,match_leng,inl]	
				# for subjct
				if ilist[1] in bdict:
					if bdict[ilist[1]][1] > score:
						bdict[ilist[1]] = [ilist[0],score,match_leng,inl]
				else:
					bdict[ilist[1]] = [ilist[0],score,match_leng,inl]
			inl = inp.readline()
		
		if fasta != "":
			sdict = fmanager.get_sizes(fasta,1)
			if sizeP > 1:
				print("Illegal size threhsold, should <= 1.")
				fasta = ""
		
		oup = open(blast+".recip","w")
		oup1= open(blast+".recip.log","w")
		oup1.write("Query\tSubjt\tLengthShort\tAlnShort\tWrite\n")
		written = {}
		for i in list(bdict.keys()):
			if i in written:
				pass
			else:
				ihit = bdict[i][0]
				if bdict[ihit][0] == i:
					write = 1
					# verify length and if from the same sp (if wself == 0)
					if fasta != "":
						L1 = float(sdict[i])        # length of i
						L2 = float(sdict[ihit])     # length of i's recip hit
						M1 = bdict[i][2]     		# match length for i-ihit
						M2 = bdict[ihit][2]  		# match length for ihit-i
						
						# find out the shorter of the two and also the match leng
						shorterL = L1
						shorterM = M1
						if L2 < L1:
							shorterL = L2
						if M2 < M1:
							shorterM = M2
						
						if float(shorterM)/float(shorterL) < sizeP:
							write = 0
						oup1.write("%s\t%s\t%f\t%f\t%i\n" % \
										(i,ihit,shorterL,shorterM,write))
					
					if(not wself and i[:i.find("_")]==ihit[:ihit.find("_")]):
						write = 0
							
					if write:
						if oflag:
							oup.write("%s\t%s\n" % (i,ihit))
						else:
							if bdict[i][0] == ihit:
								oup.write("%s\n" % bdict[i][-1])
							
				written[i] = 0
				written[ihit] = 0
		
			
	
	def match_list(self,blast,name):
	
		# read qualified names into a list
		print("Read name list...")	
		inp   = open(name,"r")
		inl   = inp.readline()
		names = {}
		while inl != "":
			if self.rmlb(inl) not in names:
				names[self.rmlb(inl)] = 0
			inl = inp.readline()

		# read blast fine
		print("Read blast...")
		inp = open(blast,"r")
		inl = inp.readline()
		bdict = {}
		# get top match for each sequence, ignore self
		while inl != "":
			L = inl.split("\t")
			if L[0] not in bdict:
				bdict[L[0]] = 1			
			inl = inp.readline()
		
		# Start matching subjects with qualified names and generate output
		print("Matching...")
		countIN = 0
		for i in list(bdict.keys()):
			if i in names:
				countIN += 1
				names[i] = 1
		
		oup = open("%s-NOT_IN-%s" % (name,blast),"w")
		for i in names:
			if names[i] == 0:
				oup.write("%s\n" % i)
		print("Total %i names, %i in query list" % (len(list(names.keys())),countIN))	
		print("Total %i queries, %i in name list" % (len(list(bdict.keys())),countIN))		
		print("\nDone!\n")
	
	#
	# This function added taxa info to the end of line. NOT TESTED.
	#
	# @param blast the blast output file in tabular format
	# @param alist a list with [taxa][gi]
	#
	def match_list2(self,blast,alist):
		
		adict = futil.file_to_dict(alist,5)
		inp   = open(blast,"r")
		oup   = open(blast+".mod","w")
		inl   = inp.readline()
		while inl != "":
			inl = self.rn_lb(inl)
			L   = inl.split("\t")
			if L[1] in adict:
				oup.write("%s\t%s\n" % (inl,adict[L[1]]))
			inl = inp.readline()
			
			

	##
	# Parse log file generated by parse_align against the list generated by  
	# match_list. For the alignments, if query-subj pairs are the same as 
	# specified in the list, the flag is [m], if not, the flag is [n]
	#
	##
	def log_vs_match(self,log,list):
		
		ldict = futil.file_to_dict(list)
		
		query  = ""
		subjt  = ""
		inp    = open(log,"r")
		oup    = open(log+".log","w")
		inline = inp.readline()
		while inline != "":
			if inline.find("Query=") != -1:
				query = inline[inline.find(" ")+1:-1]
				oup.write(inline)
			elif inline.find("[") != -1:
				subjt = inline[inline.find(":")+1:inline.find("(")]
				if subjt in ldict:
					if ldict[subjt] == query:
						oup.write("[m]%s" % inline[inline.find("]")+1:])
					else:
						oup.write("[n]%s" % inline[inline.find("]")+1:])
				else:
					oup.write("[-]%s" % inline[inline.find(" "):])
			else:
				oup.write(inline)			
			inline = inp.readline()		
				

	def delete(self,blast,glist):
		
		print("Read gene list...")
		gdict = futil.file_to_dict(glist,0)
		
		inp = open(blast,"r")
		inl = inp.readline()
		oup = open(blast+".mod","w")
		countL = 0
		print("Process line:")
		while inl != "":
			if countL % 1000000 == 0:
				print(" %i x 1m" % (countL/1000))
			L   = inl.split("\t")
			if L[0] not in gdict and L[1] not in gdict:	
				oup.write(inl)
				
			countL += 1
			inl = inp.readline()
			
		inp.close()
		oup.close()
		
		
	#
	# Simply rid of lines with -log(E value) above threshold
	#
	def threshold(self,blast,T):

		inp = open(blast,"r")
		oup = open(blast+"_T%i.out" % T,"w")
		inl = inp.readline()
		countT = countQ = 0
		while inl != "":
			if countT % 10000 == 0:
				print(" %i x 10k" % (countT/10000))
			countT += 1
			llist = inl.split("\t")
			evalue = llist[10]
			if evalue[0] == "e":
				evalue = "1"+evalue
			if evalue == "0.0":
				evalue = "1e-200"
			
			try:
				ef = -math.log10(float(evalue))
			except OverflowError:
				print("Overflow:",evalue)
				sys.exit(0)
			if ef > float(T):
				oup.write(inl)
				countQ += 1
				
			inl = inp.readline()
		
		print("Total %i scores, %i qualfied at threshold %i" % (countT,countQ,T))
		
	#
	# BLAT output can be parsed by parse_align2. Two problems:
	#  1) query sequence name has an extra space in front
	#  2) its orientation convetion is subject-based.
	#
	# So things may break down. this function change takes the output of
	# parse_align2 which originally take blat -out=blast file.
	#
	def fix_blat(self,blat):
		
		inp = open(blat)
		oup = open(blat+".fixed","w")
		inl = inp.readline()
		while inl != "":
			# 0   1   2  3   4 5 6  7  8  9  10     11
			# qid sid ID gap x x qL qR sL sR evalue score
			L = inl.split("\t")
			qL = int(L[6])
			qR = int(L[7])
			sL = int(L[8])
			sR = int(L[9])
			
			# correct query name
			if L[0][0] == " ":
				L[0] = L[0][1:]
			
			# correct order
			if qL < qR:
				oup.write(inl)
			else:
				oup.write("%s\t%i\t%i\t%i\t%i\t%s" % \
						(string.joinfields(L[:6],"\t"),qR,qL,sR,sL,
						 string.joinfields(L[10:],"\t")))
			inl = inp.readline()
			
		
	
	def chain2(self,blast,gapL=0,fragT=80,chainT=90,debug=0):
		pass
	
	#
	# Link the fragments together
	#
	# @param blast  blast output in m=8 format. If this is blat-based, make
	#               sure the output is through fix_blat().
	# @param fastaQ query fasta file
	# @param fastaS subject fasta file
	# @param gapL   the gap length allowed between fragments, default 0
	# @param fragT  the identity lower bound for fragments, those with id<fragT
	#               will not be used for chaining, default 80%
	# @param chainT overall chain identity lower bound, default 90%
	#
	def chain(self,blast,fastaQ,fastaS,gapL=0,fragT=80,chainT=90,debug=0):
		
		print("Get query sizes:",fastaQ)
		qsize = fmanager.get_sizes(fastaQ,1)
		#qsize = {"ATL8C7090 atl7C985":5000}
		#qsize = {"PUT_AT9_1244":66,"PUT_AT9_92766":105}
		#qsize = {"q1":1}
		
		print("Get subject sizes:",fastaS)
		ssize = fmanager.get_sizes(fastaS,1)
		#ssize = {"4":5000}
		#ssize = {"1":33132539,"2":19320864,"4":23328337}
		#ssize = {"s1":1}
		
		print("Read file:", blast)
		inp    = open(blast)
		inl    = inp.readline()
		bdict  = {}
		count = 0
		while inl != "":
			if count%1e5 == 0:
				print(" %i x 100k" % (count/1e5))
			count += 1
			L  = inl.split("\t")
			sL = int(L[8])
			sR = int(L[9])
			qL = int(L[6])
			qR = int(L[7])
			ID = float(L[2])
			if L[0] not in bdict:
				if sL < sR:
					bdict[L[0]] = {L[1]:{sL:[sR,ID,qL,qR,"F"]}}
				else:
					bdict[L[0]] = {L[1]:{sR:[sL,ID,qL,qR,"R"]}}
			else:
				if L[1] not in bdict[L[0]]:
					if sL < sR:
						bdict[L[0]][L[1]] = {sL:[sR,ID,qL,qR,"F"]}
					else:
						bdict[L[0]][L[1]] = {sR:[sL,ID,qL,qR,"R"]}
				else:
					if sL < sR:
						bdict[L[0]][L[1]][sL] = [sR,ID,qL,qR,"F"]
					else:
						bdict[L[0]][L[1]][sR] = [sL,ID,qL,qR,"R"]		
			inl = inp.readline()
		
		for i in bdict:
			print(i)
			for j in bdict[i]:
				print("",j)
				for k in bdict[i][j]:
					print(" ",k,bdict[i][j][k])
				
		# follow only longest chain above fragT
		print("Iterate query...")
		non    = []
		oup    = open("%s_I%iC%iG%i.chain"   %(blast,fragT,chainT,gapL),"w")
		oup1   = open("%s_I%iC%iG%i.allfrag" %(blast,fragT,chainT,gapL),"w")
		oup.write( "Query\tSubjt\t%ID\tqL\tqR\tsL\tsR\tnChains\tqSize\tsSize\tqCoord\tsCoord\tOri\n")
		oup1.write("Query\tSubjt\t%ID\tqL\tqR\tsL\tsR\tOri\n")
		
		countQ = 0  # query
		countS = 0  # subject
		countC = 0  # no qualified chain
		countF = 0  # no qualified frag
		qQ     = {} # qualified query 
		for i in bdict:               # iterate query
			if countQ%1000 == 0:
				print(" %i k" % (countQ/1000))
			countQ += 1
			for j in bdict[i]:        # iterate subj
				countS += 1
				c = list(bdict[i][j].keys())
				c.sort()

				#print c
				chains = []
				single = []
				old = 0
				for k in range(len(c)): # iterate sorted 5' coord
					# the end of the dict
					if debug:
						print("single:",single)
						print("chains:",chains)

					if debug:
						print(i,j,c[k],bdict[i][j][c[k]][1],"->",fragT)
					if k+1 == len(c):
						if single == [] and bdict[i][j][c[k]][1] >= fragT:
							# 0  1  2  3  4  5
							# sL,sR,qL,qR,ID Ori
							single = [[c[k],bdict[i][j][c[k]][0],
											bdict[i][j][c[k]][2],
											bdict[i][j][c[k]][3],
											bdict[i][j][c[k]][1],
											bdict[i][j][c[k]][4],]]				
						if single != []:
							chains.append(single)
						break
					
					# if chain has ID < fragT or
					if bdict[i][j][c[k]][1] < fragT:
						if debug:
							print(" ID < fragT")
						
					# bdict[L[0]][L[1]][sL] = [sR,ID,qL,qR,ori]
					# if the next chain has ID < fragT or
					# if 2nd qL is too far from 1st qR or
					# if adj queires are not ordered
					elif bdict[i][j][c[k+1]][1]<fragT                     or \
					  	 bdict[i][j][c[k+1]][2]-bdict[i][j][c[k]][3]>gapL or \
					  	 c[k+1]-bdict[i][j][c[k]][0]>gapL or \
					  	 bdict[i][j][c[k]][2]>bdict[i][j][c[k+1]][2]:
						case = [0,0,0,0]
						if debug:
							print(" -> next frag not qualified")
							print("   ",c[k+1],bdict[i][j][c[k+1]])
						if bdict[i][j][c[k+1]][1] < fragT:
							if debug:
								print("  case0",bdict[i][j][c[k+1]][1],fragT)
							case[0] = 1
						if bdict[i][j][c[k+1]][0]-bdict[i][j][c[k]][3]>gapL or \
						   bdict[i][j][c[k+1]][0]-bdict[i][j][c[k]][1]>gapL:
							if debug:
								print("  case1",bdict[i][j][c[k+1]][0] - c[k], gapL)
							case[1] = 1
						#if bdict[i][j][c[k]][4] != bdict[i][j][c[k+1]][0]:
						#	case[2] = 1
						if bdict[i][j][c[k]][2] > bdict[i][j][c[k+1]][2]:
							if debug:
								print("  case3",bdict[i][j][c[k]][2],bdict[i][j][c[k+1]][2])
							case[3] = 1
						# if this is empty, add the 1st element, otherwise, it is
						# already in single list
						if single == []:
							single = [[c[k],bdict[i][j][c[k]][0],
											bdict[i][j][c[k]][2],
											bdict[i][j][c[k]][3],
											bdict[i][j][c[k]][1],
											bdict[i][j][c[k]][4]]]
						chains.append(single)
						# reset
						single = []
					# both segments are good
					else:
						#print "both good:",c[k],c[k+1]
						# none added yet
						if single == []:
							single = [[c[k],bdict[i][j][c[k]][0],
											bdict[i][j][c[k]][2],
											bdict[i][j][c[k]][3],
											bdict[i][j][c[k]][1],
											bdict[i][j][c[k]][4]],
								   	  [c[k+1],bdict[i][j][c[k+1]][0],
										  	bdict[i][j][c[k+1]][2],
										  	bdict[i][j][c[k+1]][3],
										  	bdict[i][j][c[k+1]][1],
											bdict[i][j][c[k+1]][4]]]
						# only need to add second
						else:
							single.append([c[k+1],bdict[i][j][c[k+1]][0],
												  bdict[i][j][c[k+1]][2],
												  bdict[i][j][c[k+1]][3],
												  bdict[i][j][c[k+1]][1],
												  bdict[i][j][c[k+1]][4]])
					# END: if statements
				# END: for k (5' coord iteration)
				
				if debug:
					print("\nCHAINS", end=' ')
					print(chains,"\n")
				
				# BEFORE 09/19,05', output only longest chain, now output all
				# Process chains, get idx with max chain length based on 
				# bdict L and R
					
				
				# chain
				# [[sL,sR,qL,qR],[sL,sR,qL,qR],...,[sL,sR,qL,qR]]
				# more than one chain
				if len(chains) != 0:
					for k in range(len(chains)):
						
						# get the minL and maxR
						localsL = localsR = localqL = localqR = -1
						for m in chains[k]:
							if localqL == -1:
								localsL = m[0]
								localsR = m[1]
								localqL = m[2]
								localqR = m[3]
							if localqL > m[2]:
								localsL = m[0]
								localqL = m[2]
							if localqR < m[3]:
								localsR = m[1]
								localqR = m[3]						
						
						# output all frag of and generate composite ID value
						# (query based)
						cLen = 0  # composite length
						cMat = 0  # composite match
						qC   = "" # query coordinates
						sC   = "" # sbjct coordinates
						ori  = "" # orientation of all frag based on sbjct
						for m in chains[k]:
							# Query Subjt %ID qL qR sL sR Ori
							#if k == lenMax:
							oup1.write("%s\t%s\t%f\t%i\t%i\t%i\t%i\t%s\n" % \
											(i,j,m[4],m[2],m[3],m[0],m[1],m[5]))
							cLen += m[3]-m[2]+1
							cMat += (m[3]-m[2]+1)*m[4]/100
							qC   += "%i|%i," % (m[2],m[3])
							sC   += "%i|%i," % (m[0],m[1])
							ori  += m[5] 
						
						# the cLen should not be zero, just in case	
						try:
							cMat/cLen
						except:
							print("ERR: cMat/cLen:%i/%i" % (cMat,cLen))
							print(i,j,chains)
							print("EXIT, CHECK!")
							sys.exit(0)
						
						# chain composite %ID should be larger than chainT
						if cMat/cLen*100 >= chainT:
							if i not in qQ:
								qQ[i] = 1
							else:
								qQ[i] += 1
							# output chrL and chrR of the longest chain
							#if k == lenMax:
							# Query Subjt %ID qL qR sL sR nFrag qSize sSize qCoord sCoord Ori
							oup.write("%s\t%s\t%f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%s\t%s\t%s\n" % \
										(i,j,cMat/cLen*100,localqL,localqR,localsL,localsR,
										 len(chains[k]),qsize[i],ssize[j],
										 qC[:-1],sC[:-1],ori))
						else:
							countC += 1
				else:
					non.append(i)
					countF += 1
					
		print("Total %i queries"   % countQ)
		print("      %i q-s pairs" % countS)
		print("      %i q-s with no frag > fragT"   % countF)
		print("      %i q-s with no chain > chainT" % countC)
		print("      %i qualified queries"          % len(list(qQ.keys())))
		print("      %i q-s with qualified chains"  % (countS-countF-countC))
		

	def pairs(self,blast):
		
		inp = open(blast)
		oup = open(blast+".unique_pairs","w")
		inl = inp.readline()
		pairs = {}
		query = {}
		subjt = {}
		countQ = 0
		countP = 0
		countS = 0
		while inl != "":
			L = inl.split("\t")
			p = "%s\t%s\n" % (L[0],L[1])
			if p not in pairs:
				pairs[p] = 1
				oup.write(p)
				countP += 1			
			if L[0] not in query:
				query[L[0]] = 1
				countQ += 1
			if L[1] not in subjt:
				subjt[L[1]] = 1
				countS += 1
				
			inl = inp.readline()
		
		print("Queries     :",countQ)
		print("Subjects    :",countS)
		print("Unique pairs:",countP)
		
	def get_sp(self,score,sp):
		
		sp  = sp.split(",")
		inp = open(score)
		oup = open(score+".%s" % string.joinfields(sp,"_"),"w")
		inl = inp.readline()
		c = 0
		cQ = 0
		while inl != "":
			c += 1
			L   = inl.split("\t")
			sp1 = L[0][:2]
			sp2 = L[1][:2]
			if sp1 in sp and sp2 in sp:
				cQ += 1
				oup.write(inl)
			inl = inp.readline()
		print("%i scores, %i qualified" % (c,cQ))
	
	def filterout(self,blast,fasta,E,I,P,A):
		
		EV = E  # evalue threshold
		ID = I  # identity threshold
		CV = P  # alignment length coverage threshold
		AL = A  # alignment length threshold		
		
		print("Read sequences")
		fm  = FastaManager.fasta_manager()
		ssize = fm.get_sizes(fasta,1)
		
		print("\nParse blast output")
		inp = open(blast)
		oup = open(blast+".E1I40A30C0.5","w")
		inl = inp.readline()
		c = 0
		cQ = 0
		while inl != "":
			L = inl.split("\t")
			if c % 1e5 == 0:
				print(" %i x 10k" % (c/1e5))
			c += 1
			if L[0] != L[1]:
				qL = (float(L[-5])-float(L[-6])+1)/ssize[L[0]]
				sL = (float(L[-3])-float(L[-4])+1)/ssize[L[1]]
				if float(L[-2]) < EV and float(L[2]) > ID and int(L[3]) > AL:
					if (qL > sL and sL > CV) or (qL < sL and qL > CV):
						oup.write(inl)
						cQ += 1
			inl = inp.readline()		
		print("Total %i pairs, %i qualified" % (c,cQ))	
	
	# assume the first two character are sp abbreviations and only two sp.
	# in the blast output
	def xspecies(self,blast):
	
		print("ASSUMING: 1st 2 char are sp abbrv, only 2 sp")
		inp = open(blast)
		oup = open(blast+".xspecies","w")
		inl = inp.readline()
		sp = {}
		c = 0
		while inl != "":
			if c % 1e5 == 0:
				print(" %i x 10k" % (c/1e5))
			c += 1			
			L = inl.split("\t")
			sp1 = L[0][:2]
			sp2 = L[1][:2]
			if sp1 != sp2:
				oup.write(inl)
				if sp1 not in sp:
					sp[sp1] = {L[0]:1}
				else:
					sp[sp1][L[0]] = 1
				if sp2 not in sp:
					sp[sp2] = {L[1]:1}
				else:
					sp[sp2][L[1]] = 1
			inl = inp.readline()
		print("Species:",list(sp.keys()))
		for i in sp:
			print("",i,len(list(sp[i].keys())))
			
		
	
	####
	# get the top match of each query, assume tabular blast output
	####
	def gettop(self,blast):
		inp = open(blast)
		oup = open(blast+".top","w")
		inl = inp.readline()
		m   = {}
		while inl != "":
			L = inl.split("\t")
			# assume the 1st match of every query is the top match
			# also, any addition match of the same query-subject pair will not
			# be included in the output.
			if L[0] not in m:
				oup.write(inl)
				m[L[0]] = 1
			
			inl = inp.readline()

	def prepare_mcl(self,blast):
		pass
		# parse_table3
		# symmetrify
		# mcl_score
		

	def rmlb(self,astr):
		if astr[-2:] == "\r\n":
			astr = astr[:-2]
		elif astr[-1] == "\n":
			astr = astr[:-1]
		return astr

	def help(self):
		print("\npython ParseBlast.py\n")
		print(" -f  the code to excecute in ParseBlast:")
		print("     parse_db   - parse blast out stored in psql db")
		print("       REQUIRES: -c, -p, OPTIONAL: -T")
		print("     parse_file - parse blast output file")
		print("       REQUIRES: -blast")
		print("       OPTIONAL: -outbase, -T, -target")
		print("     parse_align- parse the alignment part from blast output")
		print("       REQUIRES: -blast, OPTIONAL: -T,-self,-format")
		print("       -style")
		print("     refine - parse the log file to get info on mutiple subject")
		print("       matches, REQUIRES: log")
		print("     parse_align2-parse the output file from defaut to m = 8")
		print("       REQUIRES: -blast")
		print("     parse_align2_blastPlus-parse the output file from Blast+ defaut to m = 8")
		print("       REQUIRES: -blast") 
		print("     parse_align2_blastall-parse the output file from blastall defaut to m = 8")
		print("       REQUIRES: -blast")
		print("     parse_align3-get the blast output of a particular subj.")
		print("       need -blast, -target")
		print("     parse_align4-get the squence pairs out. NEED: blast")
		print("     get_qualified - this is an extension for parse_align which")
		print("       will filter modified log file and cluster the")
		print("       qualified entries. REQUIRES:log,fasta,qtag.OPTIONAL:")
		print("       priority")
		print("     get_qualified2- this is another extension for parse_align")
		print("       this will simply get the qualified subj for one entry. No")
		print("       clustering is done. REQUIRES: log, qtag")
		print("     get_qualified3- this is similar to get_qualified. But takes")
		print("       tabular blast output and 3 threshold settings. REQUIRES:")
		print("       blast, fasta, OPTIONAL: I,L,P")
		print("     get_qualified4- similar to get_qualified3 but no clustering")
		print("       REQUIRES: blast, fasta, OPTIONAL: E,I,L,P,Q")
		print("     gettop - get the top match of every query. NEED: blast")
		print("     parse_table- parse the tabular blast output")
		print("       REQUIRES: -blast, -target, OPTIONAL: -T, verbose, wself")
		print("     parse_table2 - parse tabular output based on idenity, length")
		print("       and/or evalue and generate [Q][count][subjs]. REQUIRES: ")
		print("       -blast, OPTIONAL: E,I,L")
		print("     parse_table3 - no storage involved, simply take blast out")
		print("       then apply threshold (T) to get -log(e). REQUIRES: blast,")
		print("       OPTIONAL: E,I")
		print("     parse_table4 - parse tabular output based on idenity, length")
		print("       and/or evalue. No modification to the lines. REQUIRES: ")
		print("       -blast, OPTIONAL: E,I,L")
		print("     parse_gap  - parse blast out and get gap info. NEED: blast")
		print("     threshold  - apply threshold to a blast output, NEED: blast")
		print("        T")
		print("     symmetrify - symmetrify, normalize, and homogenize scores.")
		print("       REQUIRES: -score. OPTIONAL: -cutoff,-outbase,-homog")
		print("     check_acc  - check and correct accessions") 
		print("       REQUIRES: -acc, -score")
		print("     check_missing - check if query OR subj is missing in the")
		print("       Blast output. NEED: blast, fasta, qors")
		print("     spc_score  - get score for non-paramegnetic clustering")
		print("       REQUIRES: -score. OPTIONAL: -per_id,-homog,-o")
		print("     mcl_score  - get score matrix for MCL")
		print("       REQUIRES: -score, OPTIONAL: -cutoff,-outbase,-homog")
		print("     nei_score  - get score matrix for Neighbor in Phylip")
		print("       REQUIRES: -score (modified by symmetrify)")
		print("     mega_score - get score matrix for MEGA")
		print("       REQUIRES: -score (modified by symmetrify)")
		print("     score_matrix - plain old matrix, NEED: score (symmetrified).")
		print("     select     - get scores for a subset of the sequences used")
		print("       in the original BLAST run")
		print("       REQUIRES: -list,-score")
		print("     extract_cds- get the subject sequence out")
		print("       REQUIRES: -blast. OPTIONAL: stype, T, wu, ignorestop")
		print("     get_subj   - get a non-redundant list of matching subj")
		print("       REQUIRES: -blast, OPTIONAL: desc, style, ttype, T")
		print("     index_names- convert sequence names into indices")
		print("       REQUIRES: -score")
		print("     rename     - Change id within the fasta file")
		print("       REQUIRES: -fasta, -name")
		print("     merge_match- merge the subj coords of overlapping entries")
		print("       in a blast output. REQUIRES: -blast, -feature")
		print("     get_reciprocal - get the reciprocal best matches")
		print("       REQUIRES: blast, OPTIONAL: target, fasta, T, wself, P")
		print("     match_list - get the top match for each query sequence and")
		print("       see if the top match is in a list of names. REQUIRES:")
		print("       -blast, -name")
		print("     match_list2 - based on the passed list file, output")
		print("       qualified entries. REQUIRES: -matrix, -blast")
		print("     log_vs_match - compare the log file against the list")
		print("       generated by match_list. New flags are inserted, see")
		print("       code doc for details. REQUIRES: -log, -list")
		print("     match_fam - get subj with match to members of a specified")
		print("       group of genes. REQUIRES: blast,matrix,target. OPTIONAL:")
		print("       unk")
		print("     delete - delete all scores involving any sequence in a list")
		print("       -blast, -list")
		print("     fix_blat - fix the modified blat -out=blast output by")
		print("       parse_align2, NEED: blast (blat-based)")
		print("     chain - link fragments together. NEED: blast,fastaQ,fastaS")
		print("       OPTIONAL: gapL,fragT,chainT,debug")
		print("     pairs - get unique pairs and IDs. NEED: blast")
		print("     get_sp - get scores of particular species, NEED: score,sp")
		print("     filterout - filter blastouput. NEED: blast,fasta. OPTIONAL:")
		print("       E,I,P,A")
		print("     xspecies - get only cross species matches. NEED: blast")  
		print("     for_mega - do blast, parse score for mega. NEED: blast,")
		print("       fasta")
		print("     gettop - get top matches. NEED: blast")
		print(" -c    configuration")
		print(" -p    specify the block to excecute within config")
		print(" -log  log file generated by parse_align")
		print(" -blast blast output file, not required if -f parse_db")
		print("       for_mega, this is the blast program dir")
		print(" -blast1 for get_reciprocal, first blast output in table format")
		print(" -blast2 for get_reciprocal, second output in table format")
		print(" -cluster cluster file")
		print(" -score score file, required for mcl_score and check_acc")
		print(" -list file with a list of seq_id")
		print(" -cid  cluster id")
		print(" -o    base name for output, default=[blast]")
		print(" -acc  name for fasta file with correct accessions")	   
		print(" -cutoff cutoff value for the matrix, default = 0")
		print(" -ttype threshold type, identity[default], or evalue")
		print(" -stype sequence type, protein [1] or nucleotide [0,default]")
		print(" -T    -log(E) threshold for parsing blast scores, default 0")
		print("       or it will be compared against percentID or SIM, depend")
		print("       on which is specified in target.")
		print("       -> For extract_cds, this is the allowed gap size between")
		print("          stretches")
		print("       -> For get_subj, ID or E threshold, not transformed")
		print("       -> For get_reciprocal, length threshold")
		print(" -outbase output base name")
		print(" -H    homogenize by taking square [1],square root[2], or no")
		print("	      homogenization at all [0,default]. For e value, use 1.")
		print(" -per_id  number of top scores to get for SPC, default 10")
		print(" -stype type of subj, pep [0] or nucleotide [1,default]")
		print(" -ignorestop default no [0]")
		print(" -target for parse_file: the kind of score to parse -")
		print("	      evalue[default],percentID or percentSIM")
		print("	      for parse_table can be evalue, percentID, bit, or")
		print("	      token no. separated by ','")
		print(" -desc parse matching subj from descriptor lines [1] or from")
		print("	      alignment [0, default]")
		print(" -style the BLAST output (-m option) style, 1 (default) or 9")
		print("	      In parse_align, [0, def] means no nt involved, [1] yes")
		print(" -name for renaming id within fasta file in [new][old] format")
		print("       for match_list, it is simply a list of names, one column")
		print(" -wself include self match [1] or not [0, default]. For")
		print("       get_reciprocal, means match from the same organisms will")
		print("       kept [1] or not [0, default]. If 1, the sequence names")
		print("       HAVE to have sp abbrv followed by an underscore")
		print(" -format output align part [0,default] or only the converted")
		print("       match line [1]")
		print(" -feature The common feature for the query sequences in generat-")
		print("       ing the blast output")
		print(" -fasta For get_qualified. This file is for comparing sizes of")
		print("       cluster members, so make sure the relevant file is used")
		print(" -priority for get_qualified. Whether query sequence should take")
		print("       precedance [1] or not [0,default]")
		print(" -verbose for parse_table, get qL,qR,sL,sR")
		print(" -lenT % length treshold. If this is specified, need BOTH fasta")
		print("       and qors.")
		print(" -qors for parse_table, based on query [0, default] or subj [1]")
		print("       to apply lenT.")
		print(" -matrix file with [seq_id][group_designation] or [taxa][seq_id]")
		print("       for match_list2")
		print(" -target the target group")
		print(" -qtag for get_qualified, specify qualifying tags separated by")
		print("       ','.")
		print(" -unk  sequences of unknown group are query [1, default] or subj")
		print("       [0]")
		print(" -wu   WU format BLAST [1] or not [0, default]")
		print(" -E    evalue threshold in -log(e), default -1")
		print(" -I    identity threshold in %. Default 95")
		print(" -L    length threshold, default 150.")
		print(" -A    alignment length threshold, default 30")
		print(" -P    threshold for the proportion of match length vs sequence")
		print("       length, default 0.9")
		print(" -Q    for match length vs. sequence length comparison, use both")
		print("       query & subject [0], query only [1], or subject only [2]")
		print(" -fragT fragment %ID lower bound")
		print(" -chianT chain %ID lower bound")
		print(" -gapL max gap between fragment allowed")
		print(" -fastaQ query fasta file")
		print(" -fastaS subject fasta file")
		print(" -sp   2 letter species abbr. separated by ','")
		print(" -debug display debug string [1], default no [0]")
		print("")
		sys.exit(0)

			
#-------------------------------------------------------------------------------
# Function calls
#-------------------------------------------------------------------------------


if __name__ == '__main__':
	configFile = operation = function = blast = accfile = cluster = outbase = \
				 cid = score = list = name = log = feature = fasta = blast1 = \
				 blast2 = matrix = target = qtag = lenT = ttype = fastaQ = \
				 fastaS = sp = qors = ""

	target     = "evalue"
	homogenize = desc = wself = qors = unk = cutoff = format = priority = \
				 verbose = wu = debug = 0
	T	       = 0.0
	per_id     = 10
	stype       = 1
	style      = 1
	futil      = FileUtility.file_util()
	fmanager   = FastaManager.fasta_manager()
	link       = SingleLinkage.single_linkage()
	E          = 0
	I          = 0
	L          = 150
	P          = 0.9
	Q          = 0
	A          = 30
	fragT      = 80
	chainT     = 90
	gapL       = 0
	ignorestop = 0

	parse = parser()

	for i in range(1,len(sys.argv),2):
		if sys.argv[i] == "-c":
			configFile = sys.argv[i+1]
		elif sys.argv[i] == "-p":
			operation  = sys.argv[i+1]
		elif sys.argv[i] == "-o" or sys.argv[i] == "-outbase":
			outbase    = sys.argv[i+1]
		elif sys.argv[i] == "-blast":
			blast     = sys.argv[i+1]
		elif sys.argv[i] == "-acc":
			accfile    = sys.argv[i+1]
		elif sys.argv[i] == "-cluster":
			cluster    = sys.argv[i+1]
		elif sys.argv[i] == "-score":
			score      = sys.argv[i+1]
		elif sys.argv[i] == "-f":
			function   = sys.argv[i+1]      
		elif sys.argv[i] == "-cid":
			cid	= int(sys.argv[i+1])
		elif sys.argv[i] == "-H":
			homogenize = int(sys.argv[i+1])
		elif sys.argv[i] == "-per_id":
			per_id     = int(sys.argv[i+1])
		elif sys.argv[i] == "-cutoff":
			cutoff     = float(sys.argv[i+1])
		elif sys.argv[i] == "-T":
			T	  = float(sys.argv[i+1])
		elif sys.argv[i] == "-list":
			list       = sys.argv[i+1]
		elif sys.argv[i] == "-stype":
			stype       = int(sys.argv[i+1])
		elif sys.argv[i] == "-target":
			target     = sys.argv[i+1]
		elif sys.argv[i] == "-desc":
			desc       = int(sys.argv[i+1])
		elif sys.argv[i] == "-tokens":
			tokens     = sys.argv[i+1]
		elif sys.argv[i] == "-style":
			style      = int(sys.argv[i+1])
		elif sys.argv[i] == "-format":
			format     = int(sys.argv[i+1])
		elif sys.argv[i] == "-name":
			name       = sys.argv[i+1]
		elif sys.argv[i] == "-wself":
			wself      = int(sys.argv[i+1])
		elif sys.argv[i] == "-log":
			log	       = sys.argv[i+1]
		elif sys.argv[i] == "-feature":
			feature    = sys.argv[i+1]
		elif sys.argv[i] == "-fasta":
			fasta      = sys.argv[i+1]
		elif sys.argv[i] == "-blast1":
			blast1     = sys.argv[i+1]
		elif sys.argv[i] == "-blast2":
			blast2     = sys.argv[i+1]
		elif sys.argv[i] == "-priority":
			priority   = sys.argv[i+1]
		elif sys.argv[i] == "-verbose":
			verbose    = int(sys.argv[i+1])
		elif sys.argv[i] == "-qors":
			qors       = int(sys.argv[i+1])
		elif sys.argv[i] == "-lenT":
			lenT       = int(sys.argv[i+1])
		elif sys.argv[i] == "-matrix":
			matrix     = sys.argv[i+1]
		elif sys.argv[i] == "-target":
			target     = sys.argv[i+1]
		elif sys.argv[i] == "-qtag":
			qtag       = sys.argv[i+1]
		elif sys.argv[i] == "-unk":
			unk        = int(sys.argv[i+1])
		elif sys.argv[i] == "-wu":
			wu         = int(sys.argv[i+1])
		elif sys.argv[i] == "-E":
			E          = float(sys.argv[i+1])
		elif sys.argv[i] == "-I":
			I          = float(sys.argv[i+1])
		elif sys.argv[i] == "-L":
			L          = int(sys.argv[i+1])
		elif sys.argv[i] == "-P":
			P          = float(sys.argv[i+1])
		elif sys.argv[i] == "-A":
			A          = int(sys.argv[i+1])
		elif sys.argv[i] == "-Q":
			Q          = int(sys.argv[i+1])
		elif sys.argv[i] == "-ttype":
			ttype      = sys.argv[i+1]
		elif sys.argv[i] == "-gapL":
			gapL       = int(sys.argv[i+1])
		elif sys.argv[i] == "-fragT":
			fragT      = int(sys.argv[i+1])
		elif sys.argv[i] == "-chainT":
			chainT     = int(sys.argv[i+1])
		elif sys.argv[i] == "-fastaQ":
			fastaQ      = sys.argv[i+1]
		elif sys.argv[i] == "-fastaS":
			fastaS      = sys.argv[i+1]
		elif sys.argv[i] == "-ignorestop":
			ignorestop  = int(sys.argv[i+1])
		elif sys.argv[i] == "-sp":
			sp      = sys.argv[i+1]
		elif sys.argv[i] == "-debug":
			debug       = sys.argv[i+1]
		else:
			print("Unknown parameter:",sys.argv[i])
			parse.help()
			
	if function == "parse_db":
		if configFile == "" or operation == "":
			print("\nNeed to define config file and operation\n")
			
			print(" -help for help")
			sys.exit(0)
		dbtask = DatabaseOp
		config = dbtask.configConnect(configFile,operation)
		parse  = parser(dbtask,config)
		parse.parse_blast_db(outbase,T)
	elif function == "for_mega":
		if blast == "" or fasta == "":
			print("\nNeed to specify the blast program dir and fasta\n")	  
			print(" -help for help")
			sys.exit(0)		    
		parse.for_mega(blast,fasta)
	elif function == "xspecies":		
		if blast == "":
			print("\nNeed to specify the blast output file to be parsed\n")	  
			print(" -help for help")
			sys.exit(0)		    
		parse.xspecies(blast)
	elif function == "parse_file":		
		if blast == "":
			print("\nNeed to specify the blast output file to be parsed\n")	  
			print(" -help for help")
			sys.exit(0)		    
		parse.parse_blast_file(blast,T,target)
	elif function == "parse_align":		
		if blast == "":
			print("\nNeed to specify the blast output file to be parsed\n")	  
			print(" -help for help")
			sys.exit(0)	    
		parse.parse_align(blast,T,wself,format,style)

	elif function == "parse_align4":		
		if blast == "":
			print("\nNeed to specify the blast output\n")	  
			print(" -help for help")
			sys.exit(0)	    
		parse.parse_align4(blast)

	elif function == "refine":		
		if log == "":
			print("\nNeed log file.\n")	  
			print(" -help for help")
			sys.exit(0)	    
		parse.refine(log)
	elif function == "parse_align2":		
		if blast == "":
			print("\nNeed to specify the blast output file to be parsed\n")	  
			print(" -help for help")
			sys.exit(0)	    
		parse.parse_align2(blast)
	elif function == "parse_align2_blastPlus":
		if blast == "":
			print("\nNeed to specify the blast output file to be parsed\n")
			print(" -help for help")
			sys.exit(0)
		parse.parse_align2_blastPlus(blast)
	elif function == "parse_align2_blastall":
		if blast == "":
			print("\nNeed to specify the blast output file to be parsed\n")
			print(" -help for help")
			sys.exit(0)
		parse.parse_align2_blastall(blast)
	elif function == "parse_align3":		
		if blast == "" or target == "":
			print("\nNeed the blast output and target id\n")	  
			print(" -help for help")
			sys.exit(0)	    
		parse.parse_align3(blast,target)
	elif function == "get_qualified":		
		if log == "" or fasta == "" or qtag == "":
			print("\nNeed .log, fasta files and qtag\n")
			print(" -help for help")
			sys.exit(0)
		parse.get_qualified(log,fasta,qtag,priority)
	elif function == "get_qualified2":		
		if log == "" or qtag == "":
			print("\nNeed parse_align output and qtag\n")
			print(" -help for help")
			sys.exit(0)
		parse.get_qualified2(log,qtag)
	elif function == "get_qualified3":		
		if blast == "" or fasta == "":
			print("\nNeed blast output and fasta\n")
			print(" -help for help")
			sys.exit(0)
		parse.get_qualified3(blast,fasta,I,L,P)
	elif function == "get_qualified4":		
		if blast == "" or fasta == "":
			print("\nNeed blast output and fasta\n")
			print(" -help for help")
			sys.exit(0)
		parse.get_qualified4(blast,fasta,E,I,L,P,Q)
	elif function == "parse_table":		
		if blast == "" or target == "":
			print("\nNeed the blast output and target to be parsed\n")	       
			print(" -help for help")
			sys.exit(0)	    
		parse.parse_table(blast,target,T,verbose,wself,lenT,fasta,qors)
	elif function == "parse_table2":
		if blast == "":
			print("\nNeed the blast output\n")	       
			print(" -help for help")
			sys.exit(0)	    
		parse.parse_table2(blast,E,I,L)
	elif function == "parse_table3":		
		if blast == "":
			print("\nNeed the blast output\n")	       
			print(" -help for help")
			sys.exit(0)	    
		parse.parse_table3(blast,E,I)
	elif function == "parse_table4":		
		if "" in [blast,E,I,L]:
			print("\nNeed the blast, E, I, L\n")	       
			print(" -help for help")
			sys.exit(0)	    
		parse.parse_table4(blast,E,I,L)
	elif function == "check_acc":		
		if score == "" or accfile == "":
			print("\nNeed fasta file with correct acc and the parsed file\n")
			print(" -help for help")
			sys.exit(0)
		parse.check_acc(accfile,score)
	elif function == "check_missing":		
		if blast == "" or fasta == "" or qors == "":
			print("\nNeed fasta, blast, and qors\n")
			print(" -help for help")
			sys.exit(0)
		parse.check_missing(fasta,blast,qors)
	elif function == "symmetrify":		
		if score == "":
			print("\nNeed score file\n")
			print(" -help for help")
			sys.exit(0)
		parse.symmetrify(score,outbase,cutoff,homogenize)
	elif function == "spc_score":		
		if score == "":
			print("\nNeed blast ouput, cluster file, and cluster id\n")
			print(" -help for help")
			sys.exit(0)
		parse.get_scores_for_spc(score,outbase,per_id,homogenize)
	elif function == "mcl_score":
		if score == "":
			print("\nNeed score file\n")
			print(" -help for help")
			sys.exit(0)
		parse.get_scores_for_mcl(score,outbase,cutoff,homogenize)
	elif function == "nei_score":		
		if score == "":
			print("\nNeed symmetrified score file\n")
			print(" -help for help")
			sys.exit(0)
		parse.get_scores_for_neighbor(score)
	elif function == "mega_score":		
		if score == "":
			print("\nNeed symmetrified score file\n")
			print(" -help for help")
			sys.exit(0)
		parse.mega_score(score)
	elif function == "score_matrix":		
		if score == "":
			print("\nNeed symmetrified score file\n")
			print(" -help for help")
			sys.exit(0)
		parse.score_matrix(score)
	elif function == "select":
		if score == "" or list == "":
			print("\nNeed score file and seq_id list\n")
			print(" -help for help")
			sys.exit(0)
		parse.get_selected(list,score)
	elif function == "extract_cds":		
		if blast == "":
			print("\nNeed blast output\n")
			print(" -help for help")
			sys.exit(0)
		if T == 1.0:
			T = 0
		parse.extract_cds(blast,stype,T,wu,ignorestop)
	elif function == "get_subj":
		
		if blast == "":
			print("\nNeed blast output\n")
			print(" -help for help")
			sys.exit(0)
		parse.get_subj(blast,desc,style,ttype,T)

	elif function == "index_names":
		
		if score == "":
			print("\nNeed score file\n")
			print(" -help for help")
			sys.exit(0)
		parse.index_names(score)

	elif function == "rename":
		
		if fasta == "" or name == "":
			print("\nNeed fasta and name files\n")
			print(" -help for help")
			sys.exit(0)
		parse.rename(fasta,name)	

	elif function == "merge_match":
		
		if blast == "" or feature == "":
			print("\nNeed blast output and feature designation\n")
			print(" -help for help")
			sys.exit(0)
		parse.merge_match(blast,feature)	

	elif function == "get_reciprocal":
		
		if blast == "":
			print("\nNeed blast output\n")
			print(" -help for help")
			sys.exit(0)
		parse.get_reciprocal(blast,target,fasta,P,wself)	

	elif function == "match_list":
		
		if blast == "" or name == "":
			print("\nNeed two blast output and names\n")
			print(" -help for help")
			sys.exit(0)
		parse.match_list(blast,name)	

	elif function == "match_list2":
		
		if blast == "" or matrix == "":
			print("\nNeed two blast output and taxa-name matrix\n")
			print(" -help for help")
			sys.exit(0)
		parse.match_list2(blast,matrix)	

	elif function == "log_vs_match":
		
		if log == "" or list == "":
			print("\nNeed log file and match list\n")
			print(" -help for help")
			sys.exit(0)
		parse.log_vs_match(log,list)

	elif function == "match_fam":
		
		if blast == "" or matrix == "" or target == "":
			print("\nNeed blast output, matrix file, and target groupt\n")
			print(" -help for help")
			sys.exit(0)
		parse.match_fam(blast,matrix,target,unk)	
	elif function == "delete":
		
		if blast == "" or list == "":
			print("\nNeed blast output and a list of names\n")
			print(" -help for help")
			sys.exit(0)
		parse.delete(blast,list)	
	elif function == "threshold":
		
		if blast == "" or T == "":
			print("\nNeed blast output and threshold\n")
			print(" -help for help")
			sys.exit(0)
		parse.threshold(blast,T)	
	elif function == "parse_gap":
		
		if blast == "":
			print("\nNeed blast output\n")
			print(" -help for help")
			sys.exit(0)
		parse.parse_gap(blast)
	elif function == "fix_blat":		
		if blast == "":
			print("\nNeed blast output\n")
			print(" -help for help")
			sys.exit(0)
		parse.fix_blat(blast)
	elif function == "chain":		
		if "" in [blast,fastaQ,fastaS]:
			print("\nNeed blast output, fastaQ, fastaS\n")
			print(" -help for help")
			sys.exit(0)
		parse.chain(blast,fastaQ,fastaS,gapL,fragT,chainT,debug)
	elif function == "pairs":		
		if blast == "":
			print("\nNeed blast output\n")
			print(" -help for help")
			sys.exit(0)
		parse.pairs(blast)
	elif function == "get_sp":
		if score == "" or sp == "":
			print("\nNeed score file and species\n")
			print(" -help for help")
			sys.exit(0)
		parse.get_sp(score,sp)
	elif function == "filterout":
		if "" in [blast,fasta]:
			print("\nNeed to specify blast and fasta\n")	  
			print(" -help for help")
			sys.exit(0)		    
		parse.filterout(blast,fasta,E,I,P,A)
	# 8/10,07
	elif function == "gettop":
		if blast == "":
			print("\nNeed blast output")
			print(" -help for help")
			sys.exit(0)		    
		parse.gettop(blast)		
	else:
		print("\nNo such function defined: '%s'\n" % function)
		
		print(" -help for help")
		sys.exit(0)

'''
	# This is substantially revised on 4/14,08' for parsing the NCBI blast
	# text output generated using the web form. Not sure if the revision will
	# work on regular version or not. So here we are.
	##
	# This function gets the subject part of the alignment. But not that simple
	# though. Since most likely there will be multiple matches, the best thing
	# to do would be doing a multiple sequence alignment for each subject and
	# determine the consensus. All sequences of one subj are most likely the
	# same ori for the EST contig this module is designed for. In case it is 
	# not, some screwy situation WILL arise. This is NOT dealt with yet.
	#
	# Now, in order to make this application useful for genomic sequences. It is
	# allowed to have more than one stretch for the same subject sequence. 
	# 
	# Two outputs are generated:
	#  1. "_ext.fa" - all sequences, regardless of the length
	#  2. "_assm.fa" - the assembled one
	#
	# @param stype  sequence in protein [0] or nucleotide [1,default] coords
	# @param wu    wash u blast or not, default 0 (not).
	# @param ignorestop don't care about stop [1] or exclude it [0,defualt]
	##
	def extract_cds(self,blast,stype,T,wu,ignorestop=0):
	
		if stype:
			INC = 3
		else:
			INC = 1
		
		print "Start extract coding sequences, query:"

		#
		# this part gets all subj sequences out of the blast output file
		#
		inp = open(blast,"r")
		oup = open(blast+"_ext.fa","w")
		inline = inp.readline()
		subj = seq = L = R = ""
		hasL = 0
		cdict = {}
		
		countSubj = 0
		while inline != "":
			if inline[:7] == "Query= ":
				print "",inline[7:-1]
				pass
			# the need to find "score =" is due to the presence of multiple
			# matching areas that are disjuct.
			elif inline[0] == ">" or inline.find("Score =") != -1:
				# if sequence is not empty and ">" or "score =" is encountered
				# this signify the start of a new entry, so put the preivious
				# one into dict
				if seq != "":
					print "seq:",[seq]
					if L == '':
						print "Is this WU blast? Need to set the flag!"
						sys.exit()
					# rid of non-alpha characters
					if seq.find("-") != -1:
						segments = seq.split("-")
						seq = ""
						for i in segments:
							seq = seq + i
					if not ignorestop and seq.find("*") != -1:
						segments = seq.split("*")
						seq = ""
						for i in segments:
							seq = seq + i				   
					
					# store subj into dict and write sequence
					if not cdict.has_key(subj):
						countSubj = countSubj+1
						
						if int(L)<int(R):
							cdict[subj] = [[int(L),int(R)]]
						else:
							cdict[subj] = [[-int(L),-int(R)]]
					else:
						if int(L)<int(R):
							cdict[subj].append([int(L),int(R)])
						else:
							cdict[subj].append([-int(L),-int(R)])	   
					oup.write(">"+subj+"_"+L+"_"+R+"\n"+seq+"\n")
					
					# reset
					seq = L = R = ""					
					# only reset subj if new subj is present
					if inline[0] == ">":
						subj = ""										       
					hasL = 0
									
				if inline[0] == ">":
					subj = self.rmlb(inline)[1:]
							
			elif inline.find("Sbjct") != -1:
				llist = inline.split(" ")
				print llist
				if not hasL:
					# if wu blast, L coord is the second non-empty element.
					if wu:
						for j in llist[1:]:
							if j != "":
								L = j
								break
					else:
						L = llist[1]
					R = llist[-1][:-1]
					hasL = 1
				else:
					R = llist[-1][:-1]   # keep getting R   
				
				#print [L,R]
				
				# this is for situation where the left coord is sticked to seq
				# no space in between			   
				alphaIndex = 0
				try:
					if llist[-2][0].isdigit():
						for i in range(len(llist[-2])):
							if llist[-2][i].isalpha():
								alphaIndex = i
								break
				# there are some situations where the sbjct line contain only
				# coord and no sequence, these lines will be ignored.
				except IndexError:
					pass
									
				seq = seq + llist[-2][alphaIndex:]
								
			inline = inp.readline()	 
		
		# rid of non-alpha characters in last sequence
		if seq.find("-") != -1:
			segments = seq.split("-")
			seq = ""
			for i in segments:
				seq = seq + i
		if not ignorestop and seq.find("*") != -1:
			print "STOP"
			segments = seq.split("*")
			seq = ""
			for i in segments:
				seq = seq + i				   
		
		# store subj into dict and write sequence for last
		if not cdict.has_key(subj):
			countSubj = countSubj+1
			
			if int(L)<int(R):
				cdict[subj] = [[int(L),int(R)]]
			else:
				cdict[subj] = [[-int(L),-int(R)]]
		else:
			if int(L)<int(R):
				cdict[subj].append([int(L),int(R)])
			else:
				cdict[subj].append([-int(L),-int(R)])
		
		oup.write(">"+subj+"_"+L+"_"+R+"\n"+seq+"\n")
		
		# Separate + and - orientations
		tdict = {}
		for i in cdict.keys():
			countW = 0
			countC = 0
			wlist = []
			clist = []
			for j in range(len(cdict[i])):
				if abs(cdict[i][j][0]) < abs(cdict[i][j][1]):
					countW += 1
					wlist.append(cdict[i][j])
				elif abs(cdict[i][j][0]) > abs(cdict[i][j][1]):
					countC += 1
					clist.append(cdict[i][j])
			if countW > 0:		
				if tdict.has_key(i+"_0"):
					tdict[i+"_0"].append(wlist)
				else:
					tdict[i+"_0"] = wlist
			if countC > 0:
				if tdict.has_key(i+"_1"):		
					tdict[i+"_1"].append(clist)
				else:
					tdict[i+"_1"] = clist
		
		cdict = tdict
		#
		# consolidate the coordinates for sequences in each subj
		#
		print "\nConsolidate coordinates:"
		for i in cdict.keys():
			
			#print "Elements:",len(cdict[i])
			idxL = idxR = 0
			mdict = {} # idx of element which is eliminated		 
			
			# find entries encompassed by other sequences
			for j in range(len(cdict[i])):
				
				#print j
				if mdict.has_key(j):
					#print " skipped1, in mdict already"
					continue
					
				for k in range(j+1,len(cdict[i])):
					
					if mdict.has_key(k):
						#print " skipped2, in mdict already"
						continue		
					#print "",cdict[i][j][0],":",cdict[i][j][1],"vs",\
					#		cdict[i][k][0],":",cdict[i][k][1]
					if cdict[i][j][0] < cdict[i][k][0]:      # jL<kL
						if cdict[i][j][1] >= cdict[i][k][1]: # and jR>=kR, j(k)
							mdict[k] = 1
							#print "  j(k)"
						else:
							pass
							#print "  j/k"
					elif cdict[i][j][0] > cdict[i][k][0]:    # jL>kL
						if cdict[i][j][1] <= cdict[i][k][1]: # and jR<=kR, k(j)
							mdict[j] = 1 
							#print "  k(j)"
							break					   
						else:
							pass
							#print "  j/k"
					else:							         # jL = kL
						if cdict[i][j][1] > cdict[i][k][1]:  # and jR>kR, j(k)
							mdict[k] = 1
							#print "  j(k)"
						elif cdict[i][j][1]<cdict[i][k][1]:  # and jR<kR, k(j)
							mdict[j] = 1
							#print "  k(j)"
							break   
						else:						    # and jR=kR, k=j
							mdict[k] = 1    
							#print "  k=j"  
			
			#print cdict[i]
			#print mdict
			# clean up cdict
			tlist = []
			for j in range(len(cdict[i])):
				if not mdict.has_key(j):
					tlist.append(cdict[i][j])
			cdict[i] = tlist
			#print cdict[i]
			
			# order elements in cdict
			# L coord as key, element idx as value
			Ldict = {}
			for j in range(len(cdict[i])):
				Ldict[cdict[i][j][0]] = j
			keys = Ldict.keys()
			keys.sort()
			clist = []
			for j in keys:
				clist.append(cdict[i][Ldict[j]])
			
			# define contribution of each sequence: L,R,contributeL
			for j in range(len(clist)):
				if j == 0:
					clist[j].append(clist[j][0])
				else:
					clist[j].append(clist[j-1][1]+INC)		      
			cdict[i] = clist
			
			#print i,clist
		
		#
		# load and assemble sequences
		#
		print "\nAssemble sequences:"
		
		# load sequences into memory... this will be a problem...
		inp = open(blast+"_ext.fa","r")
		oup = open(blast+"_assm.fa","w")
		inline = inp.readline()
		sdict = {}

		countAssemblee = 0
		while inline != "":
			if inline[0] == ">":			    
				inline = inline[1:]
				#print "",inline[:-1]
				llist  = inline.split("_")
				#print llist
				subj   = ""
				L      = int(llist[-2])
				R      = int(llist[-1])
				for j in llist[:-2]:
					subj = subj + j + "_"
				subj = subj[:-1]
				
				# make sure the name has ori info
				if L < R:
					subj += "_0"
				else:
					subj += "_1"
						
				# check each cdict element and add sequence stretch
				for j in range(len(cdict[subj])):
					
					# make sure that this subj is an element in cdict[j]
					# also find out its index
					if(L!= abs(cdict[subj][j][0]) or \
					   R!= abs(cdict[subj][j][1])):
						continue
					
					countAssemblee = countAssemblee + 1
					seq = inp.readline()[:-1]
					#print "1>"
					#print cdict[subj][j][0],cdict[subj][j][2]
					#print seq
					if cdict[subj][j][0] < cdict[subj][j][2]:
						seq = seq[(cdict[subj][j][2] - \
								   cdict[subj][j][1] + 1)/INC - 1:]
						#print "High: %i-%i" % (cdict[subj][j][2],cdict[subj][j][1])
						#print seq
							
					# the contributed L coord is smaller due to the
					# presence of un-matched area within subj, add a
					# special character "^" 
					elif cdict[subj][j][0] > cdict[subj][j][2]:
						seq = "(%i-%i)" % (abs(cdict[subj][j][2]),
										   abs(cdict[subj][j][0])) + seq
						#print "Low"
						#print seq
					elif cdict[subj][j][0] == cdict[subj][j][2]:
						seq = seq
					# this is a very lame fix... I notice that some fused entry
					# has one extra residue in sequences not going through the
					# above two if-elif statements. This fix works if there
					# are sequences after this. But if this is the only segment,
					# then I will lose one residue.
					else:
						seq = seq[:-1]
						#print "Lower"
						#print seq
									
						
					cdict[subj][j].append(seq)
																
			inline = inp.readline()


		# output assembled sequence
		print "\nOutput assembled sequences..."
		
		countAssembled = 0
		for i in cdict.keys():
			#print cdict[i]
			L = abs(cdict[i][0][0])
			R = abs(cdict[i][0][1])
			seq = ""
			slist = []
			for j in range(len(cdict[i])):
				idx = cdict[i][j][3].find(")")
				disrupt = 0
				if idx != -1:
					gap = cdict[i][j][3][1:idx].split("-")
					if abs(int(gap[0])-int(gap[1])) > T:
						disrupt = 1
						
				if disrupt:					
					slist.append([i,L,R,seq])
					countAssembled = countAssembled+ 1
					# reset
					L = abs(cdict[i][j][0])
					R = abs(cdict[i][j][1])
					# rid of the connection info in front
					seq = cdict[i][j][3][idx+1:]
				else:
					seq = seq + cdict[i][j][3]
					R = abs(cdict[i][j][1])
			
			# the last sequence
			slist.append([i,L,R,seq])
			countAssembled = countAssembled+ 1
			
			for j in slist:
				oup.write(">%s_%i_%i\n%s\n" % (j[0],j[1],j[2],j[3]))
		
		print "\nTotal %i unique subjects"		   % countSubj
		print "      %i sequences to be assembled" % countAssemblee
		print "      %i assembled sequences"	% countAssembled
		print "Done!\n"

'''
