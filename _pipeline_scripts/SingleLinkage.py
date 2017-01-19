
import sys, string


class single_linkage:

	def __init__(self):
		pass

	def get_group_R(self,G,R,R1,R2):
		
		print "Group   :",G
		print "Relation:",R
		
		# entry id as key, group as value
		print "Read groups..."
		gdict = {}
		inp = open(G,"r")
		inl = inp.readline()
		while inl != "":
			L = self.rnlb(inl).split("\t")
			if gdict.has_key(L[1]):
				print "multiple assign:",L[1]
			else:
				gdict[L[1]] = L[0]
			inl = inp.readline()
		
		# group as key, a list of [id1,id2] as value
		print "Read relations..."
		edict = {}
		inp = open(R,"r")
		#first line is description
		inp.readline()
		inl = inp.readline()
		while inl != "":
			L = self.rnlb(inl).split("\t")
			# verify if this pair is from the same group
			if gdict[L[0]] == gdict[L[1]]:
				if edict.has_key(gdict[L[0]]):
					edict[gdict[L[0]]].append([L[0],L[1]]) 
				else:
					edict[gdict[L[0]]] = [[L[0],L[1]]]
			inl = inp.readline()
		
		# go through groups and get clusters
		print "Cluster each group..."
		oup1 = open(R+".clusters","w")
		oup1.write("Family\tCluster_id\tSize\tSeq_id\n")
		for i in edict:
			print ">",i
			plist = edict[i]
			oup2 = open("TMP.R","w")
			for j in plist:
				oup2.write("%s\n" % string.joinfields(j,"\t"))
			oup2.close()
			clusters = self.get_relations("TMP.R")
			c = 1
			for j in clusters:
				oup1.write("%s\t%i\t%i\t%s\n" % (i,c,len(j),
											string.joinfields(j,",")))
				c += 1
		oup1.close()
		
		print "Done!"				

	##
	# This is a private method to convert passed file into a format for single
	# linkage. What the single likage function expects is a list with:
	# [[1,2,3],
	#  [3,4],
	#  [4]
	# ]
	#
	# Well, change this scheme. This method converts passed file to a dict. Then
	# passes the dict to get_relations2() and get a cluster list back.
	#
	# @param relations  tab delimited file with a pair of whatever each line
	# @param R1	 index of the "N-term" token, default 0
	# @param R2	 index of the "C-term" token, default 1
	# @param out output file name, if not specified, return dict
	##
	def get_relations(self,relations,R1=0,R2=1,out="",isdict=0):
		
		print "Read relations..."
		if isdict:
			rdict = relations
		else:
			rdict = self.file_to_dict(relations,R1,R2)
		
		if rdict == {}:
			print "Problem: relation dict is empty"
			sys.exit(0)
		
		print "Get clusters..."
		clusters = self.dict_to_list(rdict)
		
		if out != "":
			print "Generate output..."
			oup = open(out,"w")
			oup.write("Cluster_id\tSize\tSeq_id\n")
			c = 1
			for i in clusters:
				oup.write("%i\t%i\t%s\n" % (c,len(i),
											string.joinfields(i,",")))
				c += 1
			print "Done!"
		else:
			return clusters
	
	##
	# Read relations into a dict
	##
	def file_to_dict(self,relations,R1,R2):
	
		inp = open(relations,"r")
		inline = inp.readline()	
		odict = {}  # id as key, index as value
		rdict = {}  # first token as key, second token as value
		c     = 0
		while inline != "":
			inline = self.rnlb(inline)
			llist = inline.split("\t")
			if len(llist)-1 < R2:
				print llist
				print "Insuffiecnt number of tokens, QUIT!"
				sys.exit(0)

			if not odict.has_key(llist[R1]):
				odict[llist[R1]] = c
				c += 1
			if not odict.has_key(llist[R2]):
				odict[llist[R2]] = c
				c += 1		  

			# n-term ipr (ipr1) as key, c-term ipr (ipr2) as value
			if not rdict.has_key(llist[R1]):
				rdict[llist[R1]] = [llist[R2]]
			else:
				if llist[R2] not in rdict[llist[R1]]:
					rdict[llist[R1]].append(llist[R2])
				else:
					print "Redun relation:",llist[R1],llist[R2]
			
			inline = inp.readline()
		return rdict
		
		
	##
	# Get relations from dict, the value should be a list with one or more
	# elements
	##
	def dict_to_list(self,rdict):

		idict = {}  # index as key, original id as value
		odict = {}  # id as key, index as value

		c = 0
		tlist = []
		
		print "1.indexing"
		# assign index to rdict keys first, also contruct a nonredundant list
		for i in rdict.keys():
			idict[c] = i
			odict[i] = c
			for j in rdict[i]:
				if j not in tlist:
					tlist.append(j)
			c += 1
		
		print "2.index to values"
		# assign index to rdict values
		for i in tlist:
			if not i in rdict.keys():
				idict[c] = i
				odict[i] = c
				c += 1

		print "3.generate pre-clusters"
		clusters = []
		for i in rdict.keys():
			clen = len(clusters)
			clusters.append([odict[i]])
			for j in rdict[i]:
				clusters[clen].append(odict[j])		 

		print "4.single linkage"
		clusters = self.single_linkage(clusters)
		#print "Clusters:"
		#for i in clusters:
		#	print "",i
		
		renameC = []
		for i in clusters:
			tlist = []
			for j in i:
				tlist.append(idict[j])
				
			renameC.append(tlist)

		#print "Renamed:"
		#for i in renameC:
		#	print "",i

		return renameC
		
	
	##
	# Single linkage
	##
	def single_linkage(self,cluster_list):

		# should sort the passed list according to their size so there won't
		# be things like a is merged to b, then b need to be merged to c
		sorted_dict = {}
		for i in cluster_list:
			if sorted_dict.has_key(len(i)):
				sorted_dict[len(i)].append(i)
			else:
				sorted_dict[len(i)] = [i]

		sorted_key = sorted_dict.keys()
		sorted_key.sort()
		sorted_key.reverse()
		sorted_clusters = []
		for i in sorted_key:
			alist = sorted_dict[i]
			alist.sort()
			sorted_clusters.extend(alist)

		cluster_list = sorted_clusters

		print "Preclusters:",len(cluster_list)
		for i in range(len(cluster_list)):
			cluster_list[i].sort()
			#print " ",i,cluster_list[i]
		
		merged = []  # list of merged cluster index

		# the list is sorted according to the size, so len(i) >= len(j)
		print "Go through preclusters:"
		for i in range(len(cluster_list)):
			print " ",i
			if i in merged:			    # see if i is merged
				#print "   skipped1",i
				continue

			# iterate cluster member and checked against other clusters
			j = 0
			length_i = len(cluster_list[i])
			while j < length_i:			# each cluster element
				#print " iterate element:",j			    
				linked = 0
				
				for k in range(i+1,len(cluster_list)): # iterate clusters
					#print "  iterate cluster:",k
					if k in merged:		    # see if k is merged
						#print "   skipped2",k
						continue
					
					if cluster_list[i][j] in cluster_list[k]:
						#print "   found links between clusters",i,",",k,"by",\
						#	 cluster_list[i][j]
						if k not in merged:
							merged.append(k)
						for m in cluster_list[k]:      # for addition
							if m not in cluster_list[i]:
								cluster_list[i].append(m)
						linked = 1
					
				if linked:			     # update length
					length_i = len(cluster_list[i])    
					
				j += 1
				
		#print "merged:",merged
		clusters = []
		for i in range(len(cluster_list)):
			if i not in merged:
				clusters.append(cluster_list[i])
				
		return clusters
	###
	# Sort the cluster entries so the ones that are alphabetical will be in
	# front. Also, generate a list of genes to be deleted
	###
	def sort(self,clusters):
		
		inp = open(clusters)
		oup = open(clusters+".sorted","w")
		oup2= open(clusters+".dlist","w")
		oup.write(inp.readline()) # first line is header
		inl = inp.readlines()
		
		for i in inl:
			i = i.split("\t")
			member = i[2].split(",")
			member.sort()
			member.reverse()
			oup.write("%s\t%s\t%s\n" % (i[0],i[1],member.join(",")))
			for j in member[1:]:
				oup2.write("%s\n" % j)
		
		print "Done!"
	
	def rnlb(self,astr):
		if astr[-2:] == "\r\n":
			astr = astr[:-2]
		elif astr[-1] == "\n":
			astr = astr[:-1]
		return astr

	def help(self):
		print "Usage: SingleLinkage.py "
		print "   -f  single - single linkage, NEED: R, OPTIONAL: r1,r2,out"
		print "       group  - single linkage based on groups, NEED: R, G"
		print "       sort   - Sort the cluster entries so the ones that are"
		print "                alphabetical will be in front. Also generate a"
		print "                gene list for deletion. Need: c"
		print "   -R  relations [entry_id1][entry_id2]"
		print "   -G  [group_id][entry_id]"
		print "   -r1 N-terminal token"
		print "   -r2 C-terminal token"
		print "   -o  output file name" 
		print "   -c  cluster file generated by single"                
		sys.exit(0)
		

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	f = R = G = out = c = ""
	r1 = 0
	r2 = 1
	
	single = single_linkage()
	
	for i in range(1,len(sys.argv),2):
		if sys.argv[i] == "-R":
			R  = sys.argv[i+1]
		elif sys.argv[i] == "-f":
			f	 = sys.argv[i+1]
		elif sys.argv[i] == "-G":
			G	 = sys.argv[i+1]
		elif sys.argv[i] == "-r1":
			r1	 = int(sys.argv[i+1])
		elif sys.argv[i] == "-r2":
			r2	 = int(sys.argv[i+1])
		elif sys.argv[i] == "-o":
			out  = sys.argv[i+1]
		elif sys.argv[i] == "-c":
			c    = sys.argv[i+1]
		else:
			print "Unknown option",sys.argv[i]
			single.help()

	if f == "single":
		if R == "":
			print "\nNeed to specify relations\n"
			single.help()
		single.get_relations(R,r1,r2,out)
	elif f == "group":
		if "" in [R,G]:
			print "\nNeed to specify relations\n"
			single.help()
		single.get_group_R(G,R,r1,r2)
	elif f == "sort":
		if c == "":
			print "\nNeed to specify cluster file\n"
			single.help()
		single.sort(c)
	else:
		single.help()

		
