#!/usr/bin/

##
# This module is for drawing blocks in a matrix kind of format. Written to
# display graphics of expression pattern differences.
#
# 01/10,02
#  For draw_bar function, allow negative values and draw the bar in the other
#  orientation.
##

import SVGWriter,FileUtility,sys

class draw_blocks:
	
	def __init__(self):
                pass
		#print "\n*****************"
		#print "*  DRAW BLOCKS  *"
		#print "*****************\n"
	
	##
	# @param gdom   a list of two strings with the standard gdom annotation
	# @param rdict  nt coord as key and Ka divided by Ks as values 
	# @param writer svg_writer instance with dimension set
	# @param inc    increment
	##
	def combine(self,gdom,rdict,writer,inc,step):

		Xcompress = 2
		YCompress = 2
		
		self.draw_domain_local(gdom,writer,inc)
		self.draw_matrix_local(rdict,writer,inc,step)

	#
	# This is to draw a passed matrix file and draw the values in a histogram
	# like fashion iteratively.
	#
	# @param matrix   a file with [name][value1][value2]....
	# @param step     the step size if the x-axis is not incremented by 1
	# @param dom      mark a particular domain, [seq_name][L][R]. Also provide
	#                 order information
	#
	def draw_matrix(self,matrix,step,dom):
		
		mdict  = futil.file_to_dict(matrix)
		dlist  = futil.file_to_list(dom,1,"\t")
		writer = SVGWriter.svg_writer(matrix+".svg")
		
		w = 500
		h = 200+len(mdict.keys())*300
		writer.write_header(w,h)
		inc = 0
		for i in dlist:
			
			rdict = {}
			rlist = mdict[i[0]].split("\t")
			c = 0
			for j in rlist:
				rdict[c] = float(j)*4 # assume it is passed as a value 0-1
				                      # and draw_matrix_local is drawing
				                      # from 0 to 4.
				c += step
			# now pass rdict to be drawn
			self.draw_matrix_local2(rdict,writer,inc,step,i[0],
									int(i[1]),int(i[2]))
			inc += 1
		writer.write_footer()
		
		print "Done!"
	

	def draw_matrix_local2(self,rdict,writer,inc,step,name,L,R):
		
		# set parameters
		Xoff   = 110          # x offset from the corner
		Yunit  = 300
		Yoff   = 50+Yunit*inc # y offset from the corner
		
		Yinc   = 20
		Yscale = 40
		Xscale = 3
		blockW = Xscale*3/5
		base   = Yoff+Yunit
		
		rkeys = rdict.keys()
		rkeys.sort()
		leng = rkeys[-1]
		
		writer.set_font_size(10)
		writer.draw_string(name,5,base-Yscale*4)
		
		# draw scale
		writer.set_fill_c("black")
		writer.draw_line(Xoff-5,base         ,Xoff+leng*Xscale+20,base         )
		writer.draw_line(Xoff-5,base         ,Xoff-5             ,base-Yscale*4)
		writer.draw_line(Xoff-7,base         ,Xoff-5             ,base         )
		writer.draw_line(Xoff-7,base         ,Xoff-5             ,base         )
		writer.draw_line(Xoff-7,base-Yscale*1,Xoff+leng*Xscale+20,base-Yscale*1)
		writer.draw_line(Xoff-7,base-Yscale*2,Xoff+leng*Xscale+20,base-Yscale*2)
		writer.draw_line(Xoff-7,base-Yscale*3,Xoff+leng*Xscale+20,base-Yscale*3)
		writer.draw_line(Xoff-7,base-Yscale*4,Xoff+leng*Xscale+20,base-Yscale*4)
		
		# draw the location of the passed coords
		writer.draw_rect(Xoff+L+blockW/2,base-Yscale*4.5,R-L,blockW)
		
		# now draw the histogram part
		c = 1
		for i in rkeys:
			x = i*Xscale+1

			# alternative x axis display
			if c%10 == 1:
				writer.draw_string("%i" % (int(i)+1),Xoff+x-2,base+15)
				writer.draw_line(Xoff+x+blockW/2,base,Xoff+x+blockW/2,base+2)
			c+= 1
			writer.draw_rect(Xoff+x,base-float(rdict[i])*Yscale,
									 blockW,float(rdict[i])*Yscale)
									 
	
	
	# This is for comibining both graphic domain and histograms
	def draw_matrix_local(self,rdict,writer,inc,step,name="",protein=0):
		
		# set parameters
		Xoff   = 110          # x offset from the corner
		Yunit  = 300
		Yoff   = 50+Yunit*inc # y offset from the corner
		
		blockW = 6*step/30   # for if domains are included
		#blockW = 6
		Yinc   = 20           # y increment: set 50 if domains are included	
		Yscale = 40           # y scaling factor, set 40 if domains are included
		
		base   = Yoff+Yunit
		
		rkeys = rdict.keys()
		rkeys.sort()
		if protein:
			leng = rkeys[-1]
		else:
			leng = rkeys[-1]/3
		
		if name != "":
			writer.set_font_size(10)
			writer.draw_string(name,5,base-Yscale*4)
		
		# draw scale
		writer.set_fill_c("black")
		writer.draw_line(Xoff-5,base         ,Xoff+leng+20,base         )
		writer.draw_line(Xoff-5,base         ,Xoff-5      ,base-Yscale*4)
		writer.draw_line(Xoff-7,base         ,Xoff-5      ,base         )
		writer.draw_line(Xoff-7,base         ,Xoff-5      ,base         )
		writer.draw_line(Xoff-7,base-Yscale*1,Xoff+leng+20,base-Yscale*1)
		writer.draw_line(Xoff-7,base-Yscale*2,Xoff-5      ,base-Yscale*2)
		writer.draw_line(Xoff-7,base-Yscale*3,Xoff-5      ,base-Yscale*3)
		writer.draw_line(Xoff-7,base-Yscale*4,Xoff-5      ,base-Yscale*4)
		
		writer.set_font_size(6)
		c = 1
		for i in rkeys:
			if protein:
				x = i+1
			else:
				x = i/3+1
			# alternative x axis display
			if c%2 == 1:
				writer.draw_string("%i" % (int(i)+1),Xoff+x,base+10)
				writer.draw_line(Xoff+x+blockW/2,base,Xoff+x+blockW/2,base+2)
			c+= 1
			try:
				if float(rdict[i]) > 4:
					writer.draw_rect(Xoff+x,base-4*Yscale,
									 blockW,4*Yscale)
					writer.draw_string(rdict[i][:3],Xoff+x,base-4*Yscale-10)
				else:
					writer.draw_rect(Xoff+x,base-float(rdict[i])*Yscale,
									 blockW,float(rdict[i])*Yscale)
				#writer.draw_line(Xoff+x+blockW/2,base,
				#				 Xoff+x+blockW/2,base+2)
				#if float(rdict[i]) > 1:
				#	writer.draw_line(Xoff+x+blockW/2,base,
				#				 Xoff+x+blockW/2,base+2)

			except ValueError:
				writer.draw_string(rdict[i],Xoff+x,base-5)


	def draw_domain_local(self,gdom,writer,inc):

		# set parameters
		Xoff   = 10           # x offset from the corner
		Yunit  = 300
		Yoff   = 70+Yunit*inc # y offset from the corner
		idW    = 100          # id width
		blockH = 20
		Yinc   = 50           # y increment	
		
		# list with [[seq_id,{dom1:[L,R],dom2...}],[seq_id2,{...}],...]
		# use list so the order is preserved for seqeunces
		glist = []
		d = 0
		for i in gdom:
			# skip description
			if i[0] == "#":
				continue			
			# rid of new line and empty tabs
			dlist = i.split("\t")
			c = 0
			for j in dlist:
				if j == "":
					break
				c += 1
			dlist = dlist[:c]	
			idx   = dlist[0]
			try:		
				size  = float(dlist[-1])    # scaled size
			except ValueError:
				writer.set_stroke_c("black")
				writer.set_stroke_w(0)
				writer.set_fill_c("black")
				writer.set_font_size(10)
				writer.draw_string(idx,Xoff,Yoff+(d+0.45)*Yinc)
				writer.set_font_size(8)	
				writer.draw_string(dlist[-1],idW,Yoff+(d+0.4)*Yinc)
				d += 1
				continue

			dlist = dlist[1:-1]
			
			# write seq id
			writer.set_stroke_c("black")
			writer.set_stroke_w(0)
			writer.set_fill_c("black")
			writer.set_font_size(10)
			writer.draw_string(idx,Xoff,Yoff+(d+0.45)*Yinc)
			
			# read domain into a dict, dom order as key, [mname,L,W] as value
			# notice that the coordinates are scaled in the dict
			ddict = {}
			c = 0
			for j in range(0,len(dlist),2):
				L = float(dlist[j])
				W = float(dlist[j+1][dlist[j+1].find("|")+1:])
				m = dlist[j+1][:dlist[j+1].find("|")]
				ddict[c] = [m,L,W]
				c += 1
			
			# convert coordinates based on ref domain or predefined offset
			for j in ddict.keys():
				ddict[j] = [ddict[j][0],ddict[j][1],ddict[j][2]]
				
			# draw full length backbone
			writer.set_stroke_c("black")
			writer.set_stroke_w(0.5)
			writer.set_fill_c("white")
			writer.draw_rect(idW,Yoff+(d+0.5)*Yinc,size,blockH)
			
			writer.set_font_size(8)			
			# draw domains
			for j in ddict.keys():
				
				x     = ddict[j][1]
				w     = ddict[j][2]
				
				
				if ddict[j][0] in ["TRANS","SIGNAL"]:
					writer.set_fill_c("black")
				elif ddict[j][0] in ["STYKc","TyrKc","S_TKc","pkinase"]:
					writer.set_fill_c("darkblue")
				else:
					writer.set_fill_c("beige")
				writer.draw_rect(x+idW,Yoff+(d+0.5)*Yinc,w,blockH)		
				writer.draw_string(ddict[j][0],x+idW,Yoff+(d+0.4)*Yinc)
			
			# increment seq count
			d += 1
					
	
	##
	# draw histogram-like stuff
	#
	##
	def draw_bar(self,matrix,Xzoom=100):
	
		svg_write = SVGWriter.svg_writer(matrix+".svg")
		
		# read matrix into dict
		#mdict = futil.file_to_dict(matrix,3)
		# read matrix into list, use as order
		mlist = futil.file_to_list(matrix,1,"\t")
		
		# set constants
		Xoff   = 50  # x offset
		Yoff   = 50  # y offset
		idW    = 200 # id width
		blockH = 25
		Xinc   = 300 # x increment
		Yscale = 200
		Yinc   = 30  # y increment

		W = Xoff + idW +(len(mlist[0])-1)*Xinc # identifier, each column
		H = Yoff + len(mlist)*Yinc
		svg_write.write_header(W,H)
		
		# draw identifer and rectangles
		d = 1
		for i in mlist:
			#print "",i[0]
			svg_write.draw_string(i[0], Xoff, Yoff + (d+0.5)*Yinc)	
			c = 0
			svg_write.set_stroke_c("black")
			svg_write.set_fill_c("black")
			for j in i[1:]:
				try:
					if float(j) == 0:
						pass
					# if the value is negative
					elif float(j) < 0:
						# notice that it is negative
						# if it's too negative
						if -float(j)*Xzoom > Yscale:
							svg_write.draw_rect(Xoff + idW + c*Xinc - Yscale, 
											Yoff + d*Yinc,
											Yscale, blockH, 1)
							svg_write.draw_string(j, 
											Xoff + idW + c*Xinc - Yscale - 50 , 
											Yoff + (d+0.5)*Yinc)								
						else:
							svg_write.draw_rect(Xoff + idW + c*Xinc + \
											float(j)*Xzoom, 
											Yoff + d*Yinc,
											-float(j)*Xzoom, blockH, 1)	
					else:
						if float(j)*Xzoom > Yscale:
							svg_write.draw_rect(Xoff + idW + c*Xinc, 
											Yoff + d*Yinc,
											Yscale, blockH, 1)
							svg_write.draw_string(j, 
											Xoff + idW + c*Xinc + Yscale + 10 , 
											Yoff + (d+0.5)*Yinc)								
						else:
							svg_write.draw_rect(Xoff + idW + c*Xinc, 
												Yoff + d*Yinc,
												float(j)*Xzoom, blockH, 1)
				except ValueError:
					svg_write.draw_string(j, Xoff + idW + c*Xinc, 
											 Yoff + (d+0.5)*Yinc)	
				c += 1
			
			svg_write.set_stroke_c("")
			d += 1
		
		# draw scale
		svg_write.set_stroke_c("black")
		svg_write.draw_rect(Xoff+idW, Yoff+(d+0.5)*Yinc, Yscale, blockH/5, 1)
		for i in range(5):
			svg_write.draw_line(Xoff+idW+Yscale*i/4, Yoff+(d+0.5)*Yinc -5,
								Xoff+idW+Yscale*i/4, Yoff+(d+0.5)*Yinc +5) 									
		svg_write.draw_string("%i Units (pixel:%i, xzoom:%f)" % \
							 (Yscale/Xzoom, Yscale, Xzoom), 
							  Xoff+idW, Yoff + (d+1.5)*Yinc)		
		svg_write.write_footer()
		
		print "Done!"

	
	##
	# Check the unique entries in a matrix for coloring purpose
	#
	# @param matrix
	##
	def list_unique(self,matrix):
		
		inp = open(matrix,"r")
		inl = inp.readline()
		udict = {}
		while inl != "":
			L = self.onl(inl).split("\t")
			for j in L:
				if udict.has_key(j):
					udict[j] += 1
				else:
					udict[j] = 1
			inl = inp.readline()
		oup = open(matrix+".unique","w")
		for i in udict:
			oup.write("%s\t%i\n" % (i,udict[i]))
		
		print "Done!"
		
	##
	# @param matrix  the text file specifying the matrix. First line is column
	#                header, will be used as header in SVG output. The first
	#                token in the line is the identifier.
	# @param scheme  the color scheme in the format [key1]:[color1],...
	# @param shape   the shape of objects, default rect.
	# @param cat     if the scheme file passes is categorical [0, default], or
	#                coninuous [1]. If 1, the values in the matrix will be
	#                matched to he bin range defined in the shceme file. E.g.
	#                 1 red
	#                 5 black
	#                For a value of 4, it will be black. NOT IMPLEMENTED YET!!!
	##
	def color_matrix(self,matrix,scheme,shape="rect",desc=0,cat=0):
		
		svg_write = SVGWriter.svg_writer(matrix+".svg")
		
		# if scheme don't have ":", assume it is a file
		if scheme.find(":") == -1:
			sdict = self.read_scheme(scheme)
			print "here"
			# notice that the sdict passed from read_scheme is different what is
			# anticipated, so change its format
			tdict = {}
			for i in sdict:
				#print i,sdict[i]
				tdict[i] = sdict[i][1]
			sdict = tdict
		else:
			# parse scheme into dict
			slist = scheme.split(",")
			sdict = {}
			for i in slist:
				i = i.split(":")
				sdict[i[0]] = i[1]
		
		# read matrix into dict
		#mdict = futil.file_to_dict(matrix,3)
		# read matrix into list, use as order
		mlist = futil.file_to_list(matrix,1,"\t")
		
		# set constants
		Xoff   = 50  # x offset
		Yoff   = 50  # y offset
		idW    = 100 # id width
		blockW = 25
		blockH = 25
		Xinc   = 30  # x increment
		Yinc   = 30  # y increment 

		W = Xoff + idW +(len(mlist[0])-1)*Xinc # identifier, each column
		H = Yoff + len(mlist)*Yinc
		svg_write.write_header(W,H)
		
		# if 2nd element in first line is not in the scheme, assume it is desc
		
		if desc:
			print "Assume first line is description"
			c = 0
			svg_write.draw_string(mlist[0][0],Xoff,Yoff)		
			for i in mlist[0][1:]:
				svg_write.draw_string(i, Xoff + idW + c*Xinc, Yoff)
				c += 1
			mlist = mlist[1:]
		else:
			print "Assume no description line"
		
		# based on the range the value is 
		def get_color(value):
			if sdict.has_key(value):
				return sdict[value]
			else:
				try:
					value = float(value)
					tdict = {}
					for i in sdict:
						tdict[float(i)] = sdict[i]
					tkeys = tdict.keys()
					tkeys.sort()
					color = ""
					for i in tkeys:
						if value == 0:
							break
						elif value <= i:
							color = tdict[i]
							break
					return color
				except ValueError:
					return ""
		
		# draw identifer and rectangles
		print "Draw shape:",shape
		d = 1
		not_specified = []
		for i in mlist:
			#print "",i[0]
			svg_write.draw_string(i[0], Xoff, Yoff + (d+0.5)*Yinc)	
			c = 0
			svg_write.set_stroke_c("black")
			for j in i[1:]:
				color = get_color(j)
				#print [j,color]			
				if color != "":
					svg_write.set_fill_c(color)
					if shape == "circle":
						svg_write.draw_circle(Xoff + idW + c*Xinc,
											  Yoff + (d+0.5)*Yinc, blockH/2)			
					elif shape == "rect":
						svg_write.draw_rect(Xoff + idW + c*Xinc, Yoff + d*Yinc,
											blockW, blockH, 1)
				else:
					if j not in not_specified:
						not_specified.append(j)
					#svg_write.draw_string(j,Xoff+idW+c*Xinc,Yoff+(d+0.5)*Yinc)	
					
				svg_write.set_fill_c("")
				c += 1
			
			svg_write.set_stroke_c("")
			d += 1
			
		svg_write.write_footer()
		
		print "The following keys are not specified in the scheme: between ''"
		for i in not_specified:
			print " '%s'"  % i
		print "Done!"
	
	##
	# This function will draw domains...
	##
	def draw_domain(self,gdom,scheme,R,X,Y,H,A,Yinc):
		
		# set parameters
		Xoff   = 10       # x offset from the corner
		Yoff   = 50       # y offset from the corner
		Xsize  = 2000*X   # x dimension of SVG output 
		Xref   = 600*X    # x offset for reference domain
		Xno_ref= 100*X    # x offset for those without ref or not specified
		idW    = 80*X     # id width
		blockH = H*Y
		Yinc   = Yinc*Y   # y increment
		APR    = A/X      # AA to pixel ratio 	
		
		sdict = self.read_scheme(scheme)
				
		# read gdom file and reset coords
		inp = open(gdom,"r")
		inlines = inp.readlines()

		writer = SVGWriter.svg_writer(gdom+".svg")		
		writer.write_header(Xsize,Yoff+len(inlines)*(Yinc+2))
		
		# list with [[seq_id,{dom1:[L,R],dom2...}],[seq_id2,{...}],...]
		# use list so the order is preserved for seqeunces
		glist = []
		d = 0
		missedDomain = []
		for i in inlines:
			#print [i]
			# skip description
			if i[0] == "#":
				continue
			
			# rid of new line and empty tabs
			dlist = i[:-1].split("\t")
			c = 0
			for j in dlist:
				if j == "":
					break
				c += 1
			dlist = dlist[:c]	
			idx   = dlist[0]
			
			size  = float(dlist[-1])/APR    # scaled size
			dlist = dlist[1:-1]
			# write seq id
			writer.set_stroke_c("black")
			writer.set_stroke_w(0)
			writer.set_fill_c("black")
			writer.set_font_size(12*Y)
			writer.draw_string(idx,Xoff,Yoff+(d+0.5)*Yinc)
			
			# read domain into a dict, dom order as key, [mname,L,W] as value
			# notice that the coordinates are scaled in the dict
			ddict = {}
			c = 0
			refL = ""
			for j in range(0,len(dlist),2):
				L = float(dlist[j])/APR
				W = float(dlist[j+1][dlist[j+1].find("|")+1:])/APR
				m = dlist[j+1][:dlist[j+1].find("|")]
				ddict[c] = [m,L,W]
				# only adjust according to the first reference domain
				if m in R.split(",") and refL == "":
					refL = L
				c += 1
			
			# convert coordinates based on ref domain or predefined offset
			no_ref = 0
			if R != "" and refL != "":
				shift = Xref - refL
			else:
				shift = Xno_ref
			for j in ddict.keys():
				ddict[j] = [ddict[j][0],ddict[j][1]+shift,ddict[j][2]]
				
			# draw full length backbone
			writer.set_stroke_c("black")
			writer.set_stroke_w(0.5)
			writer.set_fill_c("beige")
			writer.draw_rect(idW+shift,Yoff+(d+0.5)*Yinc+blockH/8.0*3.0,
							 size,blockH/4.0)
			
			writer.set_font_size(6*Y)
			
			# draw domains
			for j in ddict.keys():
				try:
					shape = sdict[ddict[j][0]][0]
					color = sdict[ddict[j][0]][1]
				except KeyError:
					# won't be drawn if not in the scheme
					if ddict[j][0] not in missedDomain:
						missedDomain.append(ddict[j][0])
					continue

				x     = ddict[j][1]
				w     = ddict[j][2]
				
				writer.set_fill_c(color)
				
				# deal with signal and trans right here:
				if ddict[j][0] == "SIGNAL" or ddict[j][0] == "TRANS":
					writer.draw_rect(x+idW,Yoff+(d+0.5)*Yinc+blockH/8.0*2.0
									,w,blockH/2)
				elif shape == "3":
					writer.draw_polygon(x+idW,Yoff+(d+0.5)*Yinc,w,blockH,3)
				elif shape == "4":
					writer.draw_rect(x+idW,Yoff+(d+0.5)*Yinc,w,blockH)
				elif shape == "5":
					writer.draw_polygon(x+idW,Yoff+(d+0.5)*Yinc,w,blockH,5)
				elif shape == "6":
					writer.draw_polygon(x+idW,Yoff+(d+0.5)*Yinc,w,blockH,6)
				elif shape in ["R","r","C","c"]:
					writer.draw_rounded_corners(x+idW,Yoff+(d+0.5)*Yinc,w,blockH)
					
				writer.draw_string(ddict[j][0],x+idW,Yoff+(d+0.5)*Yinc)
			
			# increment seq count
			d += 1
			
		# draw scale bar
		writer.set_fill_c("black")
		writer.set_stroke_w(2)
		writer.set_font_size(12*Y)
		writer.draw_rect(idW,Yoff+(d+0.6)*Yinc,200.0/APR,blockH/4.0)
		writer.draw_string("200 aa",idW,Yoff+(d+0.5)*Yinc)
		
		writer.write_footer()
		
		if missedDomain != []:
			print "Domains not defined in the scheme:"
			for i in missedDomain:
				print "",i
		
		print "Done!"


	##
	# This method is written to covert bootstrap values into colored circles.
	# For values not specified in the scheme file, they will not be drawn.
	# At this point, the categories are hard-wired into the function, at 10
	# per interval. Also, dimension is set at the letter page size.
	##
	def draw_obj(self,coords,scheme,size):
		
		# read scheme into dict
		sdict = self.read_scheme(scheme)
		
		writer = SVGWriter.svg_writer(coords+".svg")		
		writer.write_header(202,776)
		
		# read coords and draw obj
		inp = open(coords,"r")
		inline = inp.readline()
		
		while inline != "":
			llist = inline[:-1].split("\t")
			
			# skip those non number...
			obj = ""
			try:
				obj = str(int(llist[0])/10*10)
			except ValueError:
				inline = inp.readline()
				continue
		
			# notice that this round down the value to the nearest base 10.
			x     = float(llist[1])
			y     = float(llist[2])
			shape = sdict[obj][0]
			color = sdict[obj][1]
						
			writer.set_fill_c(color)
			
			if shape == "3":
				writer.draw_polygon(x,y,size,size,3)
			elif shape == "4":
				writer.draw_rect(x,y,size,size)
			elif shape == "5":
				writer.draw_polygon(x,y,size,size,5)
			elif shape == "6":
				writer.draw_polygon(x,y,size,size,6)
			elif shape == "R":
				writer.draw_rounded_corners(x,y,size)
			elif shape == "C":
				writer.draw_circle(x,y,size)
			
			inline = inp.readline()
		
		writer.write_footer()		
		print "Done!"
			
			
		
	##
	# A private method converting the scheme file passed into a dict. Always
	# take the last two tokens of tab-delimited lines.
	##	
	def read_scheme(self,scheme):

		# Read scheme [dom_id][count][shape][RGB] into dict
		# motif as key, [shape,color] as value
		# Shapes are: 3, 4, 5, 6, and R (rounded corner)
		sdict = {}
		inp = open(scheme,"r")		
		inline = inp.readline()
		while inline != "":
			llist = inline[:-1].split("\t")
			sdict[llist[0]] = llist[-2:]
			inline = inp.readline()
		
		print sdict
		return sdict
	
	#
	# assume the dat_table is the output of DatParser.parse_coord
	#
	# @param dat_table "CHR\tCODE\tTYPE\tORI\tL\tR\n"
	# @param names     ??
	#
	def draw_chr(self,dat_table,names):
		
		ndict = futil.file_to_dict(names,0)
			
		inp = open(dat_table,"r")
		inp.readline()
		inline = inp.readline()
		chromo = 0
		last   = ""
		
		# a nested dict with chr as key, value is a dict with one coord as key,
		# code as value.
		cdict  = {}
		Ylen   = 0
		
		Xoff   = float(100) # x offset
		Yoff   = float(50) # y offset
		W      = 2
		Yinc   = float(10) # y increment
		Yzoom  = 1e5
		textH  = 10
		countOb= 0
		
		while inline != "":
			inline = self.onl(inline)
			llist  = inline.split("\t")
			if chromo == 0:
				chromo = int(llist[0])
				cdict[chromo] = {float(llist[-2]):"N"}
			elif chromo != int(llist[0]):
				cdict[chromo][float(last[-1])] = "S"
				Ylen = Ylen + float(last[-1])/Yzoom
				countOb += 1

				chromo = int(llist[0])
				cdict[chromo] = {float(llist[-2]):"N"}

			if ndict.has_key(llist[1]):
				cdict[chromo][float(llist[-2])] = llist[1]
				
			last   = llist
			inline = inp.readline()
		
		# get the last coord for the last chr
		cdict[chromo][float(last[-1])] = "S"
		Ylen += float(last[-1])/Yzoom
		
		# draw locations		
		writer = SVGWriter.svg_writer(names+".svg")
		writer.write_header(Xoff*2,Ylen+countOb*Yinc+Yoff*2)
		
		ckeys = cdict.keys()
		ckeys.sort()
		countChr = 0
		lastS = 0
		lastT = 0
		for i in ckeys:
			
			ekeys = cdict[i].keys()
			ekeys.sort()

			# draw chr
			chrS  = ekeys[-1]/Yzoom
			writer.draw_rect(Xoff,Yoff+lastS,W,chrS)
			
			# draw features			
			for j in ekeys[1:-1]:
				coord = j/Yzoom
				writer.draw_line(Xoff-3,Yoff+lastS+coord,
								 Xoff+3,Yoff+lastS+coord)
				#print ">",Yoff+lastS+coord,lastT
				if lastT == 0 or Yoff+lastS+coord-lastT > textH:
					lastT = Yoff+lastS+coord
					writer.draw_line(Xoff-20,Yoff+lastS+coord,
									 Xoff-3,Yoff+lastS+coord)
				else:
					lastT = lastT+textH
					writer.draw_line(Xoff-20,lastT,
									 Xoff-3,Yoff+lastS+coord)
				writer.draw_string(cdict[i][j],0,lastT)
			
			lastS = lastS + chrS + Yinc
			countChr += 1
		
		writer.write_footer()
	
		print "Done!"

	#
	#
	# @param dat_table [feature][type][ori][L][R]. 
	# @param cname     chr name
	# @param cleng     chr length
	# @param S         starting position
	# @param Xzoom     x compression factor
	# @param text      draw string or not
	#
	def draw_chr3(self,dat_table,cname,cleng,S=0,Xzoom=10,text=0):
		
		inp = open(dat_table,"r")
		inline = inp.readline()
		last   = ""
		
		# a dict with L coord as key, a nested list as value
		cdict = {}
		
		# dimension preset
		Xoff   = float(100)  # x offset
		Yoff   = float(50)   # y offset
		W      = 2           # width
		Yinc   = float(10)   # y increment
		Ylen   = cleng/Xzoom # Y length 
		textH  = 10          # text height
		
		print "Read feature table..."
		g_vs_o = {}
		while inline != "":
			L = self.onl(inline).split("\t")
			try:
				float(L[-2])
				if cdict.has_key(L[-2]):
					cdict[L[-2]].append(L[:-2] + [L[-1]])
				else:
					cdict[L[-2]] = [L[:-2] + [L[-1]]]
			except ValueError:
				# title line?
				print "Format:",L
			inline = inp.readline()	
		
		# sort based on L coord
		ckeys = cdict.keys()
		ckeys.sort()
		
		# if cleng is not specified, use the last element R coord as far right
		if cleng == 0:
			cleng = int(cdict[ckeys[-1]][-1][-1]) - S
		
		# initiate svg
		writer = SVGWriter.svg_writer(dat_table+".svg")
		writer.write_header(Xoff*2+cleng/Xzoom,Ylen+Yoff*2+W*100)
		
		# draw_chr	
		writer.set_fill_c("grey")
		writer.draw_rect(Xoff,Yoff-W*0.5,cleng/Xzoom,W)
		writer.draw_string(cname,Xoff,Yoff-W*4)
		
		# legend
		writer.set_fill_c("black")
		writer.draw_rect(Xoff,Yoff+W*8,2000/Xzoom,W)
		writer.draw_string("2kb",Xoff,Yoff+W*15)
		
		# posterior probability lines
		writer.set_fill_c("")
		writer.draw_line(Xoff,Yoff+W*40 ,Xoff+cleng/Xzoom,Yoff+W*40)
		writer.draw_line(Xoff,Yoff+W*60 ,Xoff+cleng/Xzoom,Yoff+W*60)
		writer.draw_line(Xoff,Yoff+W*80 ,Xoff+cleng/Xzoom,Yoff+W*80)
		writer.draw_line(Xoff,Yoff+W*100,Xoff+cleng/Xzoom,Yoff+W*100)
		writer.draw_line(Xoff,Yoff+W*120,Xoff+cleng/Xzoom,Yoff+W*120)
		writer.draw_line(Xoff,Yoff+W*140,Xoff+cleng/Xzoom,Yoff+W*140)
		
		vC = 16
		writer.draw_line(Xoff,Yoff+W*40 ,Xoff,Yoff+W*(40-vC))
		writer.draw_line(Xoff,Yoff+W*60 ,Xoff,Yoff+W*(60-vC))
		writer.draw_line(Xoff,Yoff+W*80 ,Xoff,Yoff+W*(80-vC))
		writer.draw_line(Xoff,Yoff+W*100,Xoff,Yoff+W*(100-vC))
		writer.draw_line(Xoff,Yoff+W*120,Xoff,Yoff+W*(120-vC))
		writer.draw_line(Xoff,Yoff+W*140,Xoff,Yoff+W*(140-vC))
		
		T = 0.3285 # threshold setting 2.5% false positive rate.
		writer.set_stroke_c("red")
		writer.draw_line(Xoff,Yoff+W*40-W*T*vC ,Xoff+cleng/Xzoom,Yoff+W*40-W*T*vC)
		writer.draw_line(Xoff,Yoff+W*60-W*T*vC ,Xoff+cleng/Xzoom,Yoff+W*60-W*T*vC)
		writer.draw_line(Xoff,Yoff+W*80-W*T*vC ,Xoff+cleng/Xzoom,Yoff+W*80-W*T*vC)
		writer.draw_line(Xoff,Yoff+W*100-W*T*vC,Xoff+cleng/Xzoom,Yoff+W*100-W*T*vC)
		writer.draw_line(Xoff,Yoff+W*120-W*T*vC,Xoff+cleng/Xzoom,Yoff+W*120-W*T*vC)
		writer.draw_line(Xoff,Yoff+W*140-W*T*vC,Xoff+cleng/Xzoom,Yoff+W*140-W*T*vC)
		
		
		print "Total %i features, draw them..." % len(ckeys)
		for i in ckeys:
			cL = (float(i)-S)/Xzoom
			for j in cdict[i]:
				cR = (float(j[-1])-S)/Xzoom				
				writer.set_stroke_c("black")			
				if j[1] == "GENE":
					writer.set_fill_c("white")
					writer.draw_rect(Xoff+cL,Yoff-W*2,cR-cL,W*4)
					
					# draw orientation
					writer.set_fill_c("")
					if j[2] in ["+","F"]:
						writer.draw_polyline("%f,%f,%f,%f,%f,%f,%f,%f" % \
							(Xoff+cL      ,Yoff-W*2,
							 Xoff+cL      ,Yoff-W*3,
							 Xoff+cL+W*2  ,Yoff-W*3,
							 Xoff+cL+W*1.5,Yoff-W*3.5))
					elif j[2] in ["-","R"]:
						writer.draw_polyline("%f,%f,%f,%f,%f,%f,%f,%f" % \
							(Xoff+cR      ,Yoff-W*2,
							 Xoff+cR      ,Yoff-W*3,
							 Xoff+cR-W*2  ,Yoff-W*3,
							 Xoff+cR-W*1.5,Yoff-W*3.5))
						
						
					"""
					# draw orientation triangles
					writer.set_fill_c("black")
					if j[2] in ["+","F"]:
						writer.draw_triangle(Xoff+cR,Yoff-W*2,
											 Xoff+cR,Yoff,
											 Xoff+cL+(cR-cL)*29.0/30.0,Yoff-W*2)
						writer.draw_triangle(Xoff+cR,Yoff+W*2,
											 Xoff+cR,Yoff,
											 Xoff+cL+(cR-cL)*29.0/30.0,Yoff+W*2)
					elif j[2] in ["-","R"]:
						writer.draw_triangle(Xoff+cL,Yoff-W*2,
											 Xoff+cL,Yoff,
											 Xoff+cL+(cR-cL)/30.0,Yoff-W*2)
						writer.draw_triangle(Xoff+cL,Yoff+W*2,
											 Xoff+cL,Yoff,
											 Xoff+cL+(cR-cL)/30.0,Yoff+W*2)
					"""
					
					if text:
						writer.draw_string(j[0],Xoff+cL,Yoff+W*8)
				elif j[1] == "ORF":
					writer.set_fill_c("lightblue")
					writer.draw_rect(Xoff+cL,Yoff-W*2,cR-cL,W*4)
				elif j[1] == "intron":
					writer.set_fill_c("white")
					writer.set_stroke_c("white")
					writer.draw_rect(Xoff+cL,Yoff-W*2,cR-cL,W*4)
					writer.set_fill_c("black")
					writer.set_stroke_c("black")
					# vertical 1
					writer.draw_rect(Xoff+cL,Yoff-W*2,0.1,W*4)
					# vertocal 2
					writer.draw_rect(Xoff+cR,Yoff-W*2,0.1,W*4)
					# link 1
					writer.draw_line(Xoff+cL,Yoff,Xoff+(cL+cR)/2,Yoff-W*3)
					# link 2
					writer.draw_line(Xoff+cR,Yoff,Xoff+(cL+cR)/2,Yoff-W*3)
				elif j[1] == "PP": # assume to be posterior probability
					pp    = float(j[2])
					frame = j[0].split("_")[-1]
					writer.set_stroke_c("grey")
					writer.set_fill_c("grey")
					if pp >= T:
						writer.set_stroke_c("lightgreen")
						writer.set_fill_c("lightgreen")
					if frame == "1":
						writer.draw_rect(Xoff+cL,Yoff+W*40-W*pp*vC,cR-cL,W*pp*vC)
					elif frame == "2":
						writer.draw_rect(Xoff+cL,Yoff+W*60-W*pp*vC,cR-cL,W*pp*vC)
					elif frame == "3":
						writer.draw_rect(Xoff+cL,Yoff+W*80-W*pp*vC,cR-cL,W*pp*vC)
					elif frame == "4":
						writer.draw_rect(Xoff+cL,Yoff+W*100-W*pp*vC,cR-cL,W*pp*vC)
					elif frame == "5":
						writer.draw_rect(Xoff+cL,Yoff+W*120-W*pp*vC,cR-cL,W*pp*vC)
					elif frame == "6":
						writer.draw_rect(Xoff+cL,Yoff+W*140-W*pp*vC,cR-cL,W*pp*vC)
								
		writer.write_footer()	
		print "Done!"

	#
	# assume the dat_table is the output of DatParser.parse_coord
	#
	def draw_chr2(self,dat_table,dup):

		print "Read duplication file:",dup	
		inp    = open(dup,"r")
		inline = inp.readline()
		ddict  = {}
		while inline != "":
			inline = self.onl(inline)
			llist  = inline.split("\t")
			ddict[llist[1]] = [llist[0],"N"]
			ddict[llist[2]] = [llist[0],"S"]
			inline = inp.readline()
		
		print "Read dat table, assign coords..."
		inp = open(dat_table,"r")
		inp.readline()
		inline = inp.readline()
		chromo = 0
		last   = ""
		
		# a nested dict with chr as key, value is a dict with code as key, coord
		# as value (only one coord is saved).
		cdict  = {}
		Ylen   = 0
		
		Xoff   = float(100) # x offset
		Yoff   = float(50)  # y offset
		W      = 2
		Yinc   = float(10)  # y increment
		Yzoom  = 1e5
		textH  = 10
		countOb= 0
		triW   = 20		
		
		while inline != "":
			inline = self.onl(inline)
			llist  = inline.split("\t")
			if chromo == 0:
				chromo = int(llist[0])
				cdict[chromo] = {float(llist[-2]):"N"}
			elif chromo != int(llist[0]):
				cdict[chromo][float(last[-1])] = "S"
				Ylen = Ylen + float(last[-1])/Yzoom
				countOb += 1

				chromo = int(llist[0])
				cdict[chromo] = {float(llist[-2]):"N"}

			if ddict.has_key(llist[1]):
				cdict[chromo][float(llist[-2])] = [llist[1],
												   ddict[llist[1]][0],
												   ddict[llist[1]][1]]
				
			last   = llist
			inline = inp.readline()
		
		# get the last coord for the last chr	
		cdict[chromo][float(last[-1])] = "S"
		Ylen += float(last[-1])/Yzoom
				
		# draw locations
		print "Drawing..."
		writer = SVGWriter.svg_writer(dup+".svg")
		writer.write_header(Xoff*4,Ylen+countOb*Yinc+Yoff*2)
		
		ckeys = cdict.keys()
		ckeys.sort()
		countChr = 0
		lastS = 0
		lastT = 0
		bdict = {}
		for i in ckeys:
			
			ekeys = cdict[i].keys()
			ekeys.sort() 

			# draw chrs
			chrS  = ekeys[-1]/Yzoom
			writer.draw_rect(Xoff  ,Yoff+lastS,W,chrS)
			#writer.draw_rect(Xoff*3,Yoff+lastS,W,chrS)
			
			# draw features			
			for j in range(1,len(ekeys)-1,2):
				coordN = ekeys[j]  /Yzoom
				coordS = ekeys[j+1]/Yzoom
				
				writer.draw_rect(Xoff+5+W,Yoff+lastS+coordN,W,coordS-coordN)
				if not bdict.has_key(cdict[i][ekeys[j]][1]):
					bdict[cdict[i][ekeys[j]][1]] = [Xoff+5+W,
						Yoff+lastS+coordN+(coordS-coordN)/2,0,0]
				else:
					bdict[cdict[i][ekeys[j]][1]][2:] = [Xoff+5+W,
						Yoff+lastS+coordN+(coordS-coordN)/2]

				#writer.draw_string(cdict[i][ekeys[j]][1],
				#				   Xoff+5+W,Yoff+lastS+coordN)
			
			lastS = lastS + chrS + Yinc
			countChr += 1
		
		countBlk = 0
		linkW    = 100
		ylist    = []
		yShift   = 40
		bkeys = bdict.keys()
		bkeys.sort()
		for i in bkeys:
			"""
			# trial 1
			writer.draw_line(bdict[i][0]               ,bdict[i][1],
							 bdict[i][0]+linkW*countBlk,bdict[i][1])
			writer.draw_line(bdict[i][2]               ,bdict[i][3],
							 bdict[i][2]+linkW*countBlk,bdict[i][3])
			writer.draw_line(bdict[i][0]+linkW*countBlk,bdict[i][1],
							 bdict[i][2]+linkW*countBlk,bdict[i][3])
			
			"""
			# trial 2
			coordY = (bdict[i][1]+bdict[i][3])/2
			print i,coordY
			pushD = 0
			pushU = 0
			for j in ylist:
				if coordY != j and abs(coordY-j) < yShift:
					if pushD:
						print " pushD"
						coordY = coordY + yShift
					elif pushU:
						print " pushU"
						coordY = coordY - yShift
					elif coordY < j:
						print " pushU set"
						coordY = coordY - yShift
						pushU  = 1
					elif coordY > j:
						print " pushD set"
						coordY = coordY + yShift
						pushD  = 1
			ylist.append(coordY)
			print "",coordY
			
			writer.draw_line(bdict[i][0]      ,bdict[i][1],
							 bdict[i][0]+linkW,coordY)
			writer.draw_line(bdict[i][2]      ,bdict[i][3],
							 bdict[i][2]+linkW,coordY)
			writer.draw_string(i,bdict[i][0]+linkW,coordY)
			countBlk += 1
		
		writer.write_footer()
	
		print "Done!"

	

	##
	# Draw subfamily information based on a passed matrix
	#
	# @param matrix  a tab delimited file with id as 1st token and subfam info
	#                as the second. No header line.
	##	
	def annotate(self,matrix,shape):
		
		svg_write = SVGWriter.svg_writer(matrix+".svg")
		# read matrix into list, use as order
		mlist = futil.file_to_list(matrix,1,"\t")
		Xoff   = 50  # x offset
		Yoff   = 50  # y offset
		idW    = 100 # id width
		Xinc   = 40  # x increment
		Yinc   = 30  # y increment 
		Ygap   = 10
		Xlen   = 10		
		
		def draw_brac(Ystart,subfam,d):
			svg_write.draw_line(Xoff+idW     , Ystart+Ygap,
								Xoff+idW     , Yoff+(d+0.5)*Yinc-Ygap)
			svg_write.draw_line(Xoff+idW-Xlen, Ystart+Ygap,
								Xoff+idW     , Ystart+Ygap)			
			svg_write.draw_line(Xoff+idW-Xlen, Yoff+(d+0.5)*Yinc-Ygap,
								Xoff+idW     , Yoff+(d+0.5)*Yinc-Ygap)

			svg_write.draw_string(subfam, Xoff*2 + idW, 
								  (Ystart+Yoff+(d+0.5)*Yinc)/2)

		def draw_rect(Ystart,subfam,d):
			svg_write.draw_rect(Xoff+idW-Xlen,Ystart+Ygap,
								Xlen,Yoff+(d+0.5)*Yinc-Ygap-Ystart-Ygap)
			svg_write.draw_string(subfam, Xoff*2 + idW, 
								  (Ystart+Yoff+(d+0.5)*Yinc)/2)
		
		W = Xoff + idW +(len(mlist[0])-1)*Xinc # identifier, each column
		H = Yoff + len(mlist)*Yinc
		svg_write.write_header(W,H)
			
		d = 1
		not_specified = []
		subfam = ""
		Ystart = 0
		for i in mlist:
			#print "",i[0]
			svg_write.draw_string(i[0], Xoff, Yoff + (d+0.5)*Yinc)	
			c = 0
			
			if i[1] != "":
				if subfam != "":
					if i[1] == subfam:
						pass
					else:
						if shape == "bracket":
							draw_brac(Ystart,subfam,d)
						else:
							draw_rect(Ystart,subfam,d)
						subfam = i[1]
						Ystart = Yoff + (d+0.5)*Yinc
				else:
					subfam = i[1]
					Ystart = Yoff + (d+0.5)*Yinc
			elif subfam != "":
				if shape == "bracket":
					draw_brac(Ystart,subfam,d)
				else:
					draw_rect(Ystart,subfam,d)
				subfam = ""
			
			d += 1
		
		# draw the last one if subfam is not ""
		if subfam != "":
			if shape == "bracket":
				draw_brac(Ystart,subfam,d)
			else:
				draw_rect(Ystart,subfam,d)
			
		svg_write.write_footer()

	def get_coords(self,coords):
	
		inp   = open(coords,"r")
		cdict = {"gene":[],"exon":[]}
		for i in inp.readlines():
			i = self.onl(i)
			llist = i.split("\t")
			if llist[0] == "organism":
				cdict["org"] = llist[1]
			elif llist[0] == "chr":
				cdict["chr"] = llist[1]
			elif llist[0] == "start":
				cdict["str"] = llist[1]
			elif llist[0] == "end":
				cdict["end"] = llist[1]
			elif llist[0] == "gene":
				cdict["gene"].append(llist[1:])
			elif llist[0] == "exon":
				cdict["exon"].append(llist[1:])
			
		return cdict


	##
	# 
	# @param c1    first coord file, with the following format:
	#              organism \t mm
	#              chr      \t 3
	#              start    \t 1345524323
	#              end      \t 1356664333
	#              gene     \t L \t R \t ORI \t NAME (ori is either 1 or -1)
	#              gene     \t ...
	# @param c2    second coord file
	# @param blast the blast output of the first species against the second.
	#              NOT the other way around. At this point, the blast output
	#              is preferably generated by BlastUtility.blast_window. The
	#              query sequences are from windows and the coordinates are
	#              coded in the names. The coordinates are relative to the start
	#              of the query chromosome segment.
	# @param FLAG  for calling by another script [1], default [0]. This impact
	#              how the coordiantes are passed.
	#
	##
	def draw_match(self,c1,c2,blast,idT,lenT,Yzoom,FLAG=0):
	
		# read coords
		c1 = self.get_coords(c1)
		c2 = self.get_coords(c2)

		Xoff   = float(100)      # x offset
		Yoff   = float(100)      # y offset
		#Yzoom  = 1e3             # Y axis zoom factor
		chrW   = 10              # chr width
		
		chr1S  = int(c1["str"])  # chr1 start
		chr1E  = int(c1["end"])  # chr1 end
		chr2S  = int(c2["str"])
		chr2E  = int(c2["end"])
		
		maxH   = max(chr1E-chr1S,chr2E-chr2S)
		
		print "Draw chromosome and genes..."
		writer = SVGWriter.svg_writer("%s_idT%i_lenT%i.svg" % (blast,idT,lenT))
		writer.write_header(Xoff*4,float(maxH)/Yzoom+Yoff*2)
		
		# draw chr		
		writer.set_font_size(20)
		writer.draw_string(c1["org"],Xoff  ,20)
		writer.draw_string(c2["org"],Xoff*3,20)
		writer.set_font_size(16)
		writer.draw_string(c1["chr"],Xoff  ,50)
		writer.draw_string(c2["chr"],Xoff*3,50)
				
		S1 = str(float(chr1S)/1e6)
		S2 = str(float(chr2S)/1e6)			
		E1 = str(float(chr1E)/1e6)
		E2 = str(float(chr2E)/1e6)			
		
		writer.set_font_size(12)
		writer.draw_string("%sMb" % S1[:S1.find(".")+3],Xoff  -70,Yoff-10)
		writer.draw_string("%sMb" % S2[:S2.find(".")+3],Xoff*3+20,Yoff-10)
		writer.draw_string("%sMb" % E1[:E1.find(".")+3],Xoff  -70,
									float(chr1E-chr1S)/Yzoom+Yoff+20)
		writer.draw_string("%sMb" % E2[:E2.find(".")+3],Xoff*3+20,
									float(chr2E-chr2S)/Yzoom+Yoff+20)

		writer.set_fill_c("black")
		writer.draw_rect(Xoff  ,Yoff,chrW,float(chr1E-chr1S)/Yzoom)
		writer.draw_rect(Xoff*3,Yoff,chrW,float(chr2E-chr2S)/Yzoom)

		
		# draw genes, glist is [L,R,ori,name]
		# @param ftype 0: gene [default]
		#              1: exon 
		def draw_feature(flist,order,ftype=0):
			if order == 1:
				s = chr1S
				x = Xoff   - chrW
				r = Xoff   - 55
			else:
				s = chr2S
				x = Xoff*3 + chrW*2
				r = Xoff*3 + 30
			yL = float((int(flist[0])-s)/Yzoom)+Yoff
			yR = float((int(flist[1])-s)/Yzoom)+Yoff
			if ftype == 0:
				writer.set_fill_c("white")
				writer.set_stroke_w(0.25)
				writer.draw_rect(x-chrW/2,yL,chrW,yR-yL)
				if flist[2] == "-1":
					writer.set_fill_c("black")
					writer.draw_triangle(x-chrW/2,yL,x,yL,x-chrW/2,yL+(yR-yL)/10)
					writer.draw_triangle(x+chrW/2,yL,x,yL,x+chrW/2,yL+(yR-yL)/10)
				elif flist[2] == "1":
					writer.set_fill_c("black")
					writer.draw_triangle(x-chrW/2,yR,x,yR,x-chrW/2,yR+(yL-yR)/10)
					writer.draw_triangle(x+chrW/2,yR,x,yR,x+chrW/2,yR+(yL-yR)/10)
				writer.set_stroke_w(0)
				writer.draw_string(flist[3],r,yL+(yR-yL)/2)
			elif ftype == 1:
				writer.set_fill_c("lightgreen")
				writer.draw_rect(x-chrW/2,yL,chrW,yR-yL)
					
			

		writer.set_stroke_c("black")
		writer.set_fill_c("white")
		writer.set_font_size(5)

		# draw gene
		genes = c1["gene"]
		for j in genes:
			draw_feature(j,1)
		genes = c2["gene"]
		for j in genes:
			draw_feature(j,2)
		
		# draw exon
		exons = c1["exon"]
		for j in exons:
			draw_feature(j,1,1)
		exons = c2["exon"]
		for j in exons:
			draw_feature(j,2,1)
				
		
		# parse blast data table based on identity and length threshold
		print "Parse blast output..."
		bdict = self.parse_blast(blast,idT,lenT,chr2S)
		print " %i qualified" % len(bdict.keys())
		
		# draw matches
		print "Draw matching areas..."
		
		def draw_pair(mlist):
			# the passed list should be [id,qL,qR,sL,sR]
			if   mlist[0] >= 90:
				c = "red"
			elif mlist[0] >= 80:
				c = "orange"
			elif mlist[0] >= 70:
				c = "yellow"
			elif mlist[0] >= 60:
				c = "green"
			elif mlist[0] >= 50:
				c = "blue"
			else:
				c = "black"

			writer.set_stroke_c(c)
			writer.set_fill_c(c)
			writer.set_stroke_w(0.25)
			
			# draw left match
			yL = float(mlist[1])/Yzoom + Yoff
			yR = float(mlist[2])/Yzoom + Yoff
			if yL > yR:
				writer.draw_triangle(Xoff+chrW*3/2       ,yL,
									 Xoff+chrW*3/2-chrW/8,yR,
									 Xoff+chrW*3/2+chrW/8,yR)
				y1 = yR + (yL-yR)/2.0
			else:
				writer.draw_triangle(Xoff+chrW*3/2       ,yR,
									 Xoff+chrW*3/2-chrW/8,yL,
									 Xoff+chrW*3/2+chrW/8,yL)
				y1 = yL + (yR-yL)/2.0
			# draw right match
			yL = float(mlist[3])/Yzoom + Yoff
			yR = float(mlist[4])/Yzoom + Yoff
			if yL > yR:
				writer.draw_triangle(Xoff*3-chrW/2       ,yR,
									 Xoff*3-chrW/2-chrW/8,yL,
									 Xoff*3-chrW/2+chrW/8,yL)
				y2 = yR + (yL-yR)/2.0					 
			else:
				writer.draw_triangle(Xoff*3-chrW/2       ,yL,
									 Xoff*3-chrW/2-chrW/8,yR,
									 Xoff*3-chrW/2+chrW/8,yR)
				y2 = yL + (yR-yL)/2.0
			# link matches
			writer.draw_line(Xoff  +chrW*3/2,y1,
							 Xoff*3-chrW/2  ,y2,c)
		
		# function call
		for i in bdict.keys():
			draw_pair(bdict[i])
			
	
		writer.write_footer()	
			
	#
	# This is for drawing syntenic information and, genes of interests
	#
	# @param coords chr coord file with the following format:
	#               [sp1_chr][ori][L][R][sp2_chr][ori][L][R]
	# @param csize  file with [sp][chr][size], sp is either a or b.
	# @param name   a list of genes with [sp][chr][id][L-coord], sp
	#               is either a or b.
	#
	def draw_synteny(self,coords,csize,name,yzoom):
		
		print "Synteny coords:",coords
		print "Chr sizes     :",csize
		print "Gene coords   :",name
		
		# read chr size into a nested dict, sp as key, {chr:size} as value
		cdict = {}
		inp = open(csize,"r")
		inl = inp.readline()
		a_total = b_total = 0
		a_count = b_count = 0
		while inl != "":
			inl = self.onl(inl)
			llist = inl.split("\t")
			
			# convert chr string to number, take care of non-number situations
			try:
				llist[1] = int(llist[1])
			except ValueError:
				pass
			if cdict.has_key(llist[0]):
				cdict[llist[0]][llist[1]] = float(llist[2])
			else:
				cdict[llist[0]] = {llist[1]:float(llist[2])}
			# get total size
			if llist[0] == "a":
				a_total += float(llist[2])
				a_count += 1
			elif llist[0] == "b":
				b_total += float(llist[2])
				b_count += 1
			else:
				print "Chr name format problem, QUIT!"
				sys.exit(0)
			inl = inp.readline()
		
		# read synteny info into a list
		slist = []
		inp = open(coords,"r")
		inl = inp.readline()
		while inl != "":
			inl = self.onl(inl)
			llist = inl.split("\t")
			try:
				llist[0] = int(llist[0])
			except ValueError:
				pass
			llist[2] = float(llist[2])
			llist[3] = float(llist[3])
			try:
				llist[4] = int(llist[4])
			except ValueError:
				pass
			llist[6] = float(llist[6])
			llist[7] = float(llist[7])
						
			slist.append(llist)
			inl = inp.readline()
		
		# read genes into a dict, sp as key, a list as value with chr,id,coord
		gdict = {}
		inp = open(name,"r")
		inl = inp.readline()
		while inl != "":
			inl = self.onl(inl)
			llist = inl.split("\t")
			try:
				llist[1] = int(llist[1])
			except ValueError:
				pass
			llist[3] = float(llist[3])
			if gdict.has_key(llist[0]):
				gdict[llist[0]].append(llist[1:])
			else:
				gdict[llist[0]] = [llist[1:]]
			inl = inp.readline()
		
		# sort the gdict value first by chr then by coord
		idict = {} # intermediate dict
		gdict2= {} # sp as key, the nested tdict as value, this will be used
		           # in verifying if a block flanks a gene.
		for i in gdict:
			tlist = gdict[i]
			tdict = {}        # chr as key, a dict as value with
			                  # coord as key, a list of chr,id,coord as value
			for j in tlist:
				if tdict.has_key(j[0]):
					# assuming all coords are unique
					tdict[j[0]][j[2]] = j
				else:
					tdict[j[0]] = {j[2]:j}
			gdict2[i] = tdict
			chrs = tdict.keys()
			chrs.sort()
			idict[i] = []
			for j in chrs:
				crds = tdict[j].keys()
				crds.sort()
				for k in crds:
					idict[i].append(tdict[j][k])
		gdict = idict

		# draw chr first
		
		CHR_GAP = 30
		CHR_W   = 5
		BLOCK_W = 20
		X_OFF   = 200
		Y_OFF   = 100
		TEXT_H  = 25
		Y_ZOOM  = yzoom
		
		# see which sp has larger dimension and set SVG dimension
		if a_total + a_count*CHR_GAP > b_total + b_count*CHR_GAP:
			ysize = Y_OFF*2 + a_total/Y_ZOOM + a_count*CHR_GAP
		else:
			ysize = Y_OFF*2 + b_total/Y_ZOOM + b_count*CHR_GAP
		writer = SVGWriter.svg_writer(coords+".svg")
		writer.write_header(X_OFF*4,ysize)
		
		countSP = 1
		# store lastS in a nested list, so where chr started can be used later
		last = [{},{}]
		for i in cdict:             # iterate sp
			ikeys = cdict[i].keys() # sort chr in order
			ikeys.sort()
			countChr = 1
			lastS    = 0            # the last southmost coord
			if countSP == 2:
				F = -1
			else:
				F = 1
			for j in ikeys:         # iterate chr
				last[countSP-1][j] = lastS
				csize = cdict[i][j]/Y_ZOOM
				writer.draw_rect(X_OFF*countSP,Y_OFF+lastS,CHR_W,csize)
				writer.draw_string("chr %s"%str(j),X_OFF*countSP-CHR_W*2*F,
							       Y_OFF+lastS)
				lastS += csize + CHR_GAP
				countChr += 1
			countSP += 1
		

		# draw syntenic regions, [sp1_chr][ori][L][R][sp2_chr][ori][L][R]
		# check if the regions contain genes of interests
		for i in slist:
			
			# check if the blocks flank any gene
			sp1_flank = 0
			sp2_flank = 0
			
			# gdict2: sp:{chr1:{coord:info,...},
			#             chr2:{...           }}
			# check first a species
			try:
				clist = gdict2["a"][i[0]].keys()	
				for j in clist:
					if i[2] < j and i[3] > j:
						sp1_flank = 1
						break
			except KeyError:   # this happens when no gene is on that chr
				pass
			# check b species
			try:
				clist = gdict2["b"][i[4]].keys()	
				for j in clist:
					if i[6] < j and i[7] > j:
						sp2_flank = 1
						break
			except KeyError:   # this happens when no gene is on that chr
				pass
			color = ""
			if sp1_flank and sp2_flank:
				color = "green"
			elif sp1_flank or sp2_flank:
				color = "red"
			else:
				color = "gray"
			writer.set_fill_c(color)
			writer.set_stroke_c(color)
			
			# draw a species
			blockY = Y_OFF + i[2]/Y_ZOOM + last[0][i[0]] 
			blockL = (i[3]-i[2])/Y_ZOOM                      
			centA  = blockY + blockL/2.0                     
			writer.draw_rect(X_OFF-(BLOCK_W-CHR_W)/2,blockY,BLOCK_W,blockL)
			
			# draw b species
			blockY = Y_OFF + i[6]/Y_ZOOM + last[1][i[4]] 
			blockL = (i[7]-i[6])/Y_ZOOM                     
			centB  = blockY + blockL/2.0                    
			writer.draw_rect(X_OFF*2-(BLOCK_W-CHR_W)/2,blockY,BLOCK_W,blockL)
			
			# link
			writer.draw_line(X_OFF  -(BLOCK_W-CHR_W)/2+BLOCK_W,centA,
							 X_OFF*2-(BLOCK_W-CHR_W)/2        ,centB)
							
		# draw genes
		writer.set_fill_c("black")
		writer.set_stroke_c("black")
		# gdict: [sp]:[chr][id][L-coord]
		for i in gdict:         # iterate sp
			if i == "a":
				geneX1 = X_OFF  - (BLOCK_W-CHR_W)     # near string
				geneX2 = X_OFF  + CHR_W               # away from string
				geneX3 = geneX1 - 70                  # string x coord
				geneX4 = X_OFF  - (BLOCK_W-CHR_W)*5/2 # point to string
				idx    = 0
			else:
				geneX2 = X_OFF*2
				geneX1 = X_OFF*2 + (BLOCK_W-CHR_W)*3/2
				geneX3 = geneX1  + (BLOCK_W-CHR_W)*2 + 5
				geneX4 = X_OFF*2 + (BLOCK_W-CHR_W)*3
				idx    = 1
			lastY = 0           # remember the last y coord of the gene
			
			for j in gdict[i]:  # iterate features
				geneY1 = Y_OFF + last[idx][j[0]] + j[2]/Y_ZOOM
				geneY2 = geneY1
				writer.draw_line(geneX1,geneY1,geneX2,geneY2)
				
				geneY2 = 0
				if geneY1 < lastY + TEXT_H:
					geneY2 = lastY + TEXT_H
				else:
					geneY2 = geneY1
				lastY = geneY2
				writer.draw_line(geneX1,geneY1,geneX4,geneY2)
				writer.draw_string(j[1],geneX3,geneY2)
		
		writer.write_footer()		

	# 0         1           2           3                 4
	# Query id, Subject id, % identity, alignment length, mismatches, 
	# 5             6         7       8         9       10       11
	# gap openings, q. start, q. end, s. start, s. end, e-value, bit score
	#
	# Query id contain information, 
	#
	def parse_blast(self,blast,idT,lenT,chr2S):
		
		print " idT :",idT
		print " lenT:",lenT
		
		inp    = open(blast,"r")
		inline = inp.readline()
		c      = 0
		bdict  = {}
		while inline != "":
			if inline[0] == "#" or inline == "\n" or inline == "\r\n":
				pass
			else:
				L = inline.split("\t")
				# verify thresholds
				if float(L[2]) >= idT and int(L[3]) >= lenT:
					# the query id should have coords after underscore
					try:
						qL = int(L[0][L[0].find("_")+1:L[0].find("-")])
					except ValueError:
						print "Query ID do not contain coord info!"
						print "Assume to start from position 1."
						qL = 1
					# list into dict [id,qL,qR,sL,sR]
					bdict[c] = [float(L[2])   ,
								int(L[6])+qL-1,int(L[7])+qL-1,
								int(L[8])     ,int(L[9])]
					c += 1
			inline = inp.readline()
			
		return bdict
	

	def onl(self,line):
		if line[-2:] == "\r\n":
			line = line[:-2]
		elif line[-1] == "\n":
			line = line[:-1]
		return line
	
	def help(self):
		print "\n -f    Function to run:"
		print "         color_matrix - generate a matrix of blocks with colors"
		print "             REQUIRES: -matrix, -scheme. OPTIONAL: shape, desc"
		print "         list_unique - generate a list and count the number of"
		print "             unique entries in an matrix. REQUIRES: matrix"
		print "         draw_domain - use the output of SMARTUtility.parse_html"
		print "             to generate domain graphics. REQUIRES: -gdom,"
		print "             -scheme. OPTIONAL: -R, -X, -Y, -H, -A"
		print "         draw_obj - take a set of coords and draw objects. NEED:"
		print "             -coords, -scheme. OPTIONAL: -size"
		print "         draw_bar - draw histogram like svg. REQUIRES: -matrix"
		print "             OPTIONAL: -xzoom"
		print "         draw_chr - draw locations of sequences on chr. REQUIRES:"
		print "             -dat, -name"
		print "         draw_chr2- draw duplicated regions, REQUIRES: -dat,-dup"
		print "         annotate - draw a line between subfamily members."
		print "         draw_chr3- draw features. REQUIRES: dat (without chr"
		print "             column), cname, cleng. OPTIONAL: start, X"
		print "             REQUIRES: -matrix"
		print "         draw_match- draw matching areas between two segments"
		print "             REQUIRES: c1,c2,blast,idT,lenT. OPTIONAL: yzoom"
		print "         draw_synteny - REQUIRES: coords,csize,name, OPT: yzoom"
		print "         draw_matrix - NEED: matrix,dom, OPTIONAL: step"
		print " -matrix The text file with matrix, first line is header line."
		print " -scheme For color_matrix, means color scheme used in the format:"
		print "             [key1]:[color1],[key2]:[color2],..."
		print "         For draw_domain and draw_obj, the file with domain name,"
		print "             shape, and color specified"
		print " -gdom   the gdom output. For draw_matrix, this is a file with"
		print "         [name][L][R]"
		print " -R      reference domain"
		print " -X      X compression factor, 0 to 1, default 1. For draw_chr3,"
		print "         this is set to 10 by default"
		print " -Y      Y compression factor, 0 to 1, default 1"
		print " -H      Block height, default 20 pixels, integer only"
		print " -A      AA to pixel ratio, default 3"
		print " -Yinc   Y increment in pixels, default 30"
		print " -coords The coordinate file with [id][x,y]. For synteny, this"
		print "         is a special file with spec described in the code doc."
		print " -size   The size of the object to be drawn."
		print " -csize  chr size, a file with [sp][chr][size]"
		print " -shape  Specify the shape of objects, rect [default] or circle"
		print "         (not for annotate) or bracket (for annotate only)"
		print " -desc   1st line is description [1] or not [0, default]"
		print " -xzoom  zoom in x axis, integer only"
		print " -yzoom  compression in y axis, default 1e6"
		print " -dat    chr coord table"
		print " -name   name file, for synteny, [sp][chr][id][L_coord]"
		print " -dup    file with [block_id][L_gene_id][R_gene_id]"
		print " -c1     coord file 1"
		print " -c2     coord file 2"
		print " -blast  blast output of sequences of c1 as query, c2 as subject"
		print " -idT    match identity threshold"
		print " -lenT   match length threshold"
		print " -window window size for draw_matrix. Default is 1. This is set"
		print "         if the x axis label are not incremented by 1"
		print " -cname  chr name"
		print " -cleng  chr length"
		print " -start  starting position"
		print " -text   print sequence name [1] or not [0,default]"
		print ""
		sys.exit(0)

#-------------------------------------------------------------------------------
futil     = FileUtility.file_util()

if __name__ == '__main__':

	function = matrix  = scheme = R = gdom = coords = dat = name = dup = \
			   c1 = c2 = blast = csize = cname = cleng = "" 
	X = Y = 1.0
	A = 3.0
	H = 20
	Yinc = 30
	size = 10
	xzoom = 100
	yzoom = 1e6
	shape = "rect"
	lenT  = 100
	idT   = 80
	desc  = start = text = 0
	dblocks = draw_blocks()
	step  = 1
	
	for i in range(1,len(sys.argv),2):
		if sys.argv[i] == "-f":
			function   = sys.argv[i+1]
		elif sys.argv[i] == "-matrix":
			matrix     = sys.argv[i+1]
		elif sys.argv[i] == "-scheme":	      
			scheme     = sys.argv[i+1]
		elif sys.argv[i] == "-gdom":	      
			gdom       = sys.argv[i+1]
		elif sys.argv[i] == "-R":
			R          = sys.argv[i+1]
		elif sys.argv[i] == "-X":
			X          = float(sys.argv[i+1])
		elif sys.argv[i] == "-Y":	      
			Y          = float(sys.argv[i+1])
		elif sys.argv[i] == "-H":
			H          = int(sys.argv[i+1])
		elif sys.argv[i] == "-A":	      
			A          = float(sys.argv[i+1])
		elif sys.argv[i] == "-Yinc":
			Yinc       = int(sys.argv[i+1])
		elif sys.argv[i] == "-coords":
			coords     = sys.argv[i+1]
		elif sys.argv[i] == "-size":
			size       = int(sys.argv[i+1])
		elif sys.argv[i] == "-shape":
			shape      = sys.argv[i+1]
		elif sys.argv[i] == "-xzoom":
			xzoom      = float(sys.argv[i+1])
		elif sys.argv[i] == "-yzoom":
			yzoom      = float(sys.argv[i+1])
		elif sys.argv[i] == "-dat":
			dat     = sys.argv[i+1]
		elif sys.argv[i] == "-name":
			name       = sys.argv[i+1]
		elif sys.argv[i] == "-dup":
			dup        = sys.argv[i+1]
		elif sys.argv[i] == "-c1":
			c1     = sys.argv[i+1]
		elif sys.argv[i] == "-c2":
			c2       = sys.argv[i+1]
		elif sys.argv[i] == "-blast":
			blast        = sys.argv[i+1]
		elif sys.argv[i] == "-idT":
			idT       = int(sys.argv[i+1])
		elif sys.argv[i] == "-lenT":
			lenT      = int(sys.argv[i+1])
		elif sys.argv[i] == "-desc":
			desc      = int(sys.argv[i+1])
		elif sys.argv[i] == "-csize":
			csize     = sys.argv[i+1]
		elif sys.argv[i] == "-step":
			step      = int(sys.argv[i+1])
		elif sys.argv[i] == "-cname":
			cname     = sys.argv[i+1]
		elif sys.argv[i] == "-cleng":
			cleng     = int(sys.argv[i+1])
		elif sys.argv[i] == "-start":
			start     = int(sys.argv[i+1])
		elif sys.argv[i] == "-text":
			text      = int(sys.argv[i+1])
		else:
			print "Unknown argument:",sys.argv[i]
			print "Quit!"
			sys.exit(0)

	if function == "color_matrix":
		if matrix == "" or scheme == "":
			print "\nNeed matrix file and color scheme\n"
			dblocks.help()
		dblocks.color_matrix(matrix,scheme,shape,desc)
	elif function == "list_unique":
		if matrix == "":
			print "\nNeed matrix file"
			dblocks.help()
		dblocks.list_unique(matrix)
	elif function == "get_anything":
		if fasta == "" or coords == "":
			print "\nNeed fasta file and coordinates\n"
			dblocks.help()
		dblocks.color_matrix(matrix,scheme)
	elif function == "draw_domain":
		if gdom == "" or scheme == "":
			print "\nNeed gdom and scheme files\n"
			dblocks.help()
		dblocks.draw_domain(gdom,scheme,R,X,Y,H,A,Yinc)			
	elif function == "draw_obj":
		if coords == "" or scheme == "":
			print "\nNeed coord and scheme files\n"
			print " -h for help"
			sys.exit(0)
		dblocks.draw_obj(coords,scheme,size)			
	elif function == "draw_bar":
		if matrix == "":
			print "\nNeed matrix file\n"
			print " -h for help"
			sys.exit(0)
		dblocks.draw_bar(matrix,xzoom)			
	elif function == "draw_chr":
		if dat == "" or name == "":
			print "\nNeed coord dat file and name file\n"
			print " -h for help"
			sys.exit(0)
		dblocks.draw_chr(dat,name)			
	elif function == "draw_chr2":
		if dat == "" or dup == "":
			print "\nNeed coord dat file and duplication file\n"
			print " -h for help"
			sys.exit(0)
		dblocks.draw_chr2(dat,dup)			
	elif function == "draw_chr3":
		if dat == "" or cname == "" or cleng == "":
			print "\nNeed dat, cname, and cleng\n"
			print " -h for help"
			sys.exit(0)
		dblocks.draw_chr3(dat,cname,cleng,start,X,text)			
	elif function == "annotate":
		if matrix == "":
			print "\nNeed matrix file\n"
			print " -h for help"
			sys.exit(0)
		dblocks.annotate(matrix,shape)			
	elif function == "draw_match":
		if c1 == "" or c2 == "" or blast == "":
			print "\nNeed coord files and blast output\n"
			print " -h for help"
			sys.exit(0)
		dblocks.draw_match(c1,c2,blast,idT,lenT,yzoom)
	elif function == "draw_synteny":
		if coords == "" or csize == "" or name == "":
			print "\nNeed coords, csize, and name.\n"
			print " -h for help"
			sys.exit(0)
		dblocks.draw_synteny(coords,csize,name,yzoom)
	elif function == "draw_matrix":
		if matrix == "" or gdom == "":
			print "\nNeed matrix and gdom coord files.\n"
			print " -h for help"
			sys.exit(0)
		dblocks.draw_matrix(matrix,step,gdom)
							
	elif function == "test":
		gdom = ["At5g38990	1	SIGNAL|16	441	TRANS|23	525	STYKc|275	880",
				"At5g39000	1	SIGNAL|20	443	TRANS|23	518	STYKc|275	873"]
		matrix = "test.pairs.rates"
		dblocks.combine(gdom,matrix)			

	else:
		print "Unknown function:",function
		dblocks.help()
