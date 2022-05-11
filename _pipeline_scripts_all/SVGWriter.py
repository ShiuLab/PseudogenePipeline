

class svg_writer:

	def __init__(self,out_name):
		global string_c
		global fill_c
		global stroke_c
		global stroke_w
		global font_fam
		global font_siz
		global font_wgt
		global oup
		self.oup = open(out_name,"w")
		self.string_c = self.fill_c   = self.stroke_c = self.font_fam = \
				self.font_siz = self.font_wgt = ""
		self.stroke_w = 1
				
		
	def write_header(self,w,h):
		
		self.oup.write(\
		"<?xml version=\"1.0\" standalone=\"no\"?>"+\
		"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\"\n	 "+\
		"\"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n"+\
		"\n<svg width=\"%ipt\" height=\"%ipt\">\n" % (w,h) +\
		"\t<title>Draw blcok\n\t</title>\n");
	
	def write_footer(self):
		self.oup.write("</svg>")
		
		
	def draw_rect(self,x,y,w,h,stroke_w=0.5,G=0):
		
		if G:
			self.oup.write(
				"\t\t<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" />\n"%\
				(x,y,w,h))
		else:
			if self.fill_c == "":
				self.fill_c = "white"
			if self.stroke_c == "":
				self.stroke_c = "black"
			self.oup.write(
				"\t<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\"\n"%(x,y,w,h)+\
				"\t\tstyle=\"fill:"+self.fill_c+"; stroke:"+self.stroke_c+\
				"; stroke-width:%f\"/>\n" % self.stroke_w)
	
	def draw_rounded_corners(self,x,y,w,h,rx=0,ry=0,stroke_w=0.5,G=0):
	
 		if self.fill_c == "" or self.stroke_c == "" and not G:
			self.fill_c = "white"
			self.stroke_c = "black"   	
		
		if rx == 0:
			rx = x/4
		if ry == 0:
			ry = y/4
		
		if G:
			self.oup.write(
				"\t\t<rect x=\"%i\" y=\"%i\" rx=\"%i\" ry=\"%i\" "+\
				"width=\"%i\" height=\"%i\" />\n"%\
				(x,y,rx,ry,w,h))
		else:
			self.oup.write(
				"\t<rect x=\"%i\" y=\"%i\" " % (x,y)  +\
				"rx=\"%i\" ry=\"%i\" "       % (rx,ry)+\
				"width=\"%i\" height=\"%i\"\n"% (w,h) +\
				"\t\tstyle=\"fill:"+self.fill_c+"; stroke:"+self.stroke_c+\
				"; stroke-width:%f\"/>\n" % self.stroke_w)
	
	
	def draw_line(self,x1,y1,x2,y2,stroke_w=0.5,G=0):

		if type(x1) != type(0.1):
			x1 = float(x1)
		if type(x2) != type(0.1):
			x2 = float(x2)
		if type(y1) != type(0.1):
			y1 = float(y1)
		if type(y2) != type(0.1):
			y2 = float(y2)
 		if self.stroke_c == "" and not G:
			self.stroke_c = "black"   	

		if G:
			self.oup.write(
				"\t\t<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" />\n"%\
				(x1,y1,x2,y2))
		else:
			self.oup.write(
				"\t<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\"\n"%\
				(x1,y1,x2,y2)+\
				"\t\tstyle=\"stroke:"+self.stroke_c+\
				"; stroke-width:%f\" />\n" % self.stroke_w)	

	def draw_circle(self,cx,cy,r,stroke_w=0.5,G=0):

 		if self.fill_c == "" or self.stroke_c == "" and not G:
			self.fill_c = "white"
			self.stroke_c = "black"   	

		if G:
			self.oup.write("\t\t<circle cx=\"%i\" cy=\"%i\" r=\"%i\" />\n" % \
						  (cx,cy,r))
		else:
			self.oup.write(
				"\t<circle cx=\"%i\" cy=\"%i\" r=\"%i\"\n" % (cx,cy,r) +\
				"\t\tstyle=\"fill:"+self.fill_c+"; stroke:"+self.stroke_c+\
				"; stroke-width:%f\"/>\n" % self.stroke_w)		
		
	def draw_ellipse(self,cx,cy,rx,ry,stroke_w=0.5,G=0):
	
 		if self.fill_c == "" or self.stroke_c == "" and not G:
			self.fill_c = "white"
			self.stroke_c = "black"   	

		if G:
			self.oup.write(
				"\t\t<ellipse cx=\"%i\" cy=\"%i\" rx=\"%i\" ry=\"%i\" />\n" % \
				(cx,cy,rx,ry))
		else:
			self.oup.write(
				"\t<ellipse cx=\"%i\" cy=\"%i\" rx=\"%i\" ry=\"%i\"\n" % \
				(cx,cy,rx,ry) +\
				"\t\tstyle=\"fill:"+self.fill_c+"; stroke:"+self.stroke_c+\
				"; stroke-width:%f\"/>\n" % self.stroke_w)
	
	# clist is a string of integers that are x,y coordinate pairs, each coord
	# is separated from the others by ","
	def draw_polyline(self,cstr,stroke_w=0.5):
		
 		if self.fill_c == "" or self.stroke_c == "" and not G:
			self.fill_c = "white"
			self.stroke_c = "black"   	
		
		self.oup.write(
			"\t<polyline fill=\"%s\" stroke=\"%s\" stroke-width=\"%f\" " % \
			(self.fill_c,self.stroke_c,stroke_w) + "points=\"%s\" />\n" % cstr)
		
	
	#
	# This is written to temporarily do some jobs draw_polygon suppose to.
	#
	def draw_triangle(self,x1,y1,x2,y2,x3,y3,stroke_w=0.5):
	
		if type(x1) != type(0.1):
			x1 = float(x1)
		if type(x2) != type(0.1):
			x2 = float(x2)
		if type(x3) != type(0.1):
			x3 = float(x3)
		if type(y1) != type(0.1):
			y1 = float(y1)
		if type(y2) != type(0.1):
			y2 = float(y2)
		if type(y3) != type(0.1):
			y3 = float(y3)
			
 		if self.fill_c == "":
			self.fill_c = "white"
		if self.stroke_c == "":
			self.stroke_c = "black"   	

		coords = "%f,%f %f,%f %f,%f" % (x1,y1,x2,y2,x3,y3)
		self.oup.write(
			"\t<polygon points=\"%s\"\n" % coords +\
			"\t\tstyle=\"fill:"+self.fill_c+"; stroke:"+self.stroke_c+\
			"; stroke-width:%f\"/>\n" % self.stroke_w)
		
	
	# only deal with 3,5,or 6 sides only, the coord should be float
	#
	# @param ori   orientation, 0, left pointing; 1, right; 2, up; 3,down
	#              NOT IMPLEMENTED YET
	#
	def draw_polygon(self,x,y,w,h,sides,stroke_w=0.5,G=0,ori=0):
		
		x = float(x)
		y = float(y)
		w = float(w)
		h = float(h)
		
 		if self.fill_c == "" or self.stroke_c == "" and not G:
			self.fill_c = "white"
			self.stroke_c = "black"   	
		
		coords = ""
		if sides == 3:
			coords = "%f,%f %f,%f %f,%f" % \
					 (x,y+h/2,x+w,y,x+w,y+h)
		elif sides == 5:
			coords = "%f,%f %f,%f %f,%f %f,%f %f,%f" % \
					 (x,y+h/2,x+w/4,y,x+w,y+h/8,x+w,y+h/8*7,x+w/4,y+h)					 
		elif sides == 6:
			coords = "%f,%f %f,%f %f,%f %f,%f %f,%f %f,%f" % \
					 (x,y+h/2,x+w/5,y,x+w/5*4,y,x+w,y+h/2,x+w/5*4,y+h,x+w/5,y+h)
		
		if G:
			self.oup.write("\t\t<polygon points=\"%s\" />\n" % coords)					
		else:
			self.oup.write(
				"\t<polygon points=\"%s\"\n" % coords +\
				"\t\tstyle=\"fill:"+self.fill_c+"; stroke:"+self.stroke_c+\
				"; stroke-width:%f\"/>\n" % self.stroke_w)
			
	#
	# @param G    group or not
	# @param R    rotation angle
	#
	def draw_string(self,str,x,y,G=0,R=0): 
		
		if self.font_fam == "":
			self.font_fam = "sans-serif"
		if self.font_wgt == "":
			self.font_wgt = "normal"
		if self.font_siz == "":
			self.font_siz = 10
		if self.string_c == "":
			self.string_c = "black"
		
		# deal with trasnformation, only rotation at this moment.
		t_str = ""
		if R != 0:
			t_str = "transform=\"rotate(%i)\" " % R
		
		if G:
			self.oup.write("\t\t<text %sx=\"%i\" y=\"%i\">%s\n\t\t</text>\n" % \
						  (t_str,x,y,str))		
		else:
			self.oup.write(
				"\t<text font-family=\"%s\" " % self.font_fam + \
				"font-size=\"%ipt\" "         % self.font_siz + \
				"fill=\"%s\" "                % self.string_c + \
				"stroke=\"%s\" "              % self.stroke_c + \
				"stroke-width=\"%i\" "        % 0             + \
				"font-weight=\"%s\" "         % self.font_wgt + \
				"x=\"%i\" y=\"%i\" "          % (x,y)         + \
				"%s>%s"                       % (t_str,str)   + \
				"\n\t</text>\n")
		
	def set_string_c(self,c):
		self.string_c = c		   
	
	def set_fill_c(self,c):
		self.fill_c = c

	def set_stroke_c(self,c):
		self.stroke_c = c

	def set_stroke_w(self,w):
		self.stroke_w = w
	
	def set_font_fam(self,fam):
		self.font_fam = fam
	
	def set_font_size(self,size):
		self.font_siz = size

	def set_font_wgt(self,weight):
		self.font_wgt = weight
		
	def def_start(self):
		self.oup.write("\t<defs>\n")
	
	def def_end(self):
		self.oup.write("\t</defs>\n")
	
	def group_start(self,stroke_c,fill_c,stroke_w):
		self.oup.write("\t<g stroke=\"%s\" fill=\"%s\" stroke-width=\"%f\">\n"%
						(stroke_c,fill_c,stroke_w))
	
	def group_start_str(self,font_fam,font_wgt,font_siz,string_c):
		self.oup.write(\
			"\t<g font-family=\"%s\" " % font_fam + \
			"font-size=\"%ipt\" "      % font_siz + \
			"fill=\"%s\" "             % string_c + \
			"stroke=\"%s\" "           % "black"  + \
			"stroke-width=\"%i\" "     % 0        + \
			"font-weight=\"%s\">\n"    % font_wgt)
	
	def group_end(self):
		self.oup.write("\t</g>\n")


if __name__ == '__main__':
	
	writer = svg_writer("test.svg")
	writer.write_header(300,500)

	writer.group_start_str("sans-serif","normal",10,"black")
	writer.draw_string("TEST1",200,200,G=1)
	writer.draw_string("TEST2",200,240,G=1)
	writer.group_end()

	
	writer.group_start("red","black",2)	
	writer.draw_rect(10,10,40,20,G=1)
	writer.draw_line(100,100,200,200,G=1)
	writer.draw_circle(50,50,30,G=1)
	writer.draw_ellipse(50,100,20,40,G=1)
	writer.group_end()
	
	writer.group_start("blue","green",2)	
	writer.draw_polygon(10,200,40,20,3,G=1)
	writer.draw_polygon(10,250,40,20,5,G=1)
	writer.draw_polygon(10,300,40,20,6,G=1)
	writer.draw_rounded_corners(10,100,40,20,G=1)
	writer.group_end()
	
	
	writer.write_footer()

