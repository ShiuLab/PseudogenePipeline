#This script is designed to make an MCScanX 4col file from my 4col file.
#Gaurav's format (All Rr):
#RrC1_p1	RrC1	12152	12496
#Created by David E. Hufnagel on Nov 16, 2012
#WARNING: HIGHLY SITUATION SPECIFIC

import sys

inp = open(sys.argv[1])      #input Dave format 4col file
out = open(sys.argv[2], "w") #output MCScanX format 4col file



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#do actual processing
for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip("\n").split("\t")
        chromo = "Rr%s" % (lineLst[1].split("C")[1])
        gene = lineLst[0]
        start = lineLst[2]
        end = lineLst[3]
        newLine = "%s\t%s\t%s\t%s\n" % (chromo, gene, start, end)
        out.write(newLine)



inp.close()
out.close()

