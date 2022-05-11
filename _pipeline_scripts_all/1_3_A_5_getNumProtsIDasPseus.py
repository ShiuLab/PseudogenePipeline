#This script was designed to get a 4col file of protein fragments that were
#split and identified as pseudogenes from the pseudogene pipeline using a 4col
#of split proteins and a 4col of pseudogenes with full names in the first column
#Created by David E. Hufnagel on Aug 15, 2012

import sys

splitProts = open(sys.argv[1]) #input 4col file (new coords) of protein fragments split by fragmentation
pseus = open(sys.argv[2])      #input 4col pseudogene file (original names containing FG)
out = open(sys.argv[3], "w")   #output 4col file of protein fragments that were identified as pseudogenes



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through pseudogene file and make a list of FG names
FGlist = []
for line in pseus:
    if not line.startswith("#"):
        FG = line.split("\t")[0].split(";")[0]
        FGlist.append(FG)

#Go through splitProts and output lines with names in the FGlist
for line in splitProts:
    lineLst = line.split("\t")
    name = lineLst[0]
    if name in FGlist:
        out.write(line)



splitProts.close()
pseus.close()
out.close()
