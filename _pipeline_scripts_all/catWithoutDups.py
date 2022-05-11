#This script is designed to concatenate multiple files without duplicate lines
#Created by David E. Hufnagel on May 13, 2013
import sys

files = sys.argv[1]          #files to concatenate seperated by ","
out = open(sys.argv[2], "w") #output concatenated files



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Open all the files and extract all lines into a nonrepeating list
bigLst = []
fileLst = files.strip().split(",")
for filex in fileLst:
    fd = open(filex)
    for line in fd:
	if not line.startswith("#"):
            if not line in bigLst:
	        bigLst.append(line)

#Output all lines from the nonrepeating list
for line in bigLst:
    out.write(line)

out.close()
