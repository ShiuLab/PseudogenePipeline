#This script is designed to go through a gff and write longest=1 on every mRNA line
#Created by David E. Hufnagel on July16, 2013
import sys, os

gff = open(sys.argv[1])      #the input gff
out = open(sys.argv[2], "w") #the output modified gff

#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through gff and add longest=1 to all mRNA lines before parent
for line in gff:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[2] == "mRNA":
            newLine = line.replace(";Parent=",";longest=1;Parent=")
            out.write(newLine)
        else:
            out.write(line)



gff.close()
out.close()
