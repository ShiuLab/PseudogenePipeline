#This script was designed to take names from a 4col file and keep only fasta
#sequences with those names
#Created By David E. Hufnagel on Dec 17, 2012

import sys

fourCol = open(sys.argv[1])
fasta = open(sys.argv[2])
out = open(sys.argv[3], "w")



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

nameLst = []
for line in fourCol:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        nameLst.append(lineLst[0])

doWrite = False
for line in fasta:
    if not line.startswith("#"):
        if line.startswith(">"):
            name = line.strip("\n").strip(">") 
            if name in nameLst:
                doWrite = True
                out.write(line)
            else:
                doWrite = False
        else:
            if doWrite == True:
                out.write(line)




fourCol.close()
fasta.close()
out.close()
