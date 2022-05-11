#This script is designed to convert a GMAP output "gff" file into an m8 blast
#format.
#Created by David E. Hufnagel on July 30, 2012

import sys

inp = open(sys.argv[1])      #input GMAP "gff" format file
out = open(sys.argv[2], "w") #output m8 blast format file




#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        query = lineLst[0]
        subj = lineLst[8]
        pID = lineLst[5]
        qStart = lineLst[6].split(".")[0]
        qEnd = lineLst[6].split(".")[-1]
        length = int(qEnd) - int(qStart) + 1
        sStart = lineLst[9].split(".")[0]
        sEnd = lineLst[9].split(".")[-1]
        newLine = "%s\t%s\t%s\t%s\t-\t-\t%s\t%s\t%s\t%s\t-\t-\n" % (query, subj, pID, length, qStart, qEnd, sStart, sEnd)
        out.write(newLine)


inp.close()
out.close()
