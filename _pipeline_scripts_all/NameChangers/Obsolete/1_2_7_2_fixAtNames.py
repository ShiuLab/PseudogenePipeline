#This script is designed to change names in the first column from
#AT1G57650.1 --> AT1G57650
#Created by David E. Hufnagel on 3-16-2012
#Fixed on 5-25-2012 (was giving just the first name and not ending the line)

import sys

inp = open(sys.argv[1])      #input file
out = open(sys.argv[2], "w") #output file

#writes the users command line prompt on the first line of the output file.
#out.write("#python %s\n"%(" ".join(sys.argv)))

#inp.readline()
for line in inp:
    lineLst = line.split("\t")
    newName = lineLst[0].split(".")[0]
    newLst = [newName, ]
    for x in lineLst[1:]:
        newLst.append(x)
    newLine = "\t".join(newLst)
    out.write(newLine)

inp.close()
out.close()
