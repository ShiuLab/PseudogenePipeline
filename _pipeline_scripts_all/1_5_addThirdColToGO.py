#This script is designed to add the third column to a 2 col GO file
# geneName\tGOID --> geneName\tGOID\tGOname
#Created by David E. Hufnagel on June 28, 2012

import sys

inp = open(sys.argv[1])      #input 2col GO file
goRef = open(sys.argv[2])    #reference 2col GO terms file
out = open(sys.argv[3], "w") #output 3col GO file





#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through goRef and extract the file into a dict
goRefDict = {}    #a dict of key: GOID val: name
for line in goRef:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        goRefDict[lineLst[0]] = lineLst[1]

#Go through inp and add GO names from goRefDict
for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        GOID = lineLst[1].strip()
        if not GOID.startswith("EC"):
            GOname = goRefDict[GOID]
            newLine = "%s\t%s\t%s\n" % (lineLst[0], GOID, GOname)
            out.write(newLine)




inp.close()
goRef.close()
out.close()
