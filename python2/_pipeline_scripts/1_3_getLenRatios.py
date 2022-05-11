#This script was designed to extract pseudogene and protein coordinates and to
#output the lengths and ratios
#Created by David E. Hufnagel on May 29, 2012
#WARNING WRONGHEADED APPROACH, DOES NOT GIVE EXPECTED RESULTS

import sys

inp = open(sys.argv[1])      #input pseudogene file (...4col.true.condensed.noCM)
out = open(sys.argv[2], "w") #output file with length and ratio info



#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))
#writes the title
out.write("#pseu_len, prot_len, pseuL/protL\n")

inp.readline()
for line in inp:
    lineLst = line.split("\t")

    #get lengths
    pseuLen = abs(int(lineLst[2]) - int(lineLst[3]))
    protLen = abs(int(lineLst[0].split(";")[1].split("-")[0].split("|")[1])\
                  - int(lineLst[0].split(";")[1].split("-")[1]))

    ratio = float(pseuLen) / float(protLen)
    newLine = "%d\t%d\t%f\n" % (pseuLen, protLen, ratio)
    out.write(newLine)



inp.close()
out.close()
