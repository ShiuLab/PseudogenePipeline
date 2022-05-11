#This script is designed to remove species IDs in the format, "Ath|AT3G34353",
#from a disable_count file.
#Created by: David E. Hufnagel on June 4, 2012

import sys

inp = open(sys.argv[1])      #the input file
out = open(sys.argv[2], "w") #the output file



for line in inp:
    if line.startswith("#"):
        lineLst = line.split(" ")
        nameA = "#" + lineLst[0].split("|")[1]
        nameB = "|".join(lineLst[1].split("|")[1:])

        newLine = "%s %s %s %s %s %s %s %s" % (nameA, nameB, lineLst[2], lineLst[3], lineLst[4], lineLst[5], lineLst[6], lineLst[7])
        out.write(newLine)

    else:
        out.write(line)






inp.close()
out.close()
