#This script is designed to basically transfer the pseudogene file,
#".4col.true.condensed.noCM", into the original .disable_count format
#Created by David E. Hufnagel on May 30, 2012
"""Algorithm:
1) Go through .4col.true.condensed.noCM file and make a list of names in the file
2) Go through original unfiltered .disable_count file and output lines matching
the names in the name list.
"""

import sys

orig = open(sys.argv[1])     #the .disable_count file
true = open(sys.argv[2])     #the .4col.true.condensed.noCM file
out = open(sys.argv[3], "w") #output file (.fullyFiltered.disable_count)



#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#1) Go through .4col.true.condensed.noCM file and make a list of names in the file
nameLst = []
print true.readline()
for line in true:
    lineLst = line.split("\t")
    nameLst.append(lineLst[0])

#2) Go through original unfiltered .disable_count file and output lines matching
#the names in the name list.
write = False
for line in orig:
    lineLst = line.split(" ")
    if line.startswith("#"):
        origName = (";".join(lineLst[:3]) + ";" + "|".join(lineLst[4:])[:-1])[1:]
        if origName in nameLst:
            out.write(line)
            write = True
        else:
            write = False
    else:
        if write == True:
            out.write(line)
    




orig.close()
true.close()
out.close()
