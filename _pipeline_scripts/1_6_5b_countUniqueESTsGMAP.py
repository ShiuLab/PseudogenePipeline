#This script is designed to count the number of unique EST names in an
#unfiltered .gmap GMAP output file
#Created by David E. Hufnagel on Aug 24, 2012

import sys

inp = open(sys.argv[1])




theSet = set()
for line in inp:
    if line.startswith(">"):
        name = line.split(" ")[0].strip(">")
        theSet.add(name)

print "%s unique ESTs mapped" % (len(theSet))




inp.close()
