#This script is designed to simply print the number of unique items in a given
#column
#Created by David E. Hufnagel on June 13, 2012

import sys

inp = open(sys.argv[1])     #input file to check for unique items
col = int(sys.argv[2]) - 1  #the column to check for unique items

theSet = set()
for line in inp:
    lineLst = line.split("\t")
    theSet.add(lineLst[col])

#print theSet
print len(theSet)



inp.close()
