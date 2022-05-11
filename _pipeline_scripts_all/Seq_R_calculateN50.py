#This script is designed to calculate N50 from a .size file
#Created by David E. Hufnagel

import sys

inp = open(sys.argv[1])

#get total number of base pairs and make list from biggest to smallest
numbp = 0
bigLst = []
for line in inp:
    lineLst = line.split("\t")
    numbp += int(lineLst[1])
    bigLst.append(int(lineLst[1]))

bigLst.sort(reverse= True)

#go through list and find N50
tempTotal = 0
for n in bigLst:
    if not tempTotal < float(numbp / 2):
        break
    tempTotal += n
    
#print tempTotal
print n
#print max(bigLst)
#print numbp
    






inp.close()
