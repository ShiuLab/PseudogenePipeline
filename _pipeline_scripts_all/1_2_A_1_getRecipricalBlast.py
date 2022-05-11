#This script was designed to make a .rbm file from an all v. all protein blast
#file.
#Created by David E. Hufnagel on March 8, 2012 (header on June 7, 2012)
#Updated: added exception for A. lyrata on June 7, 2012
#WARNING: GAURAV HAS INFORMED ME THAT MY METHOD IS WRONGHEADED

import sys

inp = open(sys.argv[1])       #all vs all blast input
out = open(sys.argv[2], "w")  #.rmb output

lastKey = ""
bigDict = {}
for line in inp:
    lineLst = line.split("\t")
    currKey = lineLst[0]
          
    if currKey not in bigDict:

        #for A. lyrata 1st col
        if lineLst[0].startswith("sc") or lineLst[0].startswith("fg") or lineLst[0].startswith("Al"):
            if not (lineLst[1].startswith("sc") or lineLst[1].startswith("fg") or lineLst[1].startswith("Al")):
                key = lineLst[0]
                value = (lineLst[1], line) #a tuple with (value, line)
                bigDict[key] = value
        #for A. lyrata 2nd col
        elif lineLst[1].startswith("sc") or lineLst[1].startswith("fg") or lineLst[1].startswith("Al"):
            key = lineLst[0]
            value = (lineLst[1], line) #a tuple with (value, line)
            bigDict[key] = value
        
        #for non A. lyrata
        else:
            if lineLst[0][:2] != lineLst[1][:2]: #to protect from when the best match is a match to itself (or paralogs)
                key = lineLst[0]
                value = (lineLst[1], line) #a tuple with (value, line)
                bigDict[key] = value

trashDict = {}
for k in bigDict:
    if k not in trashDict:
        v = bigDict[k][0]
        if v in bigDict:
            if k == bigDict[v][0]:
                #out.write(k + "\t" + v + "\n")
                out.write(bigDict[k][1])
                trashDict[k] = 1
                trashDict[v] = 1

inp.close()
out.close()
