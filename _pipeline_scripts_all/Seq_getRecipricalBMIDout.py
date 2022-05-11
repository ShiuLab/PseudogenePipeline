#This script is designed to get iDs for all recipricol best matches from
#an all vs. all blast file, and output an rmb file with a 3rd column that is
#the %ID

import sys

inp = open(sys.argv[1])       #input file
out = open(sys.argv[2], "w")  #output file

lastKey = ""
bigDict = {}
for line in inp:
    lineLst = line.split("\t")
    currKey = lineLst[0]
          
    if currKey not in bigDict:
        if lineLst[0][:2] != lineLst[1][:2]: #to protect from when the best match is a match to itself
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
                ln = bigDict[k][1].split("\t")
                pair1 = ln[0]
                pair2 = ln[1]
                iD = float(ln[2])
                newline = "%s\t%s\t%s\n" % (pair1, pair2, iD)
                out.write(newline)
                trashDict[k] = 1
                trashDict[v] = 1

inp.close()
out.close()
