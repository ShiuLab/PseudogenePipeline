#This script is designed to take a file and remove duplicates without changing
#the order of the file (always removes the secondary and tertiary (etc.) instances
#Created by David E. Hufnagel on July 3, 2013
import sys

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")



def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through the file and make a dict of key: lineNum val: line
inp.seek(0)
lineDict = {} #key = line
cntDict = {}  #key = cnt
cnt = 1
for line in inp:
    if not line.startswith("#"):
        SaveIntoDict(line, cnt, lineDict)
        SaveIntoDict(cnt, line, cntDict)
        cnt += 1

#Go through lineDict and output info from cntDict
for line in lineDict:
    cnt = lineDict[line][0]
    newLine = cntDict[cnt][0]
    out.write(newLine)




inp.close()
out.close()
