#This script is designed to get % ID for all recipricol best matches from an all
#vs. all blast file.  The output should be something that can be put into MS
#Excel to make a histogram.
#Created by: David E. Hufnagel on April 2, 2012
#WARNING LIKELY DOESN'T WORK (GIVES TOO FEW RESULTS)

import sys

inp = open(sys.argv[1])       #input file
out = open(sys.argv[2], "w")  #output file
intv = float(sys.argv[3])     #the intervals for putting the sequence IDs into groups





#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#save input file into a dict
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
##        else:
##            print "same match!!: ", line

#make set of IDs
iDSet = set()
trashDict = {}
for k in bigDict:
    if k not in trashDict:
        v = bigDict[k][0]
        #print v
        #print k
        if v in bigDict:
            print "2"
            if k == bigDict[v][0]:
                print "3"
                ln = bigDict[k][1].split("\t")
                iD = float(ln[2])
                iDSet.add(iD)
                trashDict[k] = 1
                trashDict[v] = 1

#iterate through ID set and put iDs into groups based on intv
small = 0   #min value for particular group
big = intv  #max value for particular group

##make dict for each interval
x = intv
intvDict = {}
intvDict[0] = 0
while x < 100:
    intvDict[x] = 0
    x += intv

##iterate through iDset and put iDs in intvDict
print iDSet
for iden in iDSet:
    for key in intvDict:
        if iden >= key and iden < key+intv:
            if intvDict[key] == 0: #first time
                intvDict[key] = [iden]
            else:
                intvDict[key].append(iden)
        
##iterate through intvDict and output data into "histogram format"
print intvDict
for interval in intvDict:
    if intvDict[interval] != 0:
        newline = "%s\t%s\n" %(interval, len(intvDict[interval]))
        out.write(newline)



    


inp.close()
out.close()
