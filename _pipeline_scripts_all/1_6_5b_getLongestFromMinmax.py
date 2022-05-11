#this script is designed to filter a parsed .minmax GMAP output file by keeping
#only the longest EST match reason
#Created by David E. Hufnagel on Aug 31, 2012

import sys

inp = open(sys.argv[1])      #input .minmax fully parsed GMAP output file
out = open(sys.argv[2], "w") #output .longest file




def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)




#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

matchDict = {}
for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        name = lineLst[0].split("|")[0]
        SaveIntoDict(name, (lineLst[1], lineLst[2], lineLst[3].strip()), matchDict)
        
        
for key in matchDict:
    val = matchDict[key]
    if len(val) == 1:
        newLine = "%s\t%s\t%s\t%s\n" % (key, val[0][0], val[0][1], val[0][2])
        out.write(newLine)
    elif len(val) >= 2:
        #get top size match index
        ind = 0
        topSize = 0
        topSizeInd = 0
        for match in val:
            size = abs(int(val[ind][2]) - int(val[ind][1])) + 1
            if size > topSize:
                topSize = size
                topSizeInd = ind
            ind += 1

        #output info
        newLine = "%s\t%s\t%s\t%s\n" % (key, val[topSizeInd][0], \
                                        val[topSizeInd][1], val[topSizeInd][2])
        out.write(newLine)
        
    else:
        print "ERROR HERE"




inp.close()
out.close()
