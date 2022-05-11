#This script is designed to add GOIDs to .pqvalue files as an additional (9th)
#column
#Created by David E. Hufnagel on July 2, 2012

import sys

pqVal = open(sys.argv[1])    # .pqvalue file to add the GOID to       (x)
GOref = open(sys.argv[2])    # GO reference file with GOIDs and names
out = open(sys.argv[3], "w") # output .pqvalue file with GOIDs        (x.wGOID)





def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)




        
#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through ref and put it into a dict
refDict = {}
for line in GOref:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        key = lineLst[1].replace("'","_")
        key = key.replace(",","_")
        SaveIntoDict(key, lineLst[0], refDict)

#Go through pqVal and add GOIDs in the last column
for line in pqVal:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        GOID = refDict[lineLst[0]][0]
        newLine = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %\
                  (lineLst[0], lineLst[1], lineLst[2], lineLst[3],\
                   lineLst[4], lineLst[5], lineLst[6], lineLst[7].strip(), GOID)
        out.write(newLine)


        

pqVal.close()
GOref.close()
