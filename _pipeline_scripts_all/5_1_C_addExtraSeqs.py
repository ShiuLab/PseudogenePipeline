#This script is designed to look at two fastaFiles big and small, where the
#names in small are a subset of the names in big, and add all names and
#sequences present in big, but not small to small

import sys
small = open(sys.argv[1])
big = open(sys.argv[2])
out = open(sys.argv[3], "w")
   

#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through small and make a dict of key: name val: seq.  Also output all lines in this file
smallDict = {}
currName = ""   #the current name associated with a seq
seq = ""        #the current seq to be built up for each seq line
for line in small:
    if not line.startswith("#"):
        #gather info for dict
        if line.startswith(">"):
            if currName != "":
                smallDict[currName] = seq
            seq = ""
            currName = line.strip().strip(">")
        else:
            seq += line.strip()

        #output lines
        out.write(line)
            
#get the last seq on the way out
else:
    smallDict[currName] = seq

#Go through big and make a dict of key: name val: seq.
bigDict = {}
currName = ""   #the current name associated with a seq
seq = ""        #the current seq to be built up for each seq line
for line in big:
    if not line.startswith("#"):
        #gather info for dict
        if line.startswith(">"):
            if currName != "":
                bigDict[currName] = seq
            seq = ""
            currName = line.strip().strip(">")
        else:
            seq += line.strip()
            
#get the last seq on the way out
else:
    bigDict[currName] = seq

#Go through big and output names that arent in the small name dict
for bigName in bigDict:
    bigSeq = bigDict[bigName]
    if bigName not in smallDict:
        nameLine = ">%s\n" % (bigName)
        seqLine = "%s\n" % (bigSeq)
        out.write(nameLine)
        out.write(seqLine)




small.close()
big.close()
out.close()
