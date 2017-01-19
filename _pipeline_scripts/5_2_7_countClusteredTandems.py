#count tandems
import sys

inp = open(sys.argv[1])



theSet = set()
for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip("\n").split("\t")
        #name1 = lineLst[0]
        #name2 = lineLst[1]
        #theSet.add(name1)
        #theSet.add(name2)

        for item in lineLst[2:]:
            theSet.add(item)

print len(theSet)

inp.close()

        
