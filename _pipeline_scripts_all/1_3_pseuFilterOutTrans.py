#This script is designed to remove pseudogenes from the filtered set that were
#identified by repeat masker
"""Algorithm:
1) go through rep and extract names into a set
2) go through pseu and filter
3) go through out4col and filter"""


import sys

pseu = open(sys.argv[1])         #input pseudogenes in .4col format
dis = open(sys.argv[2])          #input pseudogenes in .disable count format
rep = open(sys.argv[3])          #parsed repeat masker file
outDis = open(sys.argv[4], "w")  #output pseudogenes in .disable count format
out4col = open(sys.argv[5], "w") #output pseudogenes in .4col format



nameSet = set()
for line in rep:
    lineLst = line.split("\t")
    nameSet.add(lineLst[1])

for line in pseu:
    lineLst = line.split("\t")
    if not lineLst[0] in nameLst:
        out.write(line)






pseu.close()
dis.close()
rep.close()
outDis.close()
out4col.close()
