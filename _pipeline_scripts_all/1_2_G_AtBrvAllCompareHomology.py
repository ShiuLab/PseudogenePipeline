#This script was designed to compare the MCScanX_h input homology files from
#At v. Br and At v. Br v. Al v. Rr to see if the pairwise relationships are the
#same
#Created by David E. Hufnagel on Nov 9, 2012
"""Algorithm:
1) Go through atBr and make a dict of key: at Val: br and the inverse for each
   pairwise relationship
2) Go through allSpe and make a dict of key: at Val: br and the inverse for each
   pairwise relationship
3) Make a count of when relationships are in atBr but not allSpe, when
   relationships are in allSpe but not atBr and when relationships are in both"""

import sys

atBr = open(sys.argv[1])
allSpe = open(sys.argv[2])
out = open(sys.argv[3], "w")





#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#1) Go through atBr and make a dict of key: at Val: br and the inverse for
#each pairwise relationship
atBrDict = {}
cnt = 0
for line in atBr:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        atBrDict[lineLst[0]] = lineLst[1]
        atBrDict[lineLst[1]] = lineLst[0]
        cnt += 1

print cnt
print len(atBrDict)




atBr.close()
allSpe.close()
out.close()
