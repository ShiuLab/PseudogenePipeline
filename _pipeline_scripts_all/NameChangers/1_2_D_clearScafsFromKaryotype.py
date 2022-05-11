#This script was designed to clear out all 3 lines in the karyotype file that
#came from the 1 line in the est file where there was a name with "scaffold" in
#one of the 3 in the karyotype file.  This would indicate it is part of a small
#contig, and not a real chromosome.
#WARNING: HIGHLY SITUATION SPECIFIC
"""Algorithm:
1) Go through kar and extract offending scaffold names
2) Go through est and extract all three names to be removed
3) Go through kar again and extract bad lines
4) Go through kar again write lines that aren't bad lines"""

import sys

kar = open(sys.argv[1])      #Karyotype input file
est = open(sys.argv[2])      #EST input file
out = open(sys.argv[3], "w") #new karyotype file without the superfluous population





#1) Go through kar and extract offending scaffold names
nameLst = []
for line in kar:
    if "Scaffold" in line:
        lineLst = line.split("\t")
        for name in lineLst:
            if "Scaffold" in name and ":" in name and name not in nameLst:
                nameLst.append(name)

#2) Go through est and extract all three names to be removed
badNames = []
for name in nameLst:
    est.seek(0)
    for line in est:
        if name[2:] in line:
            for estName in line.split("\t"):
                if ":" in estName and "Rr" not in estName:
                    if estName.endswith("\n"):
                        estName = estName[:-1]
                    badNames.append(estName)

#3) Go through kar again and extract bad lines
kar.seek(0)
badLines = []
for line in kar:
    #print line
    for badN in badNames:
        #if badN not in line and line not in outLst:
        if badN in line:
            badLines.append(line)

#4) Go through kar again write lines that aren't bad lines
kar.seek(0)
for line in kar:
    if line not in badLines:
        out.write(line)
        
                





kar.close()
est.close()
out.close()
