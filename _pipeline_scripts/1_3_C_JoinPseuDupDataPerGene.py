#This script was designed to take an input data file containing pseudogenization
#and speciation timing generated from 1_3_C_GetPseudogenizationTime.py and an
#input data file containing duplication time generated from
#1_3_C_GetDupTimeFromKs.py.  Creates an output info file in the format:
#pseuName  inName  outName  pseuTime  dupTime speTime dif1 dif2
#where dif1 = dupTime-pseuTime  dif2 = speTime - dupTime
#Created by David E. Hufnagel on Dec 14, 2012

import sys

pseuSpe = open(sys.argv[1])  #input data file containing pseudogenization and speciation timing generated from 1_3_C_GetPseudogenizationTime.py
dup = open(sys.argv[2])      #input data file containing duplication time generated from 1_3_C_GetDupTimeFromKs.py.
out = open(sys.argv[3], "w") #Output info file




def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)


        

#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))
out.write("#pseuName  inName  outName  pseuTime  dupTime speTime dif1 dif2\n")
out.write("#where dif1 = dupTime-pseuTime  dif2 = speTime - dupTime\n")

#Go through pseuSpe and make a dict of key: inName val: [[pseuName, outName, pseuTime, speTime],]
mainDict = {}
for line in pseuSpe:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        SaveIntoDict(lineLst[1], [lineLst[0], lineLst[2], lineLst[6], lineLst[7]], mainDict)
        
#Go through dup and export all info in the format:
#pseuName  inName  outName  pseuTime  dupTime speTime dif1 dif2
#where dif1 = dupTime-pseuTime  dif2 = speTime - dupTime
for line in dup:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        names = lineLst[0].split(";")
        for inName in names:
            if inName in mainDict:
                for group in mainDict[inName]:
                    psName = group[0]
                    outName = group[1]
                    psTime = float(group[2])
                    dupTime = float(lineLst[1])
                    speTime = float(group[3])
                    dif1 = dupTime - psTime
                    dif2 = speTime - dupTime

                    newLine = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (psName, inName, outName, psTime, dupTime, speTime, dif1, dif2)
                    out.write(newLine)




pseuSpe.close()
dup.close()
out.close()
