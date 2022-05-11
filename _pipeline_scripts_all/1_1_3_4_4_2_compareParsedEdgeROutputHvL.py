#This script was created to take in 2 egdeR output files, one for RHvRL and one
#for KHvKL, and put the relationships into 3 categories each for upregulation
#and downregulation: 1) RHvRL difference only 2) KHvKL difference only and
#3) both RHvRL and KHvKL difference.
#Created by David E. Hufnagel on Nov 21, 2012
#WARNING: HIGHLY SITUATION SPECIFIC

import sys

RHvRL = open(sys.argv[1])
KHvKL = open(sys.argv[2])
out = open(sys.argv[3], "w")



def RemoveBlancs(theList):
    newLst = []
    for item in theList:
        if item != "":
            newLst.append(item)

    return newLst



#Get the sets of genes with differential expression between RH and RL
RnegSet = set([])
RposSet = set([])
Rtitle = RHvRL.readline()
for line in RHvRL:
    lineLst = RemoveBlancs(line.strip("\n").split(" "))
    if lineLst[1].startswith("-"):
        RnegSet.add(lineLst[0])
    else:
        RposSet.add(lineLst[0])

#Get the sets of genes with differential expression between KH and KL
KnegSet = set([])
KposSet = set([])
Ktitle = KHvKL.readline()
for line in KHvKL:
    lineLst = RemoveBlancs(line.strip("\n").split(" "))
    if lineLst[1].startswith("-"):
        KnegSet.add(lineLst[0])
    else:
        KposSet.add(lineLst[0])

#Get the sets of genes with differential expression between RH and RL and differential expression between KH and KL
BothNegSet = RnegSet.intersection(KnegSet)
BothPosSet = RposSet.intersection(KposSet)

#get final R and K only gene sets
RNegOnly = RnegSet.difference(BothNegSet)
RPosOnly = RposSet.difference(BothPosSet)
KNegOnly = KnegSet.difference(BothNegSet)
KPosOnly = KposSet.difference(BothPosSet)

#Output info
out.write("##RHvRL only, overexpressed, (%s) genes\n" % (len(RPosOnly)))
for gene in RPosOnly:
    out.write("%s\n" % (gene))

out.write("##RHvRL only, underexpressed, (%s) genes\n" % (len(RNegOnly)))
for gene in RNegOnly:
    out.write("%s\n" % (gene))

out.write("##KHvKL only, overexpressed, (%s) genes\n" % (len(KPosOnly)))
for gene in KPosOnly:
    out.write("%s\n" % (gene))
    
out.write("##KHvKL only, underexpressed, (%s) genes\n" % (len(KNegOnly)))
for gene in KNegOnly:
    out.write("%s\n" % (gene))

out.write("##Both, overexpressed, (%s) genes\n" % (len(BothPosSet)))
for gene in BothPosSet:
    out.write("%s\n" % (gene))

out.write("##Both, underexpressed, (%s) genes\n" % (len(BothNegSet)))
for gene in BothNegSet:
    out.write("%s\n" % (gene))




RHvRL.close()
KHvKL.close()
out.close()
