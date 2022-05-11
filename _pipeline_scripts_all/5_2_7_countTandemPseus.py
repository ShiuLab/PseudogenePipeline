#This script is designed to take a tandem gene file created by Gaurav's script,
#find_tandemDuplicates.py in the format: Gene1  Gene2   Distance(# of genes)    Evalue
#and print the number of items, genes and pseudogenes in the input file
import sys

inp = open(sys.argv[1])  #.tandem tandem genomic features file


#Go through the tandem file and gather all, pseudos and genes into sets
allSet = set()
psSet = set()
fgSet = set()
for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")

        #Add to the allSet
        allSet.add(lineLst[0])
        allSet.add(lineLst[1])
        
        #Add the first item to the psSet or fgSet
        if lineLst[0].startswith("Ps"):
            psSet.add(lineLst[0])
        else:
            fgSet.add(lineLst[0])

        #Add the second item to the psSet or fgSet
        if lineLst[1].startswith("Ps"):
            psSet.add(lineLst[1])
        else:
            fgSet.add(lineLst[1])

#Determine the number of these genomic features by taking the lenght of the sets
allCnt = len(allSet)
psCnt = len(psSet)
fgCnt = len(fgSet)

#print the info
print
print "number of genomic features: %s" % (allCnt)
print "number of pseudogenes:      %s" % (psCnt)
print "number of genes:            %s" % (fgCnt)



inp.close()
