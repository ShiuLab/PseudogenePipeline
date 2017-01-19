#This script was designed to extract pseudogene and protein coordinates and to
#output the lengths and ratios
#Created by David E. Hufnagel on May 29, 2012

import sys, os

pseu = open(sys.argv[1])      #input pseudogene file (.4col.true.condensed.noCM)     
prot = open(sys.argv[2])      #input protein file (.size file)
out = open(sys.argv[3], "w")  #output file with length info





def SaveIntoDict(gene1, gene2, dictX):
    if gene1 not in dictX:
        dictX[gene1] = [gene2]
    else:
        print "problem!"
        dictX[gene1].append(gene2)





#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))
#writes the title
out.write("#pseu_len, prot_len, pseuL/protL, pseuL/protL*100\n")

#put protein lengths into a dict
#prot.readline()
protDict = {}
for line in prot:
    lineLst = line.split("\t")
    SaveIntoDict(lineLst[0], lineLst[1], protDict)

#get analagous pseudogene lengths and export results
pseu.readline()
pseuDict = {}
#print total
for line in pseu:
    lineLst = line.split("\t")
    pseuLen = abs(int(lineLst[0].split(";")[2].split(":")[0].split("-")[0]) -\
                  int(lineLst[0].split(";")[2].split(":")[0].split("-")[1])) + 1

    name = lineLst[0].split(";")[0]
    if name in protDict:
        protLen = protDict[name][0][:-1]
        ratio = float(pseuLen) / float(protLen)
        newLine = "%s\t%s\t%s\t%s\n" % (pseuLen, protLen, ratio, ratio * 100)#(ratio / total))
        out.write(newLine)
        
    else:
        print "problem2"





pseu.close()
out.close()
