#This script was designed to take a 4col file of pseudogenes, a 23col .dom
#protein domain file (generated from HMMscan) and a 3col file of transposable
#element related domains in the format:
#PFAMID    shortDescription    longDescription
#and output a 4col of pseudogenes without transposable element related domains
#Created by David E. Hufnagel on May 29, 2013
import sys

pseuFd = open(sys.argv[1])   #input 4col of pseudogenes
hmm = open(sys.argv[2])      #input 23col .domtblout file with pseudogenes and their PFAM domains
transp = open(sys.argv[3])   #input 3col of TE domains
out = open(sys.argv[4], "w") #output 4col of pseudogenes without TE domains



#remove all instances of an item in a list
def RemoveAll(listx, item):
    while item in listx:
        listx.remove(item)



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through transp and make a list of TE related PFAMIDs called badDom
badDom = []
for line in transp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        pfamID = lineLst[0]
        badDom.append(pfamID)

print "%s TE domains imported" % (len(badDom))

#Go through hmm and make a set of pseudogenes with PFAMIDs in badDom called badPseu
badPseu = set()
for line in hmm:
    if not line.startswith("#"):
        lineLst = line.strip().split(" ")
        RemoveAll(lineLst, "")
        pfamID = lineLst[1].split(".")[0]
        pseu = lineLst[3]
        if pfamID in badDom:
            badPseu.add(pseu)

print "%s pseudos contain TE domains" % (len(badPseu))

#Go through pseuFd and output all lines where the pseuName is not present in badPseu
for line in pseuFd:
    if not line.startswith("#"):
        lineLst = line.strip().split()
        if not lineLst[0] in badPseu:
            out.write(line)





pseuFd.close()
hmm.close()
transp.close()
out.close()
