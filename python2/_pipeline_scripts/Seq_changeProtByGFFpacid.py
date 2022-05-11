#This script was designed to take a protein fasta file with the name format:
#>name|PACid:number
#look up the PACid in the gff file and replace the name with the gene name  
#Created by David E. Hufnagel on Mar 18, 2013

import sys

prot = open(sys.argv[1])     #input prot file
gff = open(sys.argv[2])      #input gff file
out = open(sys.argv[3], "w") #output prot file with fixed names



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through the protein file and make a list of PACids
#pacLst = []
#for line in  prot:
#    if not line.startswith("#"):
#        if line.startswith(">"):
#            pacID = line.strip().split("PACid:")[1]
#            pacLst.append(pacID)

#Go through gff and make a dict of key: PACid val: geneName
nameDict = {}
for line in gff:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[2] == "mRNA":
            pacID = lineLst[8].split("pacid=")[1].split(";")[0]
            geneName = lineLst[8].split("Parent=")[1].split(";")[0]
            #geneName = geneName[1:] ### WARNING DYNAMIC LINE ###
            nameDict[pacID] = geneName

#Go through protein file again and output file with geneNames
prot.seek(0)
for line in  prot:
    if not line.startswith("#"):
        if line.startswith(">"):
            pacID = line.strip().split("PACid:")[1]
            geneName = nameDict[pacID]
            #geneName = "Alyr" + geneName### WARNING DYNAMIC LINE ###
            out.write(">%s\n" % (geneName))
        else:
            out.write(line)



prot.close()
gff.close()
out.close()
