#This script is designed to convert the MCScan collinearity file to a 5 column
#file with block number, gene1, gene2, contig1 and contig2 in that order.
#Created by David E. Hufnagel on 3-13-2012
#WARNING: IF GENE1 AND GENE2 IN THE COLL FILE ARENT ALWAYS IN THE SAME ORDER SPECIES-WISE, THIS SCIPT WILL NEED TO BE MODIFIED
#WARNING: THE ORDER OF BED1 AND BED2 MUST BE IN THE SAME ORDER AS IN THE COLLINEARITY FILE

import sys

coll = open(sys.argv[1])     #.collinearity file for syntenic block information
bed1 = open(sys.argv[2])     #first .bed file for coordinate information
bed2 = open(sys.argv[3])     #second .bed file for coordinate information
out = open(sys.argv[4], "w") #output file to contain the data




def SkipHash(fd):
    for L in fd:
        if not L.startswith("#"):
            break


   
#write the users command line prompt on the first line of the output file.
out.write('#python %s\n'%(' '.join(sys.argv)))

#skip past the title in the .collinearity file
SkipHash(coll)
#go through coll file and get Rr genes and block numbers (list of tuples.  each tuple has (block, gene1, gene2)
genes = []
for line in coll:
    if not line.startswith("##"):
        lineLst = line.split("\t")
        block = lineLst[0].split("-")[0].strip()
        theTuple = (block, lineLst[1], lineLst[2])
        genes.append(theTuple)

#skip past the title in the first .bed file
SkipHash(bed1)    
#go through the first bed file and get the contig info for the first column of genes (list of tuples.  each tuple has (block, gene1, gene2, contig1).
genes2 = []
for ln in bed1:
    lnLst = ln.split("\t")
    for tup in genes:
        if tup[1] == lnLst[1]:
            newTuple = (tup[0], tup[1], tup[2], lnLst[0])
            genes2.append(newTuple)

#skip past the title in the Second .bed file
SkipHash(bed2)
#go through the second bed file, get the contig info for the second column of genes and write the info into the output file.
for l in bed2:
    lLst = l.split("\t")
    for t in genes2:
        if t[2] == lLst[1]:
            newLine = "%s\t%s\t%s\t%s\t%s\n" %\
                      (t[0], t[1], t[2], t[3], lLst[0])
            out.write(newLine)




coll.close()
bed1.close()
bed2.close()
out.close()
