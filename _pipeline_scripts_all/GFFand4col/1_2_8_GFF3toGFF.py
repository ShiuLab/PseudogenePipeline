#This script is intended to convert GFF3 files to GFF files (for use in MCScanX)
#In this case, it also edits the gene name to remove the string, and makes sure
#The coordinates are ordered in the positive direction
#"|protein_coding_gene"
#Created by David E. Hufnagel on 2-2-2012
#WARNING: THE WHOLE SCRIPT IS SITUATION SPECIFIC

import sys

def ReverseCoords(lineLst, coor1, coor2):
    coorA = int(lineLst[2])
    coorB = int(lineLst[3][:-1])

    if coorA < coorB:
        coor1 = str(coorA)
        coor2 = str(coorB)

    elif coorB < coorA:
        coor1 = str(coorB)
        coor2 = str(coorA)

    else:
        print "***ERROR HERE***"

    return coor1, coor2


inp = open(sys.argv[1])      #The input GFF3 file
out = open(sys.argv[2], "w") #The output edited GFF file
spe = sys.argv[3]            #The two letter species identifier

for line in inp:
    lineLst = line.split("\t")

    gene = lineLst[0][:9]

    #Reverse the coordinates where necessary
    coor1 = ""
    coor2 = ""
    coor1, coor2 = ReverseCoords(lineLst, coor1, coor2)
    
    newLine = "%s\t%s\t%s\t%s\n" % (spe, gene, coor1, coor2)

    out.write(newLine)

inp.close()
out.close()
