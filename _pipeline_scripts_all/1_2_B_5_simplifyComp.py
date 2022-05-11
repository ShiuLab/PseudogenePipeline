#This script is designed to convert .comp files from the format:
#AT2G25590	['fgenesh2_kg.4_2890_AT2G46660.1']	['Bra004517', 'Bra039267']	['RrC11425_p1_0.010', 'RrC8444_p1_0.007']
#to:AT2G25590	fgenesh2_kg.4_2890_AT2G46660.1	Bra004517,Bra039267	RrC11425_p1_0.010,RrC8444_p1_0.007
#Designed by David E. Hufnagel on 5-11-2012

import sys

inp = open(sys.argv[1])      #the input
out = open(sys.argv[2], "w") #the output



for line in inp:
    lineLst = line.split("\t")
    #for first column
    new0 = ""
    if lineLst[0] != "[]":
        new0 = lineLst[0]

    #for other columns
    newLst = []
    for part in lineLst[1:]:
        partLst = eval(part)
        newLst.append(",".join(partLst))

    if new0 != "":
        newLine = new0 + "\t"
    for p in newLst:
        if p != "":
            newLine = newLine + (p + "\t")
        
    newLine = newLine[:-1] + "\n"
    out.write(newLine)
        







out.close()
