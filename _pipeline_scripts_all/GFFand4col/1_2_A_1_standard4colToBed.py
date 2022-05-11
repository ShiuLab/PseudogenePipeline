#This script is designed to take what I believe to be a standard 4-col file
#format and generate a 4 column MCScanX file which has in the past been called
#the .bed format (though the format has been changed since it was called that)
#Created by David E. Hufnagel on 3-8-2012

import sys

inp = open(sys.argv[1])      #"standard" 4-col input file (gene, chromo, start, stop)
out = open(sys.argv[2], "w") #MCScanX 4-col output file (spe, gene, start, stop)



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



for line in inp:
    lineLst = line.split("\t")
    coor1 = ""
    coor2 = ""
    coor1, coor2 = ReverseCoords(lineLst, coor1, coor2)
    newLine = "%s\t%s\t%s\t%s\n" % (lineLst[1], lineLst[0], coor1, coor2)
    out.write(newLine)




inp.close()
out.close()
