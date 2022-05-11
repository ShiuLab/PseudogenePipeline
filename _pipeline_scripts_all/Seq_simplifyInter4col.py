#This script is designed to take the 4col intergenic coordinates file outputed
#from the get_intergenic function of GFFUtil.py and make a more standard 4col
#with new names.  It also makes the smaller coordinate always first.
#Created by David E. Hufnagel on Jan 5, 2012
import sys

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")
spe = sys.argv[3]            #a short species identifier Ex: A.thaliana --> Atha



def OrderCoords(coor1, coor2):
    if type(coor1) != int or type(coor2) != int:
        print "***TYPE ERROR***"
    
    if coor1 < coor2:
        small = coor1
        big = coor2
    elif coor2 < coor1:
        small = coor2
        big = coor1
    else:
        small = coor1
        big = coor2
        print "***ERROR: SAME START AND STOP***"

    return small, big



cnt = 1
for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        num = str(cnt).zfill(6)
        name = "%sIntergenic%s" % (spe, num)
        small, big = OrderCoords(int(lineLst[1]), int(lineLst[2]))
        newLine = "%s\t%s\t%s\t%s\n" % (name, lineLst[0], small, big)
        out.write(newLine)
        cnt += 1




inp.close()
out.close()
