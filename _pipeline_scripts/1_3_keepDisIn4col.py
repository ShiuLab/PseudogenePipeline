#This script is designed to take a .disable_count file and a 4col file of
#pseudogenes, where the 4col is a subcategory of the .disable_count file and
#output a .disable_count file with only the pseudogenes from the 4col file
#Created by David E. Hufnagel on June 5, 2013
import sys

dis = open(sys.argv[1])      #input disable count file
four = open(sys.argv[2])     #input 4col file 
ref = open(sys.argv[3])      #input ref file with codeNames in the 1st col and bigName in the 2nd col
out = open(sys.argv[4], "w") #output disable count file containing only pseudogenes from the 4col file



def OrderCoords(coor1, coor2):
    coor1 = int(coor1);coor2 = int(coor2)
    
    if coor1 < coor2:
        small = coor1
        big = coor2
    elif coor2 < coor1:
        small = coor2
        big = coor1
    else:
        print "***ORDERCOORDS ERROR HERE***"
        return coor1,coor2

    return small, big



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through ref and make a dict of key: bigName val: codeName
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        refDict[lineLst[1]] = lineLst[0]

#Go through four and make a list of pseudogenes to keep, goodLst
goodLst = []
for line in four:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        goodLst.append(lineLst[0])

#Go through dis and output lines of good pseudogenes
for line in dis:
    if not line.startswith("#python"):
        if line.startswith("#"):
            lineLst = line.strip().split(" ")
            codeName = (";".join(lineLst[:3]) + ";" + "|".join(lineLst[4:]))[1:]
            if refDict[codeName] in goodLst:
                out.write(line)
                doWrite = True
            else:
                doWrite = False
        else:
            if doWrite == True:
                out.write(line)






dis.close()
four.close()
ref.close()
out.close()
