#This script was designed to make a gff of pseudogenes from a pseudogene 4col
#and a pseudogene disable_count file
#Created by David E. Hufnagel on May 14, 2013
import sys

four = open(sys.argv[1])     #input pseudogene 4col file
dis = open(sys.argv[2])      #input pseudogene disable_count file
ref = open(sys.argv[3])      #input pseudogene reference file
out = open(sys.argv[4], "w") #output pseudogene gff file



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

#Go through four and make a dict of key: psName val: (chromo, small, big)
#note: small and big is used instead of start and stop because orientation is
#acquired in the disable count file
fourDict = {}
for line in four:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        fourDict[lineLst[0]] = (lineLst[1], lineLst[2], lineLst[3])

#Go through ref and get a dict of key: bigName val: codeName for connecting
#names in four with names in dis
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        refDict[lineLst[1]] = lineLst[0]

#Go through dis, acquire information about orientation and pseudogene evidence
#code and output all info into the gff file
for line in dis:
    if not line.startswith("#python") and line.startswith("#"):
        lineLst = line.strip().split(" ")

        #determine orientation
        tempStart = lineLst[2].split(":")[1].split("-")[0]
        tempStop = lineLst[2].split(":")[1].split("-")[1]
        small, big = OrderCoords(tempStart, tempStop)
        if tempStart == str(small):
            ori = "+"
        else:
            ori = "-"

        #determine pseudogene evidence code
        evidence = ",".join(lineLst[4:])

        #gather info from refDict and fourDict
        bigName = (";".join(lineLst[:3]) + ";" + "|".join(lineLst[4:]))[1:]
        psName = refDict[bigName]
        psID = psName.split("_")[0]
        fgName = "_".join(psName.split("_")[1:])
        chromo = fourDict[psName][0]
        source = "ShiuLab"
        typex = "pseudogene"
        if ori == "+":
            start = fourDict[psName][1]
            stop = fourDict[psName][2]
        elif ori == "-":
            start = fourDict[psName][2]
            stop = fourDict[psName][1]
        score = "."
        phase = "."
        attributes = "ID=%s;Name=%s;Derives_from=%s;Note=pseudogene_evidence_code_%s"\
                     % (psID, psName, fgName, evidence) #the info spot with all kinds of info
        
        #output info
        newLine = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"\
                  % (chromo, source, typex, start, stop,\
                     score, ori, phase, attributes)
        
        out.write(newLine)





four.close()
dis.close()
ref.close()
out.close()
