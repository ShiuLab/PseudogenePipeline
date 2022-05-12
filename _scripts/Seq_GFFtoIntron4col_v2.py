#this script is designed to make an intron 4col file from a gff file with genes,
#UTRs and CDS.
#Created by David E. Hufnagel on Feb 25, 2013
import sys, os

gff = sys.argv[1]            #input gff file
out = open(sys.argv[2], "w") #output intron 4col file

def OrderCoords(coor1, coor2):
    coor1 = int(coor1);coor2 = int(coor2)
    
    if coor1 < coor2:
        small = coor1
        big = coor2
    elif coor2 < coor1:
        small = coor2
        big = coor1
    else:
        #print "WARNING: Intron start/end cooridnates are the same" #
        return coor1,coor2

    return small, big

#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through gff and for each gene...
gene = ""
intronLst = []
gff2 = open(gff + ".noAlt")
#gff2 = open(gff)
for line in gff2:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")

        #determine start coord, stop coord, gene name and chromo
        if lineLst[2] == "gene":
            
            #the very first gene
            if not gene == "":
                #make new intron coords
                nonIntronLst.sort()
                intronLst = []

                #do add intron coords to intronLst
                for coord in nonIntronLst:
                    if not (intronLst == []):
                        #if predictions are running into each other (created by UTRs), merge them
                        if intronLst[-1] == coord[0]:
                            intronLst = intronLst[:-1]
                        else:
                            intronLst.append(coord[0]-1)
                    intronLst.append(coord[1]+1)
                else:
                    intronLst = intronLst[:-1]

                #put intronLst back into tuples of 2
                newIntronLst = []  #a temporary list to put the tupled intronLst into
                cnt = 0
                for coord in intronLst:
                    if not cnt % 2:
                        newIntronLst.append((intronLst[cnt], intronLst[cnt+1]))
                    cnt += 1
                intronLst = newIntronLst              

                #output new intron coords in 4col format with the name formated as
                #gene_intron# where # is a unique number generated using a count
                num = 1
                for coord in intronLst:
                    newName = "%s_Intron%s" % (gene, str(num).zfill(5))
                    newLine = "%s\t%s\t%s\t%s\n" % (newName, chromo, coord[0], coord[1])
                    out.write(newLine)
                    num += 1

            #gather data
            nonIntronLst = []
            small, big = OrderCoords(lineLst[3], lineLst[4])
            start = small
            stop = big
            gene = lineLst[8].split("Name=")[1]
            chromo = lineLst[0]
            
        #gather all UTR and CDS info
        elif lineLst[2] == "CDS" or "UTR" in lineLst[2]:
            small, big = OrderCoords(lineLst[3], lineLst[4])
            nonIntronLst.append((small, big))

#Write info from the last gene
else:
#the very first gene
    if not gene == "":
        #make new intron coords
        nonIntronLst.sort()
        intronLst = []

        #do add intron coords to intronLst
        for coord in nonIntronLst:
            if not (intronLst == []):
                #if predictions are running into each other (created by UTRs), merge them
                if intronLst[-1] == coord[0]:
                    intronLst = intronLst[:-1]
                else:
                    intronLst.append(coord[0]-1)
            intronLst.append(coord[1]+1)
        else:
            intronLst = intronLst[:-1]

        #put intronLst back into tuples of 2
        newIntronLst = []  #a temporary list to put the tupled intronLst into
        cnt = 0
        for coord in intronLst:
            if not cnt % 2:
                newIntronLst.append((intronLst[cnt], intronLst[cnt+1]))
            cnt += 1
        intronLst = newIntronLst              

        #output new intron coords in 4col format with the name formated as
        #gene_intron# where # is a unique number generated using a count
        num = 1
        for coord in intronLst:
            newName = "%s_Intron%s" % (gene, str(num).zfill(5))
            newLine = "%s\t%s\t%s\t%s\n" % (newName, chromo, coord[0], coord[1])
            out.write(newLine)
            num += 1
            
gff2.close()
out.close()
