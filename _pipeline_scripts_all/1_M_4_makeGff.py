#This script is designed to take Guarav's pre-GFF file, a protein coding gene
#4col file and a reference name file and make a two GFF3 outputs ready for dispersion.
#Created by David E. Hufnagel on July 1, 2013
#WARNING: HIGHLY SITUATION SPECIFIC
import sys

four = open(sys.argv[1])       #input 4col radish Protein coding genes file
gff = open(sys.argv[2])        #input big pre-gff file with all sorts of info
ref = open(sys.argv[3])        #reference file for translating inputs in the big file
out = open(sys.argv[4], "w")   #output GFF3 file of maker gene predictions


def OrderCoords(coor1, coor2):
    coor1 = int(coor1);coor2 = int(coor2)
    
    if coor1 < coor2:
        small = coor1
        big = coor2
    elif coor2 < coor1:
        small = coor2
        big = coor1
    else:
        #print "***ORDERCOORDS ERROR HERE***"
        return coor1,coor2

    return small, big

def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through ref and make a dict of key: oldName val: newName
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")

        #mRNA line
        refDict[lineLst[0]] = lineLst[1]
        #gene line
        refDict[lineLst[0].split("-mRNA")[0]] = lineLst[1]

#Go through four and make a list of gene names (goodLst)
goodLst = []
for line in four:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        gnName = lineLst[0]
        goodLst.append(gnName)

#Go through gff (only "maker" results)
lineDict = {}
chromoNumLst = []
coordDict = {}
for line in gff:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[1] == "maker":
            ##gather data
            chromo = lineLst[0]
            chromoNum = int(chromo.strip("RrC"))
            start, stop = OrderCoords(lineLst[3], lineLst[4])
            
            ##Make lineDict of key: chromoNum_start_stop val: line
            key = "%s_%s_%s\n" % (chromoNum, start, stop)
            SaveIntoDict(key, line, lineDict)

            ##Make chromoNumLst (just a list of chromoNums)
            chromoNumLst.append(chromoNum)

            ##Make coordDict of key: chromoNum val: [(start, stop),()...]
            SaveIntoDict(chromoNum, (start, stop), coordDict)

#Remove dups and sort chromoNumLst
chromoNumSet = set(chromoNumLst)
chromoNumLst = list(chromoNumSet)
chromoNumLst.sort()
#print chromoNumLst[-1]

#Go through chromoNumLst, 
#print "%s lines to process:\n" % (len(chromoNumLst))
for chromoNum in chromoNumLst:
    ##Get coords from coordDict, remove duplicates and sort it
    coords = set(coordDict[chromoNum])
    coords = list(coords)
    coords.sort()
    ##Go through coords
    for pair in coords:
        ###Use lineDict to get info
        key = "%s_%s_%s\n" % (chromoNum, pair[0], pair[1])
        for line in lineDict[key]:
            #print "1: ", line.strip()
            lineLst = line.strip().split("\t")

            ###If possible Convert oldNames to newNames
            oldName = lineLst[8].split("ID=")[1].split(";")[0]
            if oldName in refDict:
                newName = refDict[oldName]
            else:
                #we need to do special work for exons and CDSs
                if lineLst[2] == "CDS" or lineLst[2] == "exon":
                    newName = refDict[oldName.split(":")[0]] + oldName.split("mRNA-")[1][1:]

                #if its just not in the dict
                else:
                    print oldName, "not in dict!"
                    newName = oldName
            
            ###If gene in goodLst Output info
            gnName = newName.split(":")[0]
            #print lineLst
            #print gnName
            #print oldName.split("-mRNA")[0]
            #print oldName
            #print newName
            #print
            #if gnName in goodLst:
                #print "2: ", line.strip()
            oldGnName = oldName.split("-mRNA")[0]
            #print line
            #print oldGnName
            #print line.replace(oldGnName,gnName)
            #print oldGnName in line.replace(oldGnName,gnName)
            #print
            out.write(line.replace(oldGnName,gnName))




four.close()
gff.close()
ref.close()
out.close()
