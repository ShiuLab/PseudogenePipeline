#Import protein files
import sys

protein = open(sys.argv[1])
proteinOut = open(sys.argv[2],"w")
ref = open(sys.argv[3],"w")



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
        print "***ERROR HERE***"

    return small, big



#Generate codenames for each name in the protein file
#Output new protein file with codenames
#Create a reference file linking codename to original name
cnt = 1
for line in protein:
    if not line.startswith("#"):
        if line.startswith(">"):
            oldName = line[1:-1]
            cntStr = str(cnt)
            codeName = cntStr.zfill(5)

            outLine = ">%s%s\n" % ("CpapProt",codeName)
            #outLine = ">%s%s\n" % ("RcomProt",codeName)
            proteinOut.write(outLine)
            refLine = "%s\t%s\n" % (outLine[1:-1],oldName)
            ref.write(refLine)
            cnt += 1
            
        else:
            proteinOut.write(line)
            

protein.close()
proteinOut.close()
ref.close()

#Import gene files
ref = open(sys.argv[3])
gene = open(sys.argv[4])
geneOut = open(sys.argv[5],"w")

#Use the reference file as a dictionary in order to output new gene file with the codenames
refDict = {}
for line in ref:
    if not line.startswith("#"):
        refLst = line.strip("\n").split("\t")
        #Key should just be the pacid numbers because they link the naming
        #conventions between the protein and gene file
        #key = refLst[1].split(":")[1].strip()
        key = refLst[1].split(" ")[1].strip()
        value = refLst[0]
        refDict[key] = value

mrnaDict = {} #a dict of mRNAs with key: geneName val: [codeName, len]
for line in gene:
    if not line.startswith("#"):
        lineLst  = line.strip().split("\t")
        
        if "mRNA" in line:
            #set the dict
            #name = refDict[line[line.index(":")+1:line.index(";")]]
            name = lineLst[-1].split(";")[0].strip("ID=")
            if name in refDict:
                outLine = "%s\t%s\t%s\t%s\n" % (refDict[name], lineLst[0], small, big)
                
                if geneName in mrnaDict:
                    if geneLen > mrnaDict[geneName][1]:
                        mrnaDict[geneName] = [name, geneLen, outLine]
                else:
                    mrnaDict[geneName] = [name, geneLen, outLine]

        elif "gene" == lineLst[2]:
            small, big = OrderCoords(int(lineLst[3]), int(lineLst[4]))
            geneLen = big - small + 1
            #geneName = lineLst[-1].split("Name=")[1]
            geneName = lineLst[-1].split(";")[0].strip("ID=")

#output info
for group in mrnaDict.values():
    geneOut.write(group[2])

ref.close()
gene.close()
geneOut.close()
