#Prepares files (for the DDF1/DDF2 project) for Shinhan's drawBlocks script in
#the format:
"""
organism  Cp12SynBlock
chr  Cpapsupercontig_74
start  688014
end  872397
gene  start  stop  ori  name(ori is either 1 or -1)
gene..."""
#Created by David E. Hufnagel on Jun 18, 2013
import sys

block = open(sys.argv[1])    #input 4col of syntenic block (yes, singular)
four = open(sys.argv[2])     #input 4col of pseudogenes
dis = open(sys.argv[3])      #input pseudogene disable_count file
ref = open(sys.argv[4])      #input reference pseudogene file
gff = open(sys.argv[5])      #input gff file for gene info
out = open(sys.argv[6], "w") #output info file for use by Shinhan's drawBlocks



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

#Go through block to get header info
cnt = 0
chromo = ""
blockStart = ""
blockStop = ""
for line in block:
    if not line.startswith("#"):
        #gather data
        lineLst = line.strip().split("\t")
        organism = lineLst[0]
        chromo = lineLst[1]
        blockStart = lineLst[2]
        blockStop = lineLst[3]

        #throw an error if the file is over 1 line long
        if cnt > 0:
            print "too many lines"
            sys.exit()

        #output the header info
        out.write("organism\t%s\n" % (organism))
        out.write("chr\t%s\n" % (chromo))
        out.write("start\t%s\n" % (blockStart))
        out.write("end\t%s\n" % (blockStop))

        #increase the count
        cnt += 1

#Go through ref and make a dict of key: bigName val: codeName
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        refDict[lineLst[1]] = lineLst[0]
        
#Go through dis and make dict of key: name val: ori (1 or -1)
oriDict = {}
for line in dis:
    if not line.startswith("#python") and line.startswith("#"):
        lineLst = line.strip().split(" ")

        #determine orientation
        tempStart = lineLst[2].split(":")[1].split("-")[0]
        tempStop = lineLst[2].split(":")[1].split("-")[1]
        small, big = OrderCoords(tempStart, tempStop)
        if tempStart == str(small):
            ori = "1"
        else:
            ori = "-1"
            
        codeName = (";".join(lineLst[:3]) + ";" + "|".join(lineLst[4:]))[1:]
        psName = refDict[codeName]
        oriDict[psName] = ori

#Go through gff and output gene info for genes within the syntenic block coordinates
for line in gff:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[2] == "gene" and lineLst[0] == chromo:
            #gather data
            gnStart = lineLst[3]
            gnStop = lineLst[4]
            #WARNING: HIGHLY SITUATION SPECIFIC
            if lineLst[0].startswith("Atha"):
                gnName = lineLst[8].split("Name=")[1].split(";")[0]
            elif lineLst[0].startswith("Alyr"):
                gnName = "Alyr" + lineLst[8].split("Name=")[1].split(";")[0]
            elif lineLst[0].startswith("Thal"):
                gnName = lineLst[8].split("Name=")[1].split(";")[0][:-3]
            elif lineLst[0].startswith("Brap"):
                gnName = lineLst[8].split("Name=")[1].split(";")[0]
            elif lineLst[0].startswith("Crub"):
                gnName = lineLst[8].split("Name=")[1].split(";")[0][:-3]
            elif lineLst[0].startswith("Cpap"):
                gnName = lineLst[8].split("Name=")[1].split(";")[0].replace(".","_")
            else:
                print "need more species functionality"
                sys.exit()
            ori = lineLst[6]
            #modify ori language
            if ori == "+":
                ori = "1"
            elif ori == "-":
                ori = "-1"
            else:
                print "ori error in Z_4_prepDataForDrawBlocks.py"

            #consider only data within coordinates and output info
            if int(gnStart) >= int(blockStart) and int(gnStop) <= int(blockStop):
                newLine = "gene\t%s\t%s\t%s\t%s\n" % (gnStart, gnStop, ori, gnName)
                out.write(newLine)

        
#Go through four and output pseudogene info for pseudogenes within the syntenic block coordinates
for line in four:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")

        #gather info
        psStart = lineLst[2]
        psStop = lineLst[3]
        psName = lineLst[0]
        ori = oriDict[psName]

        #consider only data within coordinates and output info
        if int(psStart) >= int(blockStart) and int(psStop) <= int(blockStop):
            newLine = "gene\t%s\t%s\t%s\t%s\n" % (psStart, psStop, ori, psName)
            out.write(newLine)
        







block.close()
four.close()
dis.close()
gff.close()
out.close()
