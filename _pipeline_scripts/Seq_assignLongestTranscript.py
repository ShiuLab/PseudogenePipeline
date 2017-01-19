#This script is designed to go through a gff, take in all mRNAs, determine which
#is longest for each gene, assign or reassign the longest transcript and apply
#the change to all other genetic elements
import sys, os

gff = open(sys.argv[1])      #the input gff
out = open(sys.argv[2], "w") #the output modified gff



def OrderCoords(coor1, coor2):
    coor1 = int(coor1);coor2 = int(coor2)
    
    if coor1 < coor2:
        small = coor1
        big = coor2
    elif coor2 < coor1:
        small = coor2
        big = coor1
    else:
        print "***ERROR HERE***"
        return coor1,coor2

    return small, big

def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through the gff file and make a dict of key: geneName val: [(mRNAsize, mRNApacid), ]
gffDict = {}
for line in gff:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[2] == "mRNA":
            small, big = OrderCoords(lineLst[3], lineLst[4])
            mRNAsize = big - small + 1
            mRNApacid = pacID = lineLst[8].split("pacid=")[1].split(";")[0]
            geneName = lineLst[8].split("Parent=")[1].split(";")[0]
            SaveIntoDict(geneName, (mRNAsize, mRNApacid), gffDict)

#Go through dict, sort all value lists (by mRNAsize) and make a notLongest list of mRNApacids
notLongest = []
for geneName in gffDict:
    groups = gffDict[geneName]
    groups.sort()
    if len(groups) > 1:
        lenLst = []
        for group in groups[:-1]: #In the case where the sizes are the same which is called the longest is arbitrary
            notLongest.append(group[1])

#Go through gff again and all genic features without a PACid in the notLongest list are set to longest=1 and the features in the list are set to longest=0.  All info is outputed.
gff.seek(0)
for line in gff:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[2] == "mRNA":
            #for if "longest=" is not currently a part of the gff
            if not "longest=" in lineLst[8]:
                print "\n***ADD FUNCTIONALITY HERE***\n"
            #for if "longest=" is currently a part of the gff
            else:
                mRNApacid = pacID = lineLst[8].split("pacid=")[1].split(";")[0]
                ind = lineLst[8].index("longest=") + 8
                if mRNApacid in notLongest:
                    lineLst[8] = (lineLst[8][:ind] + "0" + lineLst[8][ind+1:])
                else:
                    lineLst[8] = (lineLst[8][:ind] + "1" + lineLst[8][ind+1:])

            #output info
            newLine = "%s\t%s\n" % ("\t".join(lineLst[:8]), lineLst[8])
            out.write(newLine)
        else:
            out.write(line)



gff.close()
out.close()
