#This script is designed to get a 4col protein coordinate file of a new
#fragmented genome
#Created by David E. Hufnagel on Aug 7, 2012
#WARNING LINE 54 IS SITUATION SPECIFIC, "for tup in refDict[oldChromo.strip("At")]:"
import sys

ref = open(sys.argv[1])       #reference file with info about how the old and new contigs relate
old4 = open(sys.argv[2])      #input 4col file for the old contigs
new4 = open(sys.argv[3], "w") #output 4col file for the new contigs




def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)

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




#Write the users command line prompt on the first line of the output file.
new4.write("#python %s\n" % (" ".join(sys.argv)))

#go through ref and import it into a dict of key: oldName val: (newName:start-stop, oldName:start-stop)
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        oldName = lineLst[1].split(":")[0]
        SaveIntoDict(oldName, (lineLst[0], lineLst[1].strip()), refDict)

#print refDict

#go through old4 and translate the coords to the new contigs and output the result
cnt = 0
for line in old4:
    if not cnt % 100:
        print cnt
    if not line.startswith("#"):
        lineLst = line.split("\t")
        oldName = lineLst[0]
        oldChromo = lineLst[1]
        oldProtStart, oldProtStop = OrderCoords(int(lineLst[2]),\
                                                int(lineLst[3].strip()))
        #print oldChromo.strip("At")
        if oldChromo.strip("At") in refDict:
            #print "here"
            for tup in refDict[oldChromo.strip("At")]:
                #print "here2"
                newChromo = tup[0].split(":")[0]
                oldContStart = int(tup[1].split(":")[1].split("-")[0])
                oldContStop = int(tup[1].split(":")[1].split("-")[1])
                newContLen = int(tup[0].split(":")[1].split("-")[1])
                if oldProtStart > oldContStart and oldProtStart < oldContStop:
                    #if the whole thing is in one contig
                    if oldProtStop > oldContStart and oldProtStop < oldContStop:
                        newProtStop = oldProtStop - oldContStart + 1
                    #otherwise
                    else:
                        newProtStop = newContLen
                        
                    newProtStart = oldProtStart - oldContStart + 1
                    newLine = "%s\t%s\t%s\t%s\n" % (oldName, newChromo,\
                                                    newProtStart, newProtStop)
                    new4.write(newLine)
                        
                #if stop is in the range, but start isn't
                elif oldProtStop > oldContStart and oldProtStop < oldContStop:
                    newProtStart = 1
                    newProtStop = oldProtStop - oldContStart + 1
                    newLine = "%s\t%s\t%s\t%s\n" % (oldName, newChromo,\
                                                    newProtStart, newProtStop)
                    new4.write(newLine)
                    
                #if neither start nor stop are in the range, but the new contig is fully encompased by the protein
                elif oldContStart > oldProtStart and oldContStop < oldProtStop:
                    newProtStart = 1
                    newProtStop = newContLen
                    newLine = "%s\t%s\t%s\t%s\n" % (oldName, newChromo,\
                                                    newProtStart, newProtStop)
                    new4.write(newLine)

        else:
            print "%s not in 4col file: " % (oldChromo.strip("At"))
    cnt += 1

        
        





ref.close()
old4.close()
new4.close()
