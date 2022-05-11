#This script is designed to do all things related to handling tabular files
#except name changing (see NameChangeManager.py) like GFF -> 4col conversion or
#replacing code names with real names (I know, I know, I consider it a
#replacement and not a name change).
#Created by David E. Hufnagel on July 17, 2012

import sys, os

#Functions
class Tabular:
    #Replaces a code name in a column with the real name (or vice versa).
    def ReplaceNameWCodeName(self):
        refDict = self.ExtractFromRefFile()
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.split("\t")
                oldName = lineLst[ind]
                newName = refDict[oldName]
                lineLst.pop(ind)
                lineLst.insert(ind, newName)
                newLine = "\t".join(lineLst).strip() + "\n"
                out.write(newLine)
        sys.exit(1)

    #Takes a subset of sequences from a fasta file using the tablular file as
    #the reference (AKA all names in the tabular file are retained)
    def GetFastaFromName(self):
        #Go through the 4col "ref" file and extract names into a list
        nameLst = self.ExtractColFromFileIntoList(ref)

        #Go through the fasta file and output lines where 
        fastaDict = self.ImportFasta(inp)
        for key in fastaDict:
            if key in nameLst:
                val = fastaDict[key]
                out.write(">%s\n" % (key))
                out.write("%s\n" % (val))
        sys.exit(1)

    #Outputs the lines that are in inp, but not in ref. The comparison is done
    #on only one specified column.Outputs the lines that are in inp, but not in
    #ref.  The comparison is done on only one specified column.
    def DiffBetween2FilesAtCol(self):
        #Import the smaller file (ref) into a list
        #itemLst = self.ExtractColFromFileIntoList(ref)
	itemLst = []
        for line in ref:
            if not line.startswith("#"):
                lineLst = line.strip().split("\t")
		itemLst.append(lineLst[ind])

        #Go through the larger file and output lines where it is not in itemLst
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.strip().split("\t")
                if lineLst[ind] not in itemLst:
                    out.write(line)
        sys.exit(1)

    #makes a .size file from a tabular file with coordinates. (col is the col
    #of the first coordinate and the second is assumed to be at col + 1)
    def GetSizes(self):
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.split("\t")
                startCoord, endCoord = self.OrderCoords(int(lineLst[ind]), int(lineLst[ind + 1]))
                size = (endCoord - startCoord + 1)
                newLine = "%s\t%s\n" % (lineLst[nameInd], size) #converted them to int so it would throw an error if the wrong col was entered 
                out.write(newLine)
        sys.exit(1)

    #counts the number of unique elements in a specific column.
    def CntUnique(self):
        theSet = set()
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.split("\t")
                theSet.add(lineLst[ind])

        print "%s unique items in column %s (index %s)" % (len(theSet), col, ind)
        sys.exit(1)

    #Rearranges columns into a new specific order.
    def RearrangeColumns(self):
        #turn number list into python index values (subtract 1)
        cnt = 0
        orderLst = cols #a quick fix
        for num in orderLst:
            orderLst[cnt] = int(num)-1
            cnt += 1

        #go through the input file and do rearrangement
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.split("\t")
        
                tempLst = []
                for num in orderLst:
                    tempLst.append(lineLst[num].strip("\n"))
                newLine = "\t".join(tempLst) + "\n"
                out.write(newLine)

    #GetNotInGff - puts the names of items not in a gff file into an output file
    #with one line per item name. REQ: inp, out, cols, ref (where inp is the
    #tabular file to parse through and ref is the gff file).
    def GetNotInGff(self):

        #import gff file into a dict of key: item name val: tuple(start, stop)
        refDict = {}
        for line in ref:
            if not line.startswith("#"):
                lineLst = line.strip("\n").split("\t")
                refDict[lineLst[1]] = (lineLst[2], lineLst[3])
        
        #Go through inp ind see what isn't in the gff file
        inCnt = 0       #The count for when the inp name, and the coords are all good
        noNameCnt = 0   #The count for when the inp name is not found
        inLst = noNameLst = [] #lists to contain names in these groups
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.strip("\n").split("\t")

                #for each name in inp
                for col in cols:
                    ind = int(col) - 1
                    name = lineLst[ind]

                    #in condition
                    if name in refDict:
                        inLst.append(name)
                        inCnt += 1
                    #noName condition
                    else:
                        noNameLst.append(name)
                        noNameCnt += 1

        #output info
        out.write("##genes from the input file that were found in the gff file (%s)\n" % (inCnt))
        for item in set(inLst):
            out.write("%s\n" % (item))
        out.write("##genes from the input file that were NOT found in the gff file (%s)\n" % (noNameCnt))    
        for item in set(noNameLst):
            out.write("%s\n" % (item))

    #Makes a 4col file, in the format gene, chromo, start, end, from a protein
    #file with coordinates. Includes empty chromos and regions between chromo
    #starts/stops and protein starts/stops.  REQ: inp, out, col, ref (col is the
    #col of the first coordinate and the second is assumed to be at col + 1.
    #(NOT THE INDEX) (ref is a genomic chromosome .size file)
    def GetIntergenic(self):
        #Get chromosomal sizes from ref
        chromoDict = self.ExtractFromRefFile()
        
        #Make a list of start and end coordinates and add 0 to the beginning
        #and chromo end + 1 to the end (0 and end+1 because of a later
        #function action)
        #Note: the first coordStart will be the start of the
        #chromo (0) and the last will be the end of the chromo.
        bigDict = {} #big dict for all chromos (contains coordLsts)
        chromoLst = [] #list for one chromo
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.strip().split("\t")
                small, big = self.OrderCoords(int(lineLst[ind]), int(lineLst[ind+1]))
                if lineLst[1] not in chromoLst:                    
                    #don't end a chromo on the first line (because there's nothing to end)
                    if not chromoLst == []:
                        #end chromo
                        coordLst.append(str(int(chromoDict[chromoLst[-1]])+1))

                        bigDict[chromoLst[-1]] = coordLst
                        coordLst = []
                        
                    #start chromo
                    coordLst = ["0",]
                    coordLst.append(small)
                    coordLst.append(big)
                    chromoLst.append(lineLst[1])

                #continue chromo
                else:
                    coordLst.append(small)
                    coordLst.append(big)

        #end chromo for last chromo
        coordLst.append(str(int(chromoDict[chromoLst[-1]])+1))

        bigDict[chromoLst[-1]] = coordLst
        coordLst = []

        #Perform function that adds 1 to every odd coordinate and subtracts from every even one and puts them into groups of 2.  Ex: [0,4,6,12,16,23] --> [(1,3),(7,11),(17,22)]
        newDict = {}
        for chromo in bigDict:
            newDict[chromo] = self.GetIntergenicSubFunc(bigDict[chromo])

        #Add chromosomes without proteins on them
        for chromo in chromoDict.keys():
            if chromo not in bigDict.keys():
                newDict[chromo] = [("1", chromoDict[chromo]),]
        
        #Output info
        for chromo in newDict:
            cnt = 1
            for pair in newDict[chromo]:
                name = "%s_%s_IntergenicRegion%s" % (spe, chromo, str(cnt).zfill(5))
                newLine = "%s\t%s\t%s\t%s\n" % (name, chromo, pair[0], pair[1])
                out.write(newLine)
                cnt += 1

    #Makes a 4col file, in the format gene, chromo, start, end, from a protein
    #file with coordinates. Without empty chromos and regions between chromo
    #starts/stops and protein starts/stops. REQ: inp, out, col (col is the col
    #of the first coordinate and the second is assumed to be at col + 1. (NOT
    #THE INDEX)
    def GetIntergenicStrict(self):        
        #Make a list of start and end coordinates and add 0 to the beginning
        #and chromo end + 1 to the end (0 and end+1 because of a later
        #function action)
        #Note: the first coordStart will be the start of the
        #chromo (0) and the last will be the end of the chromo.
        bigDict = {} #big dict for all chromos (contains coordLsts)
        chromoLst = [] #list for one chromo
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.strip().split("\t")
                small, big = self.OrderCoords(int(lineLst[ind]), int(lineLst[ind+1]))
                if lineLst[1] not in chromoLst:                    
                    #don't end a chromo on the first line (because there's nothing to end)
                    if not chromoLst == []:
                        #end chromo
                        coordLst = coordLst[:-1] #cut off last value
                        if coordLst != []:
                            bigDict[chromoLst[-1]] = coordLst
                        coordLst = []
                        
                    #start chromo
                    coordLst = [] #don't add first value
                    coordLst.append(big)
                    chromoLst.append(lineLst[1])

                #continue chromo
                else:
                    coordLst.append(small)
                    coordLst.append(big)

        #end chromo for last chromo
        coordLst = coordLst[:-1] #cut off last value
        if coordLst != []:
            bigDict[chromoLst[-1]] = coordLst
        coordLst = []

        #Perform function that adds 1 to every odd coordinate and subtracts from every even one and puts them into groups of 2.  Ex: [0,4,6,12,16,23] --> [(1,3),(7,11),(17,22)]
        newDict = {}
        for chromo in bigDict:
            newDict[chromo] = self.GetIntergenicSubFunc(bigDict[chromo])
        
        #Output info
        for chromo in newDict:
            cnt = 1
            for pair in newDict[chromo]:
                name = "%s_%s_IntergenicRegion%s" % (spe, chromo, str(cnt).zfill(5))
                newLine = "%s\t%s\t%s\t%s\n" % (name, chromo, pair[0], pair[1])
                out.write(newLine)
                cnt += 1

    #Filters a file by a threshold at a specific column.  REQ: inp, out, col, thresh
    def FilterColByThresh(self):
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.split("\t")

                #set doWrite to true if the line passes the threshold
                doWrite = False
                if threshAct == "GE":
                    if float(lineLst[ind]) >= thresh:
                        doWrite = True
                elif threshAct == "G":
                    if float(lineLst[ind]) > thresh:
                        doWrite = True
                elif threshAct == "L":
                    if float(lineLst[ind]) < thresh:
                        doWrite = True
                elif threshAct == "LE":
                    if float(lineLst[ind]) <= thresh:
                        doWrite = True

                #if doWrite was set to true, output line
                if doWrite == True:
                    out.write(line)

        #remove the temp file created
        os.system("rm temp")

    #removes self-self matches in an m8 tabular BLAST output file.
    def RemoveSelfBlast(self):
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.split("\t")
                if not lineLst[0] == lineLst[1]:
                    out.write(line)

    #AddSpeToCol - Concatenates a string to the start of a col
    def AddStrToCol(self):
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.split("\t")
                lineLst[ind] = string + lineLst[ind]
                newLine = "\t".join(lineLst)
                out.write(newLine)

    #takes coordinates in a file and orders them from small to big
    def OrderCoordsF(self):
        print "in"
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.strip().split("\t")
                small, big = self.OrderCoords(int(lineLst[ind]), int(lineLst[ind + 1]))
                lineLst[ind] = str(small)
                lineLst[ind + 1] = str(big)
                newLine = "\t".join(lineLst) + "\n"
                out.write(newLine)

    #removes lines containing keywords at a given column
    def RemoveByKeyword(self):
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.split("\t")
                if not string in lineLst[ind]:
                    out.write(line)

    #takes a tabular file and concatenates features that overlap with each other
    #and separates the names with a '_'.  Does not retain features other than
    #name, chromo and coords. 
    def RemoveOverlap(self):
        #Make a list of dicts of key: name val: (startCoor, endCoor, chromoName) with all features for each chromo
        chromoLst = []        #contains nameDicts
        chromoNameSet = set() #contains chromoNames
        nameDict = {}
        lastLine = ""
        lineLen = 0           #keeps track of the length of the lines for outputting later
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.strip().split("\t")
                lineLen = len(lineLst)
                #Start a new chromo in the list
                if not lineLst[chroInd] in chromoNameSet:
                    chromoLst.append(nameDict)
                    #print lastLine[0]
                    chromoNameSet.add(lineLst[chroInd])
                    nameDict = {}
                    
                #Add line info to nameDict
                if not lineLst[nameInd] in nameDict:
                    nameDict[lineLst[nameInd]] = (lineLst[ind], lineLst[ind+1], lineLst[chroInd])
                else:
                    print "/n/n/n***Major Error in RemoveOverlap!***/n/n/n"
                lastLine = lineLst
                
        #add the last chromosome data
        else:
            chromoLst.append((nameDict))
            chromoNameSet.add(lineLst[chroInd])
            nameDict = {}

        #remove the empty dict from the chromoLst
        chromoLst = chromoLst[1:]

        #Go through list of dicts and make a list of tuples of coors for each chromo with gene info
        overlapCnt = 0
        for chromo in chromoLst:
            coorLst = []
            #import coor, name and chromo data into a list
            for gene in chromo:
                small, big = self.OrderCoords(chromo[gene][0],chromo[gene][1])
                coorLst.append((small, big, gene, chromo[gene][2]))
                
            #Sort the list
            coorLst.sort()

            #Go through the list, handle overlap and make output
            lastCoor = ""
            lastName = ""
            lastChromo = ""
            currCoor = ""
            for currInfo in coorLst:
                doSet = True
                doWrite = True
                currName = currInfo[2]
                currChromo = currInfo[3]
                #If the start and stop coordinates for currCoor are the same, set a flag to not set currCoor to lastCoor
                if currInfo[0] == currInfo[1]:
                    doSet = False
                    doWrite = False
                    print "ERROR: make sure this works! in RemoveOverlap TabularManager.py"
                
                #Look for overlap between currCoor and lastCoor, if found, bring
                #currCoor and lastCoor together and call it currCoor (soon to be lastCoor)
                if lastCoor != "" and int(currInfo[0]) >= int(lastCoor[0]) and int(currInfo[0]) <= int(lastCoor[1]):
                        #print "overlap"
                        doWrite = False
                        currName = "%s_%s" % (lastName, currName)
                        #print lastCoor
                        all4Coor = [int(lastCoor[0]), int(lastCoor[1]), int(currInfo[0]), int(currInfo[1])]
                        currCoor = (min(all4Coor), max(all4Coor), currName, currChromo)
                        overlapCnt += 1
                else:
                    currCoor = currInfo



                #if doWrite == True output info
                if lastCoor != "" and doWrite == True:
                    #print "write"
                    #print lastName
                    #print currName
                    #print lastChromo
                    #print currChromo
                    #print lastCoor
                    #print currCoor
                    #print
                    #currName = lastCoor[2]
                    #currChromo = lastCoor[3]
                    self.OutputInfo(lastName, lastChromo, lastCoor, lineLen)

                #if doSet == True set currCoor to lastCoor
                if doSet == True:
                    #print "set"
                    lastCoor = currCoor
                    lastName = currName
                    lastChromo = currChromo

            #Print the last line of last line info
            else:
                self.OutputInfo(lastName, lastChromo, lastCoor, lineLen)
        print "\n%s instances of overlap found and fixed" % (overlapCnt)
                    
    #Used in RemoveOvelap
    def OutputInfo(self, name, chromo, coords, lineLen):
        #Make an empty list of the proper length
        outLst = []
        for num in range(lineLen):
            outLst.append("-")

        #Replace name, chromo and coords in the right spots
        outLst[nameInd] = name
        outLst[chroInd] = chromo
        outLst[ind] = str(coords[0])
        outLst[col] = str(coords[1])

        newLine = "\t".join(outLst) + "\n"
        out.write(newLine)

    #Adds up all values in a column and prints the result
    def GetColSum(self):
        #add up the numbers
        numLst = []
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.strip().split("\t")
                numLst.append(float(lineLst[ind]))

        #make the output pretty by adding commas
        preDec = str(sum(numLst)).split(".")[0]
        postDec = str(sum(numLst)).split(".")[1]
        modSum = ""
        cnt = 1
        for x in preDec[::-1]:
            modSum += x
            if not cnt % 3:
                modSum += ","
            cnt += 1
        modSum = modSum[::-1].strip(",") + "." + postDec

        #output info
        newLine = "Sum: %s" % (modSum)
        print newLine

    #gets the median of all values in a column and prints the result
    def GetColMedian(self):
        #put the numbers into a list
        numLst = []
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.strip().split("\t")
                numLst.append(float(lineLst[ind]))

        #get the median of the numbers
        median = self.Median(numLst)
        
        print median

    #filters blast by removing all but the best match to the query.
    #NOT RECIPROCAL BM!
    def BestMatchBlast(self):
        queryLst = []
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.split("\t")
                query = lineLst[0]
                if query not in queryLst:
                    out.write(line)
                    queryLst.append(query)

    #makes a 4col from a file with names in it and a 4col where the
    #first file is a subset of the 4col."
    def Make4col(self):
        #go through the 4col file and make a dict of key: name val: line
        fourDict = {}
        for line in ref:
            if not line.startswith("#"):
                lineLst = line.strip().split("\t")
                fourDict[lineLst[0]] = line

        #go through the input file and output all 4col lines related to the
        #names in fourDict
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.strip().split("\t")
                newLine = fourDict[lineLst[ind]]
                out.write(newLine)
                

    def ExtractFromRefFile(self):
        refDict = {}
        for line in ref:
            if not line.startswith("#"):
                lineLst = line.split("\t")
                if FR == "F":
                    refDict[lineLst[0]] = lineLst[1].strip()
                elif FR == "R":
                    refDict[lineLst[1].strip()] = lineLst[0]
                else:
                    print "\n***Error in ExtractFromRefFile!***\n"

        return refDict

    def ExtractColFromFileIntoList(self, fd):  #fd = file directory
        itemLst = []
        for line in fd:
            if not line.startswith("#"):
                lineLst = line.strip().split("\t")
                itemLst.append(lineLst[ind])

        return itemLst

    def OrderCoords(self, coor1, coor2):
        coor1 = int(coor1);coor2 = int(coor2)

        small = coor1; big = coor2
        if coor1 < coor2:
            small = coor1
            big = coor2
        elif coor2 < coor1:
            small = coor2
            big = coor1
        else:
            small = coor1; big = coor2
            print "***ERROR HERE***"

        return small, big

    def ImportFasta(self, fasta):  #fasta is the fasta file name
        fastaDict = {}
        currName = ""   #the current name associated with a seq
        seq = ""        #the current seq to be built up for each seq line
        for line in fasta:
            if not line.startswith("#"):
                if line.startswith(">"):
                    if currName != "":
                        fastaDict[currName] = seq
                    seq = ""
                    currName = line.strip().strip(">")
                else:
                    seq += line.strip()
                    
        #get the last seq on the way out
        else:
            fastaDict[currName] = seq

        return fastaDict

    def GetIntergenicSubFunc(self, coordLst):
        cnt = 1
        chroLst = []
        pairLst = []
        for num in coordLst:
            if cnt % 2 == 1:
                num = int(num) + 1
                pairLst.append(num)
            else:
                num = int(num) - 1
                pairLst.append(num)
                chroLst.append(pairLst)
                pairLst = []
            cnt += 1

        return chroLst

    #gets the median of a list
    def Median(self, listx):
        #Make all items floats and sort the list
        newLst = []
        for item in listx:
            newLst.append(float(item))
        newLst.sort()

        #Determine if the length is an odd number
        if len(newLst) % 2:
            isOdd = True
        else:
            isOdd = False

        #Find the actual median
        if isOdd == True:
            print 1
            median = float(newLst[len(newLst)/2])
        else:
            leftMid = newLst[(len(newLst)/2)-1]
            rightMid = newLst[len(newLst)/2]
            median = float(leftMid + rightMid) / 2.0

        return median

    #Writes the users command line prompt on the first line of the output file.
    def WriteCallLine(self):
        out.write("#python %s\n" % (" ".join(sys.argv)))
        
    def Help(self):
        print "\nDescription:"
        print "This script is designed to do all things related to handling tabular files"
        print "except name changing (see NameChangeManager.py) like GFF -> 4col conversion or"
        print "replacing codes with real names (I know, I know, I consider it a"
        print "replacement and not a name change)."
        print "Created by David E. Hufnagel on July 17, 2012"
        print "    -f function:"
        print "        RemCodeName - Replaces a code name in a column with the real name"
        print "            (or vice versa). REQ: inp, out, col, ref, FR"
        print "        GetFastaFromName - Takes a subset of sequences from a"
        print "            fasta file using the tablular file as the reference."
        print "            (AKA all names in the tabular file are retained)"
        print "            REQ: inp, out, col, ref (where inp is the fasta file"
        print "                 and ref is the 4col file)"
        print "        Diff - Outputs the lines that are in inp, but not in ref."
        print "               The comparison is done on only one specified column."
        print "               REQ: inp, out, col, ref (where inp is the bigger file"
        print "                    and ref is the smaller file)"
        print "        GetSizes - makes a .size file from a tabular file with"
        print "             coordinates.  REQ: inp, out, name, col (col is the"
        print "             col of the first coordinate and the second is"
        print "             assumed to be at col + 1. (NOT THE INDEX)  Name is the"
        print "             col of the name to be used for each col"
        print "        CntUnique - counts the number of unique elements in a"
        print "             specific column. REQ: inp, col"
        print "        ShiftCols - Rearranges columns into a new"
        print "             specific order. REQ: inp, out, cols"
        print "        GetNotInGff - puts the names of items not in a gff file"
        print "             into an output file with one line per item name."
        print "             REQ: inp, out, cols, ref (where inp is the tabular"
        print "             file to parse through and ref is the gff file)"
        print "        GetIntergenic - Makes a 4col file, in the format gene,"
        print "             chromo, start, end, from a protein file with"
        print "             coordinates. Includes empty chromos and regions between"
        print "             chromo starts/stops and protein starts/stops. REQ: inp,"
        print "             out, col, ref, spe (col is the"
        print "             col of the first coordinate and the second is"
        print "             assumed to be at col + 1. (NOT THE INDEX) (ref is"
        print "             a genomic chromosome .size file)"
        print "        GetIntergenicStrict - Makes a 4col file, in the format gene,"
        print "             chromo, start, end, from a protein file with"
        print "             coordinates. Without empty chromos and regions between"
        print "             chromo starts/stops and protein starts/stops. REQ: inp,"
        print "             out, col, spe (col is the"
        print "             col of the first coordinate and the second is"
        print "             assumed to be at col + 1. (NOT THE INDEX)"
        print "        FilterColByThresh - Filters a file by a threshold at a"
        print "             specific column.  REQ: inp, out, col, thresh, threshAct"
        print "        RemoveSelfBlast - removes self-self matches in an m8 tabular BLAST"
        print "             output file. REQ: inp, out"
        print "        BestMatchBlast - filters blast by removing all but the"
        print "             best match to the query. NOT RECIPROCAL BM! REQ: inp, out"
        print "        AddStrToCol - Concatenates a string to the"
        print "             start of a col REQ: inp, out, col, string"
        print "        OrderCoords - takes coordinates and orders them from small"
        print "             to big REQ: inp, out, col (col is the"
        print "             col of the first coordinate and the second is"
        print "             assumed to be at col + 1. (NOT THE INDEX))"
        print "        RemoveByKeyword - removes lines containing keywords at a given"
        print "             column. REQ: inp, out, col, string (string is the keyword)"
        print "        RemoveOverlap - takes a tabular file and concatenates features that"
        print "             overlap with each other and separates the names with a '_'"
        print "             Does not retain features other than name, chromo and coords."
        print "             REQ: inp, out, col, name, chroCol"
        print "        GetColSum - adds up all values in a column and prints the result"
        print "             REQ: inp, col (col is in number form and NOT AN INDEX)"
        print "        GetColMedian - gets the median of all values in a column and prints the result"
        print "             REQ: inp, col (col is in number form and NOT AN INDEX)"
        print "        RemoveNotInFasta - removes lines from a tabular file that don't"
        print "             match fasta names REQ: inp, ref, out, col (inp is the tabular"
        print "             file and ref is the fasta file)"
        print "        Make4col - makes a 4col from a file with names in it and"
        print "             a 4col where the first file is a subset of the 4col."
        print "             REQ: inp, ref, out, col (where inp is the file with names,"
        print "             col is the column where the names can be found and ref is the"
        print "             greater 4col file"
        print "    -inp - the input file to be modified."
        print "    -out - the modified output file."
        print "    -col - the column to work with (number form, NOT AN INDEX)"
        print "    -cols - the columns to work with.  Comma seperated."
        print "        (number form, NOT AN INDEX)"
        print "    -ref - a reference file for remCodeName"
        print "    -FR - Forward or Reverse (for remCodeName: F = names matching the 1st col"
        print "        of the reference file are to be replaced with names matching the"
        print "        second. R = vice versa)"
        print "    -name - the column of a name (number form, NOT AN INDEX)"
        print "    -spe - A 2 character species code"
        print "    -thresh - A threshold for filtering"
        print "    -threshAct - the action to take for the threshold options:"
        print "        G,GE,L,LE which stand for greater than, greater or equal"
        print "        to, lesser than and lesser than or equal to respecively."
        print "        in other words, what to keep."
        print "    -string - a general purpose string"
        print "    -chroCol - the column where the chromosome name is contained"

        sys.exit(0)
        


FRList = ["F","f","R","r"]
TAoptions = ["G", "GE", "L", "LE"]
 
inp= FR= cols= spe= threshAct= string= ""; col= chroCol = 0; thresh= 0.0; ind= name= -1; out = "temp"
TB = Tabular()

#input handling
isEmpty = True #will be set to false if there are any parameters
out = open("temp", "w") #in case it's not changed later we need a default to not create errors
for i in range(1,len(sys.argv),2):
    if sys.argv[i] == "-f":
        function = sys.argv[i+1]
    elif sys.argv[i] == "-inp":
        inp = open(sys.argv[i+1])
    elif sys.argv[i] == "-out":
        out = open(sys.argv[i+1], "w")
    elif sys.argv[i] == "-col":
        col = int(sys.argv[i+1])
        ind = col - 1 #the index for use
    elif sys.argv[i] == "-ref":
        ref = open(sys.argv[i+1])
    elif sys.argv[i] == "-FR":
        FR = sys.argv[i+1]
        if FR not in FRList:
            print "error: improper FR!:",sys.argv[i+1]
            print "add -h to get help"
            sys.exit(0)
        elif FR == "f":
            FR == "F"
        elif FR == "r":
            FR == "R"
    elif sys.argv[i] == "-name":
        name = int(sys.argv[i+1])
        nameInd = name - 1
    elif sys.argv[i] == "-chroCol":
        chroCol = int(sys.argv[i+1])
        chroInd = chroCol - 1
    elif sys.argv[i] == "-cols":
        cols = sys.argv[i+1].split(",")
        if 0 in cols:
            print "error: improper cols!:", sys.argv[i+1]
            print "add -h to get help"
            sys.exit(0)
    elif sys.argv[i] == "-spe":
        spe = sys.argv[i+1]
        if len(spe) != 2:
            print "error: improper spe!:", sys.argv[i+1]
            print "add -h to get help"
            sys.exit(0)
    elif sys.argv[i] == "-thresh":
        thresh = float(sys.argv[i+1])
    elif sys.argv[i] == "-threshAct":
        threshAct = sys.argv[i+1]
        if threshAct not in TAoptions:
            print "error: improper threshAct!:", sys.argv[i+1]
            print "add -h to get help"
            sys.exit(0)
    elif sys.argv[i] == "-string":
        string = sys.argv[i+1]
    elif sys.argv[i] == "-h":
        TB.Help()
    else:
        print "UNKNOWN FLAG:",sys.argv[i]
        print "add -h to get help."
        sys.exit(0)
    isEmpty = False

#no parameters will bring up the help text
if isEmpty == True:
    TB.Help()

#calling funcitons
#functions without an output file
if function == "CntUnique" or function == "cntUnique":
    if inp == "" or col == 0:
        print "\n***missing parameter(s) REQ: inp, col***\n"
        TB.Help()
    TB.CntUnique()

#functions with an output file
TB.WriteCallLine()
if function == "RemCodeName" or function == "remCodeName":
    if inp == "" or out == "temp" or col == 0 or ref == "" or FR == "":
        print "\n***missing parameter(s) REQ: inp, out, col, ref, FR***\n"
        TB.Help()
    TB.ReplaceNameWCodeName()
elif function == "GetFastaFromName" or function == "getFastaFromName":
    if inp == "" or out == "temp" or col == 0 or ref == "":
        print "\n***missing parameter(s) REQ: inp, out, col, ref***\n"
        TB.Help()
    TB.GetFastaFromName()
elif function == "Diff" or function == "diff":
    if inp == "" or out == "temp" or col == 0 or ref == "":
        print "\n***missing parameter(s) REQ: inp, out, col, ref***\n"
        TB.Help()
    TB.DiffBetween2FilesAtCol()
elif function == "GetSizes" or function == "Getsizes" or function == "getSizes" or function == "getsizes":
    if inp == "" or out == "temp" or name == -1 or col == 0:
        print "\n***missing parameter(s) REQ: inp, out, name, col***\n"
        TB.Help()
    TB.GetSizes()
elif function == "ShiftCols" or function == "shiftCols":
    if inp == "" or out == "temp" or cols == "":
        print "\n***missing parameter(s) REQ: inp, out, cols***\n"
        TB.Help()
    TB.RearrangeColumns()
elif function == "GetNotInGff" or function == "getNotInGff":
    if inp == "" or out == "temp" or cols == "" or ref == "":
        print "\n***missing parameter(s) REQ: inp, out, cols, ref***\n"
        TB.Help()
    TB.GetNotInGff()
elif function == "GetIntergenic" or function == "getIntergenic":
    if inp == "" or out == "temp" or ref == "" or col == 0:
        print "\n***missing parameter(s) REQ: inp, out, ref, col***\n"
        TB.Help()
    FR = "F" #see ExtractFromRefFile()
    TB.GetIntergenic()
elif function == "GetIntergenicStrict" or function == "getIntergenicStrict":
    if inp == "" or out == "temp" or col == 0:
        print "\n***missing parameter(s) REQ: inp, out, col***\n"
        TB.Help()
    FR = "F" #see ExtractFromRefFile()
    TB.GetIntergenicStrict()
elif function == "FilterColByThresh" or function == "filterColByThresh":
    if inp == "" or out == "temp" or col == 0 or thresh == "" or threshAct == "":
        print "\n***missing parameter(s) REQ: inp, out, col, thresh, threshAct***\n"
        TB.Help()
    TB.FilterColByThresh()
elif function == "RemoveSelfBlast" or function == "removeSelfBlast":
    if inp == "" or out == "temp":
        print "\n***missing parameter(s) REQ: inp, out***\n"
        TB.Help()
    TB.RemoveSelfBlast()
elif function == "AddStrToCol" or function == "addStrToCol":
    if inp == "" or out == "temp" or col == 0 or string == "":
        print "\n***missing parameter(s) REQ: inp, out, col, string***\n"
        TB.Help()
    TB.AddStrToCol()
elif function == "OrderCoords" or function == "orderCoords":
    if inp == "" or out == "temp" or col == 0:
        print "\n***missing parameter(s) REQ: inp, out, col***\n"
        TB.Help()
    TB.OrderCoordsF()
elif function == "RemoveByKeyword" or function == "removeByKeyword":
    if inp == "" or out == "temp" or col == 0 or string == "":
        print "\n***missing parameter(s) REQ: inp, out, col, string***\n"
        TB.Help()
    TB.RemoveByKeyword()
elif function == "RemoveOverlap" or function == "removeOverlap":
    if inp == "" or out == "temp" or col == 0 or name == "" or chroCol == 0:
        print "\n***missing parameter(s) REQ: inp, out, col, name,chroCol***\n"
        TB.Help()
    TB.RemoveOverlap()
elif function == "GetColSum" or function == "getColSum":
    if inp == "" or col == 0:
        print "\n***missing parameter(s) REQ: inp, col ***\n"
        TB.Help()
    TB.GetColSum()
elif function == "GetColMedian" or function == "getColMedian":
    if inp == "" or col == 0:
        print "\n***missing parameter(s) REQ: inp, col ***\n"
        TB.Help()
    TB.GetColMedian()
elif function == "BestMatchBlast" or function == "bestMatchBlast":
    if inp == "" or out == "temp":
        print "\n***missing parameter(s) REQ: inp, out***\n"
        TB.Help()
    TB.BestMatchBlast()
elif function == "Make4col" or function == "make4col":
    if inp == "" or out == "temp" or col == 0 or ref == "":
        print "\n***missing parameter(s) REQ: inp, out, col, ref***\n"
        TB.Help()
    TB.Make4col()
else:
    print "\n***unknown function: need more functionality***\n"
