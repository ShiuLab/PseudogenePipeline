#This script is designed to take e-value files and make a data table from them
#Input format:
#forward       reverse comp.      e-val
#ATGTCGGC        GCCGACAT        0.49924\r  (that's right, not "\n"!)
#Output format:
#        TF1   TF2   TF3
# site1 e-val   -   e-val
# site2 e-val e-val   -
#  ...
#Created by David E. Hufnagel on Jan 15, 2013

import sys

inpNames = sys.argv[1]       #input e-value files.  comma seperated
out = open(sys.argv[2], "w") #output data table

#go through input files, open all of them and output info lines
out.write("#python %s\n" % (" ".join(sys.argv)))
out.write("#motif\trev_comp")
inpFiles = []
for name in inpNames.split(","):
    out.write("\t%s" % (name))
    inpFiles.append(open(name))



#get the reverse compliment of a sequence
def GetRevComp(motif):
    compDict = {"G":"C", "C":"G", "T":"A", "A":"T", ".":"."} #the dictionary of bases and their compliments
    
    revComp = []
    for ch in motif[::-1]:
        revComp.append(compDict[ch])
    revComp = "".join(revComp)
    
    return revComp



#go through opened input files and make a dict of key: site val: (TF1,TF2...).  It is assumed that each file is a different transcription factor
tableDict = {} #the big dict containing all input to be outputted
cnt = 1        #a count to keep track of how many items should be in the val after processing in the big dict for each site (either an e-val or a "-" should be added each time)
for filE in inpFiles:
    #add this TF to sites that are present in this file
    for chunk in filE:
        if not chunk.startswith("#"):
            lines = chunk.strip().split("\r")
            for line in lines:
                lineLst = line.strip().split("\t")
                motif = lineLst[0]
                revComp = lineLst[1]
                num = lineLst[2] #either e-Val or intensity
                
                if motif in tableDict or revComp in tableDict:
                    tableDict[motif].append(num)
                else:
                    #when it's showing up for the first time not in the first file
                    if cnt != 1:
                        tempCnt = 2
                        tableDict[motif] = ["-",]
                        while tempCnt < cnt:
                            tableDict[motif].append("-")
                            tempCnt += 1
                        tableDict[motif].append(num)
                    else:
                        tableDict[motif] = [num,]

    #add hashes to motifs that weren't present in this file
    for key in tableDict:
        if len(tableDict[key]) < cnt:
            tableDict[key].append("-")
    cnt += 1              

#output info
for motif in tableDict:
    out.write("\n%s\t%s" % (motif, GetRevComp(motif)))
    for tf in tableDict[motif]:  #tf = transcription factor
        out.write("\t%s" % (tf))




#close files
for filE in inpFiles:
    filE.close()
out.close()
