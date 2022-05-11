#This script was designed to compare the MCScanX_h outputs from At v. Br and
#At v. Br v. Al v. Rr to see if the groups between At and Br are the same
#in both files (they are expected to be the same)
#Created by David E. Hufnagel on Nov 9, 2012
#WARNING: HIGHLY SITUATION SPECIFIC
"""Algorithm:
1) Go through atBr and make a list of lists of tuples in the format:
   [[(At1, Br2), (At2,Br2)], (At1,Br1),(At23,Br231),(At24,Br7)]
2) Sort the list (only at the top-level)
3) Go through atBrAlRr
   a) Check if at and br are in the line.  If they are proceed.
   b) Make a list of lists of tuples in the format:
   [[(At1, Br2), (At2,Br2)], (At1,Br1),(At23,Br231),(At24,Br7)]
4) Sort the list (only at the top-level)
5) See if the lists are the same
6) Output conclusions"""

import sys

atBr = open(sys.argv[1])      #input MCScanX output .collinearity file with just At and Br
atBrAlRr = open(sys.argv[2])  #input MCScanX output .collinearity file with all species
info = open(sys.argv[3], "w") #output info file
out = open(sys.argv[4], "w")  #output file with all groups in list of list format




def ExtractPair(line):
    gene1 = "";gene2 = "";
    doExtract1 = False #for gene1
    doExtract2 = False #for gene2
    
    for ch in line:
        if ch == "\t":
            if doExtract1 == False:
                #start extraction
                if doExtract2 == False:
                    doExtract1 = True
                #end extraction
                elif doExtract2 == True:
                    doExtract2 == False
            #start extraction of second gene
            elif doExtract1 == True:
                doExtract2 = True
                doExtract1 = False
                
        #do extracting for gene1
        elif doExtract1 == True:
            gene1 += ch
        #do extracting for gene2
        elif doExtract2 == True:
            gene2 += ch

    #strip the junk off of the end
    gene2 = gene2[:-2].strip()
    
    return gene1, gene2

def ExtractPairsFromFile(bigLst, fd, filt=False):
    groupLst = []

    for line in fd.readlines()[9:]: #to skip the MCScanX header            
        if line.startswith("## "):
            if groupLst != []:
                bigLst.append(groupLst)
            groupLst = []
        else:                
            gene1, gene2 = ExtractPair(line)
            if filt == True:
                if (gene1.startswith("AT") or gene1.startswith("Bra")) and \
                   (gene2.startswith("AT") or gene2.startswith("Bra")):
                    doWrite = True
                else:
                    doWrite = False
            else:
                doWrite = True
            if doWrite == True:
                groupLst.append((gene1, gene2))



#Write the users command line prompt on the first line of the output file.
info.write("#python %s\n" % (" ".join(sys.argv)))
out.write("#python %s\n" % (" ".join(sys.argv)))

#1) Go through atBr and make a list of lists of tuples in the format:
#   [[(At1, Br2), (At2,Br2)], (At1,Br1),(At23,Br231),(At24,Br7)]
atBrBigLst = []
ExtractPairsFromFile(atBrBigLst, atBr)

#2) Sort the list (only at the top-level)
atBrBigLst.sort()

#3) Go through atBrAlRr
atBrAlRrBigLst = []
ExtractPairsFromFile(atBrAlRrBigLst, atBrAlRr, True)

#4) Sort the list (only at the top-level)
atBrAlRrBigLst.sort()

#5) See if the lists are the same
print len(atBrBigLst)
print len(atBrAlRrBigLst)

#WARNING: HIGHLY SITUATION SPECIFIC
#6) See how many pairs are in each list
atBrPairCnt = 0
for group in atBrBigLst:
    for pair in group:
        atBrPairCnt += 1

atBrAlRrPairCnt = 0
for group in atBrAlRrBigLst:
    for pair in group:
        atBrAlRrPairCnt += 1

print atBrPairCnt
print atBrAlRrPairCnt

cnt = 0
for group in atBrBigLst:
    if group not in atBrAlRrBigLst:
        cnt += 1

print "dif", cnt


atBr.close()
atBrAlRr.close()
info.close()
out.close()
