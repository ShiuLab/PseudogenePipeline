#This script is step 2 in linking Rr-Br orthologs to Br-Al-At orthologs file.
#This script will take the .comp file made in part 1 and the two column Rr-Br
#file to create the "pairwise" file

import sys

two = open(sys.argv[1])        #input file with Rr and Br
comp = open(sys.argv[2])       #comprehensive output file
pairw =open(sys.argv[3], "w")  #pairwise internal output file
pairw2 =open(sys.argv[4], "w") #pairwise final output file

#writes the users command line prompt on the first line of the output file.
pairw.write('#python %s\n'%(' '.join(sys.argv)))





def Expand(Rr, Br, Al, At):
    for r in Rr:
        for b in Br:
            for l in Al:
                for t in At:
                    newline = "%s\t%s\t%s\t%s\n"  % (r,b,l,t)
                    pairw.write(newline)

def ExpandFriends(Rr, Br, Al, At, friends, badLst):
    good = True     #says whether a line should should be written
    for b in Br:
        for r in Rr:
            #for when b is with wrong r or vice versa
            if (b,r) not in friends:
                good = False
            for l in Al:
                for t in At:
                    if good == True:
                        newL = "%s\t%s\t%s\t%s\n"  % (r,b,l,t)
                        #print newL
                        pairw.write(newL)
        #for when b has no r friend
        if good == False:
            for l in Al:
                for t in At:
                    newln = "-\t%s\t%s\t%s\n" % (b,l,t)
                    #print "***", newln
                    pairw.write(newln)
                    badLst.append(b)
    return badLst #badLst is the list of Br names where there was no Rr ortholog
            
                    

def SaveIntoDict(gene1, gene2, dictX):
    if gene1 not in dictX:
        dictX[gene1] = [gene2]
    else:
        dictX[gene1].append(gene2)




                    
#WARNING: sequence specific (to get rid of hash line)
comp.readline()

#go through two and extract data into Dict
twoDict = {}
for row in two:
    rowLst = row.split("\t")
    SaveIntoDict(rowLst[0], rowLst[1][:-1], twoDict)
#convert dict into friends list of tuples
friends = []
for k in twoDict:
    for v in twoDict[k]:
        friends.append((k,v))

#go through comp, see if there are two+ Rr's, if there are also two+ Br's, and
#if so write output accordingly (RrA-BrA, RrB-BrB only).  Otherwise, write
#output with max lines (RrA-BrA, RrA-BrB, RrA-BrC)
badLst = []
for line in comp:
    lineLst = line.split("\t")
    radLst = eval(lineLst[3])
    braLst = eval(lineLst[2])
    lyrLst = eval(lineLst[1])
    tha = lineLst[0]
    #if there are two+ Rr's and if there are also two+ Br's
    if len(radLst) > 1 and len(braLst) > 1:
        badLst = ExpandFriends(radLst, braLst, lyrLst, [tha,], friends, badLst)
        
            

    #Otherwise, write output with max lines (RrA-BrA, RrA-BrB, RrA-BrC)
    else:
        Expand(radLst, braLst, lyrLst, [tha,])

        
#Go through pairw and filter it to remove lines with "-" designated as the Rr
#ortholog when there is another line with a real Rr gene designated as the
#ortholog.  Here's how it works:  "badLst" was made already with the Br names
#where there is no Rr ortholog.  Go through the file and make a "badLst2"
#which is the Br names from "badLst" with a different line containing an Rr
#ortholog, then go through pairw again and remove lines containing a Br name
#from badLst2 AND a - line for Rr
pairw.close()
pairw = open(sys.argv[3])

#WARNING: sequence specific (to get rid of hash line)
pairw.readline()

badLst2 = []
for L in pairw:
    Llst = L.split("\t")
    if Llst[1] in badLst and Llst[0] != "-":
        #print Llst
        badLst2.append(Llst[1])


pairw.close()
pairw = open(sys.argv[3])

#WARNING: sequence specific (to get rid of hash line)
pairw.readline()

#writes the users command line prompt on the first line of the output file.
pairw2.write('#python %s\n'%(' '.join(sys.argv)))

for lin in pairw:
    linLst = lin.split("\t")
    if not (linLst[1] in badLst2 and linLst[0] == "-"):
        pairw2.write(lin)
        





two.close()
comp.close()
pairw.close()

