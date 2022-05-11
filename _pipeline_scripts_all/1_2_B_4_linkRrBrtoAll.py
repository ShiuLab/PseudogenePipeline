#This script is designed to link Rr-Br orthologs to a file with orthologs and
#paralogs for Br-At-Al and to ouput the result in a comprehensive format.
#note:  it is intended that after running this script, R_2_B_4_linkRrBrtoAll2.py
#is ran to create the "pairwise" file
#Created by David E. Hufnagel on 5-8-2012

import sys

two = open(sys.argv[1])       #input file with Rr and Br
three = open(sys.argv[2])     #input file with Br, At and Al
comp = open(sys.argv[3], "w") #comprehensive output file

#writes the users command line prompt on the first line of the output file.
comp.write('#python %s\n'%(' '.join(sys.argv)))



def SaveIntoDict(gene1, gene2, dictX):
    if gene1 not in dictX:
        dictX[gene1] = [gene2]
    else:
        dictX[gene1].append(gene2)



#WARNING: sequence specific (to get rid of hash line)
three.readline()

#go through Rr-Br file and extract pairs
twoDict = {}
for line in two:
    lineLst = line.split("\t")
    SaveIntoDict(lineLst[0], lineLst[1][:-1], twoDict)

#go through Br-At-Al file and build comp (a Br-At-Al-Rr file)
cnt1 = 0
cnt2 = 0
done = False
doneLst = [] #the list of Rr-Br orthos incorperated (just Br names)
for ln in  three:
    lnLst = ln.split("\t")
    brasLst = eval(lnLst[2]) #brassica rapa list from comp
    for p in twoDict:
        cnt1 += 1
        if not cnt1 % 1000000:
            print cnt1
        if p in brasLst:
            newline = "%s\t%s\t%s\t%s\n" % \
                      (lnLst[0], lnLst[1], lnLst[2][:-1], twoDict[p])
            comp.write(newline)
            done = True
                            #a tuple of (key, [value1, value2])
            doneLst.append(p)
    if done == False:
        newline = "%s\t%s\t%s\t[]\n" % (lnLst[0], lnLst[1], lnLst[2][:-1])
        comp.write(newline)
        done = True
    done = False

#Go through doneLst and twoDict and where the two don't overlap write
#into .comp file in the following format: [] [] [Br2123] [Rr234332]
for key in twoDict:
    if key not in doneLst:
        newline = "[]\t[]\t%s\t%s\n" % ([key], twoDict[key])
        comp.write(newline)
            
        
    

two.close()
three.close()
comp.close()
