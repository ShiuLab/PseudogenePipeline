#This script was designed to do two things 1) condense lines with Br genes with
#multiple dup times and condense lines with multiple outgroups be getting
#averages.
#Created by David E. Hufnagel on Dec 14, 2012

import sys

inp = open(sys.argv[1])     #input uncondensed info timing file
out = open(sys.argv[2], "w") #output condensed info timing file



def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)


        
#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))
info = "#pseuName  inName  outName  pseuTime  dupTime speTime dif1 dif2\n"
out.write(info)

#Make dict for the first condensation
#a dict of Br names for the first condensation with key: BrName val: [line,]
filt1Dict = {}
for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        key = "%s;%s" % (lineLst[0],lineLst[2])
        SaveIntoDict(key, lineLst, filt1Dict)
        
#Perform the first condensation and make a second dict for the second condensation
filt2Dict = {}
for key in filt1Dict:
    #get dup time values
    dupTlist = []
    for lineLst in filt1Dict[key]:
        dupTlist.append(float(lineLst[4]))

    #calculate the average
    avgDup = sum(dupTlist) / len(dupTlist)
    lineLst[4] = avgDup
    #SaveIntoDict(lineLst[0], lineLst, filt2Dict)

    #output info
    dif1 = float(lineLst[4]) - float(lineLst[3])
    dif2 = float(lineLst[5]) - float(lineLst[4])
    newLine = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
              (lineLst[0], lineLst[1], lineLst[2], lineLst[3], avgDup, \
               lineLst[5], dif1, dif2)
    out.write(newLine)    

#Perform the second condensation and output the condensed info
#It was found that there were no cases where the second condensation was necessary




inp.close()
out.close()
