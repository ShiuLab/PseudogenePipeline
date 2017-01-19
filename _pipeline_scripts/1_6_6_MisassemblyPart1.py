#This script is part of a series designed to test the rates of putative genomic
#missasembly. More specifically, it takes the distance between the EST match and
#the end of the EST (Y), the distance between the matching region on the genome
#and the end of the contig (X) and the average intron length for the species of
#interest (Z) to calculate whether Y <= X + Z  If true, add a missasemby
#evidence to the contig, with a bit more complexity.  Part I is designed to take
#in the information from the GMAP matches file and the two .size files, calculate
#whether there is a misassebly, good (not misassembly) or ambiguous and output
#the match info in the following format:
#estName   contigName   misassembled?
#with two lines per est for the upstream and downstream
#Created by David E. Hufnagel on August 20, 2012
"""Algorithm:
1) Import the contig .size file into a dict
2) Import the EST .size file into a dict
3) Go through matches file and do all processing
    a) Get info and calculate FR, Xl, Xr, Yl, Yr
    b) Determine misassembly for the left side (upstream)
        I) If Y + Z > X set isMisassembled to "ambiguous"
        II) Elif Y + Z <= X and Y < 20% est legth, set isMisassembled to"good"
        III) Elif Y + Z <= X and Y >= 20% est legth, set isMisassembled to
             "misassembled"
    c) Output match information for left side (upstream)
    d) Determine misassembly for the right side (downstream)
        I) If Y + Z > X set isMisassembled to "ambiguous"
        II) Elif Y + Z <= X and Y < 20% est legth, set isMisassembled to"good"
        III) Elif Y + Z <= X and Y >= 20% est legth, set isMisassembled to
             "misassembled"
    e) Output match information for right side (upstream)"""


import sys

match = open(sys.argv[1])      #The input blast format GMAP match file
contigSize = open(sys.argv[2]) #The input genomic fasta.size file
estSize = open(sys.argv[3])    #The input est .size file
out = open(sys.argv[4], "w")   #The output file with est match info
Z = int(sys.argv[5])           #The intron length to be used in the calculation
print "\nmatch      :", match
print "contigSize :", contigSize
print "estSize    :", estSize
print "out        :", out
print "Z          :", Z
print




def ForwardOrReverse(start, stop):
    if int(start) < int(stop):
        return "F"
    elif int(start) > int(stop):
        return "R"
    else:
        return "***FR error!!!***"




#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#1) Import the contig .size file into a dict
contigDict = {}
for line in contigSize:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        contigDict[lineLst[0]] = lineLst[1]

#2) Import the est .size file into a dict
estDict = {}
for line in estSize:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        estDict[lineLst[0]] = lineLst[1]

#3) Go through matches file and do all processing
cnt = 0
for line in match:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        
        #a) Get info and calculate FR, Xl, Xr, Yl, Yr
        estName = lineLst[0]
        
        contigName = lineLst[1]
        FR =  ForwardOrReverse(lineLst[8], lineLst[9])  #not necessary for est, because it's always "forward"
        Yl = int(lineLst[6]) - 1
        Yr = int(estDict[estName]) - int(lineLst[7])
        if FR == "F":
            Xl = int(lineLst[8]) - 1
            Xr = int(contigDict[contigName]) - int(lineLst[9])
        elif FR == "R":
            Xl = int(lineLst[9]) - 1
            Xr = int(contigDict[contigName]) - int(lineLst[8])

        #b) Determine misassembly for the left side (upstream)
        #I) If Y + Z > X set isMisassembled to "ambiguous"
        if Yl + Z > Xl:
            isMisassembled = "ambiguous"
        elif Yl + Z <= Xl:
            #II) Elif Y + Z <= X and Y < 20% est legth, set isMisassembled to "good"
            if Yl < (0.2 * float(estDict[estName])):
                isMisassembled = "good"
            #III) Elif Y + Z <= X and Y >= 20% est legth, set isMisassembled to "misassembled"
            elif Yl >= (0.2 * float(estDict[estName])):
                isMisassembled = "misassembled"

        #c) Output match information for left side (upstream)
        newLine = "%s\t%s\t%s\n" % (estName, contigName, isMisassembled)
        out.write(newLine)

        #d) Determine misassembly for the right side (downstream)
        #I) If Y + Z > X set isMisassembled to "ambiguous"
        if Yr + Z > Xr:
            isMisassembled = "ambiguous"
        elif Yr + Z <= Xr:
            #II) Elif Y + Z <= X and Y < 20% est legth, set isMisassembled to "good"
            if Yr < (0.2 * float(estDict[estName])):
                isMisassembled = "good"
            #III) Elif Y + Z <= X and Y >= 20% est legth, set isMisassembled to "misassembled"
            elif Yr >= (0.2 * float(estDict[estName])):
                isMisassembled = "misassembled"

        #e) Output match information for right side (upstream)
        newLine = "%s\t%s\t%s\n" % (estName, contigName, isMisassembled)
        out.write(newLine)
        
        cnt += 1

print cnt

match.close()
contigSize.close()
estSize.close()
out.close()
