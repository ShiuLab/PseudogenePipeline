#This script is designed to gather information about the number of pseudogenes
#with stop codons, frameshifts and neither for putting into a figure.

import sys

inp = open(sys.argv[1])      #input 4col pseudogene file (.real4col.RMfilt)
out = open(sys.argv[2], "w") #output file with info





#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

hcsCnt = 0    #high confidence stop codon count
lcsCnt = 0    #low confidence stop codon count
hcfCnt = 0    #high confidence frameshift count
lcfCnt = 0    #low confidence frameshift count
anyStop = 0   #any stop codon count
anyFrame = 0  #any frameshift count
noneCnt = 0   #none of the above count
stopOnly = 0  #any stop codon, but no frameshift count
frameOnly = 0 #any frameshift, but no stop codon count
stAndFr = 0   #any stop codon and any frameshift

for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        codeLst = lineLst[0].split(";")[3].split("|")
        #check if it is in the none group
        #print codeLst
        if codeLst[0] == "0" and codeLst[1] == "0" and codeLst[2] == "0" and codeLst[3] == "0":
            noneCnt += 1
        #if it isn't in the none group
        else:
            #hi conf stop
            if codeLst[0] != "0":
                hcsCnt += 1
            #low conf stop
            if codeLst[1] != "0":
                lcsCnt += 1
            #hi conf frameshift
            if codeLst[2] != "0":
                hcfCnt += 1
            #low conf frameshift
            if codeLst[3] != "0":
                lcfCnt += 1
            #any stop
            if codeLst[0] != "0" or codeLst[1] != "0":
                anyStop += 1
            #any frameshift
            if codeLst[2] != "0" or codeLst[3] != "0":
                anyFrame += 1
            #only stop
            if (codeLst[0] != "0" or codeLst[1] != "0") and not (codeLst[2] != "0" or codeLst[3] != "0"):
                stopOnly += 1
            #only frameshift
            if (codeLst[2] != "0" or codeLst[3] != "0") and not (codeLst[0] != "0" or codeLst[1] != "0"):
                frameOnly += 1
            #stop and frameshift
            if (codeLst[2] != "0" or codeLst[3] != "0") and (codeLst[0] != "0" or codeLst[1] != "0"):
                stAndFr += 1

out.write("High conf stop        %s\n" % (hcsCnt))
out.write("Low conf stop         %s\n" % (lcsCnt))
out.write("High conf frameshift  %s\n" % (hcfCnt))
out.write("Low conf frameshift   %s\n" % (lcfCnt))
out.write("None of the above     %s\n" % (noneCnt))
out.write("Either stop           %s\n" % (anyStop))
out.write("Either frameshift     %s\n" % (anyFrame))
out.write("Only stop             %s\n" % (stopOnly))
out.write("Only frameshift       %s\n" % (frameOnly))
out.write("stop and frameshift   %s\n" % (stAndFr))




inp.close()
out.close()
