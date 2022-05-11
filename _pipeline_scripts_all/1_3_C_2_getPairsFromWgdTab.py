#This script was designed to take the file wgdPseuBr.tab (~moghegau/projects/6_HighConfRadish/5_Pseudogenes/2_WGDpseu/2_syntenicPseu_myBRWGDpairs/wgdPseuBR.tab)
#and make a 2col pairs file from the BrPs vs. Br lines.
#Created by David E. Hufnagel on Dec 19, 2012
#WARNING: HIGHLY SITUATION SPECIFIC

import sys

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[0].startswith("PsBr"):
            brName = lineLst[0].split("_")[1]
            newLine = "%s\t%s\n" % (lineLst[0], brName)
            out.write(newLine)



inp.close()
out.close()
