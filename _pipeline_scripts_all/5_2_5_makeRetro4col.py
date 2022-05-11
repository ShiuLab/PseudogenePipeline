#This script is designed to take a .ps_class retropseudogene info file and an
#input 4col file of all pseudogenes and make two ouput 4col files.  One with
#ambiguous pseudogenes and one with retroPseudogenes
import sys

four = open(sys.argv[1])       #input 4col pseudogene file
info = open(sys.argv[2])       #input .ps_class file
ambi = open(sys.argv[3], "w")  #output ambiguous (might be retropseudogens) 4col file 
retro = open(sys.argv[4], "w") #output retropseudogenes 4col file



#Write the users command line prompt on the first line of the output file.
ambi.write("#python %s\n" % (" ".join(sys.argv)))
retro.write("#python %s\n" % (" ".join(sys.argv)))

#go through info and make a list of ambiguous pseudogenes and a list of retropseudogenes
ambiLst = []
retroLst = []
for line in info:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        codeLst = lineLst[10].strip("[]").split(",")
        if (codeLst[0] == "1" and codeLst[1] == "1") or (codeLst[0] == "0" and codeLst[1] == "0"):
            ambiLst.append(lineLst[0])
        elif codeLst[0] == "0" and codeLst[1] == "1":
            retroLst.append(lineLst[0])

#go though four and output lines based on what list they're in
for line in four:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[0] in ambiLst:
            ambi.write(line)
        elif lineLst[0] in retroLst:
            retro.write(line)




four.close()
info.close()
ambi.close()
retro.close()
