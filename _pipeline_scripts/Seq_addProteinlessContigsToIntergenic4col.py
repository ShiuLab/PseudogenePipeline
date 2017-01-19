#This script is designed to take a fresh "intergenic 4col file" with only
#contigs containing proteins and adding all the rest of the contigs
#Created by David E. Hufnagel on Feb 5, 2013
#WARNING: THE LAST LINE IN THE INTERGENIC 4COL FILE MUST BE THE ONE WITH THE HIGHEST INTERGENIC SEQUENCE NUMBER
import sys

four = open(sys.argv[1])     #the input intergenic 4col file missing proteinless contigs
genome = open(sys.argv[2])   #the input genomic .size file
out = open(sys.argv[3], "w") #the output intergenic 4col file with proteinless contigs
spe = sys.argv[4]


#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through four, make a list of chromosome names, and output all existing lines
inFourLst = []
for line in four:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        inFourLst.append(lineLst[1])
        out.write(line)
#Save the last attributed intergenic sequence number
else:
    lastNum = lineLst[0].split("Intergenic")[-1]
    lenNum = len(lastNum)
    lastNum = int(lastNum)

#Go through genome and output all lines where the name is not in the FourLst
cnt = 1
for line in genome:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        #if the name isn't in the list, add it to the 4col
        if lineLst[0] not in inFourLst:
            num = str(lastNum + cnt).zfill(lenNum)
            newLine = "%sIntergenic%s\t%s\t%s\t%s\n" % (spe, num, lineLst[0], "0", lineLst[1])
            out.write(newLine)
        cnt += 1

 


four.close()
genome.close()
out.close()
