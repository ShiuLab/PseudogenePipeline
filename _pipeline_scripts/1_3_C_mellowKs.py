#This script takes a tabular file with Ks in it and if the number is greater
#than a set threshold it is set to equal that threshold
import sys

inp = open(sys.argv[1])      #input file
out = open(sys.argv[2], "w") #output file
thresh = float(sys.argv[3])  #threshold
ind = int(sys.argv[4])       #python index of where to take the Ks from


#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        #print
        #print lineLst
        #print lineLst[ind]
        if float(lineLst[ind]) > thresh:
            #print lineLst
            lineLst[ind] = str(thresh)
            #print lineLst
            #print
        newLine = "\t".join(lineLst) + "\n"
        out.write(newLine)


inp.close()
out.close()
