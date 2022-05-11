#This script is designed to tabularize an MCSCanX output .collinearity file
#so that it can be easily parsed
#Created by David E. Hufnagel on April 8, 2013

import sys
collin = open(sys.argv[1])
out = open(sys.argv[2], "w")


#remove all instances of an item in a list
def RemoveAll(listx, item):
    while item in listx:
        listx.remove(item)


    
#Write the users command line prompt on the first line of the output file.
out.write("#python %s" % (" ".join(sys.argv)))

collin.readline();collin.readline();collin.readline();collin.readline();collin.readline();collin.readline();collin.readline();collin.readline();collin.readline()
for line in collin:
    if not line.startswith("#python"):
        #process info lines
        if line.startswith("##"):
            group = line.split(":")[0].strip("## ").replace(" ","_").strip()
            score = line.split("score=")[1].split("e_value=")[0].strip()
            e_val = line.split("e_value=")[1].split("N=")[0].strip()
            N = line.split("N=")[1].split(" ")[0].strip()
            chromos = line.strip().split(" ")[6]
            other = line.strip().split(" ")[7]
            info = "\n%s\t%s\t%s\t%s\t%s\t%s" % (group, score, e_val, N, chromos, other)
            out.write(info)

        #process pairs lines
        else:
            lineLst = line.strip().split(":")[1].split(" ")
            RemoveAll(lineLst, "")
            temp = lineLst[0].split("\t")
            RemoveAll(temp, "")
            pairs = "\t%s;%s" % (temp[0], temp[1])
            out.write(pairs)




collin.close()
out.close()
