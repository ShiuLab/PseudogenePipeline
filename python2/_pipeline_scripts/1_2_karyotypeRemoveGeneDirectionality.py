#This script is designed for going through a karyotype file and switching the
#coordinates where the start coordinate is larger than the end coordinate.
#what this does is removes the information on gene directionality, which Circos
#cannot handle

import sys

INP = sys.argv[1] #The name of the karyotype file to be filtered
OUT = sys.argv[2] #The name of the new karyotype file (this programs output)

inp = open(INP)
out = open(OUT, "w")

#print the unix command line that called this script
print out.write('#python %s\n'%(' '.join(sys.argv)))

for line in inp:
    lineLst = line.split("\t")
    if lineLst[0] == "band":  #band lines
        if int(lineLst[4]) > int(lineLst[5]):
            temp = lineLst
            lineLst = temp[:4]
            lineLst.append(temp[5])
            lineLst.append(temp[4])
            lineLst.append(temp[6])
            line = "\t".join(lineLst)

        out.write(line)

    else:
        out.write(line)

inp.close()
out.close()
