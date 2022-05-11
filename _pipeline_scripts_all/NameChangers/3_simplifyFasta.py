#simplifies fasta names from maker-Chr1-snap-gene-1.91-mRNA-1 protein AED:0.29 eAED:0.29 QI:27|1|1|1|1|1|2|39|510
#to maker-Chr1-snap-gene-1.91
#Created by David E. Hufnagel on July 5, 2012

import sys

inp = open(sys.argv[1])      #input fasta with big name
out = open(sys.argv[2], "w") #output fasta with shorter name




#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

for line in inp:
    if not line.startswith("#"):
        if line.startswith(">"):
            newLine = "-".join(line.split("\t")[0].split("-")[:-2]) + "\n"
            out.write(newLine)
        else:
            out.write(line)




inp.close()
out.close()
