#This script was designed to modify protein fasta names.  This script is
#expected to be dynamic and will likely be changed many times.
#Created by David E. Hufnagel on Dec 16, 2012

import sys

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")



for line in inp:
    if not line.startswith("#"):
        if line.startswith(">"):
            newName = "Smoe" + line[1:-1].strip("g")
            #newName = line[1:-1]#.strip("m.g")
            #newName = "_".join(newName.split("|")[:-1])[:]
            newName = newName.replace(".","_")
            #newName = newName.replace("-","_")
            #newName = newName.split("_")[0] + "_1_" + "_".join(newName.split("_")[1:])
            #if newName.startswith("A") or newName.startswith("E"):
            #    newName = "".join(newName.split("T"))
            #elif newName.startswith("GR"):
            #    newName = "".join(newName.split("_")[:-1])
            newLine = ">%s\n" % (newName)
            out.write(newLine)
        else:
            out.write(line)



inp.close()
out.close()
