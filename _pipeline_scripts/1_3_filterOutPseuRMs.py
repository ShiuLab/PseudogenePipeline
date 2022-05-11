#this script was designed to take a parsed repeat mask file
#(.real4col.mod.fa.mod.out.225Cutoff30.0divg.4col)and filter out all the repeats
#from the pseudogene files (.fullyFiltered.disable_count, .real4col, .real4col.mod.fa)
#Created by David E. Hufnagel on June 8, 2012

import sys

repM = open(sys.argv[1])          #input 4col parsed repeat masking file (.real4col.mod.fa.mod.out.225Cutoff30.0divg.4col.rlNames.filt)
dis = open(sys.argv[2])           #disable_count input file (.fullyFiltered.disable_count)
four = open(sys.argv[3])          #4col format input file    (.real4col)
fasta = open(sys.argv[4])         #fasta format input file   (.real4col.mod.fa)
disOut = open(sys.argv[5], "w")   #disable_count output file (.fullyFiltered.disable_count.RMfilt)
fourOut = open(sys.argv[6], "w")  #4col format output file   (.real4col.RMfilt)
fastaOut = open(sys.argv[7], "w") #fasta format output file  (.real4col.mod.fa.RMfilt)

print("\nyour parameters:")
print("repM: ", repM)
print("dis: ", dis)
print("four: ", four)
print("fasta: ", fasta)
print("disOut: ", disOut)
print("fourOut: ", fourOut)
print("fastaOut: ", fastaOut)





#writes the users command line prompt on the first line of the output files.
disOut.write("#python %s\n" % (" ".join(sys.argv)))
fourOut.write("#python %s\n" % (" ".join(sys.argv)))
fastaOut.write("#python %s\n" % (" ".join(sys.argv)))

#go through repM and extract unique names 
repSet = set()
for line in repM:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        repSet.add(lineLst[1])

print("set: ", len(repSet))

#print repSet

#go through dis and filter to make disOut
write = True
for line in dis:
    if line.startswith("#"):
        lineLst = line.split(" ")
        name = (";".join(lineLst[:3]) + ";" + "|".join(lineLst[4:])[:-1])[1:]
        if not name in repSet:
            write = True
            disOut.write(line)
        else:
            write = False

    else:
        if write == True:
            disOut.write(line)

#go through four and filter to make fourOut
for line in four:
    lineLst = line.split("\t")
    if not lineLst[0] in repSet:
        fourOut.write(line)

#go through fasta and filter to make fastaOut
write = True
for line in fasta:
    if line.startswith(">"):
        if not line[1:-1] in repSet:
            write = True
            fastaOut.write(line)
        else:
            write = False
    else:
        if write == True:
            fastaOut.write(line)

            
        
    
    




print("Done!")
repM.close()
dis.close()
four.close()
fasta.close()
disOut.close()
fourOut.close()
fastaOut.close()
