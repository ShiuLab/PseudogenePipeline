#This script is designed to get both nucleotide and amino acid fasta files from
#fully filtered and repeat masked pseudogenes using code names from .ref files
#in the format >PsAt000001_AT4G37685
#Created by David E. Hufnagel on June 22, 2012
#Updated on Mar 19, 2013 to be only for aa's (nts are a bit more complicated and should be done through 1_3_1_C_getPseudoCDSWrapper.py
"""Algorithm:
1) Go through ref and get a dictionary of big names(key) and code names(value)
2) Go through dis and make the amino acid fasta file [making a list of bigNames
   at the same time]
3) Go through fasta and make the nucleotide fasta file from the names in the
   list from step 2 and it's coded counterpart in the refDict"""


import sys

dis = open(sys.argv[1])        #the input pseudogene disable_count file with big names (.fullyFiltered.disable_count.RMfilt)
fasta = open(sys.argv[2])      #the input pseudogene nt fasta file with big names      (.real4col.mod.fa)
ref = open(sys.argv[3])        #the input reference file with code names and big names (.real4col.mod.fa.ref)
#outnt = open(sys.argv[4], "w") #the output nucleotide fasta file                       (.fullyFiltered.RMfilt.nt.fa)
outaa = open(sys.argv[4], "w") #the output amino acid fasta file                       (.fullyFiltered.RMfilt.aa.fa)

print
print "dis:   ", dis
print "fasta: ", fasta
print "ref:   ", ref
#print "outnt: ", outnt
print "outaa: ", outaa
print





#Write the users command line prompt on the first line of the output file.
#outnt.write("#python %s\n" % (" ".join(sys.argv)))
outaa.write("#python %s\n" % (" ".join(sys.argv)))

#1) Go through ref and get a dictionary of big names(key) and code names(value)
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        refDict[lineLst[1][:-1]] = lineLst[0]

goodBigNames = []
#2) Go through dis and make the amino acid fasta file [making a list of bigNames
#   at the same time]
cnt = 1    #a count to get rid of the query sequence line and the blank lines
for line in dis:
    if not line.startswith("#python"):
        if line.startswith("#"):
            cnt = 1
            #set variables
            lineLst = line.split(" ")
            name = (";".join(lineLst[:3]) + ";" + "|".join(lineLst[4:])[:-1])[1:]
            prot = lineLst[0][1:]
            newName = ">%s\n" % (refDict[name])

            #make outaa
            outaa.write(newName)
            
            #make goodBigNames
            goodBigNames.append(name)
        else:
            if cnt == 3: #to get rid of the query sequence line and the blank lines
                outaa.write(line)
        cnt += 1
    

###3) Go through fasta and make the nucleotide fasta file from the names in the
###   list from step 2 and it's coded counterpart in the refDict
##write = False
##for line in fasta:
##    if not line.startswith("#"):
##        if line.startswith(">"):
##            name = line[1:-1]
##            if name in goodBigNames:
##                write = True
##                prot = name.split(";")[0]
##                newName = ">Ps%s_%s\n" % (refDict[name], prot)
##                outnt.write(newName)
##            else:
##                write = False
##        else:
##            if write == True:
##                outnt.write(line)

    




dis.close()
fasta.close()
ref.close()
#outnt.close()
outaa.close()


