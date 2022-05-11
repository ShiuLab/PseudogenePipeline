#This script is designed to compare a 4col and a fasta file and state what lines
#are different in both directions (in 4col and not in fasta AND in fasta and not
#in 4col).
#Created by David E. Hufnagel on Jan 9, 2013
#WARNING: This script assumes all items in each file are unique items

import sys


fourCol = open(sys.argv[1])
fasta = open(sys.argv[2])
out = open(sys.argv[3], "w")



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#make a set of names in the 4col file
fourColSet = set()
for line in fourCol:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        fourColSet.add(lineLst[0])

#make a set of names in the fasta file
fastaSet = set()
for line in fasta:
    if line.startswith(">"):
        name = line.strip("\n").strip(">")
        fastaSet.add(name)

#compare lists
both = fourColSet.intersection(fastaSet)
fourO = fourColSet.difference(fastaSet)
fastaO = fastaSet.difference(fourColSet)   

#output info
out.write("#%s names in both\n" % (len(both)))
out.write("#%s names in 4col only\n" % (len(fourO)))
out.write("#%s names in fasta only\n" % (len(fastaO)))
#out.write("###names in both:\n")
out.write("###names in 4col only:\n")
for gene in fourO:
    newLine = "%s\n" % (gene)
    out.write(newLine)
out.write("###names in fasta only:\n")
for gene in fastaO:
    newLine = "%s\n" % (gene)
    out.write(newLine)



fourCol.close()
fasta.close()
out.close()
