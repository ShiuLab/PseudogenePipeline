#This script is designed for merging a fasta format transcript sequences and
#corresponding fasta format transcript quality values into one fastq format file

import sys

def lineMod(line):
    #this function should convert the quality values into ascii code (after adding 33)
    lineLst = line.split()
    newStr = ""
    for x in lineLst:
        x = int(x)
        x += 33
        x = chr(x)
        newStr += x
    return newStr[:-1]

FNA = sys.argv[1]
QUAL = sys.argv[2]
fna = open(FNA)

# Read the sequences in .fna file into a dict with seq name as
# key as its sequence as value.
dictx = {}
name = ""
seq = ""
for line in fna:
    if line[0] == ">":
        if seq != "":
            dictx[name] = seq
        name = line[1:].strip()
        seq = ""
    else:
        seq += line.strip()
        
dictx[name] = seq
fna.close()

#read the quality values in the .qual file that correspond with names in the
#fna based dict and add them to the sequence values as the second part of a
#tuple --> (seqence, quality values)
qual = open(QUAL)
for line in qual:
    if line[0] == ">":
        name2 = line[1:].strip()
    else:
        mod = lineMod(line)
        #form the tuple
        dictx[name2] = (dictx[name2], mod)

qual.close()

#Generate fastq file
oneLst = sys.argv[1].split(".")
twoLst = sys.argv[2].split(".")
if oneLst[0] != twoLst[0]:
    print "Warning the roots '%s' and '%s' and different!" % (oneLst, twoLst)
root = oneLst[0]
fastq = open(root + ".fastq", "w")
for i in dictx.items():
    fastq.write("@" + i[0] + "\n")
    fastq.write(i[1][0] + "\n")
    fastq.write("+\n")
    fastq.write(i[1][1] + "\n")

fastq.close()

print "Done!"


