# IMPORTS
import sys

# MAIN
print'''
Takes one input (1) a FASTA formt file. 
Note: This script is primarily meant to cut the contig/gene names found in the comment line of a FASTA file ("> ...") so
that the match what is reported by BLAST. Hence, it reduces the name to what comes before the first space.
'''

source = open(sys.argv[1],'r')
output = open(sys.argv[1] + ".truncated",'w')

for line in source:
    if line.startswith(">"):
        split_line = line.strip().split(" ")
        output.write(split_line[0]+"\n")
    else:
        output.write(line) 
output.close()
