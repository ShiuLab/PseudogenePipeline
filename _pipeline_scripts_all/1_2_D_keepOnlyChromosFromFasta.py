#This script is designed to take a genomic fasta file and remove things that
#aren't chromosomes
#Created by David E. Hufnagel on May 22, 2012
#Note: fasta files were first simplified to include only genomic chromosomes.
#      No mitcochonria, chloroplast or random unnassembled pieces

import sys

inp = open(sys.argv[1])       #Input fasta file which has everything
num = int(sys.argv[2])        #The number of chromosomes in the organism
out = open(sys.argv[3], "w")  #Output fasta file which only has chromosomes



cnt = 0
for line in inp:
    if cnt <= num:
        if line.startswith(">"):
            cnt += 1
            if cnt <= num:
                out.write(line)
        else:
            out.write(line)
    else:
        break


inp.close()
out.close()
