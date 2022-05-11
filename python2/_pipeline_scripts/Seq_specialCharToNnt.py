#This script is designed to replace special nucleotide characters to N in a
#fasta format file
#Created By David E. Hufnagel on Aug 21, 2012

import sys

inp = open(sys.argv[1])      #input unmodified fasta file with special chars
out = open(sys.argv[2], "w") #output modified fasta file with only Ns and standard nucleotides




#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Do the actual processing
spCnt = 0
cnt = 0
for line in inp:
    if not line.startswith("#"):
        if line.startswith(">"):
            out.write(line)
        else:
            seq = line.strip()
            newseq = ""
            for bp in seq:
                if bp in "ATCGatcg":
                    newseq += bp
                else:
                    newseq += "N"
                    spCnt += 1
            out.write(newseq + "\n")

            if not cnt % 100000:
                print cnt
            cnt += 1

print "%s special characters replaced" % (spCnt)
                



inp.close()
out.close()
