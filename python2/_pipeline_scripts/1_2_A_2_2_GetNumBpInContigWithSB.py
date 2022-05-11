#This script is designed to determine the number of base pairs in contigs
#containing syntenic blocks.  This was a bit tricky simply because this
#requires 3 input files to get this data.
#Created by David E. Hufnagel on 3-13-2012
#WARNING: SITUATION SPECIFIC.  To modify, consider the way that the coll file is read (currently always takes the second gene in the pair)
"""Algorithm:
1) identify contigs with syntenic blocks (3 files)
    a) go through coll file and get Rr genes involved (list of genes)
    b) go through the bed file and get the contig info for the genes
       involved (list of contigs)
2) read base pair numbers for contigs and total them
    a) go through fasta file and read lenghts of contugs in 1b list
3) output data

"""
import sys

fasta = open(sys.argv[1])    #contigs fasta file
coll = open(sys.argv[2])     #.collinearity file for syntenic block information
bed = open(sys.argv[3])      #.bed file for coordinate information
out = open(sys.argv[4], "w") #output file to contain the data




#skip past the title in the .collinearity file
for l in coll:
    if not l.startswith("#"):
        break
    
#write the users command line prompt on the first line of the output file.
out.write('#python %s\n'%(' '.join(sys.argv)))

#go through coll file and get Rr genes involved (list of genes)
genes = []
for line in coll:
    if not line.startswith("##"):
        lineLst = line.split("\t")
        genes.append(lineLst[2])

#go through the bed file and get the contig info for the genes involved (list of contigs)
contigs = []
for ln in bed:
    lnLst = ln.split("\t")
    if lnLst[1] in genes:
        contigs.append(lnLst[0])

#go through fasta file and read lenghts of contugs in 1b list and output info
readIt = False
contigCnt = 0
bigCnt = 0
for lin in fasta:
    if lin.startswith(">"):  #title lines
        if lin[1:-1] in contigs:
            readIt = True
        else:
            readIt = False

        if contigCnt != 0:   #output information
            newline = "%s\t%s\n" % (lin[1:-1], contigCnt)
            out.write(newline)

        bigCnt += contigCnt
        contigCnt = 0
        
    else:                    #sequence lines
        if readIt == True:
            contigCnt += len(lin[:-1])

newln = "%s total bp\n" % (bigCnt)    
out.write(newln)







fasta.close()
coll.close()
bed.close()
out.close()
