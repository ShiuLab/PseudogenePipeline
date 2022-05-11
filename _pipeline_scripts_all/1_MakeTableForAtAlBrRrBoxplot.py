#This script was designed to combine info from 12 .size files into a file that
#makes it easy to make boxplots of this data in R.
#Created by David E. Hufnagel on Dec 1, 2012
#WARNING: HIGHLY SITUATION SPECIFIC
import sys

AtIntr = open(sys.argv[1])  #At introns .size file
AtEx = open(sys.argv[2])    #At exons .size file
AtInte = open(sys.argv[3])  #At intergenic .size file
AlIntr = open(sys.argv[4])  #Al introns .size file
AlEx = open(sys.argv[5])    #Al exons .size file
AlInte = open(sys.argv[6])  #Al intergenic .size file
BrIntr = open(sys.argv[7])  #Br introns .size file
BrEx = open(sys.argv[8])    #Br exons .size file
BrInte = open(sys.argv[9])  #Br intergenic .size file
RrIntr = open(sys.argv[10]) #Rr introns .size file
RrEx = open(sys.argv[11])   #Rr exons .size file
RrInte = open(sys.argv[12]) #Rr intergenic .size file
out = open(sys.argv[13], "w")    #output file to be used to make the boxplot in R



def TransferFile(fd, group):
    newLst = []
    for line in fd:
        if not line.startswith("#"):
            lineLst = line.strip().split("\t")
            newLine = "%s\t%s\n" % (group, lineLst[1])
            out.write(newLine)



#write info line
out.write("group\tnum\n")

#otuput the rest of the lines
group = "AtIntrons"
TransferFile(AtIntr, group)
group = "AtExons"
TransferFile(AtEx, group)
group = "AtIntergenic"
TransferFile(AtInte, group)
group = "AlIntrons"
TransferFile(AlIntr, group)
group = "AlExons"
TransferFile(AlEx, group)
group = "AlIntergenic"
TransferFile(AlInte, group)
group = "BrIntrons"
TransferFile(BrIntr, group)
group = "BrExons"
TransferFile(BrEx, group)
group = "BrIntergenic"
TransferFile(BrInte, group)
group = "RrIntrons"
TransferFile(RrIntr, group)
group = "RrExons"
TransferFile(RrEx, group)
group = "RrIntergenic"
TransferFile(RrInte, group)





AtIntr.close()
AtEx.close()
AtInte.close()
AlIntr.close()
AlEx.close()
AlInte.close()
BrIntr.close()
BrEx.close()
BrInte.close()
RrIntr.close()
RrEx.close()
RrInte.close()
out.close()
