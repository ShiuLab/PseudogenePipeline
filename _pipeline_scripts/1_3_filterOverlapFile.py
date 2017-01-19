#This script is disigned to filter the crap out of the overlap file between a
#4 column pseudogene file and basically a GFF file.  What does this mean?  see
#below:
"""
AT5G41790;Chr5|16727709-16732388;68-194:381-1;0|0|0|0	Chr5	16727715	16732388	AT5G41790|protein_coding_gene	Chr5	16727530	16732847	Overlap
-->
AT5G41790;Chr5|16727709-16732388;59-124:1-198;0|0|0|0	Chr5	16727715	16732388        protein_coding_gene"""
#Designed by David E. Hufnagel on 5-24-2012

import sys

inp = open(sys.argv[1])      #original input file
out = open(sys.argv[2], "w") #filtered output file


for line in inp:
    lineLst = line.split("\t")
    #newFour = lineLst[4].split("|")[1]  #SPECIFICALLY FOR A. THALIANA
    newFour = "protein_coding_gene"
    newLine = "%s\t%s\t%s\t%s\t%s\n" % \
              (lineLst[0], lineLst[1], lineLst[2], lineLst[3], newFour)
    out.write(newLine)





inp.close()
out.close()
