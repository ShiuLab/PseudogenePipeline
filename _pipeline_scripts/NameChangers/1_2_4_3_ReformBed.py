# This script is intended to make my .bed file names match my .blast file names.
# In this case I will simply Change A. Thaliana gene names from this format
# (gff):
#   "ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010"
# or this format (gtf):
#   "gene_id "CUFF.1"; transcript_id "CUFF.1.1"; FPKM "3.9214100919"; frac "1.000000"; conf_lo "2.122144"; conf_hi "5.720676"; cov "3.848458";"
# to this format: "AT1G01010" or "CUFF.1".
# Created by David E. Hufnagel 10-25-2011
# Updated: 12-19-2011: made file capable of handling bed files acquired from
# gtf files (in addition to the original support for gff)

import sys

BED = sys.argv[1]      #the name of the .bed file to be processed
OUT = sys.argv[2]      #the name of the output
fileType = sys.argv[3] #the file type (gff or gtf)
bed = open(BED)
out = open(OUT, "w")


#print the unix command line that called this script
out.write('#python %s\n'%(' '.join(sys.argv)))

for line in bed:
    lineLst = line.split("\t")

    #get rid of useless header
    if line.startswith("#"):
        continue

    #for GFF files
    if fileType == "gff" or fileType == "GFF":
        if lineLst[3].startswith("ID="):
            nameLst = lineLst[3].split(";")
            if nameLst[-1].startswith("Name="):
                name = nameLst[-1].split("=")[-1]
            else:
                print "ERROR HERE"
            newLine = "%s\t%s\t%s\t%s" % (lineLst[0], lineLst[1], lineLst[2], name)
            out.write(newLine)
        else:
            out.write(line)

    #for GTF files
    elif fileType == "gtf" or fileType == "GTF":
        nameLst = lineLst[3].split(";")
        if nameLst[0].startswith("gene_id"):
            name = nameLst[0].split(" ")[-1]
            name.strip()
            name = name[1:-1]
        else:
            print "ERROR HERE"

        #### WARNING: YOU MIGHT NEED TO COMMENT THIS AND UNCOMMENT THE .BED ONE ####
        #for use in get_stretch4 (not .bed!)
        #                               gene    chromo    start coor  stop coor
        newLine = "%s\t%s\t%s\t%s\n" % (name, lineLst[0], lineLst[1], lineLst[2])

        #for use in .bed output
        #                                 chromo     start coor   stop coor  gene
        #newLine = "%s\t%s\t%s\t%s\n" % (lineLst[0], lineLst[1], lineLst[2], name)
        out.write(newLine)

    else:
        print "ERROR HERE"
            

bed.close()
out.close()
