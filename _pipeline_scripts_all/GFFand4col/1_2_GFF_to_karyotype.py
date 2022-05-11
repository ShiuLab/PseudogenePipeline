#This script is designed for extracting data from GFF files, creating a karyotype
#file and transfering the information from the GFF files to said karyotype file
#for use in the Circos software
#Warning: this script should be edited before reuse as it is specific to the situation it was created for
#Created by: David E. Hufnagel

import sys

def fixScafNums(listx):
    for x in range(len(listx)):
        listx[x] = list(listx[x]) #convert outer tuples to lists   [('name1',('1298,3964)), ('name2',('4981,6156))] --> [['name1',('1298,3964)], ['name2',('4981,6156)]]
        temp = listx[x][0].split("_")
        num = str(temp[1])
        while len(num)<4:
            num = "0" + num
        listx[x][0] =  "scaf_" + num
        #print "listx[x]: ", listx[x]
    return listx


GFF = sys.argv[1]  #GFF filename
spec = sys.argv[2] #species identifier (should be short)
BIG = sys.argv[3]  #the filename of the additional file containing the information for creating the big CDS lyrata karyotype file 
color = sys.argv[4]#the major color of the species' ideograms in Circos
gff = open(GFF)



#gather data from GFF files and store them in a list of lists
####Arabidopsis Lyrata####
big_chromo_list = []
big_gene_list = []
lineLst = []
if spec == "Al":
    ###chunks file###
    Aly_scaf_list = []
    Ldict = {}
    for line in gff:
        lineLst = line.split()
        if lineLst[2] == "CDS":
            #reverse coordinate values for "strand = -"
            if lineLst[6] == "-":
                lineLst[3], lineLst[4] = lineLst[4], lineLst[3]
            #make Ldict
            if lineLst[0] in Ldict:
                Ldict[lineLst[0]] = (Ldict[lineLst[0]][0], lineLst[4])
            else:
                Ldict[lineLst[0]] = (lineLst[3], lineLst[4])
            
            #make band lines
            ID = "_".join(lineLst[9:])
            Aly_scaf_list.append("band\tAl_%s\t%s\t%s\t%s\t%s\t%s\n" % \
                                 (lineLst[0], ID, ID, lineLst[3], lineLst[4], color))
    Llist = list(Ldict.items())
    Llist = fixScafNums(Llist)
    Llist.sort()
    #make chr lines
    Aly_chromo_list = []
    for scaf in Llist:
        Aly_chromo_list.append("chr\t-\tAl_%s\tAl_%s\t%s\t%s\tgpos25\n" % \
                               (scaf[0], scaf[0], scaf[1][0], scaf[1][1]))

    #write data into a new file
    new1 = open(GFF.split(".")[-2] + "_chunks.karyotype","w")
    for line in Aly_chromo_list:
        new1.write(line)
    for line in Aly_scaf_list:
        new1.write(line)
    new1.close()



    ###big CDS file###
    big_file = open(BIG)
    Aly_scaf_list = []
    for line in big_file:
        lineLst = line.split()
        Aly_scaf_list.append("band\tAl_%s\t%s\t%s\t%s\t%s\t%s\n" % \
                             (lineLst[1], lineLst[0], lineLst[0], lineLst[2], lineLst[3], color))
        
    #write data into a new file
    new2 = open(GFF.split(".")[-2] + "_bigCDS.karyotype","w")
    for line in Aly_chromo_list:
        new2.write(line)
    for line in Aly_scaf_list:
        new2.write(line)
    new2.close()



####Arabidopsis Thaliana####
big_chromo_list = []
big_gene_list = []
if spec == "At":
    for line in gff:
        #print "in1"
        lineLst = line.split()
        #reverse coordinate values for "strand = -"
        if lineLst[6] == "-":
            lineLst[3], lineLst[4] = lineLst[4], lineLst[3]
        #print lineLst[2]
        if lineLst[2] == "chromosome":
            #print "in2"
            ID = spec + lineLst[0][3]
            #print "\nll3", lineLst[3]
            Clen = int(lineLst[4]) - int(lineLst[3])
            Cstring = "chr\t-\t-%s\t%s\t0\t%i\t%s\n" % (ID, ID, Clen, color)
            big_chromo_list.append(Cstring)
        Gstring = ""
        if lineLst[2] == "gene":
            chromo = "At" + lineLst[0][-1]
            ID = lineLst[8].split(";")[0][3:]
            big_gene_list.append("band\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chromo, ID, ID, lineLst[3], lineLst[4], color))
            #print "name", name
            #print "chromo", chromo
            #print "ID", ID[4:]

    #write data into a new file
    new3 = open(GFF.split(".")[-2] + ".karyotype","w")
    for line in big_chromo_list:
        new3.write(line)
    for line in big_gene_list:
        new3.write(line)
    new3.close()

print "end"
