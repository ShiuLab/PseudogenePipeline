import sys
print "INP1: Filtered BLAST output"
print "INP2: 4col file"
print "INP3: Evalue threshold for calling tandem"
print "INP4: Maximum gaps between genes for calling tandem (in terms of # of genes)"

file1=open(sys.argv[1],'r')
line1=file1.readline()
E=float(sys.argv[3])
maxL=int(sys.argv[4])
#########
def add2dict(gene1,gene2,evalue,D):
    if gene1 not in D:
        D[gene1]={}
        D[gene1][gene2]=evalue
    else:
        if gene2 not in D[gene1]:
            D[gene1][gene2]=evalue
        else:
            if evalue<D[gene1][gene2]:
                D[gene1][gene2]=evalue
    return D
#########
dict1={}; m=0
while line1:
    if line1.startswith('#'):
        pass
    else:
        tab1=line1.strip().split('\t')
        g1=tab1[0]; g2=tab1[1]; ev=float(tab1[10])
        if ev<=E:
            dict1=add2dict(g1,g2,ev,dict1)
            dict1=add2dict(g2,g1,ev,dict1)
    m+=1
    if m%50000==0:
        print "BLAST lines read: ", m
    line1=file1.readline()
file1.close()
print "No. of genes in INP1: ", len(dict1.keys())

file1=open(sys.argv[2],'r') #4col file
line1=file1.readline()
dict2={}; err=0
while line1:
    if line1.startswith('#'):
        pass
    else:
        tab1=line1.strip().split('\t')
        gene=tab1[0]; chr1=tab1[1]; st1=int(tab1[2]); end1=int(tab1[3])
        if st1<end1:
            st=st1; end=end1
        elif end1<st1:
            st=end1; end=st1
        str1=('%s|%s'%(gene,end))
        if chr1 not in dict2:
            dict2[chr1]={}
            dict2[chr1][st]=[str1]
        else:
            if st not in dict2[chr1]:
                dict2[chr1][st]=[str1]
            else:
                if str1 not in dict2[chr1][st]:
                    dict2[chr1][st].append(str1)
                    #print "Strange: ", dict2[chr1][st]                    
                    err+=1
    line1=file1.readline()
file1.close()
print "Possibly erroneous cases (same start, two genes): ", err


#Now get neighborhood of each gene and ask if a good hit is in the neighborhood
n=0
out1=open(sys.argv[1]+".tandem",'w')
out1.write('#python %s\n'%(' '.join(sys.argv)))
out1.write('#Gene1\tGene2\tDistance(# of genes)\tEvalue\n')
out2=open(sys.argv[1]+".tandem.list",'w')
out2.write('#python %s\n'%(' '.join(sys.argv)))

gcount=0; dictx={}; dicty={}
for chrom in dict2:
    #print "Processing: ", chrom
    stlist=sorted(dict2[chrom].keys())
    #print "Number of genes: ", len(stlist)    
    
    for i in range(0,len(stlist)):
        basegenes=dict2[chrom][stlist[i]]
        for basegene1 in basegenes:
            basegene=basegene1.split('|')[0]
            gcount+=1
            if basegene in dict1:
                backlist=[]; frontlist=[]
                
                #Get genes before the current gene
                j=0
                while (i-j)>=0 and j<maxL:
                    #print len(stlist), i, j, i-j, "backlist"
                    k=stlist[i-j]
                    backgenes=dict2[chrom][k]
                    for backgene1 in backgenes:
                        backgene=backgene1.split('|')[0]
                        backlist.append(backgene)
                    j+=1                
                
                #Get genes after the current gene
                j=0
                while i+j<len(stlist) and j<maxL:
                    #print len(stlist), i, j, i+j, "frontlist"
                    k=stlist[i+j]
                    frontgenes=dict2[chrom][k]
                    for frontgene1 in frontgenes:
                        frontgene=frontgene1.split('|')[0]
                        frontlist.append(frontgene)
                    j+=1
                j=0

                #Find if any genes in backlist or frontlist are significant hits
                #of the basegene
                hits=dict1[basegene].keys()
                #print hits
                #print frontlist
                #print backlist
                #sys.exit()

                #FRONT
                for m in range(0,len(frontlist)):
                    fgene=frontlist[m]
                    if fgene in hits:
                        ev=dict1[basegene][fgene]
                        strx=('%s-%s'%(basegene,fgene))
                        stry=('%s-%s'%(fgene,basegene))
                        if strx not in dictx and stry not in dictx:
                            dictx[strx]=1; dictx[stry]=1
                            out1.write('%s\t%s\t%s\t%s\n'%(basegene,fgene,m,ev))
                            n+=1
                        if basegene not in dicty:
                            dicty[basegene]=1
                            out2.write('%s\n'%(basegene))
                        if fgene not in dicty:
                            dicty[fgene]=1
                            out2.write('%s\n'%(fgene))
                            
                #BACK    
                for m in range(0,len(backlist)):
                    rgene=backlist[m]
                    if rgene in hits:
                        ev=dict1[basegene][rgene]
                        strx=('%s-%s'%(basegene,rgene))
                        stry=('%s-%s'%(rgene,basegene))
                        if strx not in dictx and stry not in dictx:
                            dictx[strx]=1; dictx[stry]=1
                            out1.write('%s\t%s\t%s\t%s\n'%(basegene,rgene,m,ev))                        
                            n+=1
                        if basegene not in dicty:
                            dicty[basegene]=1
                            out2.write('%s\n'%(basegene))
                        if rgene not in dicty:
                            dicty[rgene]=1
                            out2.write('%s\n'%(rgene))

        if gcount%10000==0:
            print "Genes completed: ", gcount, " Tandem: ", n          
out1.close()
print "No. of tandem genes: ", n, len(dicty.keys())
print "Done!"
                
                
            
                
        

            
