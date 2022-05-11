import sys
print "INP1: TAIR gff file"
print "INP2: What format output (4col/5col)"
print "INP3: What features? (gene/mRNA/exon/CDS...)"
file1=open(sys.argv[1], 'r')
line1=file1.readline()
tp=sys.argv[2].strip()
feat=sys.argv[3].strip()
out1=open(sys.argv[1]+".mygff.%s"%(feat),'w')
x=0
while line1:
    x+=1
    tab1=line1.strip().split('\t')
    chr1=tab1[0]
    o=tab1[6]
    feat1=tab1[2]
    if feat1=='mRNA_TE_gene':
        pass
    else:
        if feat in feat1:
            if feat=='gene' or feat=='mRNA':
                sp=tab1[8].split(';')
                if feat=='gene':
                    note=sp[1].split('=')
                    ftp=note[1]                
            elif feat=='CDS':
                sp=tab1[8].split(',')        
            else:
                sp=tab1[8].split(';')

            #idsp1=sp[0].split('=')[1]
            #idsp2=sp[1].split('=')[1]            
            #id1=('%s|%s'%(idsp1,idsp2))

            id1=sp[0].split(' ')[1][1:-2]

            if feat=='gene':
                id1=('%s|%s'%(id1,ftp))
            if tp=='4col':        
                if o=='+':
                    st=tab1[3]
                    end=tab1[4]
                elif o=='-':
                    st=tab1[4]
                    end=tab1[3]
                out1.write('%s\t%s\t%s\t%s\n'%(id1,chr1,st,end))
            elif tp=='5col':
                st=tab1[3]
                end=tab1[4]
                out1.write('%s\t%s\t%s\t%s\t%s\n'%(id1,chr1,o,st,end))        
    line1=file1.readline()
file1.close()

    
