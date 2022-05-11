import sys, os
print "INP1: Coords file"
print "INP2: Fasta file"
print "################"
print "TYPE 1:"
print "RrC18209_p1     RrC18209        3079    3026"
print "RrC18209_p1     RrC18209        2766    1582"
print "TYPE 2:"
print "RrC18209_p1     RrC18209        2766    1582"
print "RrC18209_p1     RrC18209        3079    3026"
print "################"
file2=open(sys.argv[2],'r')
line2=file2.readline()
dict1={}
while line2:
    if line2.startswith('>')==True:
        name=line2[1:-1]
    elif line2.startswith('#')==True:
        pass
    else:
        seq=line2.strip()
        dict1[name]=seq
    line2=file2.readline()
file2.close()

#for name in dict1:
 #   print name,"Size: ",len(dict1[name])
    
file1=open(sys.argv[1], 'r')
line1=file1.readline()
m=0
n=0
list1=[]
##################################
def complement(cseq):  #taken in as a list
    cdict={'A':'T', 'T':'A', 'G':'C', 'C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}
    clist1=[]
    for item in cseq:
        if item in cdict:
            item1=cdict[item]
            clist1.append(item1)
    return ''.join(clist1)
###################################
#####################
def all_same(lst):
    return len(set(lst)) == 1
#####################

##################################
def seqextract(chrom,left,right,sdict,gene):    
    if chrom in sdict:  # direction is forward        
        seq=sdict[chrom]        
        slen=len(seq)
        myseqlist=[]
        if left<right:            
            for i in range(left-1,right):
                #print gene,chrom,left,right,i,slen,abs(left-right)
                base=seq[i].upper()                
                myseqlist.append(base)
                #print gene,left,right,i,base,myseqlist
            myseq=''.join(myseqlist)
            
        elif right<left:  # direction is reverse
            #print g1,"Reverse",right,left
            #Despite R<L, the sequence extracted is small-->large
            for i in range(right-1,left):
                #print i, len(seq)
                base=seq[i].upper()
                myseqlist.append(base)                
            myseq=''.join(myseqlist)
            #print ">>>",left, right, myseq
        else:
            #print ">>>",left, right, "NA"            
            myseq=''
    else:
        print "#####",chrom, sdict.keys()
        sys.exit()
        myseq='NNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
    coordLen=abs(left-right)+1
    seqLen=len(myseq)
    if coordLen!=seqLen:
        if seqLen!=0:
            print "CoordLen and SeqLen not matching!", gene, coordLen, seqLen        
            print chrom, left, right        
            sys.exit()
    return myseq
##################################
out1=open(sys.argv[1]+".fa",'w')
m=0
n=0; wrong=0
while line1:
    if line1.startswith('#')==True:        
        if m==0:
            pass        
        else:
            #if there is no sequence (something weird happening)
            if len(list1)==0:
                pass

            #if there is only 1 exon/CDS
            elif len(list1)==1:            
                for line in list1:
                    #print line
                    tab=line.strip().split('\t')
                    chr1=tab[1]
                    g1=tab[0]
                    if len(tab)==4:                    
                        cl=int(tab[2])
                        cr=int(tab[3])
                    elif len(tab)==5:
                        cl=int(tab[3])
                        cr=int(tab[4])
                        
                    fseq=seqextract(chr1,cl,cr,dict1,g1)
                if cl>cr:
                    orient='-'
                    nseq=fseq
                    nseq2=list(nseq)
                    nseq2.reverse()
                    fseq_rev=complement(nseq2)
                    out1.write('>%s\n%s\n'%(g1,fseq_rev))
                elif cl<cr:
                    orient='+'
                    fseq_for=fseq
                    out1.write('>%s\n%s\n'%(g1,fseq_for))

                list1=[]

            #if there are >1 exons
            else:
                n+=1
                clist=[]
                glist=[]
                #print list1
                for line in list1:
                    #print line
                    tab=line.strip().split('\t')
                    chr1=tab[1]
                    g1=tab[0]
                    if len(tab)==4:                    
                        cl=int(tab[2])
                        cr=int(tab[3])
                    elif len(tab)==5:
                        cl=int(tab[3])
                        cr=int(tab[4])     
                    tseq=seqextract(chr1,cl,cr,dict1,g1)
                    #print "#####"
                    #print tseq
                    clist.append(tseq)
                    glist.append(g1)
                    #print clist

                #If feature on negative strand and TYPE 1
                #First, REVERSE each fragment in its place to get
                #coordinates lined up in order
                #Second, COMPLEMENT
                if cl>cr:                
                    orient='-'
                    xlist=[]                    
                    for xseq in clist:
                        xseq2=list(xseq)
                        xseq2.reverse()
                        xlist.append(''.join(xseq2))
                    #print "...",xlist
                    nseq=''.join(xlist)
                    xlist=[]
                    nseq2=list(nseq)
                    #print ">>>", glist[0], nseq                    
                    fseq_rev=complement(nseq2)
                    #print ">>>", glist[0], fseq_rev
                    
                    if all_same(glist):
                        fg1=glist[0]
                    else:
                        fg1=glist[0]
                    out1.write('>%s\n%s\n'%(fg1,fseq_rev))
##                    if fseq_rev.startswith('ATG')==False:
##                        print ">",fg1, len(fseq_rev)
##                        print fseq_rev
##                        print "###########"; wrong+=1
##                        if wrong==2:                            
##                            sys.exit()
                    glist=[]
                #If feature on positive strand
                elif cl<cr:
                    orient='+'
                    fseq_for=''.join(clist)
                    if all_same(glist):
                        fg1=glist[0]
                    else:
                        fg1=glist[0]
                    out1.write('>%s\n%s\n'%(fg1,fseq_for))
                    glist=[]
                    #print ">",g1
                    #print fseq_for
                #print ">",g1,orient
                list1=[]
            #sys.exit()
    else:
        list1.append(line1)
    m+=1
    if m%2000==0:
        print m
    line1=file1.readline()
file1.close()
print "Done!"

