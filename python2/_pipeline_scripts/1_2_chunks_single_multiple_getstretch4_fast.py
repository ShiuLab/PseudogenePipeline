import sys, os
print "INP1: Coords file"
print "INP2: Fasta file"
file2=open(sys.argv[2],'r')
line2=file2.readline()
dict1={}
while line2:
    if line2.startswith('>')==True:
        name=line2[1:-1]
    elif line2.startswith('#')==True:
        pass
    else:
        seq=line2[:-1]
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
        #print chrom
        #print seq
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
            #print g1,"Reverse"
            for i in range(right-1,left):
                #print i, len(seq)
                base=seq[i].upper()
                myseqlist.append(base)                
            myseq=''.join(myseqlist)
            #print myseq
        else:
            myseq='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
    else:
        print chrom, sdict.keys()
        sys.exit()
        myseq='NNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
    return myseq
##################################
out1=open(sys.argv[1]+".fa",'w')
m=0
n=0
while line1:
    if line1.startswith('#')==True:
        if m%1000==0:
            print m
        if m==0:
            pass        
        else:
            if len(list1)==0:
                pass
            
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
                    #print tseq
                    clist.append(tseq)
                    glist.append(g1)
                    
                if cl>cr:
                    orient='-'
                    nseq=''.join(clist)
                    nseq2=list(nseq)
                    nseq2.reverse()
                    fseq_rev=complement(nseq2)
                    if all_same(glist):
                        fg1=glist[0]
                    else:
                        fg1=glist[0]
                    out1.write('>%s\n%s\n'%(fg1,fseq_rev))
                    glist=[]
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
    else:
        list1.append(line1)
    m+=1
    line1=file1.readline()
file1.close()

