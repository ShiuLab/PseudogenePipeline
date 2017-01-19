import sys
print "INP1: GMAP output file"
file1=open(sys.argv[1], 'r')
out1=open(sys.argv[1]+".breaks", 'w')
#out2=open(sys.argv[1]+".breaks.4col", 'w')
line1=file1.readline()
out1.write('#Name\tNHits\tPUTlen\tCOV\tIDT\tPSpan\tPmlen\tLocus\tChr1:St..End\tCmlen\n')
dict1={}
m=0
n=0
pdict={}
cdict={}
fdict={}
while line1:
    if line1.startswith('>')==True:
        out1.write('##\n')
        if m==0:
            tab1=line1.split(' ')
            name=tab1[0][1:]
            if name not in fdict:
                fdict[name]=1
                #out1.write('##\n')
            #get chr name
            chrlist=tab1[9].split(':')
            chrom1=chrlist[0]
            #get number of hits
            hitlist=tab1[2]
            m+=1
            #get plen
            plen1=tab1[3].split('(')
            plen=plen1[0]
            #get location
            loc=tab1[9].split(':')
            chr1=loc[0]
            span=loc[1].split('..')
            cl=span[0]
            cr=span[1]
            loc1=tab1[9]
            #get COV
            cvrg=tab1[5].split('(')
            cv=cvrg[0]
            #get IDT
            idt=tab1[6]
            
        else:
            plist=dict1[name].keys()
            stlist=[]
            stdict={}
            for i in range(0,len(plist)):
                pid=plist[i]
                sp=pid.split('-')
                st=int(sp[0])
                stlist.append(st)
                stdict[st]=pid
            #print stlist
            stlist.sort()
            #print stlist                
            for st in stlist:
                pid=stdict[st]
                #print pid
                #sys.exit()
                cid=dict1[name][pid]
                cdiff=cdict[cid]
                pdiff=pdict[pid]
                out1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s:%s\t%s\n'% \
                           (name,hitlist,plen,cv,idt,pid,pdiff+1,loc1,chrom1,cid,cdiff+1))
                csplit=cid.split('-')
                #out2.write('%s\t%s\t%s\t%s\n'%(name,chrom1,csplit[0],csplit[1]))
            
            cdict={}
            pdict={}
            dict1={}
            dict2={}
            tab1=line1.split(' ')
            name=tab1[0][1:]
            if name not in fdict:
                fdict[name]=1
                #out1.write('##\n')            
            chrlist=tab1[9].split(':')
            chrom1=chrlist[0]
            hitlist=tab1[2]
            plen1=tab1[3].split('(')
            plen=plen1[0]
            #get location
            loc=tab1[9].split(':')
            loc1=tab1[9]
            chr1=loc[0]
            span=loc[1].split('..')
            cl=span[0]
            cr=span[1]
            #get COV
            cvrg=tab1[5].split('(')
            cv=cvrg[0]
            #get IDT
            idt=tab1[6]

            m+=1
    elif line1.startswith('#')==True:
        pass
    else:
        m+=1
        msplit=line1.strip().split(' ')
        clist=[]
        plist=[]
        try:
            chr1=int(msplit[0])
            chr2=int(msplit[1])
            best=1
        except:
            print line1
            print msplit
            print "Error\n##########"
            best=0

        if best==1:
            cdiff=abs(chr1-chr2)
            clist.append(msplit[0])
            clist.append('-')
            clist.append(msplit[1])        
            cid=''.join(clist)
            cdict[cid]=cdiff
            
            put1=int(msplit[2])
            put2=int(msplit[3])
            pdiff=abs(put1-put2)
            plist.append(msplit[2])
            plist.append('-')
            plist.append(msplit[3])
            pid=''.join(plist)
            pdict[pid]=pdiff
            #print pdict

            if name not in dict1:
                dict1[name]={}
                dict1[name][pid]=cid
            else:
                if pid not in dict1[name]:
                    dict1[name][pid]=cid
        
    #m+=1
    n+=1
    line1=file1.readline()
# for last PUT
cid=dict1[name][pid]
cdiff=cdict[cid]
pdiff=pdict[pid]
out1.write('##\n')
out1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s:%s\t%s\n'% \
            (name,hitlist,plen,cv,idt,pid,pdiff+1,loc1,chrom1,cid,cdiff+1))
out1.write('##\n')
#csplit=cid.split('-')
#out2.write('%s\t%s\t%s\n'%(chrom1,csplit[0],csplit[1]))
               
                        
        
        
    
