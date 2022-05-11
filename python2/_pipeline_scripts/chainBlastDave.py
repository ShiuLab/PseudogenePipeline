from __future__ import division
import sys,os
print "INP1: Query sequence size file"
print "INP2: Parsed, Split BLAST file"
print "This script also chains ACROSS multiple contigs"

#Get the sizes of all query sequences (*.sizes file)
file1=open(sys.argv[1], 'r')
line1=file1.readline()
dict2={}
while line1:
    tab1=line1.strip().split('\t')
    if tab1[0] not in dict2:
        dict2[tab1[0]]=int(tab1[1])
    else:
        pass
    line1=file1.readline()
file1.close()


file1=open(sys.argv[2],'r')
out1=open(sys.argv[2]+".chain",'w')
out1.write('#python %s\n'%(' '.join(sys.argv)))
line1=file1.readline()
dict1={}; m=0; n=0; d=0; dictm={}
while line1:
    if line1.startswith('#'):
        if m==0:
            pass
        else:
            print "1"
            clist=[]
            for match in dict1:
                print "0"
                #Chaining the blast map
                nlist=[];total=0 
                qkeys=sorted(dict1[match].keys())
                for key in qkeys:
                    print "0000"
                    lines=dict1[match][key]
                    for line in lines:
                        print "000000"
                        sp=line.strip().split('\t')
                        #Change the following line
                        st1=int(sp[7]); end1=((int(sp[7])+int(sp[8])) - 1)
                        st2=int(sp[2]); end2=((int(sp[2])+int(sp[3])) - 1)
                        print "coords: %d\t%d\t%d\t%d\n" % (st1, end1, st2, end2)
                        
                        mlen=abs(st1-end1)+1
                        total+=mlen
                        nstr=('%s-%s:%s-%s'%(st1,end1,st2,end2))
                        nlist.append(nstr)
                nline=','.join(nlist)

                #Coverage calculations
                qsize=dict2[q]                
                cvg=(total/qsize)*100                
                clist.append(cvg)
                dictm[match]={}
                dictm[match][cvg]=('%s\t%s\t%s\t%s\t%.1f\t%s\n'%\
                           (q,qsize,match,total,cvg,nline))
                
        #Write maximum coverage contig as primary match
        maxcvg=max(clist)        
        for match in dictm:
            if maxcvg in dictm[match]:
                cline=dictm[match][maxcvg]
                out1.write(cline)
                
        for match in dictm:
            if maxcvg not in dictm[match]:
                for cvg in dictm[match]:
                    cline=dictm[match][cvg]
                    out1.write(cline)
                
        dict1={}; dictm={}; clist=[]
        out1.write('######\n')
    else:
        tab1=line1.strip().split('\t')
        #Change the following line
        q=tab1[6]; db=tab1[1]; qst=int(tab1[7]); qend=((int(tab1[7])+int(tab1[8])) - 1)
        #print "other: %s\t%s\t%d\t%d\n" % (q, db, qst, qend)
        id1=('%s-%s'%(q,db))
        if db not in dict1:
            dict1[db]={}
            dict1[db][qst]=[line1]
        else:
            if qst not in dict1[db]:
                dict1[db][qst]=[line1]
            else:
                dict1[db][qst].append(line1)
    m+=1
    line1=file1.readline()
file1.close()
out1.close()
print "Done!"
