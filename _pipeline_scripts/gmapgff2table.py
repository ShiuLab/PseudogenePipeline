import sys
print "#########"
print "INP1: GMAP output"
print "INP2: Do you care about # of hits a transcript has? (y/n)"
print "INP3: If INP2==n, do you care about quality of match of hits a transcript has? (y/n)"
print "INP4: If INP2==n, say 0; If INP2==y, provide Coverage: "
print "INP5: If INP2==n, say 0; If INP2==y, provide Identity: "
print "INP6: If INP2==n, say 0; If INP2==y, provide Minlen: "
print "#########"

#Get all values
file1=open(sys.argv[1], 'r') #the gmap output file
hitpref=sys.argv[2]  #NHITS filter (y/n)
fpref=sys.argv[3]
cov=float(sys.argv[4]); idt=float(sys.argv[5]); minlen=int(sys.argv[6])
line1=file1.readline()
m=0; r=0; d=0; q=0; f=0

out1=open(sys.argv[1]+".%icov%iidt%iminlenLOC%s.gff"%\
          (cov,idt,minlen,hitpref.upper()), 'w')
dict1={}; adict={}; bdict={}; dict2={}

out1.write('#PUT\tHits\tPLen\t??\t%LenMatch\t%ID\tPSpan\tCSpan\tChr\tCcoord1\tCcoord2\tOrient\n')
#out2.write('#PUT\tHits\tPLen\t??\t%LenMatch\t%ID\tPSpan\tCSpan\tChr\tCcoord1\tCcoord2\tOrient\n')
while line1:
    if line1.startswith('>')==True:
        tab1=line1.split()
        #print tab1
        #sys.exit()

        #Get Transcript name
        pname=tab1[0][1:]
        if pname not in bdict:
            bdict[pname]=1
            
        #Get # of hits the transcript has
        nlist=tab1[2].split('/')
        nhits=nlist[1]
        
        #Get length of match
        llist=tab1[3].split('(')
        lenm=int(llist[0])
        chrn=tab1[4]
        
        #Get coverage of match
        vlist=tab1[5].split('(')
        v1=float(vlist[0])
        
        #Get % IDT
        v2=float(tab1[6])
        
        #Get coordinates of match
        pspan=tab1[7]
        cspan1=tab1[8]
        cspan2list=tab1[9].split(':')
        chrn1=cspan2list[0]
        cspan2=cspan2list[1]
        orient=tab1[10]
        
        #print line1
        #print pname,tab1[2],lenm,chrn,v1,v2
        #sys.exit()

        #If you want to only use transcripts with a certain number of hits
        #calue must be specified in hitlist
        if hitpref=='y': #v1=coverage #v2=%idt
            hitlist=['1']
            
            if nhits in hitlist:
                if v1>=cov and lenm>=minlen and v2>=idt:
                    out1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'% \
                    (pname,tab1[2],lenm,chrn,v1,v2,pspan,cspan1,chrn1,cspan2,orient))
                    #d+=1
                    if pname not in adict:
                        adict[pname]=1
                else:
                    #out2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'% \
                     #   (pname,tab1[2],lenm,chrn,v1,v2,pspan,cspan1,chrn1,cspan2,orient))
                     pass
            else:
                #out2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'% \
                 #       (pname,tab1[2],lenm,chrn,v1,v2,pspan,cspan1,chrn1,cspan2,orient))
                pass
            
        #If you dont care about no. of matches
        else:
            #but if you care about quality of match 
            if fpref=='y':
                if v1>=cov and lenm>=minlen and v2>=idt:
                    if pname not in adict:
                        adict[pname]=1
                        out1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'% \
                            (pname,tab1[2],lenm,chrn,v1,v2,pspan,cspan1,chrn1,cspan2,orient))
                    else:
                        out1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'% \
                            (pname,tab1[2],lenm,chrn,v1,v2,pspan,cspan1,chrn1,cspan2,orient))
                        #out3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'% \
                         #   (pname,tab1[2],lenm,chrn,v1,v2,pspan,cspan1,chrn1,cspan2,orient))
                        f+=1
                else:
                    #out2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'% \
                     #   (pname,tab1[2],lenm,chrn,v1,v2,pspan,cspan1,chrn1,cspan2,orient))
                    q+=1

            #if you dont care about quality of match as well as number of hits
            #everything will be written into the output
            else:
                out1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'% \
                        (pname,tab1[2],lenm,chrn,v1,v2,pspan,cspan1,chrn1,cspan2,orient))
                #d+=1
                
        
    m+=1
    if m%200000==0:
        print m
    line1=file1.readline()
d=0
q=0
for gene in bdict:
    d+=1
    if gene not in adict:
        q+=1
        print "Gene did not satisfy criteria, not in output file: ", gene
        
print "Total entries mapped (unfiltered): ", d
print "No. of valid entries (filtered): ", d-q
print "No. of entries that mapped, but did not satisfy criteria: ", q
print "Done!"
