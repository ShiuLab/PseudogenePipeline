# take the coordinates from the gene gff directory and make a dict of gene,
# chromosome, start and end coordinates.
# then scan thorugh the map file line by line. for every pair of coords
# if atleast one of coords lies WITHIN/SPANNING the range of any coord pair in the dict,
# write
# the PUTname, chr-start-end and ATgene, chr-start-end in a new file
import sys
#print("  INP1: Whole gene gff (Generally reference file like introns/CDS/Exons etc.)")
#print("  INP2: Small gff (Generally test feature like PUT/sORF etc.)")
#print("  INP3: Whether you want to look at close features (y/n)")
#print("  INP4: If INP3=y,specify bp distance (25,50)")
file1=open(sys.argv[1], 'r')  #the WHOLE Gene gff file(Bigfile)
dict1={}
line1=file1.readline()
m=0
s=0
printed='FALSE'
x=0
############################################
# READS FILE 1 AND PUTS GENE AND ITS COORDINATES IN DICT1 IN THE FOLLOWING FORMAT
#  DICT1[GENE][CHROMOSOME][START]=[END]
############################################
while line1:  #the file sohudl be in the following format
            #Gene  -- Chrname -- XX (Optional) -- Start -- End
    if line1.startswith('#')==True:
        pass
    else:
        line1=line1[:-1]
        tab1=line1.split('\t')
        gene1=tab1[0]
        if tab1[1]=='NA':
            pass
        else:
            chrn1=tab1[1]
            if len(tab1)==4:    
                lb1=int(tab1[2])  #make changes here if start-end in diff columns
                rb1=int(tab1[3])
            elif len(tab1)==5:
                lb1=int(tab1[2])  #make changes here if start-end in diff columns
                rb1=int(tab1[3])
            elif len(tab1)>5:
                if printed=='FALSE':
                    printed='TRUE'
                    print("Assuming START as tab1[2] and END as tab1[3]")
                lb1=int(tab1[2])  #make changes here if start-end in diff columns
                rb1=int(tab1[3])
                
            if gene1 not in dict1:
                dict1[gene1]={}
                dict1[gene1][chrn1]={}
                dict1[gene1][chrn1][lb1]=rb1
                s+=1
            else:
                if chrn1 not in dict1[gene1]:
                    dict1[gene1][chrn1]={}
                    dict1[gene1][chrn1][lb1]=rb1
                    s+=1
                elif lb1 not in dict1[gene1][chrn1]:            
                    dict1[gene1][chrn1][lb1]=rb1
                    s+=1
                else:            
                    #print "Strange!!", gene1
                    pass
    x+=1
    if x%50000==0:
        print("File1: ", x)
    line1=file1.readline()
file1.close()
############################################
# FILE1 NOW IS CLOSED AND ALL ITS ENTRIES ARE IN DICT1
# PROCEED WITH ANALYZING ENTRIES IN FILE2
############################################

############################################
# READS FILE 2 AND WRITES OVERLAPPING ENTRIES TO OUTPUT
############################################
file2=open(sys.argv[2], 'r') # The second file which generally contains the test features, like PUTs/sorfs

out1_name = sys.argv[2]+".onlyoverlap"
out1=open(out1_name, 'w')
#out1=open(sys.argv[2]+".completeinternal", 'w')
ch=sys.argv[3].strip()
if ch=='y':    
    dist=int(sys.argv[4])
    #out2=open(sys.argv[2]+".%sbpclose"%str(dist), 'w')
    out3=open(sys.argv[2]+".overlap%sclosemerged"%str(dist), 'w')
dict2={}
line2=file2.readline()
mdict={}
d=0
f=0
a1=0
a2=0
a3=0
a4=0
a5=0
a6=0
a7=0
a8=0
maplist=[]
while line2:  #the file should be in this format
            #  Putname -- Chrname -- Start -- End

    # 5/11/22: the true variable defined below is missing a default, so added
    true = 0
    
    if line2.startswith('#')==True:
        pass
    else:
        #DEFINE COORDINATES OF THE GENE ON THE CURRENT LINE BEING READ
        tab2=line2.strip().split('\t')
        put=tab2[0]
        nput=put
        chr1=tab2[1]
        ################################################
        if tab2[2]=='+' or tab2[2]=='-' or tab2[2]=='.':
            orient=tab2[2]
        else:
            orient="bull"
        # IF THE FILE IS A 4 COLUMN FILE
        if len(tab2)==4:
            start=int(tab2[2])  #make changes here if start-end in diff columns
            end=int(tab2[3])
            true=1
            
        # IF THE FILE IS A 5 COLUMN FILE
        elif len(tab2)==5:
            #print "Check this now!"
            #sys.exit()
            if orient!="bull":
                start=int(tab2[3])  #make changes here if start-end in diff columns
                end=int(tab2[4])
                true=1
            else:
                start=int(tab2[2])  #make changes here if start-end in diff columns
                end=int(tab2[3])
                true=1
                
                
        # IF THE FILE IS A >5 COLUMN FILE
        elif len(tab2)>5:            
            if printed=='FALSE':
                printed='TRUE'
                print("Assuming START as tab2[2] and END as tab2[3]")
            start=int(tab2[2])  #make changes here if start-end in diff columns
            end=int(tab2[3])
            true=1
    # ABOVE THIS LINE, ONE GENE AND ITS COORDINATES ARE DEFINED
    ################################################

    ############################################

    ############################################
    # FROM HERE, OVERLAP DETERMINED BY COMPARING CURRENT ENTRY IN FILE 2 AND ALL ENTRIES IN DICT1
    # CURRENT ENTRY = THE CURRENT LINE BEING READ IN FILE2
    if true==1:
        for gene in dict1:
            for chrn in dict1[gene]:
                if chrn==chr1:
                    for lb in dict1[gene][chrn]:  #left border, right border
                    #if abs(lb-start)<=100000 or abs(lb-end)<=100000:  #this would reduce the load on the system
                        rb=dict1[gene][chrn][lb]                 #25000 is the maximum gene length
                                                            #and 3000 is the max expected lengths of all introns
                                                                # in AT. Mebbe needs to be changed
                                                                #for other plants to 100000
                        
                        #CASE1: Any of the PUT coordinate within the gene, in either orientation
                        ################################################
                        '''
                        #This loop for mapping PUTs only COMPLETELY INTERNAL to a feature
                        #                st--------->end
                        #       lb---------------------------->rb
                        #               end<---------st
                        #                       OR
                        #                st--------->end
                        #       rb<----------------------------lb
                        #               end<---------st
                        
                        if lb<=start<=rb and lb<=end<=rb:
                            out1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tInternal\n'%\
                                        (put,chr1,start,end,gene,chrn,lb,rb))
                            d+=1
                            maplist.append(put)
                        elif lb>=start>=rb and lb>=end>=rb:
                            out1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tInternal\n'%\
                                        (put,chr1,start,end,gene,chrn,lb,rb))
                            d+=1
                            maplist.append(put)
                        else:
                            pass
                        '''
                        #################################################
                        '''
                        #This loop for mapping items ONLY CLOSE TO A FEATURE
                        slb=abs(start-lb)
                        srb=abs(start-rb)
                        elb=abs(end-lb)
                        erb=abs(end-rb)
                        if slb<=dist or srb<=dist or elb<=dist or erb<=dist:
                            #out1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%sbp\n'%(put,chr1,start,end,gene,chrn,lb,rb,dist))
                            out2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%sbp\n'%(put,chr1,start,end,gene,chrn,lb,rb,dist))
                            #out3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%sbp\n'%\
                            #           (put,chr1,start,end,gene,chrn,lb,rb,dist))
                            d+=1
                            a1+=1
                            if put not in mdict:
                                mdict[put]=1
                        '''
                        #################################################
                        
                        #This loop for mapping PUTs even PARTLY INTERNAL to a feature

                        #CASE1: Any of the PUT coordinate within the gene, in either orientation
                        #   ------->  OR   ---------> OR -------->
                        #       lb-------------------------->rb
                        #   <-------  OR   <--------- OR <--------
                        #                      OR
                        #   ------->  OR   ---------> OR -------->
                        #       rb<--------------------------lb
                        #   <-------  OR   <--------- OR <--------
                        if lb<=start<=rb or lb<=end<=rb or lb>=start>=rb or lb>=end>=rb:
                            if gene==nput:
                                pass
                            else:
                                out1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tOverlap\n'%\
                                            (put,chr1,start,end,gene,chrn,lb,rb))
                                if ch=='y':
                                    out3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tOverlap\n'%\
                                                (put,chr1,start,end,gene,chrn,lb,rb))
                                d+=1
                                a7+=1
                                maplist.append(put)
                                if put not in mdict:
                                    mdict[put]=1
                                    
                        else:                        
                            #Case2: PUT in forward orientation, spanning the whole gene
                            #   st------------------------------------->end
                            #       lb---------------------------->rb
                            #                   OR
                            #   st------------------------------------->end
                            #       rb<----------------------------lb
                            if lb>=start and rb<=end:
                                if gene==nput:
                                    pass
                                else:
                                    out1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tFspan\n'%\
                                                (put,chr1,start,end,gene,chrn,lb,rb))
                                    if ch=='y':
                                        out3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tFspan\n'%\
                                                    (put,chr1,start,end,gene,chrn,lb,rb))
                                    d+=1
                                    a5+=1
                                    if put not in mdict:
                                        mdict[put]=1

                                        
                            #Case3: PUT in reverse orientation spanning the whole gene
                            #   end<-------------------------------------st
                            #        lb---------------------------->rb
                            #                   OR
                            #   end<-------------------------------------st
                            #        rb<---------------------------lb                             
                            elif lb>=end and rb<=start:
                                if gene==nput:
                                    pass
                                else:
                                    out1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tRspan\n'%\
                                                (put,chr1,start,end,gene,chrn,lb,rb))
                                    if ch=='y':
                                        out3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tRspan\n'%\
                                                    (put,chr1,start,end,gene,chrn,lb,rb))
                                    d+=1
                                    a6+=1
                                    if put not in mdict:
                                        mdict[put]=1
                            
                            else:
                                if ch=='y':
                                    # Case4: PUTs close to a gene (eg: dist=50bp)
                                    # ------>50bp            OR             50bp ------>
                                    #            lb---------------------->rb
                                    # <------50bp            OR             50bp <------
                                    
                                    #########################OR############################
                                    
                                    # ------>50bp            OR             50bp ------>
                                    #            rb<----------------------lb
                                    # <------50bp            OR             50bp <------
                                    slb=abs(start-lb)
                                    srb=abs(start-rb)
                                    elb=abs(end-lb)
                                    erb=abs(end-rb)
                                    if slb<=dist or srb<=dist or elb<=dist or erb<=dist:
                                        out1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%sbp\n'%(put,chr1,start,end,gene,chrn,lb,rb,dist))
                                        #out2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%sbp\n'%(put,chr1,start,end,gene,chrn,lb,rb,dist))
                                        out3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%sbp\n'%\
                                                    (put,chr1,start,end,gene,chrn,lb,rb,dist))
                                        d+=1
                                        a1+=1
                                        if put not in mdict:
                                            mdict[put]=1
                                else:
                                    pass
                        
                        
                #################################################
                else:
                    pass

        m+=1
        if m%100==0:
            print(m, d)            
    line2=file2.readline()
############################################
out1.close()
#out2.close()
#out3.close()
'''
print "Number of sequences completely internal: ", len(mdict.keys())           
print "Total mapping : ", d
print "Overlap: ", a7
print "Forward: ", a5
print "Reverse: ", a6
if ch=='y':
    print ("%sbp: "%dist), a1
'''

print(" output:", out1_name)


            
                                           
