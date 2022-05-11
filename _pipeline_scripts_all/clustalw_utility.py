import sys,os
print "################"
print "This script requires chunks of sequence to be aligned by clustalw"
print "INP1: Fasta file containing all sequences"
print "INP2: Type of output (clustalw/pir/phylip/fasta)"
print "INP3: Type of sequence (nt/pep)"
print "INP4: Do you want to change names?"
print "BEFORE RUNNING, TYPE THIS"
print "module load ClustalW"
print "################"

file1=open(sys.argv[1], 'r')  #the fasta file
otp=sys.argv[2]
stp=sys.argv[3]
ch=sys.argv[4]

out1=open(sys.argv[1]+".%s"%otp, 'w')
if ch=='y':
    out3=open(sys.argv[1]+".%s.newnames.2col"%otp,'w')
    out3.write('#python %s\n'%(' '.join(sys.argv)))
   
if stp=='nt':
    stp1='DNA'
elif stp=='pep':
    stp1='PROTEIN'
clustCount=0; gcount=1
dict1={}; dictn1={}

tmp1=open('%s_TMP1.out'%sys.argv[1], 'w')
line1=file1.readline()
while line1:
    if line1.startswith('#python'):
        pass    
    elif line1.startswith("#"):
        tmp1.close()
        if otp=='pir':
            os.system("clustalw2 -ALIGN -INFILE=%s_TMP1.out " \
                      "-TYPE=%s -SCORE=PERCENT -OUTPUT=PIR -OUTORDER=INPUT " \
                      ">%s_out.out"%(sys.argv[1],stp1,sys.argv[1]))
            fid=('%s_TMP1.pir'%sys.argv[1])
        elif otp=='clustalw':
            os.system("clustalw2 -ALIGN -INFILE=%s_TMP1.out " \
                      "-TYPE=%s -SCORE=PERCENT -OUTPUT=CLUSTALW -OUTORDER=INPUT " \
                      ">%s_out.out"%(sys.argv[1],stp1,sys.argv[1]))
            fid=('%s_TMP1.aln'%sys.argv[1])
        elif otp=='phylip':
            os.system("clustalw2 -ALIGN -INFILE=%s_TMP1.out " \
                      "-TYPE=%s -SCORE=PERCENT -OUTPUT=PHYLIP -OUTORDER=INPUT " \
                      ">%s_out.out"%(sys.argv[1],stp1,sys.argv[1]))
            fid=('%s_TMP1.phy'%sys.argv[1])
        elif otp=='fasta':
            os.system("clustalw2 -ALIGN -INFILE=%s_TMP1.out " \
                      "-TYPE=%s -SCORE=PERCENT -OUTPUT=FASTA -OUTORDER=INPUT " \
                      ">%s_out.out"%(sys.argv[1],stp1,sys.argv[1]))
            fid=('%s_TMP1.fasta'%sys.argv[1])

        tmp2=open(fid, 'r').readlines() #pir/aln/phy/nxs/msf/gde
        for line in tmp2:
            out1.write('%s'%line)            
        out1.write('###\n')
        tmp2=""
    
        if clustCount%200==0 and clustCount>0:
            print "Clusters completed: ", clustCount
            #sys.exit()
        
        tmp1=open('%s_TMP1.out'%sys.argv[1], 'w')
        clustCount+=1; gcount=1
    else:        
        if line1.startswith('>'):
            g1=line1[1:-1]                                    
            if ch=='y':
                nname=('%s-%s'%(clustCount,gcount))
                tmp1.write('>%s\n'%nname)
                out3.write('%s\t%s\n'%(g1,nname))
            else:
                tmp1.write('>%s\n'%g1)
            gcount+=1
        else:
            tmp1.write(line1)
    line1=file1.readline()
file1.close()
out1.close()
if ch=='y':
    out3.close()
os.system('rm -f %s_TMP1* %s_out.out'%(sys.argv[1],sys.argv[1]))
if otp=='fasta':
    os.system('python /home/moghegau/scripts/fileformat_process/join_seq_into_fasta.py %s.%s'%\
              (sys.argv[1],otp))
print "All done!"
            
        
                        
