import sys, os
print "INP1: Raxinp file containing sequences in PHYLIP format separated by ##"
print "INP2: job/commandline"
print "If the following files dont exist, please make them in CWD:"
print "mkdir _TMP_PHYFILES _FINAL_RESULTS _TMP_INFOFILES _TMP_LOGFILES _TMP_PARSIMONY _TMP_BESTTREE"

file1=open(sys.argv[1],'r')
#scount=int(sys.argv[2]); ecount=int(sys.argv[3])
ch=sys.argv[2]

if ch=='commandline':
    out1=open(sys.argv[1]+".raxout",'w')
    os.system('mkdir _TMPFILES')
elif ch=='job':
    jobfile=open('raxml_commands.txt','w')
    
line1=file1.readline()
list1=[]
pwd=os.getcwd()
ccount=1; donecount=0
tmp1=open('tmpRAX%i.phy'%ccount,'w')

while line1:        
    if line1.startswith('#python'):
        pass
    elif line1.startswith('##'):        
        tmp1.close()            
        if ch=='commandline':
            os.system('/home/moghegau/packages/RAxML-7.2.6/raxmlHPC-PTHREADS '\
                      '-f d -m PROTGAMMAJTT -T 6 -s tmpRAX%i.phy '\
                      '-n tmpRAX%i.phy.raxout' %(ccount,ccount))
        
            raxout=open('RAxML_result.tmpRAX%i.phy.raxout'%ccount,'r')
            for line in raxout:
                out1.write(line)
            out1.write('######\n')
            os.system('mv RAxML_result.tmpRAX%i.phy* tmpRAX* _TMPFILES\n'%(ccount))
            
        elif ch=='job':
            jobfile.write('python /home/moghegau/scripts/projects/2_RadishGenome/'\
                          '5_Orthology/run_raxmlJobs.py %s\n'%(ccount))
        ccount+=1
        if ccount%4000==0:
            print "Clusters written: ", ccount
        tmp1=open('tmpRAX%i.phy'%ccount,'w')            
    else:        
        tmp1.write(line1)    
    line1=file1.readline()
file1.close()
print "Total clusters written: ", ccount
print "Done!"
        
