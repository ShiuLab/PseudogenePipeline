import sys,os
print "INP1: Table file"
print "INP2: No. of lines in each fragment"
print "INP3: CDS sequence file 1"
print "INP4: CDS sequence file 2"
print "INP5: Walltime in minutes (assume 2 seconds/pair)"
print "##############"
print "INP3 and INP4 should be in the same directory as the PWD"
print "##############"
file1=open(sys.argv[1],'r')
nlines=int(sys.argv[2])
fa=sys.argv[3]; pep=sys.argv[4].strip()
#jobfile=open('prank_commands.txt','w')
line1=file1.readline()
wall=sys.argv[5]
m=0; ccount=0; n=0
pwd=os.getcwd()
os.system('mkdir Cluster%i'%ccount)
os.system('cp codeml.ctl Cluster%i'%(ccount))
tmp1=open('Cluster%i/frag%i.tab'%(ccount,ccount),'w')
while line1:
    if m<nlines:
        tmp1.write(line1)
        n+=1
    else:
        #print "Lines: ", m, nlines
        tmp1.write(line1)
        tmp1.close()        
        
        #Create job
        tmp2=open('Cluster%i/jobCode%i.sh'%(ccount,ccount),'w')
        tmp2.write('#!/bin/sh -login\n\n')
        tmp2.write('#PBS -q main\n')
        tmp2.write('#PBS -l nodes=1:ppn=1,walltime=0:%s:00,mem=1gb\n'%wall)
        tmp2.write('#PBS -d %s/Cluster%i\n'%(pwd,ccount))
        tmp2.write('module load PRANK PAML exonerate\n')        
        tmp2.write('python /home/moghegau/scripts/popgenScripts/prank_utility_pairwise.py '\
                   '%s/%s %s/%s frag%i.tab\n'%(pwd,pep,pwd,fa,ccount))
#        jobfile.write('python /home/moghegau/scripts/popgenScripts/prank_utility_pairwise.py '\
#                   '%s/%s %s/%s %s/Cluster%i/frag%i.tab\n'%\
#                   (pwd,pep,pwd,fa,pwd,ccount,ccount))
        tmp2.write('\n') #IMPORTANT
        tmp2.close()
        os.system('qsub Cluster%i/jobCode%i.sh'%(ccount,ccount))
        print ("Submitted: Cluster%i/jobCode%i.sh"%(ccount,ccount))

        #Reinitialize, create new job
        m=0; ccount+=1; n+=1        
        os.system('mkdir Cluster%i'%ccount)
        print ('Making Cluster%i...'%ccount)
        os.system('cp codeml.ctl Cluster%i'%(ccount))        
        tmp1=open('Cluster%i/frag%i.tab'%(ccount,ccount),'w')
    m+=1   
    line1=file1.readline()
file1.close()
tmp1.close()

#For last cluster
#Create job
tmp2=open('Cluster%i/jobCode%i.sh'%(ccount,ccount),'w')
tmp2.write('#!/bin/sh -login\n\n')
tmp2.write('#PBS -q main\n')
tmp2.write('#PBS -l nodes=1:ppn=1,walltime=0:%s:00,mem=1gb\n'%wall)
tmp2.write('#PBS -d %s/Cluster%i\n'%(pwd,ccount))
tmp2.write('module load PRANK PAML exonerate\n')        
tmp2.write('python /home/moghegau/scripts/popgenScripts/prank_utility_pairwise.py '\
           '%s/%s %s/%s frag%i.tab\n'%(pwd,pep,pwd,fa,ccount))
#jobfile.write('python /home/moghegau/scripts/popgenScripts/prank_utility_pairwise.py '\
#           '%s/%s %s/%s frag%i.tab\n'%(pwd,pep,pwd,fa,ccount))
tmp2.write('\n') #IMPORTANT
tmp2.close()
os.system('qsub Cluster%i/jobCode%i.sh'%(ccount,ccount))
print ("Submitted: Cluster%i/jobCode%i.sh"%(ccount,ccount))
print "Total lines written to all files: ", n

print "Done!"

    
    

