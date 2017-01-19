#Submit multiple files to calculon
#ask for base file name
#ask for number of files
#ask for database name
#generate jobs
#perform qsub
import sys
import os

def help1():
    print "REQ:"
    print "-fasta:  Base file name"
    print "-by:     How many fragments wanted?"
    print "-db:     Name of database"    
    print "-p:      Program (megablast/blastn/blastx)"
    print "-W:      Word size for BLAST (OPTIONAL)"
    print "         Default: blastn:11, blastx:3, megablast:28"
    print "-e:      Default 1 (OPTIONAL)"
    print "-a:      Number of processors/threads (Default: 4)"
    print "-out:    Base output file prefix"    
    
    sys.exit()

file1=no=db=hd=prg=ws=o=ev=time=mem=proc=""
for i in range(1,len(sys.argv)):
    if sys.argv[i]=='-fasta':
        file1=sys.argv[i+1]
    if sys.argv[i]=='-by':
        no=int(sys.argv[i+1])
    if sys.argv[i]=='-db':
        db=sys.argv[i+1]    
    if sys.argv[i]=='-p':
        prg=sys.argv[i+1]
    if sys.argv[i]=='-W':
        ws=sys.argv[i+1]
    if sys.argv[i]=='-out':
        o=sys.argv[i+1]
    if sys.argv[i]=='-e':
        ev=sys.argv[i+1]
    if sys.argv[i]=='-a':
        proc=sys.argv[i+1] 
out1=open('blast_commands.txt','w')
        
if file1=="" or db=="" or prg=="" or o=="":
    help1()
else:
    os.system('python ~/scripts/shinhanscripts/FastaManager.py -f divide -fasta %s -by %s'%(file1, str(no)))
    p=1
    while p<(no+1):        
        if ev!="":
            userev=ev
        else:
            ev=1

        if proc!="":
            nproc=proc
        else:
            nproc=4

        if ws=="":
            if prg=='blastn':
                word=11
            elif prg=='blastx':
                word=3
            elif prg=='megablast':
                word=28
        else:
            word=ws
                
        if prg=='megablast':
            out1.write('%s -i %s_%s -d %s -o %s%s -e %s -m 8 -W %s\n'%\
                       (prg,file1,p,db,o,p,userev,word))
        else:
            out1.write('blastall -p %s -i %s_%s -d %s -o %s%s -e %s -m 8 -W %s -a %s\n'%\
                       (prg,file1,p,db,o,p,userev,word,nproc))        
        
        p+=1
out1.close()
print "Done!"
        
    
       
        
    
