#made by Gaurav Moghe

import sys,os
print "FOR HPC ONLY"
print "INP1: Base ID (eg: job)"
print "INP2: # of jobs"
print "INP3: Starting job number (eg: 0/1/2)"

base=sys.argv[1]
itern=int(sys.argv[2])
startjob=int(sys.argv[3])
cwd=os.getcwd()
dirList=os.listdir(cwd)

for i in range(startjob, itern+1):    
    fn=('%s%i.sh.e'%(base,i))    
    for item in dirList:
        if item.startswith(fn):
            file1=open(item,'r').readlines()
            finish=0
            print "Now reading: ", item
            for line in file1:                
                if line.startswith('=>> PBS: job killed'):
                    print "KILLED: ", item
                    finish=0
                elif line.startswith("Overall execution time"):
                    #print "FINISHED: ", item
                    finish=1
            '''
            if finish==0:
                print "INCOMPLETE: ", item
                print "###"
            '''
            file1=""
    
'''    
m=0
filecount=0
for item in dirList:
    if item.startswith(base):
        sp=item.split('.')
        if len(sp)==3:
            sp1=sp[2]
            if sp1.startswith('e'):
                print "Now reading: ", item
                filecount+=1
                file1=open(item,'r').readlines()
                print file1
                sys.exit()
                if 'completeness_report' in file1:
                    print item, " complete"
                else:
                    print item, " INCOMPLETE"
                    m+=1                    
            print "##"

        
print "Startjob was: ", ('%s%s.sh.e%s'%(base,'0',startjob))
print "Endjob was: ", ('%s%s.sh.e%s'%(base,itern,startjob+itern))
print "No. of files read: ", filecount
print "If you see only ## above, all error files are empty"
print "Files to redo:",m

'''    
    
    
