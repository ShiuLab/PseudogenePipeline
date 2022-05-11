import sys,os
print "INP1: Base name of job (For wt1.sh.o*, say wt)"
print "INP2: Number of jobs"
prefix=sys.argv[1]
njobs=int(sys.argv[2])
pwd=os.getcwd()
dirList=os.listdir(os.getcwd())
fcount=0
for i in range(0,njobs):
    fname=('%s%i.sh.o'%(prefix,i))
    for f in dirList:
        if f.startswith(fname):
            file1=open(f,'r').readlines()
            if file1[-1]=='Done!\n':
                print "Good: ", f
            else:
                print "###"
                print "Error: ", f
                print file1
                fcount+=1
                print "###"

print "No. of jobs to redo: ", fcount
print "Done!"
