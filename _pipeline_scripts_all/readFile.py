import sys
###################################
def readFasta(fn,D):
    #Reads FASTA files
    myfile=open(fn,'r')
    line=myfile.readline()
    m=0
    while line:
        if line.startswith('#'):
            pass            
        elif line.startswith('>'):            
            if m==0:
                myname=line.strip()[1:]
            else:
                seq=''.join(seqlist)                
                if myname not in D:
                    D[myname]=seq
                else:
                    #pass
                    print "Name repeat: ",myname
                    #sys.exit()
            myname=line.strip()[1:]; seqlist=[]
            m+=1
        else:                        
            seqlist.append(line.strip())
            m+=1
        line=myfile.readline()
    myfile.close()

    #For the final sequence
    seq=''.join(seqlist)
    if myname not in D:
        D[myname]=seq
    else:
        print "Last Name repeat: ",myname
        #sys.exit()
    return D
###################################

def readIndex(fn,D,index):
    #type(index) should be integer
    myfile=open(fn,'r')
    line=myfile.readline()
    m=0
    while line:
        if line.startswith('#'):
            pass
        else:
            tab1=line.strip().split('\t')
            #g1='_'.join(tab1[index].split('_')[1:])
            g1=tab1[index]
            if g1 not in D:
                D[g1]=[line.strip()]
            else:
                if line.strip() not in D[g1]:
                    D[g1].append(line.strip())
        line=myfile.readline()
    myfile.close()
    return D

###################################

def IndexStart(fn,index,str1,outName):
    #If indexItem.startswith(str1), write line to output
    myfile=open(fn,'r')
    line=myfile.readline()
    out=open(outName,'w')
    m=0;tc=0;pc=0
    while line:
        if line.startswith('#'):
            pass
        else:
            tab1=line.strip().split('\t')            
            g1=tab1[index]
            if g1.startswith(str1):
                out.write(g1)
                pc+=1
            tc+=1
        line=myfile.readline()
    myfile.close(); out.close()
    return D,pc,tc
###################################

def readTwoIndexes(fn,D,index1,index2):
    #type(index) should be integer
    myfile=open(fn,'r')
    line=myfile.readline()
    m=0
    while line:
        if line.startswith('#'):
            pass
        else:
            tab1=line.strip().split('\t')
            #g1='_'.join(tab1[index].split('_')[1:])
            g1=tab1[index1]; g2=tab1[index2]
            if g1 not in D:
                D[g1]=[g2]
            else:
                if g2 not in D[g1]:
                    D[g1].append(g2)
        line=myfile.readline()
    myfile.close()
    return D

###################################

def writeHead(OUT):
    OUT.write('#python %s\n'%(' '.join(sys.argv)))
    return OUT

###################################


    


    
