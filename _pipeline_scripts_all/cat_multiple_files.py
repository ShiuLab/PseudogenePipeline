import sys, os
print "INP1: The base id of the files (with separator)"
print "INP2: How many files"
print "INP3: Name of output file"
print "INP4: Normal format or Genicmapped format (see below) (say normal or genicmapped)"
print "Eg: Normal format: bb1, bb2, bb3"
print "Eg: Genicmapped format: bb.frag1.onlyoverlap"
print "If INP4=genicmapped, INP5: Extension (eg: onlyoverlap,50bpclose etc.) Do not include dot"
fid=sys.argv[1]
count=int(sys.argv[2])
oid=sys.argv[3].strip()
tp=sys.argv[4].strip()

print ("You selected %s type"%tp)

list1=[]
for i in range(1,count+1):
    list1.append(fid)
    list1.append('%s'%str(i))
    if tp=='genicmapped':
        ext=sys.argv[5].strip()
        list1.append('.%s'%ext)
        
    list1.append(' ')
string1=''.join(list1)
print "The following files will be joined: "
print string1
os.system('cat %s >%s'%(string1,oid))

print "###########"
print "Output file created: ", oid
print "###########"
