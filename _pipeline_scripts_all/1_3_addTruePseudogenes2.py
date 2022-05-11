#This script is designed to add true pseudogenes (the ones in the pseudogene
#.4col file, but not the pseudogenes .4col.onlyoverlap.filtered file) to the
#.4col.onlyoverlap.filtered file to make the .4col.annotated file
#This second version was made to get rid of the anno file and work by David E.
#Hufnagel on Jan 28, 2013

import sys

fourCol = open(sys.argv[1])   #the input full genome pseudogene .4col file
filt = open(sys.argv[2])      #the input filtered overlap file
#anno = open(sys.argv[3], "w") #the output (.4col.annotated) file with both true pseudogenes and false positives
true = open(sys.argv[3], "w") #the output (.4col.true) file with only the true pseudogenes



#go through filt and make a list of names
nameLst = []
for line in filt:
    lineLst = line.split("\t")
    nameLst.append(lineLst[0])

#print nameLst

#go through fourCol make a list of lines where the name in fourCol is not in nameLst.
interLines = []
fourCol.readline()
for line in fourCol:
    lineLst = line.split("\t")
    #print lineLst[0]
    if lineLst[0] not in nameLst:
        interLines.append(line)


#writes the users command line prompt on the first line of the output file.
#anno.write('#python %s\n'%(' '.join(sys.argv)))
true.write('#python %s\n'%(' '.join(sys.argv)))

#output data        
filt.seek(0)
for line in interLines:
    #anno.write(line)
    true.write(line)
#for line in filt:
#    anno.write(line)



fourCol.close()
filt.close()
#anno.close()
true.close()
    
        
