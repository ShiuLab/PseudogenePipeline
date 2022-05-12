#This script is designed to take a fasta file with a huge name, replace the name
#with a code and make a reference file with the code and original name for the
#purpose of getting the original name back
#Updated on Jan 16, 2013

import sys

inp = open(sys.argv[1])      #input fasta file
out = open(sys.argv[2], "w") #output fasta file with shortened names
ref = open(sys.argv[3], "w") #2 col reference file with code name and original name

#writes the users command line prompt on the first line of the output file.
#out.write("#python %s\n" % (" ".join(sys.argv)))
#writes the users command line prompt on the first line of the ref file.
ref.write("#python %s\n" % (" ".join(sys.argv)))

cnt = 1
for line in inp:
    if line.startswith(">"):
        #make output line
        oldName = line[1:-1]
        num = str(cnt).zfill(6)
        fg = oldName.split(";")[0]
        outLine = ">Ps%s_%s\n" % (num, fg)
        out.write(outLine)

        #make reference line
        refLine = "%s\t%s\n" % (outLine[1:-1], oldName)
        ref.write(refLine)
        
        cnt += 1
        
    else:
        out.write(line)

inp.close()
out.close()
ref.close()

print("  shortened name file:", sys.argv[2])
print("  2 col reference file:", sys.argv[3])