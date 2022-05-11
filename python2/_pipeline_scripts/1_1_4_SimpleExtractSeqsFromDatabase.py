#This program was designed to extract a list of names from a file
#(names seperated by '\n') and use that list to extract the sequences
#from a corresponding fasta database.  Then the program will put these
#sequences into one output file.

import sys

if len(sys.argv) == 4:
    IN1 = sys.argv[1] #input file name
    DAT = sys.argv[2] #fasta database file name
    OUT = sys.argv[3] #name of output file to be made
elif sys.argv[1] == "?":
    #explanations
    print """\nThis program was designed to extract a list of names from a file
    (names seperated by '\n') and use that list to extract the sequences
    from a corresponding fasta database.  Then the program will put these
    sequences into one output file.  The script will accept names starting with
    or without a ">" symbol

    parameters: 1)\n(endline) seperated input file name 2) database file name
    3) name of output file to be made.\n\n"""

    sys.exit()

in1 = open(IN1, "rU")
dat = open(DAT)
out = open(OUT,"w")

#get names from input
nameLst = []
for row in in1:
    name = row.strip()
    if name.startswith(">"):
        name = name[1:]
    nameLst.append(name)

#search fasta database for names and write output into new file
keep = False
for x in nameLst:
    print x
    for line in dat:
        print "in0"
        if line.startswith(">"):
            name2 = line[1:]
            print name2
            print x
            print
            print name2 == x
            if name2 == x:
                keep = True
                print "in1"
                out.write(line[:-1] + "\n")
            else:
                keep = False
        elif keep == True:
            out.write(line)
            print "in2"
    
print "done!"


in1.close()
dat.close()
out.close()
