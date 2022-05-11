#This script is the wrapper for running the misassembly pipeline coming directly
#from the parseGMAP wrapper (1_7_4_5_parseGMAPWrapper.py)
#Created by David E. Hufnagel on Sep 12, 2012

import sys, os

base = sys.argv[1]       # the base for the name of all the parsed gmap files  (.gmap)
genomeSize = sys.argv[2] # the genomic fasta.size file
est = sys.argv[3]        # the final est file
Z = sys.argv[4]          # the intron length to be used in the calculation
scriptDir = sys.argv[5]  # the directory containing the scripts needed
idt = sys.argv[6]        # the %ID from GMAP parsing
cov = sys.argv[7]        # the coverage from GMAP parsing
if "/" in scriptDir:
    scriptDir.strip("/")

print "base:       ", base
print "genomeSize: ", genomeSize
print "est:        ", est
print "Z:          ", Z
print "scriptDir:  ", scriptDir
print "idt:        ", idt
print "cov:        ", cov
#sys.exit()



os.system("python %s/FastaManager.py -f get_sizes -fasta %s" % (scriptDir, est))
os.system("python %s/1_6_6_longestToBlast.py %s.breaks "\
          "%s.breaks.%sidt%scovfilt.4col.repmarked.minmax.noCM.longest "\
          "%s.breaks.%sidt%scovfilt.4col.repmarked.minmax.noCM.longest.blast" \
          % (scriptDir, base, base, idt, cov, base, idt, cov))
os.system("python %s/1_6_6_MisassemblyPart1.py "\
          "%s.breaks.%sidt%scovfilt.4col.repmarked.minmax.noCM.longest.blast "\
          "%s %s.size %s.fullyParsed.MisOut1 %s" % (scriptDir, base, idt, cov, \
                                                    genomeSize, est, base, Z))
os.system("python  %s/1_6_6_MisassemblyPart2.py %s.fullyParsed.MisOut1 "\
          "%s.fullyParsed.MisOut2" % (scriptDir, base, base))
os.system("python  %s/1_6_6_MisassemblyPart3.py %s.fullyParsed.MisOut2 "\
          "%s.fullyParsed.MisOut3" % (scriptDir, base, base))
