#This is the alternative wrapper for processing pseudogene files after making
#the 4col file
#designed by David E. Hufnagel on May 24, 2012
#updated on Jan 23, 2013

import sys, os

base = sys.argv[1]      #the base name for the files working on (w/o the .4col)
fasta = sys.argv[2]     #the fasta genome file name (of the database)
spe = sys.argv[3]       #specied identifier Ex: "Cpap"
scriptDir = sys.argv[4] #the folder containing all scripts used
scriptDir = "/" + scriptDir.strip("/") #get rid of '/' after scriptDir
intronLen = sys.argv[5] #95th %ile intron len to be used for hiConf filter
genes4col = sys.argv[6] #4col genomic genes file
repCut = sys.argv[7]    #RepeatMasker cutoff parameter
repDiv = sys.argv[8]    #RepeatMasker divergence parameter

#you will always need to load the Repeat Masker module ahead of time
#filter out pseudogene predictions overlapping with genes
os.system("python %s/1_3_addTruePseudogenes2.py %s.4col %s.4col.onlyoverlap %s.4col.true" % (scriptDir, base, base, base))
os.system("python %s/1_3_pseuOrigToTrue.py *.PS1.pairs.sw.out.disable_count %s.4col.true %s.fullyFiltered.disable_count" % (scriptDir, base, base))
#to not filter use these lines instead
#os.system("cp %s.4col %s.4col.true" % (base, base))
#os.system("cp %s.tblastn_parsed_G500.PE_I500.PS1.pairs.sw.out.disable_count %s.fullyFiltered.disable_count" % (base, base))

#Repeat masking prep
os.system("python %s/FastaManagerGaurav.py -f get_stretch4 -fasta %s -coords %s.4col.true " % (scriptDir, fasta, base))
os.system("python %s/1_3_MakeFastaCodeName.py %s.4col.true.fa %s %s.4col.true.fa.mod %s.4col.true.fa.ref" % (scriptDir, base, spe, base, base))

#Repeat Masking
os.system("grep -v '#' %s.4col.true.fa.mod > %s.4col.true.fa.mod.mod" % (base, base))
os.system("RepeatMasker -cutoff 200 -par 8 -xsmall -species viridiplantae -gff %s.4col.true.fa.mod.mod" % (base))
os.system("python %s/parse_RepeatMasker_gff.py %s.4col.true.fa.mod.mod.out %s %s" % (scriptDir, base, repCut, repDiv))
os.system("python %s/NameChangers/1_3_removeCodeNamesFromRM4col.py %s.4col.true.fa.mod.mod.out.%sCutoff%s.0divg.4col %s.4col.true.fa.ref %s.4col.true.fa.mod.mod.out.%sCutoff%s.0divg.4col.longNames RM4col 0" % (scriptDir, base, repCut, repDiv, base, base, repCut, repDiv)) #I need this just so I don't have to modify the next script (I know it's bad programming, but it was a time saver)
os.system("grep -v Simple_repeat %s.4col.true.fa.mod.mod.out.%sCutoff%s.0divg.4col.longNames | grep -v Low_complexity | grep -v Satellite > %s.4col.true.fa.mod.mod.out.%sCutoff%s.0divg.4col.longNames.filt" % (base, repCut, repDiv, base, repCut, repDiv))
os.system("python %s/1_3_filterOutPseuRMs.py %s.4col.true.fa.mod.mod.out.%sCutoff%s.0divg.4col.longNames.filt %s.fullyFiltered.disable_count %s.4col.true %s.4col.true.fa %s.fullyFiltered.disable_count.RMfilt %s.4col.true.RMfilt %s.4col.true.fa.RMfilt" % (scriptDir, base, repCut, repDiv, base, base, base, base, base, base))

#cleanup
os.system("Temp")
#os.system("mv %s.tblastn_* Temp" % (base))
os.system("mv %s.4col.true.fa.mod.mod* Temp" % (base))
os.system("mv Temp/*count .")

#do high confidence filter
os.system("python %s/FastaManagerDave.py -f get_sizes -fasta %s" % (scriptDir, fasta))
os.system("python %s/1_3_verifyPseudogenes.py %s.fullyFiltered.disable_count.RMfilt  %s.4col.true.RMfilt %s.4col.true.fa.RMfilt %s.size %s %s.fullyFiltered.disable_count.RMfilt.hiConf %s.4col.true.RMfilt.hiConf %s.4col.true.fa.RMfilt.hiConf %s" % (scriptDir, base, base, base, fasta, genes4col, base, base, base, intronLen))

#bring back to code names %s.4col.true.RMfilt %s.4col.true.fa.mod.mod.RMfilt
os.system("python %s/TabularManager.py -f RemCodeName -inp %s.4col.true.RMfilt.hiConf -out %s.4col.true.RMfilt.hiConf.cdnm -col 1 -ref %s.4col.true.fa.ref -FR R" % (scriptDir, base, base, base))
os.system("python %s/FastaManagerDave.py -f replace_names -fasta %s.4col.true.fa.RMfilt.hiConf -name %s.4col.true.fa.ref -out %s.4col.true.fa.RMfilt.hiConf.cdnm" % (scriptDir, base, base, base))

print "Pseudogenes acquired!"


