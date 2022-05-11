#This is the 2nd alternative wrapper for processing pseudogene files without
#making the 4col file. This is different from the other 2 in that it assumes
#an intergenic file was used for the BLAST, making the overlap pipeline
#before unnecessary.
#designed by David E. Hufnagel on Jan 14, 2013

import sys, os

base = sys.argv[1]      #the base name for the files working on (w/o the .4col)
fasta = sys.argv[2]     #the fasta intergenic genome file name (of the database)
spe = sys.argv[3]       #specied identifier Ex: "Cpap"
scriptDir = sys.argv[4] #the folder containing all scripts used
scriptDir = "/" + scriptDir.strip("/") #get rid of '/' after scriptDir

#you will always need to load the Repeat Masker module ahead of time
#Repeat masking prep
os.system("python %s/1_3_annotatePseudogenes.py *.PS1.pairs.sw.out.disable_count %s.4col" % (scriptDir, base))
os.system("grep -v '#' %s.4col > %s.4col.mod" % (base, base))
os.system("rm %s.4col; mv %s.4col.mod %s.4col" % (base, base, base))
os.system("python %s/FastaManagerGaurav.py -f get_stretch4 -fasta %s -coords %s.4col" % (scriptDir, fasta, base))
os.system("python %s/1_3_MakeFastaCodeName.py %s.4col.fa %s %s.4col.fa.mod %s.4col.fa.ref" % (scriptDir, base, spe, base, base))
#Repeat Masking
os.system("RepeatMasker -cutoff 200 -par 8 -xsmall -species viridiplantae -gff %s.4col.fa.mod" % (base))
os.system("python %s/parse_RepeatMasker_gff.py %s.4col.fa.mod.out 300 30" % (scriptDir, base))
os.system("python %s/NameChangers/1_3_removeCodeNamesFromRM4col.py %s.4col.fa.mod.out.300Cutoff30.0divg.4col %s.4col.fa.ref %s.4col.fa.mod.out.300Cutoff30.0divg.4col.rlNames RM4col 0" % (scriptDir, base, base, base))
os.system("grep -v Simple_repeat %s.4col.fa.mod.out.300Cutoff30.0divg.4col.rlNames | grep -v Low_complexity | grep -v Satellite > %s.4col.fa.mod.out.300Cutoff30.0divg.4col.rlNames.filt" % (base, base))
os.system("python %s/1_3_filterOutPseuRMs.py %s.4col.fa.mod.out.300Cutoff30.0divg.4col.rlNames.filt *.PS1.pairs.sw.out.disable_count %s.4col %s.4col.fa %s.disable_count.RMfilt %s.4col.RMfilt %s.4col.fa.RMfilt" % (scriptDir, base, base, base, base, base, base))

#high confidence filter

#make code names
#os.system("python %s/TabularManager.py -f RemCodeName -inp %s.4col.RMfilt -out %s.4col.RMfilt.cdnm -col 1 -ref %s.4col.fa.ref -FR R" % (scriptDir, base, base, base))
#os.system("python %s/FastaManagerDave.py -f replace_names -fasta %s.4col.fa.RMfilt -name %s.4col.fa.ref -out %s.4col.fa.RMfilt.cdnm" % (scriptDir, base, base, base))

#cleanup
os.system("mkdir Temp")
os.system("mv %s.4col.fa.mod* Temp" % (base))
os.system("mv %s.tblastn_* Temp" % (base))
os.system("mv Temp/*count .")


print "done!"
