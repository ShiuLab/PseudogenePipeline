#This is the wrapper for processing pseudogene files after performing Gaurav's
#overlap pipeline
#designed by David E. Hufnagel on May 24, 2012

import sys, os

base = sys.argv[1]     #the base name for the files working on (w/o the .4col)
fasta = sys.argv[2]    #the fasta genome file name (of the database)

#you will always need to load the Repeat Masker module ahead of time
os.system("python ~/Shiu/Scripts/1_3_filterOverlapFile.py %s.4col.onlyoverlap %s.4col.onlyoverlap.filtered" % (base, base))
os.system("python ~/Shiu/Scripts/1_3_addTruePseudogenes.py %s.4col %s.4col.onlyoverlap.filtered %s.4col.annotated %s.4col.true" % (base, base, base, base))
os.system("python ~/Shiu/Scripts/GFF\ and\ 4col/1_3_condenseAnnoPseu.py %s.4col.annotated %s.4col.annotated.condensed" % (base, base))
os.system("python ~/Shiu/Scripts/GFF\ and\ 4col/1_3_condenseAnnoPseu.py %s.4col.true %s.4col.true.condensed" % (base, base))
os.system("grep -v mitochondria %s.4col.annotated.condensed | grep -v chloroplast | grep -v ATM | grep -v ATC > %s.4col.annotated.condensed.noCM" % (base, base))
os.system("grep -v mitochondria %s.4col.true.condensed | grep -v chloroplast | grep -v ATM | grep -v ATC > %s.4col.true.condensed.noCM" % (base, base))
os.system("python ~/Shiu/Scripts/1_3_pseuOrigToTrue.py %s.tblastn_parsed_G500.PE_I500.PS1.pairs.sw.out.disable_count %s.4col.true.condensed.noCM %s.fullyFiltered.disable_count" % (base, base, base))
#Repeat masking prep
os.system("python ~/Shiu/Scripts/1_2_4_RearrangeColums.py %s.4col.true.condensed.noCM 1,2,3,4 %s.real4col" % (base, base))
os.system("grep -v '#' %s.real4col > %s.real4col.mod" % (base, base))
os.system("python ~/Shiu/Scripts/FastaManagerGaurav.py -f get_stretch4 -fasta %s -coords %s.real4col.mod " % (fasta, base))
os.system("python ~/Shiu/Scripts/1_3_MakeFastaCodeName.py %s.real4col.mod.fa At %s.real4col.mod.fa.mod %s.real4col.mod.fa.ref" % (base, base, base))
#Repeat Masking
os.system("RepeatMasker -cutoff 200 -par 10 -xsmall -species viridiplantae -gff %s.real4col.mod.fa.mod" % (base))
os.system("python ~/Shiu/Scripts/parse_RepeatMasker_gff.py %s.real4col.mod.fa.mod.out 300 30" % (base))
os.system("python ~/Shiu/Scripts/NameChangers/1_3_removeCodeNamesFromRM4col.py %s.real4col.mod.fa.mod.out.300Cutoff30.0divg.4col %s.real4col.mod.fa.ref %s.real4col.mod.fa.mod.out.300Cutoff30.0divg.4col.rlNames RM4col 0" % (base, base, base))
os.system("grep -v Simple_repeat %s.real4col.mod.fa.mod.out.300Cutoff30.0divg.4col.rlNames | grep -v Low_complexity | grep -v Satellite > %s.real4col.mod.fa.mod.out.300Cutoff30.0divg.4col.rlNames.filt" % (base, base))
os.system("python ~/Shiu/Scripts/1_3_filterOutPseuRMs.py %s.real4col.mod.fa.mod.out.300Cutoff30.0divg.4col.rlNames.filt %s.fullyFiltered.disable_count %s.real4col %s.real4col.mod.fa %s.fullyFiltered.disable_count.RMfilt %s.real4col.RMfilt %s.real4col.mod.fa.RMfilt" % (base, base, base, base, base, base, base))

print "done!"
