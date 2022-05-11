#This script is designed to make a 4col and a fasta of intergenic sequences from
#a genome file and a gff file
#Created by David E. Hufnagel on Jan 5, 2012

import sys, os


genome = sys.argv[1]
gff = sys.argv[2]
scriptDir = sys.argv[3]
scriptDir = "/" + scriptDir.strip("/") #get rid of '/' after scriptDir
spe = sys.argv[4]   #a short species identifier Ex: A.thaliana --> Atha


os.system("python %s/FastaManagerGaurav.py -f get_sizes -fasta %s" % (scriptDir, genome))
os.system("python %s/GFFUtil.py -f getfeat -gff %s -col 2 -feat gene" % (scriptDir, gff))
os.system("python %s/TabularManager.py -f removeOverlap -inp %s_col2_gene -out %s_col2_gene.noOverlap -name 9 -col 4 -chroCol 1" % (scriptDir, gff, gff))
os.system("python %s/GFFUtil.py -f get_intergenic -gff %s_col2_gene.noOverlap -size %s.size" % (scriptDir, gff, genome))
os.system("python %s/Seq_simplifyInter4col.py %s_col2_gene.noOverlap.intergenic_coord %s.intergenic.4col %s" % (scriptDir, gff, gff, spe))
os.system("python %s/Seq_addProteinlessContigsToIntergenic4col.py %s.intergenic.4col %s.size %s.intergenic.4col.mod %s" % (scriptDir, gff, genome, gff, spe))
os.system("python %s/FastaManagerGaurav.py -f get_stretch4 -fasta %s -coords %s.intergenic.4col.mod" % (scriptDir, genome, gff))
os.system("rm %s.intergenic.4col" % (gff))
os.system("mv %s.intergenic.4col.mod %s.intergenic.4col" % (gff, gff))
os.system("mv %s.intergenic.4col.mod.fa %s.intergenic.fa" % (gff, gff))
os.system("rm %s_col2_gene*" % (gff))
os.system("rm %s_col2_gene.noOverlap.intergenic_coord" % (gff))
